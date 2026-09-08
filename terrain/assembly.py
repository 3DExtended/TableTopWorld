"""Whole-flower assembly: terrain surface + base plate, welded into one solid.

Written fresh for the explicit-mesh redesign - the old terrain/assembly.py
(CSG-based) was deleted in Phase 4 along with the rest of the discarded
pipeline (see the project plan's Phase 4 note).
"""

from __future__ import annotations

import numpy as np
import trimesh

from terrain.base_plate import build_floor_cap
from terrain.boundary_noise import hash_to_unit_interval
from terrain.constants import BASE_PLATE_DEPTH_MM
from terrain.layout import FlowerLayout
from terrain.magnets import DEFAULT_MAGNET_BORES, MagnetBores
from terrain.standability import check_min_standable_hexes
from terrain.surface_mesh import build_flower_open_solid
from terrain.tileset import Tileset


def pick_standable_hexes(
    seed: int,
    *,
    min_standable: int,
    forced_count: int | None = None,
    hex_count: int = FlowerLayout.HEX_CELL_COUNT,
) -> list[int]:
    """Which of a flower's hex cells get a flat, standable plateau.

    How MANY varies per flower rather than being all 7 every time, so a
    map has a mix of open, buildable flowers and broken/rugged ones. The
    count is drawn from the flower's own seed in
    [min_standable, hex_count], so it stays a pure function of authored
    data - two independent builds of the same flower agree, same as every
    other seeded choice in this codebase. `forced_count` pins it instead
    (0 = no plateaus at all, i.e. fully organic terrain).

    Uses boundary_noise's splitmix hash, never Python's built-in hash(),
    which is salted per process by PYTHONHASHSEED and would silently
    break that reproducibility.
    """
    if forced_count is None:
        span = hex_count - min_standable
        count = min_standable + (
            int(hash_to_unit_interval(seed, 0xF1A7) * (span + 1)) if span > 0 else 0
        )
    else:
        count = forced_count
    count = max(0, min(hex_count, count))
    # Deterministic shuffle: order the cells by a per-cell hash, take the
    # first `count`. Sorted back into index order so callers see a stable,
    # readable set rather than hash order.
    ranked = sorted(
        range(hex_count), key=lambda h: hash_to_unit_interval(seed, 0x5A5A, h)
    )
    return sorted(ranked[:count])


def build_flower_mesh(
    tileset: Tileset,
    flower_id: str,
    *,
    subdivisions_per_edge: int = 8,
    jitter_amplitude: float = 0.05,
    xy_jitter_mm: float = 1.0,
    interior_relief_mm: float = 1.0,
    include_hex_grooves: bool = True,
    groove_depth_mm: float = 1.5,
    groove_width_mm: float = 6.0,
    groove_profile: list[tuple[float, float]] | None = None,
    standable_hexes: int | None = None,
    plate_depth_mm: float = BASE_PLATE_DEPTH_MM,
    magnet_bores: MagnetBores | None = DEFAULT_MAGNET_BORES,
) -> trimesh.Trimesh:
    """Build one flower's complete printable solid: the terrain surface
    (decisions #1-#10) whose walls run straight down to the print bed at
    z = -plate_depth_mm (the flat "basement" of decision #13, physical
    millimetres - never scaled), closed by a floor cap, with one blind
    magnet bore per silhouette edge cut into those walls (terrain/magnets.py).
    One Trimesh, no boolean union anywhere: the floor shares the walls'
    own bottom-rim vertex indices, and the bore patches share the wall
    quads' rim vertices (see terrain/base_plate.py for why that matters).

    A flower's declared roads (tileset.py's FlowerDef.roads) are passed
    straight through as road_water_side_pairs; only one road per flower is
    currently supported (surface_mesh.py raises NotImplementedError for
    more than one), and combining a road with include_hex_grooves on the
    same flower isn't wired up yet either (same limitation, unchanged from
    Phase 5).
    """
    flower = tileset.flowers[flower_id]
    layout = tileset.layout()
    level_z = tileset.meta.heights.z
    if magnet_bores is not None and plate_depth_mm < magnet_bores.min_plate_depth_mm:
        raise ValueError(
            f"plate_depth_mm={plate_depth_mm} is too shallow for the magnet bores: "
            f"need at least {magnet_bores.min_plate_depth_mm:.2f} mm (bore top at "
            f"{magnet_bores.top_above_bed_mm:.2f} mm plus a {magnet_bores.roof_mm} mm "
            "roof up to the level-0 surface)"
        )
    bottom_z = -plate_depth_mm

    road_water_side_pairs = flower.roads + flower.water

    solid = build_flower_open_solid(
        flower.side_corner_heights,
        layout,
        level_z,
        flower.seed,
        bottom_z=bottom_z,
        subdivisions_per_edge=subdivisions_per_edge,
        jitter_amplitude=jitter_amplitude,
        xy_jitter_mm=xy_jitter_mm,
        interior_relief_mm=interior_relief_mm,
        include_hex_grooves=include_hex_grooves and not road_water_side_pairs,
        groove_depth_mm=groove_depth_mm,
        groove_width_mm=groove_width_mm,
        groove_profile=groove_profile,
        hex_height_levels={
            hex_idx: flower.hexes[str(hex_idx)].height_level
            for hex_idx in pick_standable_hexes(
                flower.seed,
                min_standable=tileset.meta.min_standable_hexes,
                forced_count=standable_hexes,
            )
        },
        road_water_side_pairs=road_water_side_pairs,
        magnet_bores=magnet_bores,
    )

    all_faces = solid.faces + build_floor_cap(solid.top_triangles, solid.top_count)

    return trimesh.Trimesh(
        vertices=np.array(solid.vertices), faces=np.array(all_faces), process=True
    )


def standability_report(
    mesh: trimesh.Trimesh, tileset: Tileset
) -> tuple[bool, int, int]:
    """(meets_minimum, actual_standable_count, min_required)."""
    layout = tileset.layout()
    ok, count = check_min_standable_hexes(
        mesh, layout, tileset.meta.min_standable_hexes
    )
    return ok, count, tileset.meta.min_standable_hexes


def build_preview_mesh(tileset: Tileset, **build_kwargs) -> trimesh.Trimesh:
    """Every flower in tileset.meta's preview_map, each built standalone
    via build_flower_mesh() and placed at its correct world position via
    FlowerLayout.flower_grid_to_xy() (the fix for the placement bug
    flagged since Phase 1 - the old axial_to_xy()/flower_center_spacing
    did not correspond to true edge-sharing adjacency).

    This is a preview SCENE, not one continuously-welded solid: each
    flower is its own separate, independently-watertight print (they mate
    physically via magnets, not shared mesh topology), simply translated
    into the position it would occupy on the table. A placement's `rot`
    is applied as a rigid rotation about the flower's own center for
    visual variety - tileset.py's validate_tileset() only verifies the
    reversed-declaration side contract for rot=0 placements, so a rotated
    flower's fit against its neighbors isn't guaranteed correct here.
    """
    layout = tileset.layout()
    parts: list[trimesh.Trimesh] = []
    for placement in tileset.preview_map:
        mesh = build_flower_mesh(tileset, placement.id, **build_kwargs)
        if placement.rot:
            mesh = mesh.copy()
            mesh.apply_transform(
                trimesh.transformations.rotation_matrix(
                    placement.rot * (np.pi / 3), [0, 0, 1]
                )
            )
        x, y = layout.flower_grid_to_xy(*placement.at)
        mesh.apply_translation((x, y, 0.0))
        parts.append(mesh)
    return trimesh.util.concatenate(parts)
