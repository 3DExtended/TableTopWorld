"""Whole-flower assembly: terrain surface + base plate, welded into one solid.

Written fresh for the explicit-mesh redesign - the old terrain/assembly.py
(CSG-based) was deleted in Phase 4 along with the rest of the discarded
pipeline (see the project plan's Phase 4 note).
"""

from __future__ import annotations

import numpy as np
import trimesh

from terrain.base_plate import build_base_plate_parts
from terrain.constants import (
    BASE_PLATE_DEPTH,
    MAGNET_CENTER_Z,
    MAGNET_RADIUS,
)
from terrain.standability import check_min_standable_hexes
from terrain.surface_mesh import build_flower_open_solid
from terrain.tileset import Tileset


def build_flower_mesh(
    tileset: Tileset,
    flower_id: str,
    *,
    subdivisions_per_edge: int = 8,
    jitter_amplitude: float = 0.05,
    interior_relief_mm: float = 1.0,
    include_hex_grooves: bool = True,
    plate_depth: float = BASE_PLATE_DEPTH,
    magnet_center_z: float = MAGNET_CENTER_Z,
    magnet_radius: float = MAGNET_RADIUS,
) -> trimesh.Trimesh:
    """Build one flower's complete printable solid: the terrain surface
    (decisions #1-#10) welded directly to a flat base plate carrying the
    magnet bores (decision #13), as one Trimesh - no boolean union (see
    terrain/base_plate.py's own docstring for why, and how the weld
    actually works: shared vertex indices at the flat bottom rim, not two
    independently-capped solids glued together).

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
    bottom_z = -plate_depth

    road_water_side_pairs = flower.roads + flower.water

    vertices, faces, bottom_rim_indices = build_flower_open_solid(
        flower.side_corner_heights,
        layout,
        level_z,
        flower.seed,
        bottom_z=bottom_z,
        subdivisions_per_edge=subdivisions_per_edge,
        jitter_amplitude=jitter_amplitude,
        interior_relief_mm=interior_relief_mm,
        include_hex_grooves=include_hex_grooves and not road_water_side_pairs,
        road_water_side_pairs=road_water_side_pairs,
    )

    extra_vertices, plate_faces = build_base_plate_parts(
        vertices,
        bottom_rim_indices,
        layout,
        subdivisions_per_edge=subdivisions_per_edge,
        bottom_z=bottom_z,
        plate_depth=plate_depth,
        magnet_center_z=magnet_center_z,
        magnet_radius=magnet_radius,
    )

    all_vertices = vertices + extra_vertices
    all_faces = faces + plate_faces

    return trimesh.Trimesh(
        vertices=np.array(all_vertices), faces=np.array(all_faces), process=True
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
