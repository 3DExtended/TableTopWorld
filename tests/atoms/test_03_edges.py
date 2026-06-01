"""Layer 3 — edge atoms: bevels and magnet holes on simple bases."""

from __future__ import annotations

from terrain.catalog import EdgeProfileCatalog
from terrain.constants import HEXAGON_BEVEL_SIZE, MAGNET_CENTER_Z, TERRAIN_Z
from terrain.edges import EdgeGeometry, add_bevel, add_magnet_hole_on_side
from terrain.mesh import FlowerMeshBuilder
from terrain.render.planner import plan_flower
from tests.helpers.atom_builders import atom_all_ground, atom_height_step


def _ground_slab(tileset, flower, hex_idx: int = 1):
    builder = FlowerMeshBuilder(tileset)
    return builder.build_hex_solid(flower, hex_idx)


def test_13_single_magnet_on_ground_hex(visual) -> None:
    """One magnet bore through ground hex 1 side wall (fixed MAGNET_CENTER_Z)."""
    tileset = atom_all_ground()
    flower = tileset.flowers["atom_all_ground"]
    layout = tileset.layout()
    builder = FlowerMeshBuilder(tileset)
    base = builder.build_hex_solid(flower, 1)
    edge = next(e for e in layout.exterior_edges() if e.key == "1-0")
    solid = add_magnet_hole_on_side(base, edge.line_2d, MAGNET_CENTER_Z)
    spec_magnet = plan_flower(tileset, flower).magnet_holes[0]
    assert spec_magnet.edge_key == "1-0"
    assert spec_magnet.center_z == MAGNET_CENTER_Z
    assert "intersection()" not in str(base)
    assert "difference()" in str(solid)
    visual(
        "13_single_magnet",
        solid,
        "Ground hex slab from z=0; horizontal magnet bore on edge 1-0 at fixed height.",
        subdir="03_edges",
    )


def test_14_single_bevel_on_hex(visual) -> None:
    """One chamfer on a tall ring hex (no magnet brim) at hex top Z."""
    tileset = atom_height_step()
    flower = tileset.flowers["atom_height_step"]
    layout = tileset.layout()
    builder = FlowerMeshBuilder(tileset)
    hex_idx = 3
    base = builder.build_hex_solid(flower, hex_idx)
    ext = next(e for e in layout.exterior_edges() if e.key == f"{hex_idx}-0")
    edge = next(
        e
        for e in layout.hex_bevel_edges()
        if e.hex_idx == hex_idx
        and (
            e.line_2d == ext.line_2d
            or e.line_2d == (ext.line_2d[1], ext.line_2d[0])
        )
    )
    z_top = builder.hex_top_z(flower, hex_idx)
    z_anchor = z_top
    line_3d = (
        (edge.line_2d[0][0], edge.line_2d[0][1], z_anchor),
        (edge.line_2d[1][0], edge.line_2d[1][1], z_anchor),
    )
    solid = add_bevel(
        base,
        line_3d,
        HEXAGON_BEVEL_SIZE,
        z_anchor,
        toward_xy=layout.cell_center(hex_idx),
    )
    assert z_top == TERRAIN_Z["high"]
    assert z_anchor == TERRAIN_Z["high"]
    scad = str(solid)
    assert "difference()" in scad
    assert "polyhedron" in scad
    assert "intersection()" not in scad.split("difference()")[0]
    visual(
        "14_single_bevel",
        solid,
        "High hex 3 — chamfer on one outward side clips the top corner (no lip above the bevel).",
        subdir="03_edges",
    )


def test_15_all_magnets_flat_plateau(visual) -> None:
    """All 18 magnet holes on uniform ground flower (mesh + magnets only)."""
    tileset = atom_all_ground()
    flower = tileset.flowers["atom_all_ground"]
    builder = FlowerMeshBuilder(tileset)
    solid = builder.build_flower(flower)
    edge_geom = EdgeGeometry(tileset.layout())
    catalog = EdgeProfileCatalog()
    solid = edge_geom.apply_magnets(solid, flower, catalog)
    spec = plan_flower(tileset, flower)
    assert len(spec.magnet_holes) == 18
    visual(
        "15_all_magnets_ground_plateau",
        solid,
        "Flat ground 7-hex plate with 18 magnet holes — ring should be perforated.",
        subdir="03_edges",
    )


def test_16_bevels_on_height_step(visual) -> None:
    """Bevels applied to seven-hex height-step mesh (no magnets)."""
    from tests.helpers.atom_builders import atom_height_step

    tileset = atom_height_step()
    flower = tileset.flowers["atom_height_step"]
    builder = FlowerMeshBuilder(tileset)
    solid = builder.build_flower(flower)
    max_z = builder.max_flower_z(flower)
    hex_top_z = lambda idx: builder.hex_mesh_top_z(flower, idx)
    solid = EdgeGeometry(tileset.layout()).apply_bevels(solid, hex_top_z)
    spec = plan_flower(tileset, flower)
    assert spec.bevel.z_anchor == max_z
    visual(
        "16_bevels_on_height_step",
        solid,
        "Height-step union with chamfer on every hex side (all seven cells).",
        subdir="03_edges",
    )
