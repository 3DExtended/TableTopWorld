"""Bevel cutter must overlap the hex solid in Z (no double z_anchor offset)."""

from __future__ import annotations

import pytest

from tests.helpers.atom_builders import atom_all_ground
from terrain.constants import (
    BASE_PLATE_DEPTH,
    BEVEL_TOP_SLAB_DEPTH,
    BEVEL_Z_INSET,
    HEXAGON_BEVEL_SIZE,
)
from terrain.edges import add_bevel
from terrain.layout import FlowerLayout
from terrain.mesh import FlowerMeshBuilder


def test_bevel_cutter_side_below_top_at_z_anchor() -> None:
    """Side wedge below z_anchor; top tools centered on z_anchor."""
    tileset = atom_all_ground()
    flower = tileset.flowers["atom_all_ground"]
    layout = tileset.layout()
    builder = FlowerMeshBuilder(tileset)
    builder.build_hex_solid(flower, 1)
    edge = next(e for e in layout.exterior_edges() if e.key == "1-0")
    z_anchor = BASE_PLATE_DEPTH
    line_3d = (
        (edge.line_2d[0][0], edge.line_2d[0][1], z_anchor),
        (edge.line_2d[1][0], edge.line_2d[1][1], z_anchor),
    )
    tool_scad = str(add_bevel(None, line_3d, HEXAGON_BEVEL_SIZE, z_anchor))
    z_side = z_anchor - BEVEL_Z_INSET
    assert tool_scad.count(f"translate(v = [0, 0, {z_side}])") == 1
    assert tool_scad.count(f"translate(v = [0, 0, {z_anchor}])") == 1
    depth = HEXAGON_BEVEL_SIZE
    assert f", {-depth}]" in tool_scad or f",-{depth}]" in tool_scad.replace(" ", "")
    assert f", {depth}]" in tool_scad or f",{depth}]" in tool_scad.replace(" ", "")
    half = BEVEL_TOP_SLAB_DEPTH / 2
    assert f", {-half}]" in tool_scad or f",{-half}]" in tool_scad.replace(" ", "")
    assert f", {half}]" in tool_scad or f",{half}]" in tool_scad.replace(" ", "")


def test_bevel_top_wedge_inward_offset() -> None:
    """Top shelf inward line is offset 0.75× bevel depth from the exterior edge."""
    import math

    from terrain.edges import _bevel_line_xy, _bevel_top_shelf_wedge

    layout = FlowerLayout()
    edge = next(e for e in layout.exterior_edges() if e.key == "1-0")
    line_xy = _bevel_line_xy(
        (
            (edge.line_2d[0][0], edge.line_2d[0][1], BASE_PLATE_DEPTH),
            (edge.line_2d[1][0], edge.line_2d[1][1], BASE_PLATE_DEPTH),
        )
    )
    scad = str(_bevel_top_shelf_wedge(line_xy, HEXAGON_BEVEL_SIZE))
    import re

    nums = [float(x) for x in re.findall(r"[-+]?\d*\.\d+|\d+", scad.split("points = ")[1].split("]")[0])]
    in_x, in_y = nums[0], nums[1]
    offset = math.hypot(in_x - line_xy[0][0], in_y - line_xy[0][1])
    assert offset == pytest.approx(HEXAGON_BEVEL_SIZE * 0.75, rel=1e-4)


def test_path_feature_tools_exclude_bevel_cutters() -> None:
    """Road/water sweeps must not embed chamfer polyhedrons (bevels are edge-only)."""
    from solid2 import union

    from tests.atoms.test_04_features import _road_subtraction_tools, _water_subtraction_tools
    from tests.helpers.atom_builders import atom_road_only, atom_water_only

    for tileset, tools_fn in (
        (atom_road_only(), _road_subtraction_tools),
        (atom_water_only(), _water_subtraction_tools),
    ):
        flower = next(iter(tileset.flowers.values()))
        scad = str(union()(*tools_fn(tileset, flower)))
        assert "polyhedron" not in scad


def test_bevel_subtracts_from_ground_hex() -> None:
    tileset = atom_all_ground()
    flower = tileset.flowers["atom_all_ground"]
    layout = tileset.layout()
    builder = FlowerMeshBuilder(tileset)
    base = builder.build_hex_solid(flower, 1)
    edge = next(e for e in layout.exterior_edges() if e.key == "1-0")
    z_anchor = builder.hex_top_z(flower, 1)
    line_3d = (
        (edge.line_2d[0][0], edge.line_2d[0][1], z_anchor),
        (edge.line_2d[1][0], edge.line_2d[1][1], z_anchor),
    )
    solid = add_bevel(base, line_3d, HEXAGON_BEVEL_SIZE, z_anchor)
    assert z_anchor == BASE_PLATE_DEPTH
    scad = str(solid)
    assert "difference()" in scad
    assert "polyhedron" in scad
