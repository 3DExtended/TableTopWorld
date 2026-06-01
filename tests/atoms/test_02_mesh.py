"""Layer 2 — mesh atoms: single hex prisms, slopes, topping holes."""

from __future__ import annotations

import re

from terrain.constants import BASE_PLATE_DEPTH, FLOWER_BOTTOM_Z, TERRAIN_Z
from terrain.mesh import FlowerMeshBuilder
from terrain.render.planner import plan_flower
from terrain.tileset import load_tileset
from tests.helpers.atom_builders import atom_all_ground, atom_height_step, atom_slope_pad

ROOT = __import__("pathlib").Path(__file__).resolve().parent.parent.parent
DEFAULT = ROOT / "tilesets" / "default.yaml"


def test_05_hex_prism_ground_standable(visual) -> None:
    """Single ground standable hex (cell 1): 2 mm base plate from flower bottom to Z=0."""
    tileset = atom_all_ground()
    flower = tileset.flowers["atom_all_ground"]
    builder = FlowerMeshBuilder(tileset)
    solid = builder.build_hex_solid(flower, 1)
    spec = plan_flower(tileset, flower).hexes[1]
    assert spec.top_z == 0.0
    assert spec.prism_height == BASE_PLATE_DEPTH
    scad = str(solid)
    assert f"linear_extrude(height = {BASE_PLATE_DEPTH})" in scad
    visual(
        "05_hex_prism_ground",
        solid,
        "One ground hex — short slab from Z=-2 to Z=0 only.",
        subdir="02_mesh",
    )


def test_06_hex_prism_middle_standable(visual) -> None:
    """Center hex at middle terrain: extrusion height 6 (from -2 to +4)."""
    tileset = atom_height_step()
    flower = tileset.flowers["atom_height_step"]
    builder = FlowerMeshBuilder(tileset)
    solid = builder.build_hex_solid(flower, 0)
    spec = plan_flower(tileset, flower).hexes[0]
    assert spec.terrain == "middle"
    assert spec.top_z == TERRAIN_Z["middle"]
    assert spec.prism_height == TERRAIN_Z["middle"] - FLOWER_BOTTOM_Z
    visual(
        "06_hex_prism_middle",
        solid,
        "Center middle hex — taller block topping at Z=4.",
        subdir="02_mesh",
    )


def test_07_hex_prism_high_standable(visual) -> None:
    """Hex 3 at high terrain: extrusion to Z=8."""
    tileset = atom_height_step()
    flower = tileset.flowers["atom_height_step"]
    builder = FlowerMeshBuilder(tileset)
    solid = builder.build_hex_solid(flower, 3)
    spec = plan_flower(tileset, flower).hexes[3]
    assert spec.top_z == TERRAIN_Z["high"]
    visual(
        "07_hex_prism_high",
        solid,
        "High hex 3 — tallest single cell in height-step atom.",
        subdir="02_mesh",
    )


def test_08_slope_ramp_ground_to_middle(visual) -> None:
    """Slope hex 5: ground base plus ramp up to middle neighbor height (Z=4)."""
    tileset = atom_slope_pad()
    flower = tileset.flowers["atom_slope_pad"]
    builder = FlowerMeshBuilder(tileset)
    solid = builder.build_hex_solid(flower, 5)
    spec = plan_flower(tileset, flower).hexes[5]
    assert spec.slope_ramp_height == TERRAIN_Z["middle"]
    assert spec.slope_ramp_base_z == 0.0
    heights = [float(x) for x in re.findall(r"linear_extrude\(height = ([\d.]+)\)", str(solid))]
    assert BASE_PLATE_DEPTH in heights
    assert TERRAIN_Z["middle"] in heights
    visual(
        "08_slope_ramp",
        solid,
        "Ground slope hex with wedge ramp — top should meet middle plateau.",
        subdir="02_mesh",
    )


def test_09_topping_hole_tool(visual) -> None:
    """Cylinders at topping positions (hex 1 & 2) — subtraction tools only."""
    from solid2 import cylinder, union

    tileset = atom_all_ground()
    flower = tileset.flowers["atom_all_ground"]
    builder = FlowerMeshBuilder(tileset)
    layout = tileset.layout()
    tools = []
    for hex_idx in flower.topping_hexes:
        z = builder.hex_height_at(flower, hex_idx)
        cx, cy = layout.cell_center(hex_idx)
        tools.append(
            cylinder(h=2.0, r=0.53, center=True).translate([cx, cy, z])
        )
    solid = union()(*tools)
    visual(
        "09_topping_hole_tools",
        solid,
        "Two vertical topping cylinders at hex 1 and 2 centers — holes for mini bases.",
        subdir="02_mesh",
    )


def test_10_topping_subtracted_on_hex(visual) -> None:
    """Hex 1 with topping hole cut — should show cylindrical pocket on top."""
    tileset = atom_all_ground()
    flower = tileset.flowers["atom_all_ground"]
    builder = FlowerMeshBuilder(tileset)
    part = builder.build_hex_solid(flower, 1)
    parts = builder.subtract_topping_holes([part], flower)
    assert len(parts) == 1
    visual(
        "10_topping_subtracted",
        parts[0],
        "Ground hex 1 after topping subtraction — dimple on top face.",
        subdir="02_mesh",
    )


def test_11_seven_hex_union_no_features(visual) -> None:
    """Full 7-hex union (mesh only): height-step atom without edges or cuts."""
    tileset = atom_height_step()
    flower = tileset.flowers["atom_height_step"]
    builder = FlowerMeshBuilder(tileset)
    solid = builder.build_flower(flower)
    spec = plan_flower(tileset, flower)
    assert spec.max_z == TERRAIN_Z["high"]
    assert len(spec.hexes) == 7
    visual(
        "11_seven_hex_union",
        solid,
        "Seven hexes unioned — center middle, hex 3 high, rest ground.",
        subdir="02_mesh",
    )


def test_12_hill_north_single_hex_from_default(visual) -> None:
    """Reference: hill_north hex 5 slope from production tileset."""
    tileset = load_tileset(DEFAULT)
    flower = tileset.flowers["hill_north"]
    builder = FlowerMeshBuilder(tileset)
    solid = builder.build_hex_solid(flower, 5)
    visual(
        "12_hill_north_slope_hex5",
        solid,
        "Production hill_north slope hex — compare with 08_slope_ramp.",
        subdir="02_mesh",
    )
