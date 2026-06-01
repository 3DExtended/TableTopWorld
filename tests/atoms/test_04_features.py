"""Layer 4 — feature atoms: road and water channel cutters."""

from __future__ import annotations

from solid2 import cube, union

from terrain.constants import STREET_INDENT_HEIGHT, WATER_INDENT_HEIGHT
from terrain.features import FeatureCutters
from terrain.mesh import FlowerMeshBuilder
from terrain.render.planner import plan_flower
from tests.helpers.atom_builders import atom_road_only, atom_water_only


def _road_subtraction_tools(tileset, flower) -> list:
    layout = tileset.layout()
    cutters = FeatureCutters(layout, resolution=32)
    max_z = max(tileset.terrain_z(h.terrain) for h in flower.hexes.values())
    tools: list = []
    counter = 0
    for entry_j, exit_j in flower.roads:
        counter += 1
        tools.extend(
            cutters._path_pair_tools(
                entry_j,
                exit_j,
                n_gon=4,
                spin=45,
                width_scalar=1.0,
                indent_height=STREET_INDENT_HEIGHT,
                host_z=max_z,
                counter=counter,
            )
        )
    return tools


def _water_subtraction_tools(tileset, flower) -> list:
    layout = tileset.layout()
    cutters = FeatureCutters(layout, resolution=32)
    water_host_z = tileset.water_z("ground") + tileset.meta.model_step
    tools: list = []
    counter = 0
    for entry_j, exit_j in flower.water:
        counter += 1
        tools.extend(
            cutters._path_pair_tools(
                entry_j,
                exit_j,
                n_gon=6,
                spin=60,
                width_scalar=0.75,
                indent_height=WATER_INDENT_HEIGHT,
                host_z=water_host_z,
                counter=counter,
            )
        )
    return tools


def test_17_road_cut_tools_only(visual) -> None:
    """Road path_sweep subtraction tools (junction 0 → 4) shown alone."""
    tileset = atom_road_only()
    flower = tileset.flowers["atom_road_only"]
    tools = _road_subtraction_tools(tileset, flower)
    solid = union()(*tools)
    spec = plan_flower(tileset, flower).path_cuts[0]
    assert spec.kind == "road"
    assert spec.entry_junction == 0
    assert spec.exit_junction == 4
    visual(
        "17_road_cut_tools",
        solid,
        "Road cutter geometry only — curved channel from junction 0 to 4.",
        subdir="04_features",
    )


def test_18_road_cut_into_ground_slab(visual) -> None:
    """Road subtracted from thick ground block — visible groove."""
    tileset = atom_road_only()
    flower = tileset.flowers["atom_road_only"]
    block = cube([30, 30, 4]).translate([0, 0, -2])
    cutters = FeatureCutters(tileset.layout(), resolution=32)
    solid = cutters.apply_features(block, flower, tileset)
    visual(
        "18_road_on_block",
        solid,
        "Ground slab minus road — groove should match 17_road_cut_tools.",
        subdir="04_features",
    )


def test_19_road_on_seven_hex_plateau(visual) -> None:
    """Road cut on full flat 7-hex mesh (features only, no magnets/bevels)."""
    tileset = atom_road_only()
    flower = tileset.flowers["atom_road_only"]
    builder = FlowerMeshBuilder(tileset)
    solid = builder.build_flower(flower)
    solid = FeatureCutters(tileset.layout(), 32).apply_features(solid, flower, tileset)
    visual(
        "19_road_on_plateau",
        solid,
        "Flat flower mesh with road channel — path crosses the plate.",
        subdir="04_features",
    )


def test_20_water_cut_tools_only(visual) -> None:
    """Water channel cutters (hexagon profile) shown alone."""
    tileset = atom_water_only()
    flower = tileset.flowers["atom_water_only"]
    tools = _water_subtraction_tools(tileset, flower)
    solid = union()(*tools)
    spec = plan_flower(tileset, flower).path_cuts[0]
    assert spec.kind == "water"
    visual(
        "20_water_cut_tools",
        solid,
        "Water cutter tools only — wider hex profile than road.",
        subdir="04_features",
    )


def test_21_water_on_seven_hex_plateau(visual) -> None:
    """Water channel on flat 7-hex mesh."""
    tileset = atom_water_only()
    flower = tileset.flowers["atom_water_only"]
    builder = FlowerMeshBuilder(tileset)
    solid = builder.build_flower(flower)
    solid = FeatureCutters(tileset.layout(), 32).apply_features(solid, flower, tileset)
    visual(
        "21_water_on_plateau",
        solid,
        "Flat flower with water groove — deeper/wider than road.",
        subdir="04_features",
    )
