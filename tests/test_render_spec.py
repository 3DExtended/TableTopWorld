"""Unit tests against text render specs (no SCAD string parsing)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from terrain.assembly import AssemblyExporter
from terrain.catalog import EdgeProfileCatalog
from terrain.constants import (
    BASE_PLATE_DEPTH,
    FLOWER_BOTTOM_Z,
    MAGNET_CENTER_Z,
    MAGNET_DEPTH,
    MAGNET_RADIUS,
    TERRAIN_Z,
)
from terrain.layout import FlowerLayout
from terrain.render.format import format_render_spec
from terrain.render.planner import plan_flower
from terrain.tileset import load_tileset

ROOT = Path(__file__).resolve().parent.parent
DEFAULT = ROOT / "tilesets" / "default.yaml"
EXTERIOR_EDGE_COUNT = len(FlowerLayout(5.1961525).exterior_edge_keys())


@pytest.fixture
def tileset():
    return load_tileset(DEFAULT)


def test_flat_plains_all_ground_hexes(tileset) -> None:
    spec = plan_flower(tileset, tileset.flowers["flat_plains"])
    assert spec.flower_bottom_z == FLOWER_BOTTOM_Z
    assert spec.max_z == BASE_PLATE_DEPTH
    for cell in spec.hexes:
        assert cell.terrain == "ground"
        assert cell.top_z == BASE_PLATE_DEPTH
        assert cell.prism_height == BASE_PLATE_DEPTH
        assert cell.slope_ramp_height is None


def test_hill_north_mixed_heights_and_slope_ramp(tileset) -> None:
    spec = plan_flower(tileset, tileset.flowers["hill_north"])
    by_idx = {h.hex_idx: h for h in spec.hexes}
    assert by_idx[1].terrain == "ground" and by_idx[1].top_z == BASE_PLATE_DEPTH
    assert by_idx[0].terrain == "middle" and by_idx[0].top_z == TERRAIN_Z["middle"]
    assert by_idx[3].terrain == "high" and by_idx[3].top_z == TERRAIN_Z["high"]
    slope = by_idx[5]
    assert slope.role == "slope"
    assert slope.slope_ramp_height == TERRAIN_Z["middle"] - BASE_PLATE_DEPTH
    assert slope.slope_ramp_base_z == BASE_PLATE_DEPTH


def test_all_hexes_share_flower_bottom_z(tileset) -> None:
    spec = plan_flower(tileset, tileset.flowers["hill_north"])
    for cell in spec.hexes:
        assert cell.bottom_z == FLOWER_BOTTOM_Z


@pytest.mark.parametrize("flower_id", ["flat_plains", "hill_north", "river_grove"])
def test_mating_edges_all_get_magnet_holes(tileset, flower_id: str) -> None:
    catalog = EdgeProfileCatalog()
    flower = tileset.flowers[flower_id]
    spec = plan_flower(tileset, flower)
    mating = sum(1 for e in spec.exterior_edges if catalog.is_mating_profile(e.profile))
    assert mating == EXTERIOR_EDGE_COUNT
    assert len(spec.magnet_holes) == EXTERIOR_EDGE_COUNT
    assert {m.edge_key for m in spec.magnet_holes} == {e.edge_key for e in spec.exterior_edges}


def test_magnet_z_uniform(tileset) -> None:
    for flower_id in ("flat_plains", "hill_north"):
        spec = plan_flower(tileset, tileset.flowers[flower_id])
        assert all(m.center_z == MAGNET_CENTER_Z for m in spec.magnet_holes)


def test_flat_plains_has_road_cut(tileset) -> None:
    spec = plan_flower(tileset, tileset.flowers["flat_plains"])
    roads = [p for p in spec.path_cuts if p.kind == "road"]
    assert len(roads) == 1
    assert roads[0].entry_junction == 0
    assert roads[0].exit_junction == 4


def test_river_grove_has_water_cut(tileset) -> None:
    spec = plan_flower(tileset, tileset.flowers["river_grove"])
    water = [p for p in spec.path_cuts if p.kind == "water"]
    assert len(water) >= 1


def test_format_render_spec_is_stable_json(tileset) -> None:
    spec = plan_flower(tileset, tileset.flowers["flat_plains"])
    text = format_render_spec(spec)
    parsed = json.loads(text)
    assert parsed["flower_id"] == "flat_plains"
    assert len(parsed["hexes"]) == 7


def test_assembly_exporter_describe_flower(tileset) -> None:
    exporter = AssemblyExporter(tileset, resolution=32)
    spec = exporter.describe_flower("hill_north")
    assert spec.flower_id == "hill_north"
    assert spec.resolution == 32


def test_topping_holes_on_default_hexes(tileset) -> None:
    spec = plan_flower(tileset, tileset.flowers["flat_plains"])
    assert len(spec.topping_holes) == 2
    indices = {h.hex_idx for h in spec.topping_holes}
    assert indices == {1, 2}
    for hole in spec.topping_holes:
        assert hole.radius == MAGNET_RADIUS
        assert hole.depth == MAGNET_DEPTH
