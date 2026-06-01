"""Terrain height step and mesh extrusion tests."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from terrain.assembly import AssemblyExporter
from terrain.constants import (
    BASE_PLATE_DEPTH,
    FLOWER_BOTTOM_Z,
    TERRAIN_LEVELS,
    TERRAIN_Z,
)
from terrain.mesh import FlowerMeshBuilder
from terrain.tileset import load_tileset

ROOT = Path(__file__).resolve().parent.parent
DEFAULT = ROOT / "tilesets" / "default.yaml"


def test_terrain_z_equal_step_deltas() -> None:
    levels = list(TERRAIN_LEVELS)
    for a, b in zip(levels, levels[1:]):
        assert TERRAIN_Z[b] - TERRAIN_Z[a] == 4.0


def test_terrain_z_matches_meta_model_step() -> None:
    tileset = load_tileset(DEFAULT)
    step = tileset.meta.model_step
    for i, level in enumerate(TERRAIN_LEVELS):
        assert tileset.terrain_z(level) == i * step


def test_terrain_z_absolute_levels() -> None:
    assert TERRAIN_Z["ground"] == 0.0
    assert TERRAIN_Z["middle"] == 4.0
    assert TERRAIN_Z["high"] == 8.0


@pytest.mark.parametrize(
    ("hex_idx", "terrain", "expected_top"),
    [
        (1, "ground", 0.0),
        (0, "middle", 4.0),
        (3, "high", 8.0),
    ],
)
def test_hill_north_hex_top_z(hex_idx: int, terrain: str, expected_top: float) -> None:
    tileset = load_tileset(DEFAULT)
    flower = tileset.flowers["hill_north"]
    builder = FlowerMeshBuilder(tileset)
    assert flower.hexes[str(hex_idx)].terrain == terrain
    assert builder.hex_top_z(flower, hex_idx) == expected_top


def test_ground_standable_uses_base_plate_below_zero() -> None:
    tileset = load_tileset(DEFAULT)
    flower = tileset.flowers["hill_north"]
    builder = FlowerMeshBuilder(tileset)
    scad = str(builder.build_hex_solid(flower, 1))
    assert f"linear_extrude(height = {BASE_PLATE_DEPTH})" in scad
    assert f"translate(v = [0, 0, {FLOWER_BOTTOM_Z}])" in scad
    assert "linear_extrude(height = 8.0)" not in scad


def test_middle_standable_extrudes_from_flower_bottom() -> None:
    tileset = load_tileset(DEFAULT)
    flower = tileset.flowers["hill_north"]
    builder = FlowerMeshBuilder(tileset)
    scad = str(builder.build_hex_solid(flower, 0))
    assert "linear_extrude(height = 6.0)" in scad
    assert f"translate(v = [0, 0, {FLOWER_BOTTOM_Z}])" in scad


def test_high_standable_extrudes_from_flower_bottom() -> None:
    tileset = load_tileset(DEFAULT)
    flower = tileset.flowers["hill_north"]
    builder = FlowerMeshBuilder(tileset)
    scad = str(builder.build_hex_solid(flower, 3))
    assert "linear_extrude(height = 10.0)" in scad
    assert f"translate(v = [0, 0, {FLOWER_BOTTOM_Z}])" in scad


def test_flower_hexes_share_common_bottom_z() -> None:
    tileset = load_tileset(DEFAULT)
    flower = tileset.flowers["hill_north"]
    builder = FlowerMeshBuilder(tileset)
    bottom = f"translate(v = [0, 0, {FLOWER_BOTTOM_Z}])"
    for hex_idx in range(7):
        scad = str(builder.build_hex_solid(flower, hex_idx))
        assert bottom in scad


def test_slope_ramp_span_equals_level_delta() -> None:
    tileset = load_tileset(DEFAULT)
    flower = tileset.flowers["hill_north"]
    builder = FlowerMeshBuilder(tileset)
    # hex 5: ground slope; neighbors top out at middle (z=4), not high (hex 3)
    scad = str(builder.build_hex_solid(flower, 5))
    heights = [float(x) for x in re.findall(r"linear_extrude\(height = ([\d.]+)\)", scad)]
    assert BASE_PLATE_DEPTH in heights
    assert 4.0 in heights


def test_flat_plains_ground_tops_at_zero() -> None:
    tileset = load_tileset(DEFAULT)
    flower = tileset.flowers["flat_plains"]
    builder = FlowerMeshBuilder(tileset)
    for i in range(7):
        assert builder.hex_top_z(flower, i) == 0.0


def test_flat_plains_export_magnet_at_ground_top() -> None:
    tileset = load_tileset(DEFAULT)
    exporter = AssemblyExporter(tileset, resolution=32)
    scad = str(exporter.build_flower("flat_plains"))
    assert "translate(v = [0, 0, 2.0])" not in scad
    assert "translate(v = [0, 0, 0.0])" in scad
