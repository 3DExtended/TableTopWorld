"""Exterior edge magnets and bevels."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from terrain.assembly import AssemblyExporter
from terrain.catalog import EdgeProfileCatalog
from terrain.constants import MAGNET_CENTER_Z
from terrain.layout import FlowerLayout
from terrain.tileset import load_tileset

ROOT = Path(__file__).resolve().parent.parent
DEFAULT = ROOT / "tilesets" / "default.yaml"

EXTERIOR_EDGE_COUNT = len(FlowerLayout(5.1961525).exterior_edge_keys())


def _magnet_hole_count(scad: str) -> int:
    return scad.count("rotate(a = [90, 0, 0])")


def _magnet_z_values(scad: str) -> list[float]:
    """Z translate immediately before each magnet cylinder chain."""
    pattern = (
        r"translate\(v = \[0, 0, ([-\d.]+)\]\)\s*\{"
        r"\s*translate\(v = \[[^\]]+\]\)\s*\{"
        r"\s*translate\(v = \[[^\]]+\]\)\s*\{"
        r"\s*rotate\(a = \[0, 0, [^\]]+\]\)\s*\{"
        r"\s*rotate\(a = \[90, 0, 0\]\)"
    )
    return [float(m) for m in re.findall(pattern, scad)]


def test_exterior_edge_count_is_18() -> None:
    assert EXTERIOR_EDGE_COUNT == 18


@pytest.mark.parametrize("flower_id", ["flat_plains", "hill_north", "river_grove"])
def test_all_mating_exterior_edges_get_magnets(flower_id: str) -> None:
    tileset = load_tileset(DEFAULT)
    catalog = EdgeProfileCatalog()
    flower = tileset.flowers[flower_id]
    mating_count = sum(
        1 for p in flower.edges.values() if catalog.is_mating_profile(p)
    )
    assert mating_count == EXTERIOR_EDGE_COUNT

    exporter = AssemblyExporter(tileset, resolution=32)
    scad = str(exporter.build_flower(flower_id))
    assert _magnet_hole_count(scad) == EXTERIOR_EDGE_COUNT


def test_magnet_z_uniform_at_legacy_height() -> None:
    tileset = load_tileset(DEFAULT)
    exporter = AssemblyExporter(tileset, resolution=32)
    for flower_id in ("flat_plains", "hill_north"):
        scad = str(exporter.build_flower(flower_id))
        z_values = _magnet_z_values(scad)
        assert len(z_values) == EXTERIOR_EDGE_COUNT
        assert all(z == pytest.approx(MAGNET_CENTER_Z) for z in z_values)


def test_flat_plains_magnets_applied_after_road_cuts() -> None:
    tileset = load_tileset(DEFAULT)
    exporter = AssemblyExporter(tileset, resolution=32)
    scad = str(exporter.build_flower("flat_plains"))
    assert _magnet_hole_count(scad) == EXTERIOR_EDGE_COUNT
    last_road = scad.rfind("path_sweep")
    last_magnet = scad.rfind("rotate(a = [90, 0, 0])")
    assert last_road >= 0
    assert last_magnet > last_road
