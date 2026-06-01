"""Magnet collar is part of ring-hex mesh so bores intersect without separate fins."""

from __future__ import annotations

from terrain.assembly import AssemblyExporter
from terrain.constants import FLOWER_BOTTOM_Z, MAGNET_CENTER_Z, MAGNET_WALL_TOP_Z
from terrain.mesh import FlowerMeshBuilder
from terrain.tileset import load_tileset
from tests.helpers.atom_builders import atom_all_ground
from tests.helpers.visual_scad import ROOT

DEFAULT = ROOT / "tilesets" / "default.yaml"
COLLAR_HEIGHT = MAGNET_WALL_TOP_Z - FLOWER_BOTTOM_Z


def test_ground_ring_hex_includes_magnet_collar() -> None:
    tileset = atom_all_ground()
    flower = tileset.flowers["atom_all_ground"]
    builder = FlowerMeshBuilder(tileset)
    scad = str(builder.build_hex_solid(flower, 1))
    assert "intersection()" in scad
    assert f"linear_extrude(height = {COLLAR_HEIGHT})" in scad


def test_ground_hex_magnet_bore_subtracts() -> None:
    tileset = atom_all_ground()
    exporter = AssemblyExporter(tileset, resolution=32)
    scad = str(exporter.build_flower("atom_all_ground"))
    assert "difference()" in scad
    assert f"translate(v = [0, 0, {MAGNET_CENTER_Z}])" in scad


def test_flat_plains_ring_hexes_have_collars() -> None:
    tileset = load_tileset(DEFAULT)
    builder = FlowerMeshBuilder(tileset)
    flower = tileset.flowers["flat_plains"]
    scad = str(builder.build_hex_solid(flower, 1))
    assert "intersection()" in scad
