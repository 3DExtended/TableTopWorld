"""Magnets use fixed Z and ground wall height — no separate magnet collar mesh."""

from __future__ import annotations

from terrain.assembly import AssemblyExporter
from terrain.constants import BASE_PLATE_DEPTH, MAGNET_CENTER_Z
from terrain.mesh import FlowerMeshBuilder
from tests.helpers.atom_builders import atom_all_ground


def test_ground_ring_hex_has_no_magnet_collar() -> None:
    tileset = atom_all_ground()
    flower = tileset.flowers["atom_all_ground"]
    builder = FlowerMeshBuilder(tileset)
    scad = str(builder.build_hex_solid(flower, 1))
    assert "intersection()" not in scad
    assert f"linear_extrude(height = {BASE_PLATE_DEPTH})" in scad
    assert "translate(v = [0, 0, -" not in scad


def test_ground_hex_magnet_bore_subtracts() -> None:
    tileset = atom_all_ground()
    exporter = AssemblyExporter(tileset, resolution=32)
    scad = str(exporter.build_flower("atom_all_ground"))
    assert "difference()" in scad
    assert f"translate(v = [0, 0, {MAGNET_CENTER_Z}])" in scad
