#!/usr/bin/env python3
"""
Legacy flat 7-hex flower export.

Delegates to the terrain package for flat_plains parity (road, magnets, bevels).
Run from repo root or printableFiles/; writes hexagon.scad beside this file.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from solid2 import set_global_fa, set_global_fn, set_global_fs

from terrain.assembly import AssemblyExporter
from terrain.tileset import load_tileset

RESOLUTION = 100
SCALE = 5
TILESET = ROOT / "tilesets" / "default.yaml"
OUTPUT = HERE / "hexagon.scad"


def main() -> None:
    set_global_fn(RESOLUTION)
    set_global_fa(RESOLUTION)
    set_global_fs(RESOLUTION)
    tileset = load_tileset(TILESET)
    exporter = AssemblyExporter(tileset, resolution=RESOLUTION)
    solid = exporter.build_flower("flat_plains")
    solid.scale(SCALE).save_as_scad(str(OUTPUT))
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
