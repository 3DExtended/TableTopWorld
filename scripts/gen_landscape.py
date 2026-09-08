"""Generate tilesets/landscape.yaml: a 19-flower landscape.

The flowers fill a hexagon of radius 2 on the flower grid (every (q, r)
with max(|q|, |r|, |q+r|) <= 2). Every number in the tileset is derived
from ONE continuous elevation field over world millimetres - hex levels
at the hex centres, side corner heights at the corner positions - so the
two flowers that share a side always declare the same corners and the
landscape reads as one piece instead of 19 tiles. A river meanders from
the west edge to the east edge through a valley cut into that field, and
a road runs south to north and fords the river in the centre column.

Run from the repo root:

    .venv/bin/python scripts/gen_landscape.py

and preview with

    .venv/bin/python -m terrain.cli render preview --tileset tilesets/landscape.yaml --output output/landscape.stl
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from terrain.layout import FlowerLayout  # noqa: E402

HEX_OUTER_WIDTH = 5.1961525
SCALE = 5
LEVEL_COUNT = 4
GRID_DELTAS = ((1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1))  # side k -> neighbour delta
RADIUS = 2

# Paths as sequences of flower-grid cells, in flow direction. Consecutive
# cells must be grid neighbours, and a path may not leave a flower through
# a side adjacent to the one it entered (that is a 120-degree bend, tighter
# than a road or river is wide - the validator rejects it).
RIVER = [(-2, 0), (-1, 0), (-1, 1), (0, 1), (1, 0), (2, 0)]
ROAD = [(0, -2), (0, -1), (0, 0), (0, 1), (0, 2)]

# Elevation field ingredients, in world mm: (flower cell, height in levels, radius)
HILLS = [((1, 1), 2.4, 150.0), ((-1, -1), 1.7, 110.0), ((2, -2), 1.3, 95.0), ((-2, 2), 0.9, 90.0)]
VALLEY_DEPTH = 2.6
VALLEY_HALF_WIDTH = 38.0
BASE_LEVEL = 1.0


def column_name(q: int) -> str:
    return "abcde"[q + RADIUS]


def flower_id(q: int, r: int) -> str:
    return f"land_{column_name(q)}{r + RADIUS + 1}"


def path_segments(cells: list[tuple[int, int]]) -> dict[tuple[int, int], tuple[int, int]]:
    """{cell: (entry_side, exit_side)} for a path through `cells`; the
    first cell enters from outside the map opposite its first step and the
    last cell leaves the map continuing its entry direction."""
    steps = []
    for a, b in zip(cells, cells[1:]):
        delta = (b[0] - a[0], b[1] - a[1])
        if delta not in GRID_DELTAS:
            raise ValueError(f"{a} -> {b} are not grid neighbours")
        steps.append(GRID_DELTAS.index(delta))
    segments = {}
    entry = (steps[0] + 3) % 6
    for i, cell in enumerate(cells):
        exit_side = steps[i] if i < len(steps) else (entry + 3) % 6
        if (exit_side - entry) % 6 not in (2, 3, 4):
            raise ValueError(f"{cell}: exit side {exit_side} is adjacent to entry side {entry}")
        segments[cell] = (entry, exit_side)
        entry = (exit_side + 3) % 6
    return segments


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--output", type=Path, default=ROOT / "tilesets" / "landscape.yaml")
    args = parser.parse_args(argv)
    layout = FlowerLayout(HEX_OUTER_WIDTH * SCALE)
    cells = [
        (q, r)
        for q in range(-RADIUS, RADIUS + 1)
        for r in range(-RADIUS, RADIUS + 1)
        if max(abs(q), abs(r), abs(q + r)) <= RADIUS
    ]
    origin = {cell: layout.flower_grid_to_xy(*cell) for cell in cells}
    river = path_segments(RIVER)
    road = path_segments(ROAD)

    # The river's world polyline: crossing -> ring hex -> centre hex -> ring hex -> crossing per flower.
    river_poly: list[tuple[float, float]] = []
    for cell in RIVER:
        ox, oy = origin[cell]
        entry, exit_side = river[cell]
        # exits are never adjacent to the entry, so the route always runs through the centre hex
        pts = [
            layout.junction_center(entry),
            layout.cell_center(entry + 1),
            layout.cell_center(0),
            layout.cell_center(exit_side + 1),
            layout.junction_center(exit_side),
        ]
        for x, y in pts:
            world = (ox + x, oy + y)
            if not river_poly or math.hypot(world[0] - river_poly[-1][0], world[1] - river_poly[-1][1]) > 1e-6:
                river_poly.append(world)  # consecutive flowers share the crossing point

    def dist_to_river(x: float, y: float) -> float:
        best = float("inf")
        for (ax, ay), (bx, by) in zip(river_poly, river_poly[1:]):
            dx, dy = bx - ax, by - ay
            t = ((x - ax) * dx + (y - ay) * dy) / (dx * dx + dy * dy)
            t = min(max(t, 0.0), 1.0)
            best = min(best, math.hypot(x - (ax + t * dx), y - (ay + t * dy)))
        return best

    def elevation(x: float, y: float) -> float:
        e = BASE_LEVEL
        for cell, height, radius in HILLS:
            hx, hy = origin[cell]
            d2 = (x - hx) ** 2 + (y - hy) ** 2
            e += height * math.exp(-d2 / (radius * radius))
        d = dist_to_river(x, y)
        e -= VALLEY_DEPTH * math.exp(-(d * d) / (VALLEY_HALF_WIDTH * VALLEY_HALF_WIDTH))
        return e

    def level_at(x: float, y: float) -> int:
        return int(min(max(round(elevation(x, y)), 0), LEVEL_COUNT - 1))

    lines = [
        "# GENERATED by scripts/gen_landscape.py - edit that script, not this file.",
        "#",
        "# A 19-flower landscape: a hexagon of flowers (radius 2 on the flower",
        "# grid) whose hex levels and side corner heights all come from one",
        "# continuous elevation field, so every shared side matches by",
        "# construction. Plains at level 1, a valley at level 0 along the river",
        "# (west edge to east edge), hills up to level 3 in the north-east and",
        "# west, and a road from the south edge to the north edge that fords the",
        "# river in flower land_c4.",
        "",
        "meta:",
        "  height_step_mm: 15",
        f"  level_count: {LEVEL_COUNT}",
        f"  scale: {SCALE}",
        f"  hex_outer_width: {HEX_OUTER_WIDTH}",
        "  min_standable_hexes: 1",
        "",
        "preview_map:",
    ]
    for cell in cells:
        lines.append(f"  - {{ id: {flower_id(*cell)}, at: [{cell[0]}, {cell[1]}], rot: 0 }}")
    lines += ["", "flowers:"]
    for index, cell in enumerate(cells):
        ox, oy = origin[cell]
        fid = flower_id(*cell)
        levels = []
        for h in range(FlowerLayout.HEX_CELL_COUNT):
            cx, cy = layout.cell_center(h)
            levels.append(level_at(ox + cx, oy + cy))
        sides = []
        for k in range(FlowerLayout.SIDE_COUNT):
            corners = layout.side_corners(k)
            sides.append([level_at(ox + x, oy + y) for x, y in corners])
        for k in range(FlowerLayout.SIDE_COUNT):
            assert sides[k][3] == sides[(k + 1) % 6][0], (fid, k)
        features = []
        if cell in river:
            features.append("river")
        if cell in road:
            features.append("road")
        lines.append(f"  # {fid}: flower grid {cell}, levels {levels}" + (f", {' + '.join(features)}" if features else ""))
        lines.append(f"  {fid}:")
        lines.append("    hexes:")
        for h, lvl in enumerate(levels):
            lines.append(f'      "{h}": {{ height_level: {lvl} }}')
        lines.append("    side_corner_heights:")
        for k, side in enumerate(sides):
            lines.append(f"      {k}: [{', '.join(str(v) for v in side)}]")
        lines.append(f"    seed: {101 + index}")
        if cell in road:
            entry, exit_side = road[cell]
            lines.append(f"    roads: [[{entry}, {exit_side}]]")
        else:
            lines.append("    roads: []")
        if cell in river:
            entry, exit_side = river[cell]
            lines.append(f"    water: [[{entry}, {exit_side}]]")
        else:
            lines.append("    water: []")
        lines.append("")
    out = args.output
    out.write_text("\n".join(lines))
    counts = {}
    for cell in cells:
        ox, oy = origin[cell]
        for h in range(FlowerLayout.HEX_CELL_COUNT):
            cx, cy = layout.cell_center(h)
            lvl = level_at(ox + cx, oy + cy)
            counts[lvl] = counts.get(lvl, 0) + 1
    print(f"wrote {out} ({len(cells)} flowers); hexes per level: {dict(sorted(counts.items()))}")


if __name__ == "__main__":
    main()
