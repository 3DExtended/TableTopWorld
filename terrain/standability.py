"""Post-generation standability check - design decision #11.

"'>=1 standable hex' is checked after generation, not guaranteed by
construction, and is itself configurable (a flower may validly have 0)."
This deliberately checks the REAL generated mesh (sampling actual Z at
each hex cell's own corners/center), not the authored tileset data - a
hex could be nominally "flat" in intent but still end up sloped in the
real output if its corners sit at different boundary height levels, and
a flower with mismatched corner heights is exactly the interesting case
this check needs to catch.
"""

from __future__ import annotations

import trimesh

from terrain.layout import FlowerLayout


def hex_cell_z_range(mesh: trimesh.Trimesh, layout: FlowerLayout, hex_idx: int) -> float:
    """Max - min Z among the mesh's own top-surface vertices that fall
    within hex cell `hex_idx`'s 2D footprint - 0.0 for a perfectly flat
    cell, larger for a sloped or cliffed one."""
    polygon = layout.ring_vertices(hex_idx)
    cx, cy = layout.cell_center(hex_idx)
    # A hex cell's own corners are themselves the most informative points
    # (shared boundary between this cell and its neighbors) but sit
    # exactly ON the polygon boundary, where ray-casting point_in_polygon
    # is not reliable (an exact vertex/edge hit can go either way).
    # Inflate slightly around the cell center before testing containment
    # so the cell's own corners always land unambiguously inside.
    inflated = [(cx + (x - cx) * 1.001, cy + (y - cy) * 1.001) for x, y in polygon]
    xs = [p[0] for p in inflated]
    ys = [p[1] for p in inflated]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    from terrain.triangulate import point_in_polygon

    # The base plate's floor/walls share XY footprint with the terrain's
    # own boundary rim (by construction - see terrain/base_plate.py) but
    # sit at lower Z, extruded straight down to bottom_z. A plain scan of
    # mesh.vertices would pick those up too, at exactly the rim's XY, and
    # report the plate's own depth as if it were terrain slope. Grouping
    # by (rounded) XY and keeping only the highest Z per position isolates
    # the actual top surface, regardless of plate_depth (which this module
    # has no other way to know from the mesh alone).
    top_z_by_xy: dict[tuple[float, float], float] = {}
    for x, y, z in mesh.vertices:
        if (
            min_x - 1e-6 <= x <= max_x + 1e-6
            and min_y - 1e-6 <= y <= max_y + 1e-6
            and point_in_polygon((float(x), float(y)), inflated)
        ):
            key = (round(float(x), 3), round(float(y), 3))
            z = float(z)
            if z > top_z_by_xy.get(key, float("-inf")):
                top_z_by_xy[key] = z
    if not top_z_by_xy:
        return float("inf")
    zs = list(top_z_by_xy.values())
    return max(zs) - min(zs)


def standable_hex_count(
    mesh: trimesh.Trimesh, layout: FlowerLayout, *, flatness_tolerance_mm: float = 1.0
) -> int:
    """How many of the flower's 7 hex cells are flat within tolerance in
    the REAL generated mesh."""
    return sum(
        1
        for hex_idx in range(FlowerLayout.HEX_CELL_COUNT)
        if hex_cell_z_range(mesh, layout, hex_idx) <= flatness_tolerance_mm
    )


def check_min_standable_hexes(
    mesh: trimesh.Trimesh,
    layout: FlowerLayout,
    min_standable_hexes: int,
    *,
    flatness_tolerance_mm: float = 1.0,
) -> tuple[bool, int]:
    """Returns (meets_minimum, actual_count). Does not raise - decision
    #11 explicitly permits a flower with 0 standable hexes; callers
    decide what to do with a failing result (e.g. cli.py warns)."""
    count = standable_hex_count(mesh, layout, flatness_tolerance_mm=flatness_tolerance_mm)
    return count >= min_standable_hexes, count
