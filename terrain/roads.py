"""Roads and rivers as embedded surface constraints - decisions #8-#10.

A road or river is declared as (entry_side, exit_side): two of a flower's
6 side indices. Its centerline runs from the junction corner where
entry_side begins, straight through the interior, to the junction corner
where exit_side begins, at a height linearly interpolated between those
two corners' declared height levels (reusing decision #9's "same corner
positions", now for height too - a plain, level, "engineered path" look
rather than adding independent relief). A junction corner is exactly the
point two neighbor flowers already declare must match at a shared
boundary (decisions #2-#4), so a road/river threaded through that same
point continues seamlessly into whichever flower sits next door
(decision #8).

The road is embedded as a REAL mesh edge chain via mapbox_earcut, not
scipy's Delaunay: the flower's boundary loop is split into two arcs at
the road's entry/exit points, each arc is closed off with the road
polyline to form a simple (possibly non-convex) polygon, and each is
triangulated independently. Ear-clipping triangulation (unlike Delaunay)
always preserves every edge of its input polygon - a mathematical
guarantee, not a probabilistic one. This matters here because this
flower's hex-grid geometry is highly symmetric: scipy Delaunay plus
iterative segment subdivision (tried first) left a stubborn, non-
vanishing fraction of constraint segments unresolved no matter how much
the polyline was subdivided or perturbed - see git history for that
attempt. Adjacent entry/exit sides are a genuine geometric degenerate
case (the short arc's enclosed area is ~0, since a road between
immediately-adjacent junctions nearly retraces the boundary itself) and
are rejected rather than silently producing a sliver.

Every hex cell also gets a groove tracing its own 6-edge outline, so hex
boundaries stay visible for distance counting (decision #10) - but that
uses a different, simpler guarantee (terrain.heightfield.build_flower_cells:
each hex is triangulated independently as its own convex polygon, whose
hull edges are always part of any Delaunay triangulation of its own
points) and doesn't need earcut.
"""

from __future__ import annotations

from typing import Callable, Sequence

from terrain.layout import FlowerLayout, Point2D
from terrain.triangulate import earcut_triangulate_polygon

Vertex3D = tuple[float, float, float]
Triangle = tuple[int, int, int]


def side_entry_point(layout: FlowerLayout, side_idx: int) -> Point2D:
    """The junction corner where side `side_idx` begins."""
    return layout.side_corners(side_idx % FlowerLayout.SIDE_COUNT)[0]


def side_entry_boundary_loop_index(side_idx: int, subdivisions_per_edge: int) -> int:
    """Index into build_flower_boundary_loop()'s output of that same
    point: side k's fine contour occupies loop indices
    [k*3*subdivisions_per_edge : (k+1)*3*subdivisions_per_edge), and its
    first point is exactly side_entry_point(layout, k) (verified against
    the real boundary loop output)."""
    return side_idx * FlowerLayout.EDGES_PER_SIDE * subdivisions_per_edge


def subdivide_line(p1: Point2D, p2: Point2D, steps: int) -> list[Point2D]:
    return [
        (p1[0] + (p2[0] - p1[0]) * i / steps, p1[1] + (p2[1] - p1[1]) * i / steps)
        for i in range(steps + 1)
    ]


def _points_close(a: Point2D, b: Point2D, tol: float) -> bool:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 < tol * tol


def internal_hex_edges(
    layout: FlowerLayout, *, tol: float = 1e-4
) -> list[tuple[Point2D, Point2D]]:
    """The hex-to-hex boundary edges strictly inside the flower (not part
    of the flower's own 18-edge exterior boundary) - the ones that need an
    engraved groove, since only they aren't already visually marked by the
    flower's own true boundary contour.

    Deliberately does NOT filter by comparing hex_bevel_edges()'s keys
    against exterior_edge_keys(): the two use the SAME "{hex}-{n}" string
    format for two DIFFERENT numbering schemes - exterior_edge_keys()'s
    `n` is a side-position 0..2 (sorted by angle), while hex_bevel_edges()'s
    `n` is the real local edge index 0..5 - and they only coincidentally
    agree for a ring hex whose exterior_vertex_indices() happens to equal
    exactly (0, 1, 2). For any other rotation (which most ring hexes have),
    that comparison silently keeps the wrong 3 of 6 edges as "exterior".
    Comparing against exterior_vertex_indices() directly (the same local
    edge-index space hex_edge_line() uses) avoids the collision entirely.

    Each internal edge is also computed twice, once from each of the two
    hexes it borders, via independent coordinate paths - so the two
    computed copies can differ by float noise well past 1e-9 and must be
    matched by distance tolerance, not by rounding-and-hashing (tried
    first: out of the 12 true shared edges, rounding to 6 decimals merged
    only 1).
    """
    exterior_local_indices = {
        h: set(layout.exterior_vertex_indices(h)) for h in layout.RING_HEX_INDICES
    }
    kept: list[tuple[Point2D, Point2D]] = []
    for hex_idx in range(layout.HEX_CELL_COUNT):
        skip = exterior_local_indices.get(hex_idx, set())
        for edge_idx in range(layout.EDGES_PER_HEX):
            if edge_idx in skip:
                continue
            p1, p2 = layout.hex_edge_line(hex_idx, edge_idx)
            if any(
                (_points_close(p1, q1, tol) and _points_close(p2, q2, tol))
                or (_points_close(p1, q2, tol) and _points_close(p2, q1, tol))
                for q1, q2 in kept
            ):
                continue
            kept.append((p1, p2))
    return kept


class DegenerateRoadError(ValueError):
    """Raised when entry_side/exit_side are adjacent - see module docstring."""


def build_flower_mesh_with_road(
    boundary_loop_3d: list[Vertex3D],
    entry_side: int,
    exit_side: int,
    entry_height_level: int,
    exit_height_level: int,
    level_z: Callable[[int], float],
    *,
    subdivisions_per_edge: int = 8,
    road_subdivisions: int = 20,
    min_side_gap: int = 2,
) -> tuple[list[Vertex3D], list[Triangle]]:
    """Split the flower into two regions at the road/river's entry and
    exit points and triangulate each via earcut, guaranteeing every road
    segment is a real mesh edge. Returns (vertices_3d, triangles) for the
    ENTIRE flower's top surface (both regions combined) - this replaces
    the whole-flower interior triangulation for a flower with a road
    (no freeform interior noise in this region; the surface is a plain
    interpolation between the boundary contour and the road).
    """
    side_gap = (exit_side - entry_side) % FlowerLayout.SIDE_COUNT
    if side_gap < min_side_gap or side_gap > FlowerLayout.SIDE_COUNT - min_side_gap:
        raise DegenerateRoadError(
            f"entry_side {entry_side} and exit_side {exit_side} are too close "
            f"({side_gap} apart) - the short arc between them is nearly "
            "zero-area and cannot be triangulated meaningfully"
        )

    boundary_2d = [(v[0], v[1]) for v in boundary_loop_3d]
    n = len(boundary_2d)
    i_entry = side_entry_boundary_loop_index(entry_side, subdivisions_per_edge)
    i_exit = side_entry_boundary_loop_index(exit_side, subdivisions_per_edge)
    entry_pt, exit_pt = boundary_2d[i_entry], boundary_2d[i_exit]

    road_2d = subdivide_line(entry_pt, exit_pt, road_subdivisions)
    z0, z1 = level_z(entry_height_level), level_z(exit_height_level)
    road_3d: list[Vertex3D] = [
        (x, y, z0 + (z1 - z0) * i / road_subdivisions)
        for i, (x, y) in enumerate(road_2d)
    ]

    def arc(start: int, end: int) -> list[int]:
        return list(range(start, end + 1)) if start <= end else (
            list(range(start, n)) + list(range(0, end + 1))
        )

    arc_forward = arc(i_entry, i_exit)  # entry -> exit
    arc_backward = arc(i_exit, i_entry)  # exit -> entry

    vertices: list[Vertex3D] = list(boundary_loop_3d)
    road_start_idx = len(vertices)
    vertices.extend(road_3d[1:-1])  # interior road points only; endpoints reuse boundary indices
    road_indices = [i_entry] + list(range(road_start_idx, len(vertices))) + [i_exit]

    def region_polygon(arc_indices: list[int], road_indices_order: list[int]) -> tuple[list[int], list[Point2D]]:
        idx_sequence = arc_indices + road_indices_order[1:-1]
        points = [
            (vertices[i][0], vertices[i][1]) if i < len(vertices) else None
            for i in idx_sequence
        ]
        return idx_sequence, [(vertices[i][0], vertices[i][1]) for i in idx_sequence]

    idx_a, points_a = region_polygon(arc_forward, list(reversed(road_indices)))
    idx_b, points_b = region_polygon(arc_backward, road_indices)

    triangles: list[Triangle] = []
    for idx_sequence, points in ((idx_a, points_a), (idx_b, points_b)):
        local_triangles = earcut_triangulate_polygon(points)
        for i, j, k in local_triangles:
            triangles.append((idx_sequence[i], idx_sequence[j], idx_sequence[k]))

    return vertices, triangles
