"""Polygon triangulation, isolating the third-party dependency.

`triangle` (true constrained Delaunay triangulation) does not build on
this repo's Python 3.14.6 - its C extension depends on `longintrepr.h`,
removed from CPython's internal headers in 3.13+, and it has no
maintained fix (see project plan). This wraps `scipy.spatial.Delaunay`
(unconstrained) instead: triangulate the convex hull of all supplied
points, then discard any triangle whose centroid falls outside the
source polygon boundary - this correctly handles concave boundaries
(the flower's zigzag silhouette) as long as no triangle's centroid
happens to land inside a boundary concavity.

Callers needing a triangle edge to fall exactly on a constraint segment
(grooves, road rails) must supply points densely enough along that
segment that unconstrained Delaunay cannot skip it with a longer chord -
scoped, testable per-shape rather than a general constrained-Delaunay
guarantee.
"""

from __future__ import annotations

from typing import Sequence

import mapbox_earcut as earcut
import numpy as np
from scipy.spatial import Delaunay

Point2D = tuple[float, float]
Triangle = tuple[int, int, int]


def point_in_polygon(point: Point2D, polygon: Sequence[Point2D]) -> bool:
    """Standard ray-casting point-in-polygon test."""
    x, y = point
    inside = False
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        if (y1 > y) != (y2 > y):
            x_intersect = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
            if x < x_intersect:
                inside = not inside
    return inside


def _signed_area(p1: Point2D, p2: Point2D, p3: Point2D) -> float:
    return (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1])


def triangulate_polygon(
    boundary: Sequence[Point2D], interior_points: Sequence[Point2D] = ()
) -> tuple[list[Point2D], list[Triangle]]:
    """Triangulate a simple polygon (boundary, CCW, no self-intersection),
    optionally with extra interior Steiner points.

    Returns (all_points, triangles): all_points is boundary points in
    order followed by interior_points in order, and each triangle is a
    tuple of indices into all_points, wound CCW (matching a +Z-up
    "looking down" viewpoint) regardless of scipy's own simplex order.
    """
    all_points = list(boundary) + list(interior_points)
    if len(all_points) < 3:
        raise ValueError("need at least 3 points to triangulate")

    arr = np.array(all_points)
    delaunay = Delaunay(arr)

    triangles: list[Triangle] = []
    for simplex in delaunay.simplices:
        i, j, k = (int(x) for x in simplex)
        p1, p2, p3 = all_points[i], all_points[j], all_points[k]
        centroid = ((p1[0] + p2[0] + p3[0]) / 3, (p1[1] + p2[1] + p3[1]) / 3)
        if not point_in_polygon(centroid, boundary):
            continue
        if _signed_area(p1, p2, p3) < 0:
            i, j = j, i
        triangles.append((i, j, k))
    return all_points, triangles


def earcut_triangulate_with_holes(
    outer: Sequence[Point2D], holes: Sequence[Sequence[Point2D]] = ()
) -> tuple[list[Point2D], list[Triangle]]:
    """Ear-clipping triangulation of a polygon with holes (e.g. a magnet
    bore cut through a flat base-plate wall panel). Returns (all_points,
    triangles): all_points is `outer` followed by each hole's points in
    order, and each triangle is wound CCW (matching `outer`'s own
    orientation) regardless of earcut's internal simplex order.

    mapbox_earcut's `rings` argument is CUMULATIVE END indices into the
    flattened point array (not lengths, and not start indices) - its last
    entry must equal the total point count, verified empirically since
    the library raises otherwise.
    """
    all_points = list(outer)
    ring_ends = [len(all_points)]
    for hole in holes:
        all_points.extend(hole)
        ring_ends.append(len(all_points))

    arr = np.array(all_points, dtype=np.float64)
    rings = np.array(ring_ends, dtype=np.uint32)
    flat = earcut.triangulate_float64(arr, rings)
    triangles: list[Triangle] = []
    for t in range(0, len(flat), 3):
        i, j, k = (int(x) for x in flat[t : t + 3])
        p1, p2, p3 = all_points[i], all_points[j], all_points[k]
        if _signed_area(p1, p2, p3) < 0:
            i, j = j, i
        triangles.append((i, j, k))
    return all_points, triangles


def earcut_triangulate_polygon(polygon_2d: Sequence[Point2D]) -> list[Triangle]:
    """Ear-clipping triangulation of one simple polygon (no holes),
    wound CCW - always preserves every edge of `polygon_2d` as a real
    triangle edge, unlike scipy Delaunay above.

    Needed whenever a polygon has extra points that lie exactly on (not
    just near) its own boundary - collinear points strictly between two
    real corners, e.g. a hex cell's fine-subdivided exterior edges, or a
    road/river polyline splitting a region. Delaunay triangulation of a
    point set with collinear runs is a degenerate case with no unique
    answer, and scipy/Qhull resolves it by producing triangles whose
    centroid can fall fractionally outside the source polygon - a real,
    observed failure (see git history for build_flower_cells and
    terrain/roads.py), not a hypothetical one. Ear-clipping has no such
    ambiguity: it operates on the polygon's given vertex order, not on
    point geometry, so collinear runs triangulate cleanly.
    """
    arr = np.array(polygon_2d, dtype=np.float64)
    rings = np.array([len(polygon_2d)], dtype=np.uint32)
    flat = earcut.triangulate_float64(arr, rings)
    triangles: list[Triangle] = []
    for t in range(0, len(flat), 3):
        i, j, k = (int(x) for x in flat[t : t + 3])
        p1, p2, p3 = polygon_2d[i], polygon_2d[j], polygon_2d[k]
        if _signed_area(p1, p2, p3) < 0:
            i, j = j, i
        triangles.append((i, j, k))
    return triangles
