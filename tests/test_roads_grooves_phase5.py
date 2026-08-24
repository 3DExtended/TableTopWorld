"""Phase 5: grooves and roads/rivers as embedded surface constraints.

Verifies real embedding (edges actually present in the exported mesh's
edge list), not mere proximity - the gap the plan specifically calls out
("groove/channel edges queried directly from the exported mesh's edge
list, proving real embedding, not proximity").
"""

from __future__ import annotations

import math

import pytest

from terrain.heights import HeightLevels
from terrain.layout import FlowerLayout
from terrain.roads import internal_hex_edges, side_entry_point, subdivide_line
from terrain.surface_mesh import build_flower_surface_mesh

LEVELS = HeightLevels()
FLAT_SIDES = {k: (0, 0, 0, 0) for k in range(FlowerLayout.SIDE_COUNT)}


@pytest.fixture
def layout() -> FlowerLayout:
    return FlowerLayout()


def _find_vertex(mesh, point3d, tol: float = 1e-6) -> int:
    d = (
        (mesh.vertices[:, 0] - point3d[0]) ** 2
        + (mesh.vertices[:, 1] - point3d[1]) ** 2
        + (mesh.vertices[:, 2] - point3d[2]) ** 2
    ) ** 0.5
    idx = int(d.argmin())
    assert d[idx] < tol, f"no vertex within {tol} of {point3d} (closest dist {d[idx]})"
    return idx


def _has_edge(mesh, idx_a: int, idx_b: int) -> bool:
    edges = {tuple(sorted(e)) for e in mesh.edges_unique.tolist()}
    return tuple(sorted((idx_a, idx_b))) in edges


def test_flower_with_grooves_is_still_watertight(layout: FlowerLayout) -> None:
    mesh = build_flower_surface_mesh(
        FLAT_SIDES, layout, LEVELS.z, seed=3, include_hex_grooves=True
    )
    assert mesh.is_watertight
    assert mesh.is_winding_consistent
    assert mesh.volume > 0


def test_hex_grooves_are_embedded_as_direct_mesh_edges(layout: FlowerLayout) -> None:
    """Every internal hex-to-hex edge must be embedded as a chain of REAL
    direct mesh edges corner-to-corner, not merely a run of vertices that
    happen to be close together.

    Subdivided at the same subdivisions_per_edge granularity as the
    flower's own exterior boundary (unlike the exterior boundary, an
    internal groove has no cross-flower contract to satisfy - decision
    #5's interior noise is a pure function of (x, y), so the neighboring
    cell across the same edge queries the exact same points and gets a
    bit-identical Z back with no shared-vertex-index bookkeeping needed).
    Subdividing internal edges too (not just exterior ones) is what fixed
    a real "low poly" look: before this, every internal edge - including
    all 6 of hex 0's own edges - was a single raw corner-to-corner
    segment, so a uniform-height flower's center cell rendered as a
    literal 6-triangle fan.
    """
    subdivisions_per_edge = 8
    mesh = build_flower_surface_mesh(
        FLAT_SIDES,
        layout,
        LEVELS.z,
        seed=3,
        include_hex_grooves=True,
        subdivisions_per_edge=subdivisions_per_edge,
    )
    internal_edges = internal_hex_edges(layout)
    assert len(internal_edges) == 12  # 7 hexes, 42 edges total, 18 exterior

    for p1, p2 in internal_edges:
        chain = subdivide_line(p1, p2, subdivisions_per_edge)
        vertex_indices = []
        for x, y in chain:
            d_xy = (
                (mesh.vertices[:, 0] - x) ** 2 + (mesh.vertices[:, 1] - y) ** 2
            ) ** 0.5
            idx = int(d_xy.argmin())
            assert d_xy[idx] < 1e-6, f"groove point ({x}, {y}) missing from mesh"
            vertex_indices.append(idx)
        for a, b in zip(vertex_indices, vertex_indices[1:]):
            assert _has_edge(mesh, a, b), (
                f"groove segment {p1}->{p2}: no direct mesh edge between "
                f"consecutive subdivided points at indices {a}, {b}"
            )


def test_road_endpoints_equal_side_corner_positions(layout: FlowerLayout) -> None:
    mesh = build_flower_surface_mesh(
        FLAT_SIDES,
        layout,
        LEVELS.z,
        seed=3,
        road_water_side_pairs=[(0, 3)],
    )
    p_start = side_entry_point(layout, 0)
    p_end = side_entry_point(layout, 3)
    for x, y in (p_start, p_end):
        d_xy = ((mesh.vertices[:, 0] - x) ** 2 + (mesh.vertices[:, 1] - y) ** 2) ** 0.5
        assert d_xy.min() < 1e-6


def test_road_is_embedded_as_a_direct_edge_chain(layout: FlowerLayout) -> None:
    mesh = build_flower_surface_mesh(
        FLAT_SIDES,
        layout,
        LEVELS.z,
        seed=3,
        road_water_side_pairs=[(0, 3)],
        road_subdivisions=10,
    )
    p1 = side_entry_point(layout, 0)
    p2 = side_entry_point(layout, 3)
    points = subdivide_line(p1, p2, 10)
    vertex_indices = []
    for x, y in points:
        d_xy = ((mesh.vertices[:, 0] - x) ** 2 + (mesh.vertices[:, 1] - y) ** 2) ** 0.5
        idx = int(d_xy.argmin())
        assert d_xy[idx] < 1e-6
        vertex_indices.append(idx)
    for a, b in zip(vertex_indices, vertex_indices[1:]):
        assert _has_edge(mesh, a, b)


def test_road_continues_across_a_shared_flower_boundary() -> None:
    """Flower A's road exits at side_entry_point(A, 0); its true neighbor
    in direction 0 physically shares that exact point as
    side_entry_point(B, 4) (verified: side k's chain means B's side3-last
    == B's side4-first, and the reversed cross-flower contract maps that
    to A's side0-first). A road entering B there must land at the same
    XY and the same Z as where A's road leaves off.
    """
    layout = FlowerLayout()
    dx, dy = layout.neighbor_flower_offset(0)

    mesh_a = build_flower_surface_mesh(
        FLAT_SIDES, layout, LEVELS.z, seed=11, road_water_side_pairs=[(2, 0)]
    )
    mesh_b = build_flower_surface_mesh(
        FLAT_SIDES, layout, LEVELS.z, seed=12, road_water_side_pairs=[(4, 1)]
    )
    mesh_b.apply_translation((dx, dy, 0.0))

    shared_xy = side_entry_point(layout, 0)
    idx_a = _find_vertex(
        mesh_a,
        (shared_xy[0], shared_xy[1], mesh_a.vertices[:, 2].mean()),
        tol=1e9,  # XY-driven match; Z unknown ahead of time
    )
    # Re-locate by XY only for an exact XY match, then read its real Z.
    d_xy = (
        (mesh_a.vertices[:, 0] - shared_xy[0]) ** 2
        + (mesh_a.vertices[:, 1] - shared_xy[1]) ** 2
    ) ** 0.5
    idx_a = int(d_xy.argmin())
    assert d_xy[idx_a] < 1e-6
    a_point = mesh_a.vertices[idx_a]

    d_xy_b = (
        (mesh_b.vertices[:, 0] - a_point[0]) ** 2
        + (mesh_b.vertices[:, 1] - a_point[1]) ** 2
    ) ** 0.5
    idx_b = int(d_xy_b.argmin())
    assert d_xy_b[idx_b] < 1e-6
    b_point = mesh_b.vertices[idx_b]

    assert b_point[2] == pytest.approx(a_point[2], abs=1e-6)
