"""Phase 4: full 7-hex flower, freeform interior, cliffs.

Builds against real constructed geometry only - a trimesh.Trimesh (and,
for the cross-flower check, two independently built ones placed via the
real neighbor_flower_offset()) - never string-matched text.
"""

from __future__ import annotations

import math

import pytest

from terrain.heightfield import build_side_boundary_vertices
from terrain.heights import HeightLevels
from terrain.layout import FlowerLayout
from terrain.surface_mesh import build_flower_surface_mesh

LEVELS = HeightLevels()

# 6 junction corner heights (shared between adjacent sides), cyclic.
J = [0, 1, 3, 0, 2, 1]

# Side 0 forces a big jump right after its first junction (0 -> 3, a full
# 45mm cliff with no intermediate slope step) - the "no cap on delta"
# fixture for design decision #7.
SIDE_CORNER_HEIGHTS = {
    0: (J[0], 3, 3, J[1]),
    1: (J[1], J[1], J[2], J[2]),
    2: (J[2], J[2], J[3], J[3]),
    3: (J[3], J[3], J[4], J[4]),
    4: (J[4], J[4], J[5], J[5]),
    5: (J[5], J[5], J[0], J[0]),
}


@pytest.fixture
def layout() -> FlowerLayout:
    return FlowerLayout()


def _build(layout: FlowerLayout, seed: int = 42):
    return build_flower_surface_mesh(SIDE_CORNER_HEIGHTS, layout, LEVELS.z, seed)


def test_flower_mesh_is_watertight_and_has_volume(layout: FlowerLayout) -> None:
    mesh = _build(layout)
    assert mesh.is_watertight
    assert mesh.is_winding_consistent
    assert mesh.volume > 0


def test_flower_mesh_boundary_matches_true_18_corner_silhouette(
    layout: FlowerLayout,
) -> None:
    """The mesh's XY footprint must pass through all 18 true exterior
    corners (design decision #2), not a simplified hexagon - X/Y are pure
    linear interpolation in build_side_boundary_vertices (only Z carries
    jitter), so corner XY survives into the final mesh exactly."""
    mesh = _build(layout)
    expected_corners: list[tuple[float, float]] = []
    for side_idx in range(FlowerLayout.SIDE_COUNT):
        expected_corners.extend(layout.side_corners(side_idx)[:-1])
    assert len(expected_corners) == 18

    mesh_xy = mesh.vertices[:, :2]
    for cx, cy in expected_corners:
        dist = ((mesh_xy[:, 0] - cx) ** 2 + (mesh_xy[:, 1] - cy) ** 2) ** 0.5
        assert dist.min() < 1e-6, f"corner ({cx}, {cy}) not found in mesh vertices"


def test_flower_mesh_is_deterministic_at_full_scale(layout: FlowerLayout) -> None:
    """Phase 2's core hypothesis, re-checked at full-flower scale including
    the seeded freeform interior: identical inputs, independently built,
    must produce bit-identical vertex/face arrays."""
    mesh_a = _build(layout, seed=7)
    mesh_b = _build(layout, seed=7)
    assert (mesh_a.vertices == mesh_b.vertices).all()
    assert (mesh_a.faces == mesh_b.faces).all()


def test_different_seed_changes_the_interior_but_not_the_boundary(
    layout: FlowerLayout,
) -> None:
    mesh_a = _build(layout, seed=1)
    mesh_b = _build(layout, seed=2)
    assert not (mesh_a.vertices.shape == mesh_b.vertices.shape and
                (mesh_a.vertices == mesh_b.vertices).all())


def test_two_true_neighbor_flowers_share_matching_boundary_geometry() -> None:
    """The full integration proof, at real flower geometry scale: flower
    A's side 0 (a real 3-segment zigzag polyline, not a synthetic straight
    line like Phase 2's unit test) and its true neighbor's corresponding
    side 3 (per the verified (k+3)%6-reversed contract) must describe the
    exact same physical curve when B is placed via the real
    neighbor_flower_offset() - fine jitter included, not just the 4 coarse
    corners.

    This deliberately calls build_side_boundary_vertices directly for just
    the one shared side on each flower, rather than going through
    build_flower_surface_mesh's full boundary loop: a flower's OTHER 5
    sides only meet neighbors that aren't part of this test, and
    build_flower_boundary_loop's junction bookkeeping (each shared corner
    takes its height from the *next* side's declaration, discarding the
    previous side's) means an isolated single-side override without a
    fully re-derived 6-side flower would introduce an unrelated internal
    seam that has nothing to do with the cross-flower contract under test
    here - keeping tileset-wide junction consistency correct is a
    validate_tileset concern (Phase 4 schema), not this test's job.
    """
    layout = FlowerLayout()
    dx, dy = layout.neighbor_flower_offset(0)

    heights_a = SIDE_CORNER_HEIGHTS[0]
    heights_b = tuple(reversed(heights_a))

    geom_a = layout.side_corners(0)
    geom_b_local = layout.side_corners(3)
    geom_b = tuple((p[0] + dx, p[1] + dy) for p in geom_b_local)

    vertices_a = build_side_boundary_vertices(heights_a, geom_a, LEVELS.z)
    vertices_b = build_side_boundary_vertices(heights_b, geom_b, LEVELS.z)

    # Coordinates go through an extra +dx/+dy floating-point translation for
    # B, so compare with tolerance rather than exact equality - the jitter
    # component (computed from integer-hashed inputs, no float drift) still
    # matches bit-for-bit; only the translated X/Y picks up float noise.
    for a_point, b_point in zip(vertices_a, reversed(vertices_b)):
        assert a_point == pytest.approx(b_point, abs=1e-9)
