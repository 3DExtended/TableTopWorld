"""Phase 1: boundary-side grouping and true flower-to-flower adjacency.

These turn the grill session's brute-force-verified claims into permanent
regression tests against real FlowerLayout geometry, not assumptions:

- a flower's 18 exterior edges group into exactly 6 sides of 3 consecutive
  edges / 4 corners each (design decision #2)
- two flowers placed via a true neighbor offset actually share one side's
  4 corners exactly (the pre-mesh check for design decision #3/#4)
"""

import math

import pytest

from terrain.layout import FlowerLayout


@pytest.fixture
def layout() -> FlowerLayout:
    return FlowerLayout()


def test_exactly_six_sides_of_four_corners(layout: FlowerLayout) -> None:
    for side_idx in range(FlowerLayout.SIDE_COUNT):
        corners = layout.side_corners(side_idx)
        assert len(corners) == FlowerLayout.EDGES_PER_SIDE + 1


def test_side_groups_cover_all_18_exterior_edges_exactly_once(
    layout: FlowerLayout,
) -> None:
    groups = layout._side_groups()
    assert len(groups) == FlowerLayout.SIDE_COUNT
    all_keys = [e.key for group in groups for e in group]
    assert len(all_keys) == 18
    assert len(set(all_keys)) == 18  # no edge counted twice
    assert set(all_keys) == set(layout.exterior_edge_keys())


def test_invalid_side_index_raises(layout: FlowerLayout) -> None:
    with pytest.raises(ValueError):
        layout.side_corners(6)


def _points_close(a: tuple[float, float], b: tuple[float, float], tol: float = 1e-6) -> bool:
    return math.hypot(a[0] - b[0], a[1] - b[1]) < tol


def _corner_sets_match(
    a: tuple, b: tuple, tol: float = 1e-6
) -> bool:
    """True if two 4-corner chains describe the same boundary run, forward or reversed."""
    if all(_points_close(x, y, tol) for x, y in zip(a, b)):
        return True
    if all(_points_close(x, y, tol) for x, y in zip(a, reversed(b))):
        return True
    return False


@pytest.mark.parametrize("direction_idx", range(6))
def test_true_neighbor_flower_shares_exactly_one_side(
    layout: FlowerLayout, direction_idx: int
) -> None:
    """Two flowers placed via neighbor_flower_offset() must share exactly
    one side's 4 corners exactly - the geometric precondition for the whole
    corner-height boundary-matching contract (decisions #2-#4)."""
    dx, dy = layout.neighbor_flower_offset(direction_idx)

    def world_corners(side_idx: int, offset: tuple[float, float]) -> tuple:
        local = layout.side_corners(side_idx)
        return tuple((p[0] + offset[0], p[1] + offset[1]) for p in local)

    a_sides = {i: world_corners(i, (0.0, 0.0)) for i in range(FlowerLayout.SIDE_COUNT)}
    b_sides = {i: world_corners(i, (dx, dy)) for i in range(FlowerLayout.SIDE_COUNT)}

    matches = [
        (ai, bi)
        for ai, a_corners in a_sides.items()
        for bi, b_corners in b_sides.items()
        if _corner_sets_match(a_corners, b_corners)
    ]

    assert len(matches) == 1, (
        f"direction {direction_idx}: expected exactly one matching side pair, "
        f"found {matches}"
    )


def test_all_six_neighbor_directions_are_distinct_offsets(layout: FlowerLayout) -> None:
    offsets = [layout.neighbor_flower_offset(i) for i in range(6)]
    for i in range(6):
        for j in range(i + 1, 6):
            assert not _points_close(offsets[i], offsets[j], tol=1e-3), (i, j)
