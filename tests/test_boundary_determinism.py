"""Phase 2: the core hypothesis - deterministic boundary contours.

Design decision #4 is that two independently-generated flowers whose
shared side declares matching corner heights produce bit-identical
physical geometry along the ENTIRE fine jagged contour, not just at the
4 coarse corners, with nothing exchanged between them but those 4
numbers. If this doesn't hold exactly, nothing built on top of it
(Phase 3 onward) is worth building.

A flower's side k always meets a neighbor's side (k+3)%6 in REVERSED
corner order (verified computationally against real FlowerLayout
geometry in Phase 1) - so the neighbor must declare the reversed
corner-height sequence for its own local walk direction. The trickiest
part of this contract is that plain per-point noise keyed on
(corner_heights, position) would satisfy only the 4 corner heights and
produce uncorrelated jitter in between; terrain/boundary_noise.py
canonicalizes (sequence, position) before hashing specifically to rule
that out, and several tests below exist only to catch a regression of
that canonicalization.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap

from terrain.heightfield import build_side_boundary_vertices

LEVEL_MM = 15.0


def level_z(level: int) -> float:
    return level * LEVEL_MM


SIDE_GEOM = [(0.0, 0.0), (10.0, 0.0), (20.0, 0.0), (30.0, 0.0)]


def test_independent_builds_are_bit_identical() -> None:
    corner_heights = (1, 2, 0, 3)
    a = build_side_boundary_vertices(corner_heights, SIDE_GEOM, level_z)
    b = build_side_boundary_vertices(corner_heights, list(SIDE_GEOM), level_z)
    assert a == b


def test_changing_a_corner_height_changes_the_output() -> None:
    base = build_side_boundary_vertices((1, 2, 0, 3), SIDE_GEOM, level_z)
    changed = build_side_boundary_vertices((1, 2, 1, 3), SIDE_GEOM, level_z)
    assert base != changed


def test_reversed_neighbor_declaration_matches_physical_contour() -> None:
    """The actual cross-flower scenario: flower A walks a side p0->p3 with
    corner_heights (1,2,3,4); its neighbor shares that physical boundary
    but walks it p3->p0 locally, so per the verified (k+3)%6-reversed
    contract it must declare (4,3,2,1) for its own local walk. The two
    independently-built contours must describe the exact same physical
    curve - same points, same heights, same jitter - just enumerated in
    opposite order.

    Heights are deliberately non-palindromic here - see
    test_palindromic_sequence_also_matches_across_the_shared_boundary
    below for the palindrome case, which is NOT just a hypothetical
    variant of this one: it was broken in production (found by visually
    inspecting an assembled two-flower preview whose matching side
    happened to be the very common all-flat (0,0,0,0)) despite this test
    passing throughout, because a palindrome exercises a genuinely
    different code path in canonicalize_sequence_position (the
    sequence == reversed_seq branch, never hit by a non-palindrome).
    """
    forward = build_side_boundary_vertices((1, 2, 3, 4), SIDE_GEOM, level_z)
    reversed_geom = list(reversed(SIDE_GEOM))
    backward = build_side_boundary_vertices((4, 3, 2, 1), reversed_geom, level_z)
    assert forward == list(reversed(backward))


def test_palindromic_sequence_also_matches_across_the_shared_boundary() -> None:
    """Regression test for a real bug: canonicalize_sequence_position's
    `sequence <= reversed_seq` check could not distinguish "forward" from
    "backward" for a palindromic sequence like (0,0,0,0) or (1,2,2,1),
    since a palindrome's reverse is the literal same tuple - so it never
    flipped position for either flower, and two flowers sharing a
    palindromic side only matched at the exact midpoint of the run,
    mismatched everywhere else. A flat_plains-style all-zero side is the
    most common real case this affects.

    Unlike the non-palindrome case above, BOTH sides of a shared
    palindromic boundary declare the identical sequence (there is no
    separate "reversed" declaration to write, since reversing a
    palindrome doesn't change it) - so this test passes the SAME sequence
    to both calls, only reversing the geometry (side_geom) to model the
    neighbor's opposite walk direction.
    """
    for seq in [(0, 0, 0, 0), (1, 2, 2, 1), (3, 3, 3, 3)]:
        forward = build_side_boundary_vertices(seq, SIDE_GEOM, level_z)
        reversed_geom = list(reversed(SIDE_GEOM))
        backward = build_side_boundary_vertices(seq, reversed_geom, level_z)
        assert forward == list(reversed(backward)), f"mismatch for palindrome {seq}"


def test_xy_jitter_also_matches_across_the_shared_boundary() -> None:
    """Regression test for the sideways (in-plane) jitter added on top of
    the existing Z-only jitter: a naive "rotate my own p0->p1 direction by
    90 degrees" perpendicular would put the wiggle on physically OPPOSITE
    sides for the two flowers sharing this edge, since their local walk
    directions along the same physical segment are always opposite (see
    build_side_boundary_vertices' own docstring) - true even for a
    palindrome, where the sequence carries no directional information at
    all. Covers both the non-palindrome and palindrome cases, since they
    take different paths through canonicalize_sequence_position."""
    forward = build_side_boundary_vertices(
        (1, 2, 3, 4), SIDE_GEOM, level_z, xy_jitter_mm=2.0
    )
    reversed_geom = list(reversed(SIDE_GEOM))
    backward = build_side_boundary_vertices(
        (4, 3, 2, 1), reversed_geom, level_z, xy_jitter_mm=2.0
    )
    assert forward == list(reversed(backward))

    for seq in [(0, 0, 0, 0), (1, 2, 2, 1), (3, 3, 3, 3)]:
        fwd = build_side_boundary_vertices(seq, SIDE_GEOM, level_z, xy_jitter_mm=2.0)
        bwd = build_side_boundary_vertices(
            seq, reversed_geom, level_z, xy_jitter_mm=2.0
        )
        assert fwd == list(reversed(bwd)), f"mismatch for palindrome {seq}"


def test_xy_jitter_leaves_declared_corners_exactly_in_place() -> None:
    """The sideways offset must taper to 0 at the true corners (t=0, t=1)
    - many other code paths (hex_edge_line, side_corners, road entry
    points, the canonical-vertex registry) key off the exact declared
    corner XY, so moving it would break far more than just this
    function."""
    verts = build_side_boundary_vertices(
        (1, 2, 3, 4), SIDE_GEOM, level_z, xy_jitter_mm=5.0
    )
    assert verts[0][:2] == SIDE_GEOM[0]
    assert verts[8][:2] == SIDE_GEOM[1]
    assert verts[16][:2] == SIDE_GEOM[2]
    assert verts[-1][:2] == SIDE_GEOM[3]


def test_forward_declaration_on_both_sides_does_not_match() -> None:
    """Sanity check on the contract itself: declaring the SAME (not
    reversed) sequence on both sides of a shared boundary must NOT
    produce a matching contour - otherwise the reversed-declaration test
    above would be proving nothing distinctive."""
    forward = build_side_boundary_vertices((1, 2, 3, 4), SIDE_GEOM, level_z)
    reversed_geom = list(reversed(SIDE_GEOM))
    also_forward = build_side_boundary_vertices((1, 2, 3, 4), reversed_geom, level_z)
    assert forward != list(reversed(also_forward))


def test_wrong_edge_count_raises() -> None:
    import pytest

    with pytest.raises(ValueError):
        build_side_boundary_vertices((1, 2, 3), SIDE_GEOM, level_z)


_SUBPROCESS_SCRIPT = textwrap.dedent(
    """
    import json
    import sys
    sys.path.insert(0, {repo_root!r})
    from terrain.heightfield import build_side_boundary_vertices

    def level_z(level):
        return level * 15.0

    side_geom = [(0.0, 0.0), (10.0, 0.0), (20.0, 0.0), (30.0, 0.0)]
    vertices = build_side_boundary_vertices((1, 2, 0, 3), side_geom, level_z)
    print(json.dumps(vertices))
    """
)


def test_determinism_is_independent_of_pythonhashseed(tmp_path) -> None:
    """Guards against accidentally relying on Python's salted built-in
    hash() instead of the fully-specified splitmix64 mixer - a real risk
    since hash() would look identical to a passing test within one
    process but silently break across separate print runs."""
    import pathlib

    repo_root = str(pathlib.Path(__file__).resolve().parent.parent)
    script = tmp_path / "run_once.py"
    script.write_text(_SUBPROCESS_SCRIPT.format(repo_root=repo_root))

    import os

    outputs = []
    for seed in ("0", "1", "12345"):
        env = dict(os.environ)
        env["PYTHONHASHSEED"] = seed
        result = subprocess.run(
            [sys.executable, str(script)],
            capture_output=True,
            text=True,
            env=env,
            check=True,
        )
        outputs.append(result.stdout.strip())

    assert len(set(outputs)) == 1, "output differs across PYTHONHASHSEED values"
