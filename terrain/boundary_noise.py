"""Deterministic hash-based value noise for jagged boundary contours.

Bit-exact reproducibility (design decision #4) rules out Python's built-in
hash() (salted per-process by PYTHONHASHSEED) and any compiled/libm-based
noise library (not guaranteed bit-identical across platforms, machines, or
Python versions). This module uses only integer bit-mixing (a
splitmix64-style hash) and polynomial interpolation - no trig, no external
RNG state, no salted hashing.
"""

from __future__ import annotations

import math

_MASK64 = (1 << 64) - 1
_SEED_CONSTANT = 0x2545F4914F6CDD1D  # arbitrary fixed constant, never changes


def _splitmix64(x: int) -> int:
    """A fully-specified 64-bit bit-mixing hash: same input -> same output
    on any machine, any process, any Python version - unlike Python's
    built-in hash(), which is salted per-process by PYTHONHASHSEED."""
    x = (x + 0x9E3779B97F4A7C15) & _MASK64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & _MASK64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & _MASK64
    return (x ^ (x >> 31)) & _MASK64


def hash_to_unit_interval(*ints: int) -> float:
    """Deterministically hash a tuple of integers to a float in [0, 1)."""
    acc = _SEED_CONSTANT
    for value in ints:
        acc = _splitmix64(acc ^ (value & _MASK64))
    # Top 53 bits -> exact float64 mantissa precision, no rounding surprises.
    return (acc >> 11) / float(1 << 53)


def _smootherstep(t: float) -> float:
    """Ken Perlin's quintic smootherstep - polynomial only, no trig."""
    return t * t * t * (t * (t * 6 - 15) + 10)


def _lattice_jitter(seed_ints: tuple[int, ...], lattice_index: int) -> float:
    """Deterministic jitter in [-1, 1] for one coarse lattice point."""
    return hash_to_unit_interval(*seed_ints, lattice_index) * 2.0 - 1.0


def canonicalize_sequence_position(
    sequence: tuple[int, ...], position: float, total: float
) -> tuple[tuple[int, ...], float]:
    """Canonicalize a (corner-height sequence, position-along-run) pair.

    A boundary side declared forward by one flower is geometrically required
    to be declared as the REVERSE sequence by the neighbor on the other side
    (verified: a flower's side k always meets a neighbor's side (k+3)%6 in
    reversed corner order). Without canonicalizing before hashing, two
    flowers whose sides match at the 4 coarse corners would still produce
    uncorrelated jitter in between them - matching numbers alone would not
    actually guarantee matching physical geometry, only matching endpoints.

    A PALINDROMIC sequence (e.g. the very common all-flat (0,0,0,0)) is a
    real edge case here, not just a hypothetical one: sequence ==
    reversed_seq as VALUES, so comparing them carries no information about
    which of the two flowers is walking the shared edge "forward" - both
    sides declare the literal same tuple. The two calls (S, t) and (S,
    total - t) must still land on the same canonical point, so position
    itself has to be mirrored around the run's own midpoint in this case
    (min(position, total - position)) rather than left as-is. Getting this
    wrong doesn't fail loudly - it silently produces a boundary that only
    matches its neighbor at the exact midpoint, mismatched everywhere else
    (found by visually inspecting an assembled two-flower preview whose
    matching side happened to be all-zero, i.e. exactly this case).
    """
    reversed_seq = tuple(reversed(sequence))
    if sequence < reversed_seq:
        return sequence, position
    if sequence > reversed_seq:
        return reversed_seq, total - position
    return sequence, min(position, total - position)


def sample_noise_1d(
    seed_ints: tuple[int, ...],
    position: float,
    total: float,
    *,
    lattice_step: float = 4.0,
) -> float:
    """Deterministic, smoothly-varying noise in [-1, 1] at a continuous 1D
    position along a run of length `total`, built by interpolating between
    coarse hashed lattice points - correlated ruggedness, not per-point
    static. Canonicalized first (see canonicalize_sequence_position) so the
    result is identical regardless of which direction the caller happens to
    be walking the shared boundary from.
    """
    seed_ints, position = canonicalize_sequence_position(seed_ints, position, total)
    lattice_pos = position / lattice_step
    i0 = math.floor(lattice_pos)
    i1 = i0 + 1
    t = lattice_pos - i0
    v0 = _lattice_jitter(seed_ints, i0)
    v1 = _lattice_jitter(seed_ints, i1)
    return v0 + (v1 - v0) * _smootherstep(t)


def sample_noise_2d(seed: int, x: float, y: float, *, lattice_step: float = 6.0) -> float:
    """Deterministic, smoothly-varying 2D value noise in [-1, 1], seeded by
    a per-flower integer seed - used only for freeform interior sculpting
    (design decision #5), never for the boundary contract (decision #4
    uses sample_noise_1d instead, with its own canonicalization)."""
    lx, ly = x / lattice_step, y / lattice_step
    ix0, iy0 = math.floor(lx), math.floor(ly)
    ix1, iy1 = ix0 + 1, iy0 + 1
    tx, ty = _smootherstep(lx - ix0), _smootherstep(ly - iy0)

    def corner(ix: int, iy: int) -> float:
        return hash_to_unit_interval(seed, ix, iy) * 2.0 - 1.0

    v00, v10 = corner(ix0, iy0), corner(ix1, iy0)
    v01, v11 = corner(ix0, iy1), corner(ix1, iy1)
    vx0 = v00 + (v10 - v00) * tx
    vx1 = v01 + (v11 - v01) * tx
    return vx0 + (vx1 - vx0) * ty
