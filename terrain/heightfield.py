"""Deterministic boundary-contour construction - design decisions #2-#4.

build_side_boundary_vertices is the single most important function in the
redesign: a pure function of (corner heights, side geometry, subdivision
count) whose output is bit-identical across independently-built flowers
whenever those inputs match, so matching numbers alone guarantee a physical
match with nothing else exchanged or stored.
"""

from __future__ import annotations

import math
from typing import Callable, Sequence

from terrain.boundary_noise import sample_noise_1d, sample_noise_2d
from terrain.layout import FlowerLayout
from terrain.triangulate import point_in_polygon

Point2D = tuple[float, float]
Vertex3D = tuple[float, float, float]
Triangle = tuple[int, int, int]

# Fraction of a hex cell's AREA held perfectly flat and noise-free at the
# cell's declared height_level - its "plateau", the part a miniature can
# actually stand on. The remaining area is an outer band that blends the
# plateau into whatever the cell's shared edges are doing (a neighbouring
# cell's own plateau, or the flower's declared silhouette contour), so
# height differences between cells still read as real steps/cliffs.
#
# Peter's requirement: "the hex does not need to be flat everywhere, but at
# least 2/3 should be flat (and without noise too)". Scaling a hexagon
# about its centre by k scales its area by k^2, so a 2/3-AREA plateau
# reaches sqrt(2/3) ~= 0.816 of the way from centre to edge.
PLATEAU_AREA_FRAC = 2.0 / 3.0


def build_side_boundary_vertices(
    corner_heights: tuple[int, int, int, int],
    side_geom: Sequence[Point2D],
    level_z: Callable[[int], float],
    *,
    subdivisions_per_edge: int = 8,
    jitter_amplitude: float = 0.3,
    xy_jitter_mm: float = 0.0,
) -> list[Vertex3D]:
    """Build the fine, jagged 3D contour along one flower side.

    Pure deterministic function of (corner_heights, side_geom, subdivision
    count) only - same inputs always produce bit-identical output (design
    decision #4). corner_heights and side_geom must describe the boundary
    walked in the same direction (4 points, 3 edges - see
    FlowerLayout.side_corners); a neighbor sharing this side must declare
    the reversed corner-height sequence for its own local walk direction
    (verified: a flower's side k always meets a neighbor's side (k+3)%6 in
    reversed corner order) - build_side_boundary_vertices itself does not
    need to know about this, since terrain/boundary_noise.py canonicalizes
    (sequence, position) before hashing so the reversed declaration still
    produces a bit-identical physical contour.

    `xy_jitter_mm` additionally displaces each point SIDEWAYS (in-plane,
    perpendicular to the local edge), so the boundary's XY path itself
    wiggles instead of only its Z height - matching decision #4's
    determinism the same way Z-jitter does, but a signed 2D displacement
    has a handedness problem Z (an unsigned scalar) never had: naively
    using "rotate MY OWN p0->p1 direction by 90 degrees" would put the
    bump on physically OPPOSITE sides for the two flowers sharing this
    edge, because their local walk directions along the same physical
    segment are always opposite (this holds even for a palindromic
    corner-height sequence - the walk-direction mismatch is a pure
    geometry fact, unrelated to what heights are declared, so the
    sequence-based tie-break canonicalize_sequence_position uses for
    position can't resolve it). Fixed by deriving the perpendicular's
    sign from the sub-edge's own two endpoint coordinates instead of
    "which one is p0 vs p1": always point from the lexicographically
    smaller local point to the larger one. Lexicographic point comparison
    is translation-invariant (subtracting either flower's own placement
    offset from both endpoints before comparing doesn't change which one
    sorts first), so both flowers derive the identical world-space
    perpendicular for the physically shared sub-edge regardless of which
    one of them calls p0 vs p1, sequence, or reversal. The corners
    themselves (t=0 and t=1) are never displaced - only interior points
    move - since many other code paths key off the exact declared corner
    XY (hex_edge_line, side_corners, road entry points, the canonical-
    vertex registry).
    """
    n_edges = len(side_geom) - 1
    if len(corner_heights) != n_edges + 1:
        raise ValueError(
            f"corner_heights must have {n_edges + 1} values for "
            f"{n_edges} edges, got {len(corner_heights)}"
        )
    one_level_z = level_z(1) - level_z(0)
    total_positions = float(n_edges * subdivisions_per_edge)

    vertices: list[Vertex3D] = []
    for edge_idx in range(n_edges):
        p0, p1 = side_geom[edge_idx], side_geom[edge_idx + 1]
        edge_dx, edge_dy = p1[0] - p0[0], p1[1] - p0[1]
        edge_len = math.hypot(edge_dx, edge_dy)
        canon_dx, canon_dy = edge_dx, edge_dy
        if edge_len > 1e-9 and (round(p1[0], 6), round(p1[1], 6)) < (
            round(p0[0], 6),
            round(p0[1], 6),
        ):
            canon_dx, canon_dy = -edge_dx, -edge_dy
        perp_x, perp_y = (
            (-canon_dy / edge_len, canon_dx / edge_len) if edge_len > 1e-9 else (0.0, 0.0)
        )
        z0 = level_z(corner_heights[edge_idx])
        z1 = level_z(corner_heights[edge_idx + 1])
        is_last_edge = edge_idx == n_edges - 1
        step_count = subdivisions_per_edge + (1 if is_last_edge else 0)
        for step in range(step_count):
            t = step / subdivisions_per_edge
            base_x = p0[0] + edge_dx * t
            base_y = p0[1] + edge_dy * t
            base_z = z0 + (z1 - z0) * t
            global_pos = edge_idx * subdivisions_per_edge + step
            jitter = sample_noise_1d(
                corner_heights, float(global_pos), total_positions, channel=0
            )
            z = base_z + jitter * jitter_amplitude * one_level_z
            window = 4.0 * t * (1.0 - t)  # 0 at the true corners, 1 at t=0.5
            xy_jitter = sample_noise_1d(
                corner_heights, float(global_pos), total_positions, channel=1
            )
            offset = xy_jitter * window * xy_jitter_mm
            x = base_x + perp_x * offset
            y = base_y + perp_y * offset
            vertices.append((x, y, z))
    return vertices


def build_flower_boundary_loop(
    side_corner_heights: dict[int, tuple[int, int, int, int]],
    layout: FlowerLayout,
    level_z: Callable[[int], float],
    *,
    subdivisions_per_edge: int = 8,
    jitter_amplitude: float = 0.3,
    xy_jitter_mm: float = 0.0,
) -> list[Vertex3D]:
    """The full closed flower boundary contour: one build_side_boundary_vertices
    run per side, concatenated. FlowerLayout.side_corners()'s 6 sides chain
    corner-to-corner around the flower (side k's 4th corner is exactly side
    (k+1)%6's 1st corner - verified against the real geometry), so each
    side's last point is dropped before appending the next side's run to
    avoid a duplicate vertex at the shared corner.
    """
    loop: list[Vertex3D] = []
    for side_idx in range(FlowerLayout.SIDE_COUNT):
        corners = layout.side_corners(side_idx)
        heights = side_corner_heights[side_idx]
        verts = build_side_boundary_vertices(
            heights,
            corners,
            level_z,
            subdivisions_per_edge=subdivisions_per_edge,
            jitter_amplitude=jitter_amplitude,
            xy_jitter_mm=xy_jitter_mm,
        )
        loop.extend(verts[:-1])
    return loop


def _interior_grid_points(
    boundary_2d: Sequence[Point2D], step: float
) -> list[Point2D]:
    """A regular grid of points strictly inside the boundary polygon, used
    as extra Steiner points for the freeform interior (decision #5). Points
    outside the boundary are simply dropped - triangulate_polygon already
    needs a point-in-polygon filter for its own triangles, so a grid point
    too close to the boundary just yields a smaller/uneven triangle there,
    not an error."""
    xs = [p[0] for p in boundary_2d]
    ys = [p[1] for p in boundary_2d]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    margin = step * 0.5
    points: list[Point2D] = []
    x = min_x + margin
    while x < max_x - margin:
        y = min_y + margin
        while y < max_y - margin:
            if point_in_polygon((x, y), boundary_2d):
                points.append((x, y))
            y += step
        x += step
    return points


def build_flower_pslg(
    side_corner_heights: dict[int, tuple[int, int, int, int]],
    layout: FlowerLayout,
    level_z: Callable[[int], float],
    seed: int,
    *,
    subdivisions_per_edge: int = 8,
    jitter_amplitude: float = 0.3,
    xy_jitter_mm: float = 0.0,
    interior_grid_step: float | None = None,
    interior_relief_mm: float = 6.0,
) -> tuple[list[Point2D], list[Vertex3D], int]:
    """The full per-flower planar-straight-line-graph: the deterministic
    18-edge boundary loop (decision #4) plus a freeform seeded interior
    grid (decision #5, "never affecting the boundary contract" - interior
    height is independent noise, only loosely anchored to the boundary's
    own average height so it doesn't wander off to an arbitrary absolute
    level).

    Returns (all_points_2d, all_vertices_3d, boundary_point_count): the
    first `boundary_point_count` entries are the boundary loop in order,
    the remainder are interior grid points in grid-scan order.
    """
    boundary_3d = build_flower_boundary_loop(
        side_corner_heights,
        layout,
        level_z,
        subdivisions_per_edge=subdivisions_per_edge,
        jitter_amplitude=jitter_amplitude,
        xy_jitter_mm=xy_jitter_mm,
    )
    boundary_2d = [(v[0], v[1]) for v in boundary_3d]
    # Corner samples (one per declared side corner, at every
    # subdivisions_per_edge-th boundary point) drive an inverse-distance-
    # weighted height blend below, so interior points near a high side of
    # the boundary get pulled up with it and points near a low side get
    # pulled down - see build_flower_cells' interior_base_z for why a flat
    # single average was wrong.
    corner_samples = boundary_3d[::subdivisions_per_edge]

    def interior_base_z(x: float, y: float) -> float:
        total_w = 0.0
        total_wz = 0.0
        for bx, by, bz in corner_samples:
            w = 1.0 / ((x - bx) ** 2 + (y - by) ** 2 + 1.0)
            total_w += w
            total_wz += w * bz
        return total_wz / total_w

    step = interior_grid_step or (layout.hex_outer_width * 0.6)
    interior_2d = _interior_grid_points(boundary_2d, step)
    interior_3d = [
        (x, y, interior_base_z(x, y) + sample_noise_2d(seed, x, y) * interior_relief_mm)
        for x, y in interior_2d
    ]

    all_points_2d = boundary_2d + interior_2d
    all_vertices_3d = boundary_3d + interior_3d
    return all_points_2d, all_vertices_3d, len(boundary_2d)


def build_flower_cells(
    side_corner_heights: dict[int, tuple[int, int, int, int]],
    layout: FlowerLayout,
    level_z: Callable[[int], float],
    seed: int,
    *,
    subdivisions_per_edge: int = 8,
    jitter_amplitude: float = 0.3,
    xy_jitter_mm: float = 0.0,
    interior_relief_mm: float = 6.0,
    groove_depth_mm: float = 0.0,
    groove_width_mm: float = 2.0,
    groove_profile: list[tuple[float, float]] | None = None,
    hex_height_levels: dict[int, int] | None = None,
    plateau_area_frac: float = PLATEAU_AREA_FRAC,
) -> tuple[list[Vertex3D], list[Vertex3D], list[Triangle]]:
    """Triangulate each of the flower's 7 hex cells INDEPENDENTLY, as its
    own simple hexagon, instead of one Delaunay pass over the whole
    flower's interior.

    `groove_depth_mm` (decision #10) carves a real, visible/tactile
    recessed channel tracing each cell's own 6-edge outline, so hex-cell
    boundaries are distinguishable at a glance (and by touch) for
    counting distances - not just a guaranteed mesh edge along the seam
    (which is all "groove" meant before: a topological guarantee, not an
    actual depression in the surface). It's built as a dedicated thin
    strip of NEW geometry per wedge, between the true edge (untouched -
    still the literal shared vertex the boundary contract/neighboring
    cell rely on) and an "inset" row pulled exactly `groove_width_mm`
    toward the apex at full undisturbed height, which becomes the coarse
    interior grid's new outer boundary. This keeps the depression confined
    to a strip of the requested width regardless of subdivisions_per_edge
    - an earlier version depressed whichever existing interior grid rows
    happened to fall within groove_width_mm of the edge, which could
    (and, at the shipped default resolution, did) vanish to a completely
    invisible 0mm if a single grid step was already wider than
    groove_width_mm. The offset/depth taper to exactly 0 at each wedge
    corner (a `min(t, 1-t)` window over the strip's own j=0..n parameter),
    so the strip's corner points degenerate to the corner itself - already
    shared via spokes - with no separate stitching needed where 3 wedges'
    worth of strips would otherwise meet.

    This exists specifically so hex-to-hex boundaries can be engraved as
    real grooves (decision #10): triangulating per-cell, rather than one
    whole-flower pass, sidesteps the highly symmetric shared vertices
    (hex0's 6 corners, each touching 3 cells at once) that plain
    unconstrained Delaunay, even after heavy subdivision and several
    attempted symmetry-breaking jitter schemes, could not be made to
    reliably include as direct edges - see git history for that approach.

    Each cell is split into 6 wedges (its center, plus one of its 6
    edges), and each wedge is triangulated as a proper triangular
    lattice - the standard "subdivide a triangle into n^2 smaller
    triangles" grid, at subdivisions_per_edge granularity along all 3 of
    the wedge's sides (2 spokes from the center, 1 boundary edge) - not a
    plain fan from a single vertex. A fan was tried first: simpler, and
    correct in the sense of being watertight, but it only adds detail
    radiating from one point, which reads as visibly low-poly/spiky
    rather than genuinely subdivided (the exact complaint a human raised
    looking at the actual STL, with a reference image of the intended
    triangular-grid look). Two triangulation strategies were tried and
    rejected before the fan, for the SAME cell-boundary geometry (still
    relevant background for why a fan, not Delaunay/ear-clipping, was the
    fan-era fallback, and why the lattice below reuses the exact same
    "shared vertex object" discipline rather than re-litigating either
    problem): scipy Delaunay produced triangles whose centroid fell
    outside the hexagon on the exactly-collinear fine subdivision points
    every boundary edge carries (jitter only perturbs Z, not X/Y); ear-
    clipping (mapbox_earcut) had no trouble with that WITHIN one cell's
    own triangulation, but once internal (hex-to-hex) edges also carry
    fine subdivision points, a shared edge is triangulated INDEPENDENTLY
    by the two cells on either side of it, and ear-clipping is free to
    legitimately bridge straight across an exactly-collinear intermediate
    point with a direct chord instead of routing through it - the two
    independent cells generally made different choices there, leaving the
    shared chord double-counted (an edge used by 4 triangles instead of
    2). The triangular-lattice grid below sidesteps both failure modes
    the same way the fan did: every grid point on a wedge's boundary
    (either spoke, or the far edge) is either the cell's own already-
    computed apex/spoke vertex or an already-computed edge_pts entry -
    literally the same shared vertex object two neighboring wedges (or
    cells) resolve to, never two independently-computed values hoped to
    be equal - so there is no ambiguous choice for any triangulation
    method to make differently on either side of a shared edge.

    Each hex's exterior-facing edges reuse build_side_boundary_vertices
    with the exact same inputs the flower's own boundary loop uses (same
    corner heights, same slice of the same side's geometry), so those
    points are bit-identical to build_flower_boundary_loop's output and
    weld seamlessly into the same wall construction. Interior-facing
    edges and any interior Steiner points get their height from
    sample_noise_2d, a pure function of (seed, x, y) - so two hexes
    independently querying the SAME shared vertex always agree, and
    trimesh's vertex welding (merge_vertices) stitches the 7 patches into
    one continuous surface.

    Returns (all_vertices_3d, boundary_indices, triangles): boundary_indices
    is a list of indices into all_vertices_3d, in the same corner-to-corner
    order as build_flower_boundary_loop's own output, for wall/cap
    construction. Deliberately NOT a second call to
    build_flower_boundary_loop returning its own fresh vertex list: even
    though that would be a pure function of the same inputs and so
    mathematically bit-identical, in practice trimesh's vertex-merge during
    Trimesh(process=True) did not reliably weld the two independently-
    computed copies of each exterior boundary point (the same float-noise
    phenomenon documented in terrain/roads.py's internal_hex_edges), which
    left the top surface's boundary and the wall's top edge as two
    disconnected sets of vertices at the same position - a non-watertight
    mesh. Returning indices into the ONE already-built vertex list side-
    steps the welding question entirely: there is only ever one vertex at
    each boundary position, referenced by both the top surface and the
    walls.
    """
    exterior_local_indices = {
        h: layout.exterior_vertex_indices(h) for h in layout.RING_HEX_INDICES
    }

    # Precompute each side's full fine contour once (bit-identical to
    # build_flower_boundary_loop's per-side calls).
    side_points: dict[int, list[Vertex3D]] = {}
    for side_idx in range(FlowerLayout.SIDE_COUNT):
        side_points[side_idx] = build_side_boundary_vertices(
            side_corner_heights[side_idx],
            layout.side_corners(side_idx),
            level_z,
            subdivisions_per_edge=subdivisions_per_edge,
            jitter_amplitude=jitter_amplitude,
            xy_jitter_mm=xy_jitter_mm,
        )

    # local_edge_idx -> (fine 3D points, (side_idx, abs start position in
    # side_points[side_idx])), per ring hex. The abs position lets us later
    # recover, for each fine point, its place in the flower-wide boundary
    # order without a second geometry computation.
    exterior_fine_points: dict[int, dict[int, list[Vertex3D]]] = {
        h: {} for h in layout.RING_HEX_INDICES
    }
    exterior_fine_meta: dict[int, dict[int, tuple[int, int]]] = {
        h: {} for h in layout.RING_HEX_INDICES
    }
    for side_idx in range(FlowerLayout.SIDE_COUNT):
        group = layout._side_groups()[side_idx]
        for pos, e in enumerate(group):
            start = pos * subdivisions_per_edge
            end = start + subdivisions_per_edge + 1
            local_edge_idx = exterior_local_indices[e.hex_idx][e.side_idx]
            exterior_fine_points[e.hex_idx][local_edge_idx] = side_points[side_idx][
                start:end
            ]
            exterior_fine_meta[e.hex_idx][local_edge_idx] = (side_idx, start)

    # The flower's 18 true silhouette corners (4 per side, positions
    # 0/sd/2sd/3sd - dedup happens naturally via the rounded-key dict since
    # consecutive sides share their junction corner exactly). A ring hex's
    # own local edge numbering doesn't line up with "exterior" vs "interior"
    # at these points: the edge immediately following (or preceding) that
    # hex's own 3-edge exterior arc is classified interior for THIS hex, but
    # still starts exactly at one of these 18 corners - a real flower
    # boundary corner with an EXACT declared height (level_z, no noise), not
    # an interior Steiner point. Looking it up here (rather than always
    # calling interior_z) is what keeps that corner's height consistent
    # between "the hex whose exterior arc ends there" and "the hex whose
    # interior groove edge starts there" - without this, the two would
    # independently disagree on that corner's Z, leaving the mesh with a
    # real hole (found and fixed via this exact symptom: a duplicate vertex
    # position with two very different Z values, at each of the flower's 6
    # main side-to-side junctions).
    # A single canonical-vertex registry, keyed by rounded (x, y), used for
    # EVERY point along an internal (hex-to-hex) edge - not just the 18
    # true exterior corners. Seeded from the exterior boundary's own
    # points so an internal edge that happens to touch one of those 18
    # corners reuses its exact height (no noise) rather than
    # recomputing one independently and disagreeing.
    #
    # This also covers a subtler case only the internal-edge subdivision
    # above introduces: a "triple point" where hex0's edge to hexA, and
    # a DIFFERENT edge (hexA's own edge to hexB), both end at the same
    # physical corner. Each hex computes that corner's (x, y) via its own
    # independent layout.hex_edge_line() call - mathematically the same
    # point, but not always bit-identical (the same non-associative-float
    # issue documented on internal_edge_points below, one level up).
    # Resolving every point through this ONE registry (first caller wins,
    # everyone else reuses that exact stored (x, y, z)) guarantees any two
    # edges that touch "the same" rounded position end up with the
    # literal same vertex, however many cells/edges converge there.
    canonical_registry: dict[tuple[float, float], Vertex3D] = {}
    for pts in side_points.values():
        for k in range(FlowerLayout.EDGES_PER_SIDE + 1):
            v = pts[k * subdivisions_per_edge]
            canonical_registry[(round(v[0], 3), round(v[1], 3))] = v
    corner_samples = list(canonical_registry.values())

    def interior_base_z(x: float, y: float) -> float:
        """Inverse-distance-weighted blend of the flower's 18 declared
        corner heights. Interior points near a high side of the boundary
        (e.g. hill_peak's cliff side) get pulled up with it and points
        near a low side get pulled down, so the terrain ramps smoothly
        from the boundary's own actual height into the interior. A human
        caught the previous behavior (every interior point anchored to
        one flat flower-wide average height, regardless of how high or
        low its own local boundary actually was) as "a hard cliff but no
        gradient slope stretching the border height into the flower" -
        that flat anchor was correct at keeping the interior from
        wandering to an arbitrary absolute level (decision #5) but wrong
        in throwing away all spatial relationship to the boundary it's
        supposed to blend into."""
        total_w = 0.0
        total_wz = 0.0
        for bx, by, bz in corner_samples:
            w = 1.0 / ((x - bx) ** 2 + (y - by) ** 2 + 1.0)
            total_w += w
            total_wz += w * bz
        return total_wz / total_w

    def interior_z(x: float, y: float) -> float:
        return interior_base_z(x, y) + sample_noise_2d(seed, x, y) * interior_relief_mm

    def canonical_point(x: float, y: float) -> Vertex3D:
        key = (round(x, 3), round(y, 3))
        existing = canonical_registry.get(key)
        if existing is not None:
            return existing
        v = (x, y, interior_z(x, y))
        canonical_registry[key] = v
        return v

    # An earlier round tried lowering a hex-to-hex corner's own shared
    # vertex (used by up to 3 cells' grooves AND the surrounding regular
    # terrain) to avoid each groove climbing back to a fixed peak there.
    # That fix numerically resolved the peak but visually produced a worse
    # artifact: dipping ONE vertex shared by 3 independently-triangulated
    # cells' ordinary (non-groove) terrain created a small sinkhole, each
    # cell interpolating toward it from a different direction. The groove
    # ring is now built as an exact scale-toward-apex transform (see the
    # groove-ring comment below) that reaches full depth at every ring
    # point, corners included, without needing the true corner vertex
    # itself touched - so that hack is gone; the true corner vertex stays
    # exactly as declared, exterior and interior alike.

    # Cache of internal (hex-to-hex) edges' own subdivided point lists,
    # keyed by the edge's two rounded endpoint positions (order-
    # independent), so both cells that border a given edge get the exact
    # same list (just possibly reversed) instead of each independently
    # resolving every point through canonical_point - equivalent, but
    # avoids the registry-lookup cost per shared edge twice over.
    internal_edge_cache: dict[frozenset[tuple[float, float]], list[Vertex3D]] = {}

    # The skirt's sunken shared-edge vertices, keyed by the ORIGINAL point's
    # rounded (x, y). A hex-to-hex edge is walked once by each of the two
    # cells that border it (in opposite directions), and both must end up
    # referencing the literal same lowered vertex - first caller computes
    # and stores it, the other reuses it - rather than two independently
    # computed values that merely ought to be equal. Same discipline, and
    # the same reason, as canonical_registry above.
    sunken_edge_registry: dict[tuple[float, float], Vertex3D] = {}

    def internal_edge_points(p1: Point2D, p2: Point2D) -> list[Vertex3D]:
        key1 = (round(p1[0], 3), round(p1[1], 3))
        key2 = (round(p2[0], 3), round(p2[1], 3))
        cache_key = frozenset((key1, key2))
        cached = internal_edge_cache.get(cache_key)
        if cached is None:
            dx, dy = p2[0] - p1[0], p2[1] - p1[1]
            pts: list[Vertex3D] = [
                canonical_point(p1[0] + dx * (step / subdivisions_per_edge),
                                p1[1] + dy * (step / subdivisions_per_edge))
                for step in range(subdivisions_per_edge + 1)
            ]
            internal_edge_cache[cache_key] = pts
            return pts
        return pts if (round(cached[0][0], 3), round(cached[0][1], 3)) == key1 else list(
            reversed(cached)
        )

    n = subdivisions_per_edge
    all_vertices: list[Vertex3D] = []
    all_triangles: list[Triangle] = []
    side_position_to_global_index: dict[tuple[int, int], int] = {}

    for hex_idx in range(FlowerLayout.HEX_CELL_COUNT):
        exterior_edges = exterior_local_indices.get(hex_idx, ())
        offset = len(all_vertices)

        # edge_pts[edge_idx]: n+1 points from corner_i (index 0) to
        # corner_{i+1} (index n), matching hex_edge_line(hex_idx,
        # edge_idx)'s own (p1, p2) direction - the SAME point objects
        # exterior_fine_points/internal_edge_points already produced
        # (never recomputed), so this cell's own boundary literally
        # shares vertex objects with its neighbors, not just equal values.
        edge_pts: list[list[Vertex3D]] = []
        for edge_idx in range(FlowerLayout.EDGES_PER_HEX):
            if edge_idx in exterior_edges:
                edge_pts.append(exterior_fine_points[hex_idx][edge_idx])
            else:
                p1, p2 = layout.hex_edge_line(hex_idx, edge_idx)
                edge_pts.append(internal_edge_points(p1, p2))

        # Split the cell into 6 wedges (center -> corner_i -> corner_{i+1})
        # and triangulate each as a proper triangular lattice - not a fan
        # from one vertex - so relief detail shows evenly across the
        # whole cell face, not just along its boundary (the visual "too
        # low poly" complaint a human immediately noticed: a fan only
        # adds triangles radiating from a single point, giving a spiky
        # rather than a genuinely subdivided look).
        cx, cy = layout.cell_center(hex_idx)

        # Perpendicular distance from this cell's own center to one of its
        # (ideal, unjittered) edges - a regular hexagon's apothem. Computed
        # once per cell, from layout's raw geometry rather than any
        # particular wedge's own (possibly xy-jittered) far_edge, so every
        # wedge scales its groove ring by the exact same factor - see the
        # groove-ring comment below for why that single shared factor is
        # what keeps ring corners exact instead of needing a taper.
        _ideal_p1, _ideal_p2 = layout.hex_edge_line(hex_idx, 0)
        _adx, _ady = _ideal_p2[0] - _ideal_p1[0], _ideal_p2[1] - _ideal_p1[1]
        _alen = math.hypot(_adx, _ady)
        apothem = (
            abs(_adx * (cy - _ideal_p1[1]) - _ady * (cx - _ideal_p1[0])) / _alen
            if _alen > 1e-9
            else 0.0
        )

        # --- Plateau shaping (decision #11 / height_level, finally wired up).
        # hex_radial(p) is p's distance from this cell's centre measured in
        # units of the cell's own apothem, using the hexagon's own metric
        # (max over the 6 outward edge normals) rather than a circle - so
        # hex_radial == 1 exactly on the cell's edges and ~1.155 at its
        # corners, and a level set of it is a smaller concentric hexagon.
        _edge_normals: list[tuple[float, float]] = []
        for _e in range(FlowerLayout.EDGES_PER_HEX):
            _p1, _p2 = layout.hex_edge_line(hex_idx, _e)
            _ex, _ey = _p2[0] - _p1[0], _p2[1] - _p1[1]
            _el = math.hypot(_ex, _ey)
            if _el <= 1e-9:
                continue
            # Outward normal: flip whichever perpendicular points away from
            # the centre, so this works regardless of edge winding.
            _nx, _ny = -_ey / _el, _ex / _el
            if (_p1[0] - cx) * _nx + (_p1[1] - cy) * _ny < 0.0:
                _nx, _ny = -_nx, -_ny
            _edge_normals.append((_nx, _ny))

        _plateau_r = math.sqrt(plateau_area_frac) if plateau_area_frac > 0.0 else 0.0
        _plateau_z: float | None = None
        if hex_height_levels is not None and hex_idx in hex_height_levels:
            _plateau_z = level_z(hex_height_levels[hex_idx])

        def hex_radial(x: float, y: float) -> float:
            if apothem <= 1e-9 or not _edge_normals:
                return 0.0
            return max(
                (x - cx) * nx + (y - cy) * ny for nx, ny in _edge_normals
            ) / apothem

        def cell_z(x: float, y: float) -> float:
            """This cell's own surface height: dead flat and noise-free at
            the declared height_level across the plateau, then blended out
            to the ordinary interior_z field by the time it reaches the
            cell's edges. Only ever used for points strictly INSIDE one
            cell - points shared with a neighbouring cell (edge points,
            corners) keep going through canonical_point/interior_z, which
            is what lets two cells at different height_levels still agree
            exactly on the vertices they share (the disagreement that made
            an earlier attempt at consuming height_level get reverted)."""
            if _plateau_z is None:
                return interior_z(x, y)
            h = hex_radial(x, y)
            if h <= _plateau_r:
                return _plateau_z
            t = (h - _plateau_r) / max(1.0 - _plateau_r, 1e-9)
            t = min(1.0, max(0.0, t))
            # smoothstep: no crease where the blend leaves the plateau
            t = t * t * (3.0 - 2.0 * t)
            return _plateau_z + (interior_z(x, y) - _plateau_z) * t

        apex = (cx, cy, cell_z(cx, cy))

        # Skirt sizing, computed once per cell. _skirt_scale is the central-
        # scaling factor mapping this cell's true edge onto the inset "rim"
        # line where the skirt has climbed back to full terrain height;
        # _skirt_depth is how far the shared edge itself is sunk below it.
        _effective_profile = groove_profile
        if _effective_profile is None and groove_depth_mm > 0.0:
            _effective_profile = [(1.0, -groove_depth_mm)]
        _has_groove = bool(_effective_profile) and groove_width_mm > 0.0 and apothem > 1e-9
        _skirt_scale = 1.0
        _skirt_depth = 0.0
        if _has_groove:
            _skirt_width_frac, _skirt_z_off = _effective_profile[-1]
            _skirt_scale = 1.0 - (groove_width_mm * _skirt_width_frac) / apothem
            _skirt_depth = abs(_skirt_z_off)

        # spokes[edge_idx]: n+1 points from apex (index 0) to corner_i
        # (index n, = edge_pts[edge_idx][0], the SAME object). Shared
        # between the two wedges that meet at that corner (wedge edge_idx
        # and wedge (edge_idx-1)%6) - computed once per cell, not once
        # per wedge, the same shared-vertex-object discipline as
        # everywhere else in this module.
        #
        # rim_spokes[edge_idx]: the same, but stopping at the skirt's inner
        # rim (apex + (corner-apex) * _skirt_scale) rather than the true
        # corner. The coarse interior grid uses THESE as its boundary, so
        # it hands off to the skirt exactly where the skirt begins; the
        # skirt strip alone owns everything between the rim and the sunken
        # edge. Two adjacent wedges of the same cell derive this from the
        # same apex/scale/corner, so they land on the identical point with
        # no stitching needed. Without a skirt _skirt_scale is 1.0 and this
        # is bit-identical to spokes.
        spokes: list[list[Vertex3D]] = []
        rim_spokes: list[list[Vertex3D]] = []
        for edge_idx in range(FlowerLayout.EDGES_PER_HEX):
            corner = edge_pts[edge_idx][0]
            dx, dy = corner[0] - apex[0], corner[1] - apex[1]
            spoke = [apex]
            for step in range(1, n):
                t = step / n
                x = apex[0] + dx * t
                y = apex[1] + dy * t
                spoke.append((x, y, cell_z(x, y)))
            spoke.append(corner)
            spokes.append(spoke)

            if _has_groove:
                rim_spoke = [apex]
                for step in range(1, n):
                    t = step / n
                    x = apex[0] + dx * _skirt_scale * t
                    y = apex[1] + dy * _skirt_scale * t
                    rim_spoke.append((x, y, cell_z(x, y)))
                rx = apex[0] + dx * _skirt_scale
                ry = apex[1] + dy * _skirt_scale
                # Full, undisturbed terrain height: the rim is where the
                # skirt has finished climbing, so it must agree with the
                # coarse grid it borders (which is interior_z everywhere).
                # Must stay bit-identical to what the strip's own rim row
                # computes at j=0/n - both scale the same `corner` by the
                # same factor and call interior_z on the result.
                rim_spoke.append((rx, ry, cell_z(rx, ry)))
                rim_spokes.append(rim_spoke)
            else:
                rim_spokes.append(spoke)

        vertex_index: dict[Vertex3D, int] = {}

        def local_index(v: Vertex3D) -> int:
            idx = vertex_index.get(v)
            if idx is None:
                idx = len(vertex_index)
                vertex_index[v] = idx
                all_vertices.append(v)
            return idx

        for edge_idx in range(FlowerLayout.EDGES_PER_HEX):
            # Barycentric grid over the wedge (apex, corner_i,
            # corner_{i+1}): I counts steps from apex toward corner_i (the
            # I axis), J counts steps from apex toward corner_{i+1} (the J
            # axis), I + J <= n. This is the standard "subdivide a
            # triangle into n^2 smaller triangles" lattice - exactly the
            # picture asked for, applied to each of the hex's 6 wedges,
            # not a fan from a single vertex.
            far_edge = edge_pts[edge_idx]
            spoke_i = spokes[edge_idx]
            spoke_j = spokes[(edge_idx + 1) % FlowerLayout.EDGES_PER_HEX]
            rim_spoke_i = rim_spokes[edge_idx]
            rim_spoke_j = rim_spokes[(edge_idx + 1) % FlowerLayout.EDGES_PER_HEX]
            grid_cache: dict[tuple[int, int], Vertex3D] = {}

            # far_edge[0] and far_edge[n] are this wedge's OWN raw corner
            # computations (this side's independent noise sample at that
            # position) - for an internal, same-hex wedge corner they're
            # bit-identical to spoke_i[n]/spoke_j[n] (same list slice), but
            # at one of the flower's 6 side-to-side silhouette junctions
            # they are NOT: the coarse grid (grid_point's J==0/I==0
            # branches) always resolves that corner via spoke_i[n]/
            # spoke_j[n], which for an internal-edge-adjacent corner chains
            # through canonical_point's dedup registry to whichever side
            # was processed last - a different, non-bit-identical value
            # from this side's own far_edge[n]. Using far_edge[j] directly
            # for the corner-adjacent strip points reproduced exactly that
            # discrepancy (found via a real hill_peak build: two vertices
            # at the same rounded XY, Z differing by ~0.07mm, a genuine
            # non-watertight hole) - so the strip must reuse the SAME
            # spoke corner objects the coarse grid does, not far_edge's own.
            def far_edge_at(j: int) -> Vertex3D:
                if j == 0:
                    return spoke_i[n]
                if j == n:
                    return spoke_j[n]
                return far_edge[j]

            inset_line = far_edge
            if _has_groove:
                # THE SKIRT. Each hex cell's surface slopes DOWN as it
                # approaches a shared hex-to-hex edge, and the shared edge
                # itself is the low line - so two neighbouring cells form
                # one single V between them.
                #
                # The previous design left the shared edge at full terrain
                # height and sank a channel just inside it, on each side
                # independently. Because BOTH cells did that, the untouched
                # shared edge was left standing as a thin wall between two
                # channels: measured on a real flat_plains build, a
                # cross-section straight across one hex-to-hex edge read
                # 0.0 / -0.6 / 0.0 / -0.6 / 0.0 over about 4mm - a 0.6mm
                # tall knife-edge blade running the full length of every
                # hex boundary. A human reported it repeatedly as a "spike"
                # (and, viewed end-on, that is exactly what it looks like)
                # and asked for a skirt instead: "i dont want a groove that
                # goes up down up again."
                #
                # Two rows per wedge, then, not a channel:
                #   edge row - the true edge's own shared vertices, sunk by
                #              _skirt_depth (see sunken_edge_registry: the
                #              neighbouring cell reuses these exact lowered
                #              vertices, so the two skirts meet in one V
                #              rather than each owning a separate floor).
                #   rim row  - _skirt_scale inward, back at full terrain
                #              height, where the coarse interior grid takes
                #              over (grid_point uses rim_spoke_i/j to match).
                #
                # EVERY hex edge is sunk by the same _skirt_depth, with no
                # taper anywhere - exterior (flower silhouette) edges
                # included. A single hex therefore reads as "/" on all six
                # of its sides (flat across the face, dropping to -depth
                # exactly at the edge), and the V only appears once two
                # hexes are put together, each contributing one half. That
                # is the shape asked for: "i want the V to form only after
                # connecting two hexagons together. so a single hex should
                # have a '/' instead of the 'V'."
                #
                # Sinking the silhouette does NOT break the cross-flower
                # boundary contract (decision #4). The contract requires
                # that two flowers meeting along a side compute the same
                # physical contour, and both derive it from the same
                # deterministic build_side_boundary_vertices output; both
                # then subtract the same constant _skirt_depth from it, so
                # the two contours stay bit-identical. It is also what
                # makes joined flowers behave like joined hexes: their two
                # half-skirts meet to form the same V any interior hex
                # boundary has, instead of the seam being the one hex
                # boundary on the whole piece with no engraved line.
                # A taper back to full height at the 18 declared corners
                # was tried and dropped: it matches between flowers too,
                # but it pinches the V shut at every corner, scalloping
                # each edge instead of leaving one clean constant-depth
                # line.
                def sunken_edge_at(j: int) -> Vertex3D:
                    v = far_edge_at(j)
                    key = (round(v[0], 3), round(v[1], 3))
                    existing = sunken_edge_registry.get(key)
                    if existing is not None:
                        return existing
                    sunk = (v[0], v[1], v[2] - _skirt_depth)
                    sunken_edge_registry[key] = sunk
                    return sunk

                def rim_at(j: int) -> Vertex3D:
                    # Central scaling toward the apex: a constant factor
                    # maps the straight true edge onto an exactly parallel
                    # line a constant distance inside it, endpoints
                    # included - so the two wedges meeting at a corner
                    # scale that shared corner by the same factor and land
                    # on the same rim point with no stitching. At j=0/n
                    # this reproduces rim_spoke_i[n]/rim_spoke_j[n]
                    # bit-for-bit (same expression, same inputs).
                    ex, ey, _ = far_edge_at(j)
                    ix = apex[0] + (ex - apex[0]) * _skirt_scale
                    iy = apex[1] + (ey - apex[1]) * _skirt_scale
                    return (ix, iy, cell_z(ix, iy))

                edge_row = [sunken_edge_at(j) for j in range(n + 1)]
                rim_row = [rim_at(j) for j in range(n + 1)]
                for j in range(n):
                    a = local_index(edge_row[j])
                    b = local_index(edge_row[j + 1])
                    c = local_index(rim_row[j])
                    d = local_index(rim_row[j + 1])
                    for tri in ((a, b, c), (b, d, c)):
                        if len(set(tri)) == 3:
                            all_triangles.append(
                                (tri[0] + offset, tri[1] + offset, tri[2] + offset)
                            )
                inset_line = rim_row

            def grid_point(I: int, J: int) -> Vertex3D:
                # Uses rim_spoke_i/rim_spoke_j, not spoke_i/spoke_j: with a
                # skirt, the coarse interior fill's own outer boundary is
                # the skirt's inner RIM, and the skirt strip alone owns
                # everything from there out to the sunken edge. This is
                # also what keeps the sunken shared edge from dragging the
                # ordinary terrain down with it - the true corner where
                # three cells meet is a skirt vertex only, never a corner
                # of any cell's interior grid, so lowering it cannot pull
                # three independently-triangulated interiors into a
                # sinkhole (an artifact an earlier round hit and a human
                # rejected). Without a skirt, rim_spoke_i/j ARE spoke_i/j.
                cached = grid_cache.get((I, J))
                if cached is not None:
                    return cached
                if J == 0:
                    v = rim_spoke_i[I]
                elif I == 0:
                    v = rim_spoke_j[J]
                elif I + J == n:
                    v = inset_line[J]
                else:
                    x = apex[0] + (rim_spoke_i[n][0] - apex[0]) * (I / n) + (
                        rim_spoke_j[n][0] - apex[0]
                    ) * (J / n)
                    y = apex[1] + (rim_spoke_i[n][1] - apex[1]) * (I / n) + (
                        rim_spoke_j[n][1] - apex[1]
                    ) * (J / n)
                    v = (x, y, cell_z(x, y))
                grid_cache[(I, J)] = v
                return v

            for I in range(n):
                for J in range(n - I):
                    p1 = local_index(grid_point(I, J))
                    p2 = local_index(grid_point(I + 1, J))
                    p3 = local_index(grid_point(I, J + 1))
                    all_triangles.append((p1 + offset, p2 + offset, p3 + offset))
                    if I + J < n - 1:
                        p4 = local_index(grid_point(I + 1, J + 1))
                        all_triangles.append((p2 + offset, p4 + offset, p3 + offset))

        for edge_idx in exterior_edges:
            side_idx, start = exterior_fine_meta[hex_idx][edge_idx]
            for i, v in enumerate(edge_pts[edge_idx][:-1]):
                # The vertex actually IN the mesh at this boundary position
                # is the sunken one once a skirt exists (the silhouette is
                # skirted like any other hex edge, so joined flowers form
                # the same V two joined hexes do). The walls and base plate
                # consume these indices, so they must reference the sunken
                # vertex - not the pre-sink value, which no cell ever
                # emitted. Falls back to v when there is no skirt.
                key = (round(v[0], 3), round(v[1], 3))
                mesh_v = sunken_edge_registry.get(key, v)
                side_position_to_global_index[(side_idx, start + i)] = (
                    offset + vertex_index[mesh_v]
                )

    boundary_indices = [
        side_position_to_global_index[(side_idx, abs_pos)]
        for side_idx in range(FlowerLayout.SIDE_COUNT)
        for abs_pos in range(len(side_points[side_idx]) - 1)
    ]

    return all_vertices, boundary_indices, all_triangles


def _barycentric(
    p: Point2D, a: Point2D, b: Point2D, c: Point2D
) -> tuple[float, float, float] | None:
    (x, y), (ax, ay), (bx, by), (cx, cy) = p, a, b, c
    d = (by - cy) * (ax - cx) + (cx - bx) * (ay - cy)
    if abs(d) < 1e-12:
        return None
    w_a = ((by - cy) * (x - cx) + (cx - bx) * (y - cy)) / d
    w_b = ((cy - ay) * (x - cx) + (ax - cx) * (y - cy)) / d
    return w_a, w_b, 1.0 - w_a - w_b


def sample_height(
    x: float,
    y: float,
    points_2d: Sequence[Point2D],
    vertices_3d: Sequence[Vertex3D],
    triangles: Sequence[tuple[int, int, int]],
    *,
    tol: float = 1e-7,
) -> float:
    """Barycentric-interpolate the Z of an already-triangulated, already
    heightfield-lifted surface at an arbitrary (x, y) inside it. Used to
    give groove/road constraint points a height that follows the real
    surface they're being engraved into, rather than a flat or
    independently-noised value (design decision #10: grooves must
    "follow the actual terrain surface however steep")."""
    for i, j, k in triangles:
        bary = _barycentric((x, y), points_2d[i], points_2d[j], points_2d[k])
        if bary is None:
            continue
        w_a, w_b, w_c = bary
        if w_a >= -tol and w_b >= -tol and w_c >= -tol:
            return (
                w_a * vertices_3d[i][2]
                + w_b * vertices_3d[j][2]
                + w_c * vertices_3d[k][2]
            )
    raise ValueError(f"point ({x}, {y}) is not inside any triangle")
