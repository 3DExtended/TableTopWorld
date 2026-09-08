"""Flower-to-flower magnet bores: 18 blind holes in the base plate walls.

One bore per silhouette edge, at the edge's midpoint (Peter's June 2026
design, measured from printableFiles/hexagonWithRoad.stl). Two
neighbouring flowers share three physical edges, so their bores line up
1:1 - no polarity bookkeeping is needed in the geometry, magnets are
glued in by hand.

A bore is a BLIND recess: a cylinder `depth_mm` deep into the solid plate,
open on the outer wall face. Cutting it into the wall without a CSG
boolean - and without handing a library any choice it could get wrong -
works like this, per edge:

1. The wall is a strip of quads, one per fine boundary segment, between
   the terrain's silhouette rim (top, T_j) and the same points dropped to
   the print bed (bottom, B_j). The fine points within `flat_window_mm`
   of the edge midpoint carry NO XY jitter (terrain.heightfield.
   build_side_boundary_vertices fades it out there), so the columns
   around the midpoint all lie in the edge's own vertical plane.
2. The columns L..R-1 spanning the bore (chosen so their outer points
   clear the circle by `side_margin_mm`) are replaced. A rectangular
   COLLAR [u_L, u_R] x [v_floor, v_cap] is drawn around the mouth, with a
   new row Q_j just above the bed and a new row P_j halfway between the
   circle's top and the lowest rim point above it (the roof check keeps
   that gap open). The ray from the bore centre through every one of the
   mouth's vertices is extended to the collar's boundary and the hit
   point inserted there, so the annulus between mouth and collar splits
   into radial sectors, each fanned from its mouth vertex - every
   triangle has a vertex strictly inside a convex polygon and an edge on
   its boundary, or is a radial sliver, so none can be degenerate or
   inverted no matter how coarse the columns are. (A plain angular merge
   of the two rings was tried first and inverts triangles whenever a
   collar corner has to fan across circle chords it cannot see.)
3. Trapezoids join the collar's rows to the bed (B_j) and to the rim
   (T_j), fanned from a corner so the inserted ray hits are real
   vertices on both sides; the two neighbouring columns L-1 and R get
   the same treatment for the collar's side edges, so nothing is left as
   a T-junction. Those two columns may carry jitter again.
4. The mouth is joined to an inner circle `depth_mm` inside the solid by
   a cylinder wall, closed by a fan disc.

Everything reuses the existing rim vertex indices, so the patch welds to
the untouched quads on either side and to the floor below by shared
vertices, never by coordinate matching. The neighbour's copy of the same
edge fades its jitter over the identical window (a function of distance
from the shared midpoint only), so the boundary contract of decision #4
is untouched.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from terrain.constants import (
    MAGNET_BORE_DEPTH_MM,
    MAGNET_BORE_RADIUS_MM,
    MAGNET_CENTER_Z_MM,
)
from terrain.layout import FlowerLayout

Vertex3D = tuple[float, float, float]
Triangle = tuple[int, int, int]
Point2D = tuple[float, float]


@dataclass(frozen=True)
class MagnetBores:
    radius_mm: float = MAGNET_BORE_RADIUS_MM
    depth_mm: float = MAGNET_BORE_DEPTH_MM
    center_above_bed_mm: float = MAGNET_CENTER_Z_MM
    segments: int = 24
    # Material that must remain between the bore's top and the lowest
    # terrain point above it (a level-0 hex with its skirt and jitter dips
    # ~1.4 mm below the level-0 surface).
    roof_mm: float = 2.0
    # Material kept beside the bore, between the hole and the columns
    # where the jittered wall resumes.
    side_margin_mm: float = 0.75

    @property
    def top_above_bed_mm(self) -> float:
        return self.center_above_bed_mm + self.radius_mm

    @property
    def min_plate_depth_mm(self) -> float:
        return self.top_above_bed_mm + self.roof_mm

    def flat_window_mm(self, edge_len: float, subdivisions_per_edge: int) -> float:
        """Half-length of the jitter-free run around each exterior edge's
        midpoint: the collar's outermost columns start at the first fine
        point beyond radius + side margin, i.e. within one fine segment
        of it, so this covers the whole collar at every resolution."""
        return self.radius_mm + self.side_margin_mm + edge_len / subdivisions_per_edge


DEFAULT_MAGNET_BORES = MagnetBores()


def _orient(tri: Triangle, vertices: Sequence[Vertex3D], want: tuple[float, float, float]) -> Triangle:
    """Return `tri` wound so its normal points along `want` (dot > 0).
    Only for triangles that are non-degenerate by construction."""
    a, b, c = (vertices[i] for i in tri)
    ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    ac = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    nx = ab[1] * ac[2] - ab[2] * ac[1]
    ny = ab[2] * ac[0] - ab[0] * ac[2]
    nz = ab[0] * ac[1] - ab[1] * ac[0]
    dot = nx * want[0] + ny * want[1] + nz * want[2]
    if abs(dot) < 1e-12:
        raise RuntimeError("degenerate triangle in a magnet bore patch")
    return (tri[0], tri[2], tri[1]) if dot < 0 else tri


def _signed_area(a: Point2D, b: Point2D, c: Point2D) -> float:
    return 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))


def cut_magnet_bores(
    vertices: list[Vertex3D],
    top_count: int,
    boundary_indices: Sequence[int],
    layout: FlowerLayout,
    *,
    subdivisions_per_edge: int,
    bottom_z: float,
    bores: MagnetBores,
) -> tuple[list[Vertex3D], list[Triangle], set[int]]:
    """Build the bore patches for all 18 silhouette edges.

    `vertices` is the top surface followed by its flat-bottom copy (bottom
    vertex of top index i is top_count + i), `boundary_indices` the top
    rim in boundary order (side k, edge e, fine point j at position
    (3k + e) * subdivisions_per_edge + j). Returns (new_vertices, faces,
    skipped_positions): faces index into vertices + new_vertices; the
    caller must NOT emit its plain wall quad for any boundary position in
    skipped_positions (the patch replaces them).
    """
    sd = subdivisions_per_edge
    total = len(boundary_indices)
    new_vertices: list[Vertex3D] = []
    faces: list[Triangle] = []
    skipped: set[int] = set()
    z_center = bottom_z + bores.center_above_bed_mm
    r = bores.radius_mm
    m = bores.segments
    two_pi = 2.0 * math.pi

    def add(v: Vertex3D) -> int:
        new_vertices.append(v)
        return len(vertices) + len(new_vertices) - 1

    def resolve(i: int) -> Vertex3D:
        return vertices[i] if i < len(vertices) else new_vertices[i - len(vertices)]

    for side_idx in range(FlowerLayout.SIDE_COUNT):
        corners = layout.side_corners(side_idx)
        for e in range(FlowerLayout.EDGES_PER_SIDE):
            p0, p1 = corners[e], corners[e + 1]
            edge_len = math.hypot(p1[0] - p0[0], p1[1] - p0[1])
            tx, ty = (p1[0] - p0[0]) / edge_len, (p1[1] - p0[1]) / edge_len
            mid = ((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2)
            ox, oy = ty, -tx
            if ox * mid[0] + oy * mid[1] < 0:  # outward = away from the flower centre
                ox, oy = -ox, -oy
            # CCW in the (u = along edge, v = z) frame means normal t x z-hat
            # = (ty, -tx, 0); flip every planar triangle if that points inward.
            plane_flip = (ty * ox - tx * oy) < 0
            outward = (ox, oy, 0.0)

            base = (side_idx * FlowerLayout.EDGES_PER_SIDE + e) * sd

            def u_of(j: int) -> float:
                return (j / sd - 0.5) * edge_len

            def top(j: int) -> int:
                return boundary_indices[(base + j) % total]

            def bot(j: int) -> int:
                return top_count + top(j)

            def uv(i: int) -> Point2D:
                x, y, z = resolve(i)
                return ((x - mid[0]) * tx + (y - mid[1]) * ty, z)

            clear = r + bores.side_margin_mm
            left = [j for j in range(sd + 1) if u_of(j) <= -clear]
            right = [j for j in range(sd + 1) if u_of(j) >= clear]
            L, R = (left[-1] if left else -1), (right[0] if right else sd + 1)
            if L < 1 or R > sd - 1:
                raise ValueError(
                    f"subdivisions_per_edge={sd} is too coarse to cut a magnet bore "
                    f"(need a fine point beyond {clear:.2f} mm from the edge midpoint "
                    "on both sides, with a column to spare)"
                )
            for j in range(L - 1, R + 1):
                skipped.add(base + j)

            lowest_top = min(vertices[top(j)][2] for j in range(L, R + 1))
            if lowest_top < z_center + r + bores.roof_mm * 0.5:
                raise ValueError(
                    f"side {side_idx} edge {e}: terrain surface ({lowest_top:.2f}) leaves "
                    f"less than {bores.roof_mm * 0.5:.1f} mm of roof above the magnet bore "
                    f"(top at {z_center + r:.2f})"
                )
            v_cap = 0.5 * ((z_center + r) + lowest_top)
            v_floor = 0.5 * (bottom_z + (z_center - r))
            u_L, u_R = u_of(L), u_of(R)

            def emit_planar(tri: Triangle) -> None:
                area = _signed_area(uv(tri[0]), uv(tri[1]), uv(tri[2]))
                if abs(area) < 1e-9:
                    raise RuntimeError("degenerate triangle in a magnet bore collar")
                if (area < 0) != plane_flip:
                    tri = (tri[0], tri[2], tri[1])
                faces.append(tri)

            def on_plane(u: float, v: float) -> Vertex3D:
                return (mid[0] + tx * u, mid[1] + ty * u, v)

            # New vertices: the collar's two rows (vertically above B_j),
            # the bore mouth, the inner circle and the blind end's centre.
            Q, P = {}, {}
            for j in range(L, R + 1):
                x, y, _ = vertices[bot(j)]
                Q[j] = add((x, y, v_floor))
                P[j] = add((x, y, v_cap))
            angles = [two_pi * k / m for k in range(m)]
            hole_uv = [(r * math.cos(a), z_center + r * math.sin(a)) for a in angles]
            mouth = [add(on_plane(u, v)) for u, v in hole_uv]
            inner = [
                add((mid[0] + tx * u - ox * bores.depth_mm,
                     mid[1] + ty * u - oy * bores.depth_mm, v))
                for u, v in hole_uv
            ]
            center_inner = add((mid[0] - ox * bores.depth_mm, mid[1] - oy * bores.depth_mm, z_center))

            # Collar boundary, one sorted (key, vertex) list per edge; the
            # ray through each mouth vertex is extended to the boundary and
            # its hit inserted (or an existing vertex reused if it lands on one).
            edges: dict[str, list[tuple[float, int]]] = {
                "bottom": [(u_of(j), Q[j]) for j in range(L, R + 1)],
                "top": [(u_of(j), P[j]) for j in range(L, R + 1)],
                "left": [(v_floor, Q[L]), (v_cap, P[L])],
                "right": [(v_floor, Q[R]), (v_cap, P[R])],
            }
            hits: list[int] = []
            for a in angles:
                c, sn = math.cos(a), math.sin(a)
                cands = []
                if c > 1e-12:
                    cands.append((u_R / c, "right"))
                if c < -1e-12:
                    cands.append((u_L / c, "left"))
                if sn > 1e-12:
                    cands.append(((v_cap - z_center) / sn, "top"))
                if sn < -1e-12:
                    cands.append(((v_floor - z_center) / sn, "bottom"))
                t, edge = min(cands)
                hu, hv = t * c, z_center + t * sn
                key = hu if edge in ("top", "bottom") else hv
                lst = edges[edge]
                found = next((idx for kk, idx in lst if abs(kk - key) < 1e-7), None)
                if found is None:
                    found = add(on_plane(hu, hv))
                    lst.append((key, found))
                    lst.sort()
                hits.append(found)

            bottom_row = [idx for _, idx in edges["bottom"]]      # u ascending
            top_row = [idx for _, idx in edges["top"]]            # u ascending
            left_col = [idx for _, idx in edges["left"]]          # v ascending
            right_col = [idx for _, idx in edges["right"]]        # v ascending
            ring = bottom_row + right_col[1:] + top_row[::-1][1:] + left_col[::-1][1:-1]
            pos = {idx: t for t, idx in enumerate(ring)}
            n_ring = len(ring)

            # Annulus: one radial sector per mouth segment, fanned from C_k.
            for k in range(m):
                k1 = (k + 1) % m
                t = pos[hits[k]]
                t_end = pos[hits[k1]]
                while t != t_end:
                    emit_planar((mouth[k], ring[t], ring[(t + 1) % n_ring]))
                    t = (t + 1) % n_ring
                emit_planar((mouth[k], hits[k1], mouth[k1]))

            # Strips between the collar's rows and the bed / the rim, per
            # column, fanned from a corner off the row so inserted hits are
            # real vertices of these triangles too.
            def between(row: list[int], first: int, last: int) -> list[int]:
                return row[row.index(first) : row.index(last) + 1]

            for j in range(L, R):
                top_poly = between(top_row, P[j], P[j + 1]) + [top(j + 1), top(j)]
                for a_idx, b_idx in zip(top_poly[:-1], top_poly[1:]):
                    if top(j) in (a_idx, b_idx):
                        continue
                    emit_planar((top(j), a_idx, b_idx))
                bot_poly = [bot(j), bot(j + 1)] + between(bottom_row, Q[j], Q[j + 1])[::-1]
                for a_idx, b_idx in zip(bot_poly[1:-1], bot_poly[2:]):
                    emit_planar((bot(j), a_idx, b_idx))

            # Neighbouring columns (possibly jittered, hence 3D orientation),
            # fanned from their outer bottom corner.
            all_v = vertices + new_vertices
            left_poly = [bot(L), *left_col, top(L), top(L - 1)]
            for a_idx, b_idx in zip(left_poly[:-1], left_poly[1:]):
                faces.append(_orient((bot(L - 1), a_idx, b_idx), all_v, outward))
            right_poly = [top(R + 1), top(R), *right_col[::-1], bot(R)]
            for a_idx, b_idx in zip(right_poly[:-1], right_poly[1:]):
                faces.append(_orient((bot(R + 1), a_idx, b_idx), all_v, outward))

            # Cylinder wall (normal at the bore axis) and blind end (normal
            # back out of the opening).
            for k in range(m):
                k1 = (k + 1) % m
                o0, o1 = mouth[k], mouth[k1]
                i0, i1 = inner[k], inner[k1]
                u, v = hole_uv[k]
                want = (-(tx * u), -(ty * u), -(v - z_center))
                faces.append(_orient((o0, o1, i1), all_v, want))
                faces.append(_orient((o0, i1, i0), all_v, want))
                faces.append(_orient((center_inner, i1, i0), all_v, outward))

    return new_vertices, faces, skipped
