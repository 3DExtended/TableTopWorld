"""The terrain height field: flat pads, S-curve bands, roads and rivers.

This is the interior model Peter approved from docs/mockups on 2026-09-08
(docs/mockups/scripts/mock_common.py is the idealised original). Every hex
cell has a declared height_level; the field between them is:

* Level blend. Inside hex i the level is l_i. Towards each of its 6 edges
  a smootherstep band (from `organic_r` of the apothem out to the edge)
  blends in the neighbour's level, so a shared edge sits at the mean of
  the two hexes and a corner at the mean of the three that meet there -
  the same value whichever hex evaluates it, which is what keeps the
  independently triangulated cells continuous. A low-frequency wobble
  of the band's start makes the step fronts sinuous rather than
  concentric hexagons.
* Plateau (standable) hexes mask a dead-flat pad on top of that blend:
  their own level inside `pad_r` of the apothem (55% of the area), a
  smootherstep to the blend over the last 6 mm to the edge. Organic
  hexes get the broad blend as is, plus a rolling relief.
* The silhouette is sacred (decision #4): along a flower's 18 exterior
  edges the field is driven to the declared contour - the lerp between
  the side's corner heights - by an override that reaches exactly 1 on
  the edge line and, on the line, listens to that edge alone, so every
  interior point next to the boundary meets the same contour both
  neighbouring flowers build.
* Noise: +/-noise_mm fine texture in the bands and over organic hexes,
  never on a pad.
* Roads and rivers are polylines with rounded bends from one side's
  middle-edge midpoint (layout.junction_center, shared by exactly two
  flowers) through hex centres to another. A road is a smooth bed 1 mm below the (smoothed)
  terrain with near-vertical edges, allowed to stand proud as an
  embankment; a river is a flat-bottomed channel whose smoothstep banks
  blend from the local terrain down to the floor and only ever cut
  down. Both cross the silhouette perpendicular to
  it, and their bed is pinned to the declared contour height at the
  crossing, so the neighbour's copy of the same road meets it exactly;
  boundary_z applies the identical rule to the silhouette points.

Nothing here is bit-exact across platforms by contract - only the
silhouette contour (terrain/boundary_noise.py) needs that - but it is a
pure function of the declared data, the seed and the position.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Collection, Sequence

import numpy as np

from terrain.boundary_noise import _smootherstep, sample_noise_2d
from terrain.layout import FlowerLayout

Point2D = tuple[float, float]
PathSpec = tuple[int, int, tuple[int, ...]]  # (entry_side, exit_side, via hex indices)


@dataclass(frozen=True)
class FieldParams:
    # Plateau hexes: flat inside pad_r of the apothem (55% of the area),
    # a smootherstep band from there to the edge.
    pad_r: float = math.sqrt(0.55)
    # Organic hexes: the band starts almost at the centre - a broad roll.
    organic_r: float = 0.25
    # Sinuous fronts: the band start moves by up to this fraction of the
    # apothem, driven by a slow noise.
    wobble: float = 0.12
    # Within this fraction of the apothem of a silhouette edge, the other
    # silhouette edges' say in the contour override fades to nothing, so
    # ON an edge the field is exactly that edge's declared contour.
    edge_silence_r: float = 0.15
    wobble_cell_mm: float = 30.0
    # Fine texture in the bands / on organic hexes (never on a pad).
    noise_mm: float = 0.5
    noise_cell_mm: float = 8.0
    # Rolling relief of organic hexes, fading out at their edges.
    organic_relief_mm: float = 2.0
    organic_relief_cell_mm: float = 20.0
    # Roads: 16 mm wide, bed 1 mm below the smoothed terrain, 0.6 mm edge.
    road_half_width_mm: float = 8.0
    road_depth_mm: float = 1.0
    road_edge_mm: float = 0.6
    # Rivers: 22 mm wide, 2.5 mm deep flat bed, 5 mm smoothstep banks.
    river_half_width_mm: float = 11.0
    river_depth_mm: float = 2.5
    river_bank_mm: float = 5.0
    # A river bed never goes below this (level 0 is z = 0): keeps
    # >= 1.05 mm of roof over the wall bores, whose top is at -3.45.
    river_bed_min_z: float = -2.4
    # Roads/rivers are straight legs between hex centres joined by circular
    # fillets; the radius must exceed every half-width so the nearest point
    # on the centreline (hence the bed height) never jumps inside the bed.
    bend_radius_mm: float = 15.0
    # Bed smoothing window along the path, and the distance over which
    # the bed is pinned to the contour height at a silhouette crossing.
    bed_smooth_mm: float = 12.0
    crossing_pin_mm: float = 12.0
    # Spline control point pushed outward past each crossing, so the path
    # crosses the silhouette perpendicular to it (both flowers agree).
    path_outward_mm: float = 20.0

    @property
    def pad_min_r(self) -> float:
        """Smallest band start a plateau can have once the wobble has
        pushed it inward (the wobble is scaled by 1 - n^2 at n = pad_r)."""
        return self.pad_r - self.wobble * (1.0 - self.pad_r**2)


DEFAULT_FIELD_PARAMS = FieldParams()


def fillet_polyline(
    points: Sequence[Point2D], radius: float, step_mm: float = 0.5, min_radius: float = 0.0
) -> np.ndarray:
    """Straight legs through `points` with every bend rounded by a circular
    arc of `radius` (reduced where a leg is too short for it). A bend's
    inner medial axis then lies `radius` away from the centreline, so any
    band narrower than that has a continuous nearest point everywhere -
    which a Catmull-Rom spline through the same points did not give: its
    60-degree bends at the ring hex centres were tighter than a river's
    half-width, and the bed height jumped by millimetres on the inside of
    every bend."""
    pts = [np.asarray(q, dtype=float) for q in points]
    out: list[np.ndarray] = [pts[0]]
    for i in range(1, len(pts) - 1):
        prev, cur, nxt = pts[i - 1], pts[i], pts[i + 1]
        a, b = cur - prev, nxt - cur
        la, lb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
        if la < 1e-9 or lb < 1e-9:
            raise ValueError("path control points must not repeat")
        a, b = a / la, b / lb
        turn = math.acos(min(max(float(a @ b), -1.0), 1.0))
        if turn < 1e-6:
            continue  # straight through: nothing to round
        r = radius
        t = r * math.tan(turn / 2.0)
        cap = 0.5 * min(la, lb)
        if t > cap:
            t, r = cap, cap / math.tan(turn / 2.0)
        if r < min_radius:
            raise ValueError(
                f"a {math.degrees(turn):.0f}-degree bend at {tuple(np.round(cur, 2))} can only be "
                f"rounded to {r:.1f} mm, below the {min_radius:.1f} mm the band's width needs"
            )
        p_in, p_out = cur - a * t, cur + b * t
        bis = b - a
        bis /= np.linalg.norm(bis)
        centre = cur + bis * (r / math.cos(turn / 2.0))
        ang0 = math.atan2(*(p_in - centre)[::-1])
        sweep = math.copysign(turn, float(a[0] * b[1] - a[1] * b[0]))
        n = max(2, int(math.ceil(r * abs(sweep) / step_mm)))
        for k in range(n + 1):
            ang = ang0 + sweep * k / n
            out.append(centre + r * np.array([math.cos(ang), math.sin(ang)]))
        assert np.linalg.norm(out[-1] - p_out) < 1e-6
    out.append(pts[-1])
    return np.array(out)


def _resample(poly: np.ndarray, step_mm: float) -> np.ndarray:
    """Re-sample a polyline at (close to) uniform arc-length spacing, so a
    moving average over it is a moving average in millimetres."""
    seg = np.linalg.norm(np.diff(poly, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    n = max(2, int(math.ceil(s[-1] / step_mm)) + 1)
    t = np.linspace(0.0, s[-1], n)
    return np.column_stack([np.interp(t, s, poly[:, 0]), np.interp(t, s, poly[:, 1])])


def _insert_point(poly: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Insert `p`, which lies on the polyline, as a vertex (no-op if it
    already is one)."""
    a, b = poly[:-1], poly[1:]
    ab = b - a
    t = np.clip(((p - a) * ab).sum(axis=1) / np.maximum((ab * ab).sum(axis=1), 1e-12), 0.0, 1.0)
    d = np.linalg.norm(p - (a + t[:, None] * ab), axis=1)
    i = int(d.argmin())
    if d[i] > 1e-6:
        raise ValueError(f"point {p} is {d[i]:.3g} mm off the path")
    if np.linalg.norm(poly[i] - p) < 1e-9 or np.linalg.norm(poly[i + 1] - p) < 1e-9:
        return poly
    return np.vstack([poly[: i + 1], p[None, :], poly[i + 1 :]])


@dataclass
class Crossing:
    """Where a path crosses the silhouette: the middle-edge midpoint of a
    side, the unit direction along that edge, the outward unit normal,
    the declared contour height there, and the arc length along the path."""

    point: Point2D
    edge_dir: Point2D
    outward: Point2D
    z: float
    s: float = 0.0


class Path:
    """A road or river centreline: dense polyline, arc length, bed height."""

    def __init__(
        self,
        kind: str,
        control_points: Sequence[Point2D],
        crossings: list[Crossing],
        bend_radius_mm: float,
        min_radius_mm: float = 0.0,
    ):
        self.kind = kind
        poly = _resample(fillet_polyline(control_points, bend_radius_mm, min_radius=min_radius_mm), 0.5)
        # Make each crossing an exact sample: its bed is then pinned to the
        # declared height exactly (not interpolated between two samples a
        # hair away from it), which is what the neighbour computes too.
        for c in crossings:
            poly = _insert_point(poly, np.asarray(c.point, dtype=float))
        self.poly = poly
        seg = np.diff(self.poly, axis=0)
        seg_len = np.linalg.norm(seg, axis=1)
        self.s = np.concatenate([[0.0], np.cumsum(seg_len)])
        self._a = self.poly[:-1]
        self._ab = seg
        self._len2 = np.maximum((seg * seg).sum(axis=1), 1e-12)
        self._seg_len = seg_len
        self.bed_z: np.ndarray | None = None
        self.crossings = crossings
        for c in self.crossings:
            _, c.s = self.dist_s(*c.point)

    def dist_s(self, x: float, y: float) -> tuple[float, float]:
        """(distance to the path, arc length of the closest point)."""
        p = np.array([x, y])
        ap = p - self._a
        t = np.clip((ap * self._ab).sum(axis=1) / self._len2, 0.0, 1.0)
        closest = self._a + t[:, None] * self._ab
        d = np.linalg.norm(p - closest, axis=1)
        i = int(d.argmin())
        return float(d[i]), float(self.s[i] + t[i] * self._seg_len[i])

    def set_bed(self, raw_z: np.ndarray, smooth_mm: float, pin_mm: float) -> None:
        """Smooth the terrain height along the path and pin it to each
        crossing's contour height. Samples beyond a crossing (outside the
        flower) are replaced by that crossing's height first, so nothing
        the neighbour would compute differently leaks into the average."""
        z = np.array(raw_z, dtype=float)
        for c in self.crossings:
            if c.s <= self.s[len(self.s) // 2]:
                z[self.s < c.s] = c.z
            else:
                z[self.s > c.s] = c.z
        ds = float(np.mean(np.diff(self.s))) if len(self.s) > 1 else 1.0
        k = max(1, int(round(smooth_mm / 2.0 / ds)))
        kernel = np.ones(2 * k + 1) / (2 * k + 1)
        smooth = np.convolve(np.pad(z, k, mode="edge"), kernel, mode="valid")
        for c in self.crossings:
            u = np.clip(1.0 - np.abs(self.s - c.s) / pin_mm, 0.0, 1.0)
            u = u * u * u * (u * (u * 6 - 15) + 10)
            smooth = smooth * (1.0 - u) + c.z * u
        self.bed_z = smooth

    def bed(self, s: float) -> float:
        assert self.bed_z is not None, "set_bed() first"
        return float(np.interp(s, self.s, self.bed_z))


@dataclass(frozen=True)
class _EdgeInfo:
    nx: float
    ny: float
    neighbour: int | None
    # exterior edges only: endpoints and their declared contour heights
    p0: Point2D | None
    p1: Point2D | None
    z0: float
    z1: float


class TerrainField:
    def __init__(
        self,
        layout: FlowerLayout,
        hex_levels: dict[int, int],
        plateau_hexes: Collection[int],
        side_corner_heights: dict[int, tuple[int, int, int, int]],
        level_z: Callable[[int], float],
        seed: int,
        *,
        params: FieldParams = DEFAULT_FIELD_PARAMS,
        roads: Sequence[PathSpec] = (),
        rivers: Sequence[PathSpec] = (),
    ) -> None:
        self.layout = layout
        self.params = params
        self.seed = seed
        self.level_z = level_z
        self.hex_levels = dict(hex_levels)
        self.plateau_hexes = set(plateau_hexes)
        self.centres = [layout.cell_center(h) for h in range(FlowerLayout.HEX_CELL_COUNT)]
        self.apothem = layout.hex_outer_width * math.sqrt(3.0) / 2.0
        self._level_mm = {h: level_z(lvl) for h, lvl in self.hex_levels.items()}

        # declared contour heights per exterior edge, keyed by its endpoints
        exterior: dict[frozenset, tuple[Point2D, Point2D, float, float]] = {}
        for side_idx in range(FlowerLayout.SIDE_COUNT):
            corners = layout.side_corners(side_idx)
            heights = side_corner_heights[side_idx]
            for e in range(FlowerLayout.EDGES_PER_SIDE):
                key = frozenset((_rkey(corners[e]), _rkey(corners[e + 1])))
                exterior[key] = (corners[e], corners[e + 1], level_z(heights[e]), level_z(heights[e + 1]))

        self._edges: list[list[_EdgeInfo]] = []
        for h in range(FlowerLayout.HEX_CELL_COUNT):
            cx, cy = self.centres[h]
            infos = []
            for k in range(FlowerLayout.EDGES_PER_HEX):
                p1, p2 = layout.hex_edge_line(h, k)
                mx, my = (p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2
                nx, ny = mx - cx, my - cy
                norm = math.hypot(nx, ny)
                nx, ny = nx / norm, ny / norm
                tx, ty = cx + 2 * self.apothem * nx, cy + 2 * self.apothem * ny
                neighbour = None
                for j, (jx, jy) in enumerate(self.centres):
                    if math.hypot(jx - tx, jy - ty) < 1e-6:
                        neighbour = j
                        break
                if neighbour is not None:
                    infos.append(_EdgeInfo(nx, ny, neighbour, None, None, 0.0, 0.0))
                    continue
                ext = exterior.get(frozenset((_rkey(p1), _rkey(p2))))
                if ext is None:
                    raise RuntimeError(f"hex {h} edge {k}: neither internal nor a declared exterior edge")
                q0, q1, z0, z1 = ext
                infos.append(_EdgeInfo(nx, ny, None, q0, q1, z0, z1))
            self._edges.append(infos)

        self._side_mid_z = {}
        for side_idx in range(FlowerLayout.SIDE_COUNT):
            h1, h2 = side_corner_heights[side_idx][1], side_corner_heights[side_idx][2]
            self._side_mid_z[side_idx] = 0.5 * (level_z(h1) + level_z(h2))

        self.roads = [self._make_path("road", spec) for spec in roads]
        self.rivers = [self._make_path("river", spec) for spec in rivers]
        for path in self.roads + self.rivers:
            raw = np.array([self.terrain_z_at(x, y) for x, y in path.poly])
            path.set_bed(raw, params.bed_smooth_mm, params.crossing_pin_mm)

    # ---------------------------------------------------------------- paths
    def _crossing(self, side_idx: int) -> Crossing:
        """A path crosses side k through the midpoint of its MIDDLE edge
        (layout.junction_center), perpendicular to that edge. The outward
        direction is the edge's own normal - NOT the neighbour-flower
        offset, which is ~19 degrees off it: leaving along the offset put
        a kink just inside every crossing, and next to a kink the nearest
        point on the path (so the bed height) jumps - an 11 mm step in
        the road bed of river_bend."""
        corners = self.layout.side_corners(side_idx)
        (ax, ay), (bx, by) = corners[1], corners[2]
        q = ((ax + bx) / 2.0, (ay + by) / 2.0)
        jx, jy = self.layout.junction_center(side_idx)
        assert math.hypot(q[0] - jx, q[1] - jy) < 1e-6
        ex, ey = bx - ax, by - ay
        norm = math.hypot(ex, ey)
        ex, ey = ex / norm, ey / norm
        nx, ny = -ey, ex
        # side k's middle edge belongs to ring hex k+1: outward points away from its centre
        cx, cy = self.centres[side_idx + 1]
        if (q[0] - cx) * nx + (q[1] - cy) * ny < 0.0:
            nx, ny = -nx, -ny
        assert abs(abs((q[0] - cx) * nx + (q[1] - cy) * ny) - self.apothem) < 1e-6
        return Crossing(point=q, edge_dir=(ex, ey), outward=(nx, ny), z=self._side_mid_z[side_idx])

    def _make_path(self, kind: str, spec: PathSpec) -> Path:
        entry, exit_, via = spec
        if entry == exit_:
            raise ValueError("a road/river must enter and leave on different sides")
        c_in, c_out = self._crossing(entry), self._crossing(exit_)
        ring_in = entry + 1  # side k's middle edge belongs to ring hex k+1
        ring_out = exit_ + 1
        route = FlowerLayout.path_route(entry, exit_, via)
        # The crossing sits on the ring hex's own edge, so crossing -> ring
        # centre is a straight leg along the edge normal; with the fillet
        # radius below an apothem the path stays exactly perpendicular to
        # the silhouette for more than any road/river half-width. That is
        # what lets boundary_z() apply the crossing rule from declared data
        # alone: a silhouette point's distance to the path is its offset
        # along the edge. The outward extension keeps that true for the
        # last silhouette points too.
        o = self.params.path_outward_mm
        pts: list[Point2D] = [
            (c_in.point[0] + c_in.outward[0] * o, c_in.point[1] + c_in.outward[1] * o),
            c_in.point,
            *[self.centres[h] for h in route],
            c_out.point,
            (c_out.point[0] + c_out.outward[0] * o, c_out.point[1] + c_out.outward[1] * o),
        ]
        r = self.params.bend_radius_mm
        reach = max(self.params.road_half_width_mm + self.params.road_edge_mm, self.params.river_half_width_mm)
        if r <= reach:
            raise ValueError(f"bend_radius_mm {r} must exceed the widest half-width {reach}")
        # 60-degree bends at the ring centres: the arc's tangent length must leave the crossing leg straight
        assert r * math.tan(math.pi / 6) <= self.apothem - reach
        half = (
            self.params.road_half_width_mm + self.params.road_edge_mm
            if kind == "road"
            else self.params.river_half_width_mm
        )
        return Path(kind, pts, [c_in, c_out], r, min_radius_mm=half + 0.5)

    def socket_allowed(self, hex_idx: int, clearance_mm: float) -> bool:
        """A top socket needs its whole lattice ring on the flat pad: the
        hex centre must be `clearance_mm` clear of every road/river."""
        if hex_idx not in self.plateau_hexes:
            return False
        cx, cy = self.centres[hex_idx]
        for path in self.roads + self.rivers:
            reach = (
                self.params.road_half_width_mm + self.params.road_edge_mm
                if path.kind == "road"
                else self.params.river_half_width_mm
            )
            d, _ = path.dist_s(cx, cy)
            if d < reach + clearance_mm:
                return False
        return True

    # ---------------------------------------------------------------- field
    def plateau_z(self, hex_idx: int) -> float | None:
        return self._level_mm[hex_idx] if hex_idx in self.plateau_hexes else None

    def hex_of(self, x: float, y: float) -> int:
        best, best_d = 0, float("inf")
        for h, (cx, cy) in enumerate(self.centres):
            d = (x - cx) ** 2 + (y - cy) ** 2
            if d < best_d:
                best, best_d = h, d
        return best

    def _noise(self, x: float, y: float) -> float:
        p = self.params
        a = sample_noise_2d(self.seed, x, y, lattice_step=p.noise_cell_mm)
        b = sample_noise_2d(self.seed + 7, x, y, lattice_step=p.noise_cell_mm / 2.0)
        return (a + 0.5 * b) / 1.5

    def level_mm(self, hex_idx: int, x: float, y: float) -> tuple[float, float]:
        """(blended level height in mm, band factor 0..1).

        The band factor is 0 where the hex is purely at its own level and
        exactly 1 on every one of its edges; it scales the fine noise on a
        plateau and fades the rolling relief of an organic hex.

        Every hex, plateau or not, blends its neighbours in with the SAME
        broad `organic_r` band, so on a shared edge both hexes evaluate
        the same numbers (own level + neighbour at weight 1 + the corner
        hexes at weights that only depend on the point's position). A
        plateau then masks its flat pad on top: inside `pad_r` the level
        is its own, and the mask fades to 0 over the band to the edge, so
        it never disturbs what the neighbour computes there. Giving the
        plateau its own narrower band instead (an earlier version) made
        the two hexes of a plateau/organic pair disagree by several mm on
        their shared edge near each corner - a sawtooth along the seam.
        """
        p = self.params
        cx, cy = self.centres[hex_idx]
        A = self.apothem
        r0 = p.organic_r
        wob = p.wobble * sample_noise_2d(self.seed + 11, x, y, lattice_step=p.wobble_cell_mm)
        own = self._level_mm[hex_idx]
        num = own
        den = 1.0
        wmax = 0.0
        hex_dist = 0.0
        ext_w: list[float] = []
        ext_t: list[float] = []
        ext_n: list[float] = []
        for e in self._edges[hex_idx]:
            n = ((x - cx) * e.nx + (y - cy) * e.ny) / A
            nc = min(max(n, 0.0), 1.0)
            n += wob * (1.0 - nc * nc)
            hex_dist = max(hex_dist, n)
            t = (n - r0) / (1.0 - r0)
            t = min(max(t, 0.0), 1.0)
            w = _smootherstep(t)
            if w > wmax:
                wmax = w
            if e.neighbour is not None:
                num += w * self._level_mm[e.neighbour]
                den += w
            elif w > 0.0:
                assert e.p0 is not None and e.p1 is not None
                ex, ey = e.p1[0] - e.p0[0], e.p1[1] - e.p0[1]
                u = ((x - e.p0[0]) * ex + (y - e.p0[1]) * ey) / (ex * ex + ey * ey)
                u = min(max(u, 0.0), 1.0)
                ext_w.append(w)
                ext_t.append(e.z0 + (e.z1 - e.z0) * u)
                ext_n.append(min(n, 1.0))
        level = num / den
        if ext_w:
            keep = 1.0
            for w in ext_w:
                keep *= 1.0 - w
            W = 1.0 - keep
            # The contour to drive towards: each exterior edge's lerp at the
            # nearest point of its segment, weighted by its band - but an
            # edge's say is silenced next to ANOTHER exterior edge of this
            # hex, so on an edge the field is exactly that edge's contour and
            # at a corner exactly the corner height. Without that, the
            # adjacent edge's corner-clamped height pulled the field up to a
            # millimetre off the contour along the edge, on whichever of the
            # two flowers sharing it owns the adjacent edge with the same ring
            # hex: a crease along the seam.
            silence = []
            for i in range(len(ext_w)):
                factor = 1.0
                for j, nj in enumerate(ext_n):
                    if j != i:
                        factor *= _smootherstep(min(max((1.0 - nj) / p.edge_silence_r, 0.0), 1.0))
                silence.append(factor)
            weight_sum = sum(w * f for w, f in zip(ext_w, silence))
            if weight_sum > 1e-12:
                T = sum(w * f * t for w, f, t in zip(ext_w, silence, ext_t)) / weight_sum
            else:  # at a corner: every edge meeting there says the corner height
                at_corner = [t for t, nn in zip(ext_t, ext_n) if nn > 1.0 - 1e-9]
                T = sum(at_corner) / len(at_corner) if at_corner else ext_t[max(range(len(ext_w)), key=ext_w.__getitem__)]
            level = (1.0 - W) * level + W * T
        if hex_idx not in self.plateau_hexes:
            return level, wmax
        t = (hex_dist - p.pad_r) / (1.0 - p.pad_r)
        band = _smootherstep(min(max(t, 0.0), 1.0))
        return (1.0 - band) * own + band * level, band

    def terrain_z(self, hex_idx: int, x: float, y: float) -> float:
        """Height without roads/rivers, as seen from hex `hex_idx`."""
        p = self.params
        level, band = self.level_mm(hex_idx, x, y)
        noise = p.noise_mm * self._noise(x, y)
        if hex_idx in self.plateau_hexes:
            return level + band * noise
        relief = p.organic_relief_mm * sample_noise_2d(
            self.seed + 23, x, y, lattice_step=p.organic_relief_cell_mm
        )
        return level + noise + (1.0 - band) * relief

    def terrain_z_at(self, x: float, y: float) -> float:
        return self.terrain_z(self.hex_of(x, y), x, y)

    def _apply_paths(self, x: float, y: float, z: float) -> float:
        p = self.params
        for road in self.roads:
            d, s = road.dist_s(x, y)
            if d <= p.road_half_width_mm + p.road_edge_mm:
                bed = road.bed(s) - p.road_depth_mm
                if d <= p.road_half_width_mm:
                    z = bed
                else:
                    t = (d - p.road_half_width_mm) / p.road_edge_mm
                    z = bed * (1.0 - t) + z * t
        for river in self.rivers:
            d, s = river.dist_s(x, y)
            if d <= p.river_half_width_mm:
                z = self._river_profile(z, d, river.bed(s))
        return z

    def _river_profile(self, z: float, d: float, bed: float) -> float:
        """Carve a river into terrain height `z` at distance `d` from the
        centreline whose (smoothed) bed height is `bed`: a flat floor
        river_depth_mm below the bed, banks that blend from the LOCAL
        terrain down to it over river_bank_mm, never raising anything.
        Blending from the local terrain (not from the bed) is what keeps
        the channel edge continuous where the ground slopes across the
        river - measured from the bed it was a 5 mm wall on the uphill
        side. The floor is clamped so a level-0 river keeps its roof over
        the wall bores."""
        p = self.params
        bank = min(max((p.river_half_width_mm - d) / p.river_bank_mm, 0.0), 1.0)
        bank = bank * bank * (3.0 - 2.0 * bank)
        floor = max(bed - p.river_depth_mm, p.river_bed_min_z)
        return min(z, (1.0 - bank) * z + bank * floor)

    def z(self, hex_idx: int, x: float, y: float) -> float:
        """Final height as seen from hex `hex_idx` (roads/rivers applied)."""
        return self._apply_paths(x, y, self.terrain_z(hex_idx, x, y))

    def z_at(self, x: float, y: float) -> float:
        """Final height, hex chosen by nearest centre. On a shared edge or
        corner every hex computes the same value by construction."""
        return self.z(self.hex_of(x, y), x, y)

    def boundary_z(self, x: float, y: float, z: float) -> float:
        """Roads/rivers applied to a silhouette point (x, y, z): the same
        crossing rule the neighbouring flower applies, expressed purely in
        declared data - the crossing point, the edge direction, the
        contour height at the crossing and the widths."""
        p = self.params
        for path in self.roads + self.rivers:
            for c in path.crossings:
                dx, dy = x - c.point[0], y - c.point[1]
                if abs(dx * c.outward[0] + dy * c.outward[1]) > 3.0:
                    continue  # not this crossing's edge
                u = abs(dx * c.edge_dir[0] + dy * c.edge_dir[1])
                if path.kind == "road":
                    if u <= p.road_half_width_mm + p.road_edge_mm:
                        bed = c.z - p.road_depth_mm
                        if u <= p.road_half_width_mm:
                            z = bed
                        else:
                            t = (u - p.road_half_width_mm) / p.road_edge_mm
                            z = bed * (1.0 - t) + z * t
                elif u <= p.river_half_width_mm:
                    z = self._river_profile(z, u, c.z)
        return z

    def river_bank_factor(self, x: float, y: float) -> float:
        """0 outside every river, rising over the bank to 1 on the flat
        bed. Used to fade the hex-line skirt out under water (a V-line in
        a river bed is invisible and would eat the wall bores' roof at
        level 0). Near a silhouette crossing the path runs straight and
        perpendicular through it, so both flowers get the same value."""
        p = self.params
        best = 0.0
        for river in self.rivers:
            d, _ = river.dist_s(x, y)
            bank = min(max((p.river_half_width_mm - d) / p.river_bank_mm, 0.0), 1.0)
            bank = bank * bank * (3.0 - 2.0 * bank)
            best = max(best, bank)
        return best

    def path_crosses_hex(self, hex_idx: int) -> bool:
        """Does a road or river reach into the hex's flat pad (the region
        the standability check measures)? A plateau hex it crosses cannot
        be standable; the road bed itself is flat but sits a step lower."""
        cx, cy = self.centres[hex_idx]
        pad = self.params.pad_min_r * self.apothem
        for path in self.roads + self.rivers:
            reach = (
                self.params.road_half_width_mm + self.params.road_edge_mm
                if path.kind == "road"
                else self.params.river_half_width_mm
            )
            d, _ = path.dist_s(cx, cy)
            if d < reach + pad:
                return True
        return False


def _rkey(p: Point2D) -> tuple[float, float]:
    return (round(p[0], 3), round(p[1], 3))
