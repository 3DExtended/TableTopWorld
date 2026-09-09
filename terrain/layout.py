"""7-hex flower layout: positions, polygons, exterior edges, junction lines."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence, Tuple

Line2D = Tuple[Tuple[float, float], Tuple[float, float]]
Line3D = Tuple[Tuple[float, float, float], Tuple[float, float, float]]
Point2D = Tuple[float, float]


def _generate_hexagon_vertices(radius: float) -> list[Point2D]:
    vertices: list[Point2D] = []
    for i in range(6):
        angle_rad = math.radians(60 * i)
        vertices.append((radius * math.cos(angle_rad), radius * math.sin(angle_rad)))
    return vertices


def subdivide_hexagon(hexagon_vertices: Sequence[Point2D]) -> list[Point2D]:
    center = (0.0, 0.0)
    triangles: list[Point2D] = []
    n = len(hexagon_vertices)
    for i in range(n):
        v1 = hexagon_vertices[i]
        v2 = hexagon_vertices[(i + 1) % n]
        triangles.extend([center, v1, v2])
    return triangles


def convert_outer_hexagon_size_to_inner(outer_hexagon_size: float) -> float:
    angle_rad = math.radians(60 * 0)
    temp_x_1 = outer_hexagon_size * math.cos(angle_rad)
    temp_y_1 = outer_hexagon_size * math.sin(angle_rad)
    angle_rad = math.radians(60 * 1)
    temp_x_2 = outer_hexagon_size * math.cos(angle_rad)
    temp_y_2 = outer_hexagon_size * math.sin(angle_rad)
    inner_x = temp_x_1 + 0.5 * (temp_x_2 - temp_x_1)
    inner_y = temp_y_1 + 0.5 * (temp_y_2 - temp_y_1)
    return math.sqrt(inner_x**2 + inner_y**2)


def rotate_point_2d(point: Point2D, degrees: float) -> Point2D:
    rad = math.radians(degrees)
    c, s = math.cos(rad), math.sin(rad)
    x, y = point
    return (x * c - y * s, x * s + y * c)


@dataclass(frozen=True)
class ExteriorEdge:
    key: str
    hex_idx: int
    side_idx: int
    line_2d: Line2D


@dataclass(frozen=True)
class HexBevelEdge:
    """One of six sides on a flower cell; chamfer targets the hex top rim."""

    key: str
    hex_idx: int
    edge_idx: int
    line_2d: Line2D


class FlowerLayout:
    """Flower graph: 7 cells, 18 exterior edges, 6 road/water junctions."""

    # Ring hex 1–6; sides 0–2 are outward-facing (legacy keys "h-s").
    RING_HEX_INDICES = tuple(range(1, 7))
    HEX_CELL_COUNT = 7
    EDGES_PER_HEX = 6
    EXTERIOR_SIDES_PER_HEX = 3
    JUNCTION_COUNT = 6
    # A "side" is the boundary run shared with one potential neighbor flower:
    # 3 consecutive exterior edges spanning 2 ring hexes, 4 corner vertices.
    SIDE_COUNT = 6
    EDGES_PER_SIDE = 3

    def __init__(self, hex_outer_width: float = 5.1961525) -> None:
        self.hex_outer_width = hex_outer_width
        self.inner_hexagon_size = convert_outer_hexagon_size_to_inner(hex_outer_width)
        self._cell_centers: dict[int, Point2D] = {0: (0.0, 0.0)}
        self._cell_polygons: dict[int, list[Point2D]] = {}
        self._ring_vertices: list[list[Point2D]] = []
        self._build_cells()

    def _ring_center(self, ring_index: int) -> Point2D:
        """ring_index 0..5 maps to hex cells 1..6."""
        angle_rad = math.radians(60 * ring_index + 30)
        x = self.hex_outer_width * math.cos(angle_rad)
        y = self.hex_outer_width * math.sin(angle_rad)
        magnitude = 1.0 / math.sqrt(x**2 + y**2)
        return (
            magnitude * x * self.inner_hexagon_size * 2,
            magnitude * y * self.inner_hexagon_size * 2,
        )

    def _build_cells(self) -> None:
        center_verts = _generate_hexagon_vertices(self.hex_outer_width)
        self._cell_polygons[0] = subdivide_hexagon(center_verts)
        for i in range(6):
            hex_idx = i + 1
            center = self._ring_center(i)
            self._cell_centers[hex_idx] = center
            verts = [
                (v[0] + center[0], v[1] + center[1]) for v in center_verts
            ]
            self._ring_vertices.append(verts)
            self._cell_polygons[hex_idx] = subdivide_hexagon(center_verts)
            # Translate polygon points to cell center
            self._cell_polygons[hex_idx] = [
                (p[0] + center[0], p[1] + center[1])
                for p in subdivide_hexagon(center_verts)
            ]

    @staticmethod
    def edge_id(hex_idx: int, side_idx: int) -> str:
        return f"{hex_idx}-{side_idx}"

    @classmethod
    def parse_edge_id(cls, key: str) -> tuple[int, int]:
        parts = key.split("-")
        if len(parts) != 2:
            raise ValueError(f"invalid edge key {key!r}")
        return int(parts[0]), int(parts[1])

    def exterior_edge_keys(self) -> list[str]:
        keys: list[str] = []
        for h in self.RING_HEX_INDICES:
            for s in range(self.EXTERIOR_SIDES_PER_HEX):
                keys.append(self.edge_id(h, s))
        return keys

    def cell_center(self, hex_idx: int) -> Point2D:
        return self._cell_centers[hex_idx]

    def cell_polygon(self, hex_idx: int) -> list[Point2D]:
        return self._cell_polygons[hex_idx]

    def ring_vertices(self, hex_idx: int) -> list[Point2D]:
        if hex_idx == 0:
            return [
                (v[0], v[1]) for v in _generate_hexagon_vertices(self.hex_outer_width)
            ]
        return self._ring_vertices[hex_idx - 1]

    def exterior_vertex_indices(self, hex_idx: int) -> tuple[int, int, int]:
        """Vertex indices i for the three edges (i, i+1) on the flower exterior.

        Ring hexes share the same local vertex winding but sit at different angles
        around the center, so outward edges are not the same index triple on every
        hex. side_idx 0..2 maps to these indices in CCW order around the cell.
        """
        if hex_idx not in self.RING_HEX_INDICES:
            raise ValueError(f"hex {hex_idx} has no exterior edges")
        verts = self.ring_vertices(hex_idx)
        center = self.cell_center(hex_idx)
        dist_center = math.hypot(center[0], center[1])
        candidates: list[tuple[float, int]] = []
        for i in range(6):
            p1, p2 = verts[i], verts[(i + 1) % 6]
            mid = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
            if math.hypot(mid[0], mid[1]) <= dist_center + 1e-9:
                continue
            angle = math.atan2(mid[1] - center[1], mid[0] - center[0])
            candidates.append((angle, i))
        if len(candidates) != self.EXTERIOR_SIDES_PER_HEX:
            found = [i for _, i in candidates]
            raise RuntimeError(
                f"hex {hex_idx}: expected {self.EXTERIOR_SIDES_PER_HEX} exterior edges, "
                f"found {len(candidates)} ({found})"
            )
        candidates.sort(key=lambda item: item[0])
        return (candidates[0][1], candidates[1][1], candidates[2][1])

    def exterior_edge_line(self, hex_idx: int, side_idx: int) -> Line2D:
        """Outward-facing edge segment for a ring hex side (0..2)."""
        if side_idx not in range(self.EXTERIOR_SIDES_PER_HEX):
            raise ValueError(f"side {side_idx} out of range for hex {hex_idx}")
        verts = self.ring_vertices(hex_idx)
        vi = self.exterior_vertex_indices(hex_idx)[side_idx]
        vj = (vi + 1) % 6
        return (verts[vi], verts[vj])

    def hex_edge_line(self, hex_idx: int, edge_idx: int) -> Line2D:
        """CCW edge segment on a cell (edge_idx 0..5 joins vertex i to i+1)."""
        if hex_idx not in range(self.HEX_CELL_COUNT):
            raise ValueError(f"hex {hex_idx} out of range")
        if edge_idx not in range(self.EDGES_PER_HEX):
            raise ValueError(f"edge {edge_idx} out of range for hex {hex_idx}")
        verts = self.ring_vertices(hex_idx)
        return (verts[edge_idx], verts[(edge_idx + 1) % self.EDGES_PER_HEX])

    def hex_bevel_edges(self) -> list[HexBevelEdge]:
        """All six sides of every flower cell (magnets use exterior_edges only)."""
        edges: list[HexBevelEdge] = []
        for hex_idx in range(self.HEX_CELL_COUNT):
            for edge_idx in range(self.EDGES_PER_HEX):
                edges.append(
                    HexBevelEdge(
                        key=self.edge_id(hex_idx, edge_idx),
                        hex_idx=hex_idx,
                        edge_idx=edge_idx,
                        line_2d=self.hex_edge_line(hex_idx, edge_idx),
                    )
                )
        return edges

    def exterior_edges(self) -> list[ExteriorEdge]:
        edges: list[ExteriorEdge] = []
        for h in self.RING_HEX_INDICES:
            for s in range(self.EXTERIOR_SIDES_PER_HEX):
                edges.append(
                    ExteriorEdge(
                        key=self.edge_id(h, s),
                        hex_idx=h,
                        side_idx=s,
                        line_2d=self.exterior_edge_line(h, s),
                    )
                )
        return edges

    # The 6 canonical flower-to-flower neighbor directions, expressed as
    # small-hex axial offsets (q, s). Derived from one verified 3-shared-edge
    # offset, (3, -1), rotated by the axial 60-degree rotation formula
    # (q, s) -> (-s, q + s). Brute-force verified (grill session) to be the
    # maximum possible shared-edge count between two validly-tiled neighbor
    # flowers: 3 edges / 4 corners, never 1 or more than 3.
    #
    # An earlier, wrong placement convention (a generic pointy-top
    # axial_to_xy()/flower_center_spacing formula) did NOT correspond to
    # true edge-sharing adjacency (confirmed numerically at the time) -
    # removed once flower_grid_to_xy() (below) was verified as the correct
    # replacement, used by terrain/assembly.py for preview_map placement.
    NEIGHBOR_SMALL_HEX_OFFSETS: tuple[tuple[int, int], ...] = (
        (3, -1), (1, 2), (-2, 3), (-3, 1), (-1, -2), (2, -3),
    )

    def _small_hex_basis_vector(self, k: int) -> Point2D:
        """World direction of one single-hex-to-hex step k (0..5), matching
        _ring_center's own angle convention (60*k + 30 degrees) so this stays
        consistent with the rest of this class's actual hex arrangement."""
        angle_rad = math.radians(60 * k + 30)
        r = self.hex_outer_width
        return (r * math.sqrt(3) * math.cos(angle_rad), r * math.sqrt(3) * math.sin(angle_rad))

    def neighbor_flower_offset(self, direction_idx: int) -> Point2D:
        """World XY center offset to place a TRUE edge-sharing neighbor flower.

        direction_idx 0..5 selects one of the 6 canonical neighbor directions
        (60 degrees apart). No rotation is needed for the placed flower - see
        design decision #2. Expressed as q*e0 + s*e1 where e0/e1 are the
        single-hex-step basis vectors for steps 0 and 1 (this class's own
        angle convention, NOT a generic pointy-top axial formula).
        """
        if direction_idx not in range(6):
            raise ValueError(f"direction {direction_idx} must be 0..5")
        q, s = self.NEIGHBOR_SMALL_HEX_OFFSETS[direction_idx]
        e0 = self._small_hex_basis_vector(0)
        e1 = self._small_hex_basis_vector(1)
        return (q * e0[0] + s * e1[0], q * e0[1] + s * e1[1])

    def flower_grid_to_xy(self, q: int, r: int) -> Point2D:
        """World XY center for a flower placed at flower-grid axial
        coordinates (q, r) - the correct replacement for the old, wrong
        axial_to_xy()/flower_center_spacing (see the NOTE above
        NEIGHBOR_SMALL_HEX_OFFSETS: that pairing does not correspond to
        true edge-sharing adjacency at all).

        Standard axial hex coordinates: q*e0 + r*e1, where e0/e1 are
        neighbor_flower_offset(0)/(1) - verified numerically that this
        exactly reproduces all 6 neighbor_flower_offset(k) directions at
        the expected small-integer (q, r) combinations
        ((1,0), (0,1), (-1,1), (-1,0), (0,-1), (1,-1) for k=0..5), so it
        generalizes correctly to any flower placement, not just direct
        neighbors.
        """
        e0 = self.neighbor_flower_offset(0)
        e1 = self.neighbor_flower_offset(1)
        return (q * e0[0] + r * e1[0], q * e0[1] + r * e1[1])

    def _side_groups(self) -> list[list[ExteriorEdge]]:
        """Group the 18 exterior edges into 6 boundary sides.

        Each side is 3 edges (spanning 2 ring hexes: 2 edges of one + 1 of
        its ring-neighbor, or vice versa) shared with one potential neighbor
        flower. A naive "sort 18 edges by angle, chunk into runs of 3"
        actually groups each *ring hex's own* 3 edges together instead
        (verified against neighbor_flower_offset() - that grouping produces
        zero true edge-for-edge matches with a real adjacent flower). The
        correct grouping instead buckets each edge by which of the 6 real
        neighbor_flower_offset() directions its midpoint angle is closest
        to, since that's the direction a physically-adjacent flower is
        actually placed in.
        """
        edges = self.exterior_edges()

        def edge_angle(e: ExteriorEdge) -> float:
            (x1, y1), (x2, y2) = e.line_2d
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2
            return math.atan2(my, mx)

        def angular_diff(a: float, b: float) -> float:
            d = abs(a - b) % (2 * math.pi)
            return min(d, 2 * math.pi - d)

        direction_angles = [
            math.atan2(*reversed(self.neighbor_flower_offset(k)))
            for k in range(self.SIDE_COUNT)
        ]

        buckets: list[list[ExteriorEdge]] = [[] for _ in range(self.SIDE_COUNT)]
        for e in edges:
            ea = edge_angle(e)
            best_k = min(
                range(self.SIDE_COUNT),
                key=lambda k: angular_diff(ea, direction_angles[k]),
            )
            buckets[best_k].append(e)

        def signed_diff(a: float, b: float) -> float:
            """a - b, wrapped into (-pi, pi] - avoids +-180 degree sort bugs
            that a plain angle sort would hit near the wraparound point."""
            return (a - b + math.pi) % (2 * math.pi) - math.pi

        for k, bucket in enumerate(buckets):
            if len(bucket) != self.EDGES_PER_SIDE:
                raise RuntimeError(
                    f"side {k}: expected {self.EDGES_PER_SIDE} edges, "
                    f"got {len(bucket)} ({[e.key for e in bucket]})"
                )
            center = direction_angles[k]
            bucket.sort(key=lambda e: signed_diff(edge_angle(e), center))
        return buckets

    @classmethod
    def path_route(cls, entry: int, exit_side: int, via: Sequence[int] = ()) -> list[int]:
        """Hex cells a road or river visits between the middle edges of
        sides `entry` and `exit_side`: the ring hex behind each side, any
        `via` hexes between them, and the centre hex when the two sides are
        not adjacent and no via is given. Consecutive repeats are dropped."""
        ring_in, ring_out = entry + 1, exit_side + 1
        route = [ring_in, *via, ring_out]
        if not via and ring_in != ring_out and (exit_side - entry) % cls.SIDE_COUNT not in (1, 5):
            route = [ring_in, 0, ring_out]
        return [h for i, h in enumerate(route) if i == 0 or h != route[i - 1]]

    @staticmethod
    def side_id(side_idx: int) -> str:
        return f"side-{side_idx}"

    def side_corners(
        self, side_idx: int
    ) -> tuple[Point2D, Point2D, Point2D, Point2D]:
        """The 4 corner points of one of the flower's 6 sides, in boundary order.

        Two flowers may only share a side if their declared corner heights
        match in this same order (see design decision #3).
        """
        if side_idx not in range(self.SIDE_COUNT):
            raise ValueError(f"side {side_idx} must be 0..{self.SIDE_COUNT - 1}")
        chain_edges = self._side_groups()[side_idx]

        def close(a: Point2D, b: Point2D, tol: float = 1e-6) -> bool:
            return math.hypot(a[0] - b[0], a[1] - b[1]) < tol

        points = [chain_edges[0].line_2d[0], chain_edges[0].line_2d[1]]
        for e in chain_edges[1:]:
            p1, p2 = e.line_2d
            if close(p1, points[-1]):
                points.append(p2)
            elif close(p2, points[-1]):
                points.append(p1)
            else:
                raise RuntimeError(
                    f"side {side_idx}: edge {e.key} does not chain from "
                    f"previous corner {points[-1]}"
                )
        if len(points) != self.EDGES_PER_SIDE + 1:
            raise RuntimeError(
                f"side {side_idx}: expected {self.EDGES_PER_SIDE + 1} corners, "
                f"got {len(points)}"
            )
        return (points[0], points[1], points[2], points[3])

    def junction_lines(self, junction_idx: int) -> tuple[Line2D, Line2D, Line2D]:
        """Three hex-side lines meeting at junction (port of getOuterHexFlowerLines)."""
        if junction_idx not in range(self.JUNCTION_COUNT):
            raise ValueError(f"junction {junction_idx} must be 0..5")
        index = junction_idx
        top_hexagon_index = index % 6
        hex_right_index = (index - 1) % 6
        top_hex = self._ring_vertices[top_hexagon_index]
        hex_right = self._ring_vertices[hex_right_index]
        top_line = (top_hex[(2 + index - 1) % 6], top_hex[(1 + index - 1) % 6])
        top_right_line = (
            top_hex[(2 + index - 1 - 1) % 6],
            top_hex[(2 + index - 1 - 2) % 6],
        )
        adjacent_line = (
            hex_right[(2 + index - 1) % 6],
            hex_right[(1 + index - 1) % 6],
        )
        return (top_line, top_right_line, adjacent_line)

    @staticmethod
    def center_of_three_lines(
        line_a: Line2D, line_b: Line2D, line_c: Line2D
    ) -> Point2D:
        """Port of getCenterOfThreeLines."""
        cx = line_a[0][0] + 0.5 * (line_c[1][0] - line_a[0][0])
        cy = line_a[0][1] + 0.5 * (line_c[1][1] - line_a[0][1])
        return (cx, cy)

    def junction_center(self, junction_idx: int) -> Point2D:
        a, b, c = self.junction_lines(junction_idx)
        return self.center_of_three_lines(a, b, c)

    def transform_edge_to_world(
        self,
        edge: ExteriorEdge,
        origin: Point2D,
        rot_steps: int,
    ) -> tuple[Point2D, Point2D, Point2D]:
        """Return midpoint, outward normal (unit), and second endpoint in world space."""
        deg = rot_steps * 60
        p1 = rotate_point_2d(edge.line_2d[0], deg)
        p2 = rotate_point_2d(edge.line_2d[1], deg)
        p1 = (p1[0] + origin[0], p1[1] + origin[1])
        p2 = (p2[0] + origin[0], p2[1] + origin[1])
        mid = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
        dx, dy = p2[0] - p1[0], p2[1] - p1[1]
        length = math.hypot(dx, dy) or 1.0
        tangent = (dx / length, dy / length)
        # Outward normal (perpendicular, away from flower center)
        center = rotate_point_2d((0.0, 0.0), deg)
        center = (center[0] + origin[0], center[1] + origin[1])
        to_mid = (mid[0] - center[0], mid[1] - center[1])
        n1 = (-tangent[1], tangent[0])
        n2 = (tangent[1], -tangent[0])
        if n1[0] * to_mid[0] + n1[1] * to_mid[1] > 0:
            normal = n1
        else:
            normal = n2
        return mid, normal, p2
