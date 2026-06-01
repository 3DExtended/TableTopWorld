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


def axial_to_xy(q: int, r: int, spacing: float) -> Point2D:
    """Axial hex coords to world XY (pointy-top)."""
    x = spacing * (math.sqrt(3) * q + math.sqrt(3) / 2 * r)
    y = spacing * (1.5 * r)
    return (x, y)


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


class FlowerLayout:
    """Flower graph: 7 cells, 18 exterior edges, 6 road/water junctions."""

    # Ring hex 1–6; sides 0–2 are outward-facing (legacy keys "h-s").
    RING_HEX_INDICES = tuple(range(1, 7))
    EXTERIOR_SIDES_PER_HEX = 3
    JUNCTION_COUNT = 6

    def __init__(self, hex_outer_width: float = 5.1961525) -> None:
        self.hex_outer_width = hex_outer_width
        self.inner_hexagon_size = convert_outer_hexagon_size_to_inner(hex_outer_width)
        self.flower_center_spacing = self.inner_hexagon_size * 4.0
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

    def exterior_edge_line(self, hex_idx: int, side_idx: int) -> Line2D:
        """Outward edge segment for ring hex (matches legacy vertex indexing)."""
        if hex_idx not in self.RING_HEX_INDICES:
            raise ValueError(f"hex {hex_idx} has no exterior edges")
        if side_idx not in range(self.EXTERIOR_SIDES_PER_HEX):
            raise ValueError(f"side {side_idx} out of range for hex {hex_idx}")
        verts = self.ring_vertices(hex_idx)
        # Legacy getOuterHexFlowerLines uses indices (2+index-1)%6 etc. for junction 0
        # on hex 1: sides map to vertex pairs used in bevel/magnet code
        vi = (2 + side_idx) % 6
        vj = (1 + side_idx) % 6
        return (verts[vi], verts[vj])

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
