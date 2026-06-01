"""Exterior edge bevels and magnet holes."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Tuple

import numpy as np
from solid2 import cylinder, difference, linear_extrude, polygon, polyhedron, union

from terrain.catalog import EdgeProfileCatalog
from terrain.constants import (
    BEVEL_TOP_SLAB_DEPTH,
    BEVEL_Z_INSET,
    HEXAGON_BEVEL_SIZE,
    MAGNET_CENTER_Z,
    MAGNET_DEPTH,
    MAGNET_RADIUS,
)
from terrain.layout import FlowerLayout, Line3D

if TYPE_CHECKING:
    from solid2 import OpenSCADObject

    from terrain.tileset import FlowerDef

Line2D = Tuple[Tuple[float, float], Tuple[float, float]]


def shift_line_3d(line: Line3D, x: float, y: float, z: float) -> Line3D:
    p1, p2 = line
    direction_vector = np.array([p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]])
    magnitude = np.linalg.norm(direction_vector)
    if magnitude == 0:
        return line
    unit_direction = direction_vector / magnitude
    if unit_direction[0] == 0 and unit_direction[1] == 0:
        perpendicular_vector1 = np.array([1.0, 0.0, 0.0])
    else:
        perpendicular_vector1 = np.array(
            [-unit_direction[1], unit_direction[0], 0.0]
        )
    perpendicular_vector1 /= np.linalg.norm(perpendicular_vector1)
    perpendicular_vector2 = np.cross(unit_direction, perpendicular_vector1)
    translated_p1 = p1 + x * perpendicular_vector1 + y * perpendicular_vector2 + z * unit_direction
    translated_p2 = p2 + x * perpendicular_vector1 + y * perpendicular_vector2 + z * unit_direction
    return (tuple(translated_p1), tuple(translated_p2))


_BEVEL_POLYHEDRON_FACES = [
    (0, 2, 4),
    (1, 5, 3),
    (0, 4, 1),
    (4, 5, 1),
    (4, 2, 3),
    (5, 4, 3),
    (2, 0, 1),
    (3, 2, 1),
]


def _bevel_line_xy(line: Line3D) -> Line3D:
    return (
        (
            line[0][0] + 0.0001 * (line[1][0] - line[0][0]),
            line[0][1] + 0.0001 * (line[1][1] - line[0][1]),
            0.0,
        ),
        (
            line[1][0] + 0.0001 * (line[0][0] - line[1][0]),
            line[1][1] + 0.0001 * (line[0][1] - line[1][1]),
            0.0,
        ),
    )


def _bevel_polyhedron(
    line_in: Line3D,
    line_out: Line3D,
    line_apex: Line3D,
) -> OpenSCADObject:
    """Six-point prism along edge; indices 0–1 in, 2–3 out, 4–5 apex (legacy winding)."""
    return polyhedron(
        points=[
            line_in[0],
            line_in[1],
            line_out[0],
            line_out[1],
            line_apex[0],
            line_apex[1],
        ],
        faces=_BEVEL_POLYHEDRON_FACES,
    )


def _bevel_side_wedge(line: Line3D, depth: float) -> OpenSCADObject:
    """Chamfer on the exterior vertical wall (legacy orientation: +0.75d is inward)."""
    line_in = shift_line_3d(line, depth * 0.75, 0, 0)
    line_out = shift_line_3d(line, -depth * 0.75, 0, 0)
    line_apex = shift_line_3d(line, 0, -depth, 0)
    return _bevel_polyhedron(line_in, line_out, line_apex)


def _bevel_top_shelf_wedge(line: Line3D, depth: float) -> OpenSCADObject:
    """Sloped cut from the exterior edge toward hex interior (45° top chamfer)."""
    line_in = shift_line_3d(line, depth * 0.75, 0, 0)
    line_apex = shift_line_3d(line, depth * 0.75, depth, 0)
    return _bevel_polyhedron(line_in, line, line_apex)


def _bevel_top_flat_box(line: Line3D, depth: float) -> OpenSCADObject:
    """Axis-aligned box along the edge; bites through the flat top face at z_anchor."""
    inward = depth * 0.75
    half = BEVEL_TOP_SLAB_DEPTH / 2
    edge_lo = shift_line_3d(line, 0, -half, 0)
    edge_hi = shift_line_3d(line, 0, half, 0)
    in_lo = shift_line_3d(line, inward, -half, 0)
    in_hi = shift_line_3d(line, inward, half, 0)
    points = [
        edge_lo[0],
        edge_lo[1],
        in_lo[0],
        in_lo[1],
        edge_hi[0],
        edge_hi[1],
        in_hi[0],
        in_hi[1],
    ]
    faces = [
        (0, 1, 3),
        (0, 3, 2),
        (4, 6, 7),
        (4, 7, 5),
        (0, 4, 5),
        (0, 5, 1),
        (2, 3, 7),
        (2, 7, 6),
        (0, 2, 6),
        (0, 6, 4),
        (1, 5, 7),
        (1, 7, 3),
    ]
    return polyhedron(points=points, faces=faces)


def add_bevel(
    obj: OpenSCADObject | None,
    line: Line3D,
    depth: float,
    z_anchor: float,
) -> OpenSCADObject:
    """Subtract chamfer at z_anchor through the top exterior corner."""
    line_xy = _bevel_line_xy(line)
    z_side = z_anchor - BEVEL_Z_INSET
    top = (
        _bevel_top_shelf_wedge(line_xy, depth) + _bevel_top_flat_box(line_xy, depth)
    ).translateZ(z_anchor)
    tool = _bevel_side_wedge(line_xy, depth).translateZ(z_side) + top
    if obj is None:
        return tool
    return obj - tool


def angle_with_x_axis(line_for_hole: Line2D) -> tuple[float, float]:
    point1, point2 = line_for_hole
    dx = point2[0] - point1[0]
    dy = point2[1] - point1[1]
    angle_rad = math.atan2(dy, dx)
    return angle_rad, math.degrees(angle_rad)


def magnet_wall_slab(
    line: Line2D,
    z_bottom: float,
    z_top: float,
    cell_center: tuple[float, float],
    *,
    thickness: float = 1.2,
    inset: float = 0.45,
) -> OpenSCADObject:
    """Volume along an exterior edge, shifted inward so it merges with the hex body."""
    p1, p2 = line
    dx, dy = p2[0] - p1[0], p2[1] - p1[1]
    length = math.hypot(dx, dy) or 1.0
    angle = math.degrees(math.atan2(dy, dx))
    mid = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
    cx, cy = cell_center
    to_center = (cx - mid[0], cy - mid[1])
    dist = math.hypot(to_center[0], to_center[1]) or 1.0
    place = (
        mid[0] + inset * to_center[0] / dist,
        mid[1] + inset * to_center[1] / dist,
    )
    height = z_top - z_bottom
    bar = linear_extrude(height=height)(
        polygon(
            points=[
                (-length / 2, -thickness / 2),
                (length / 2, -thickness / 2),
                (length / 2, thickness / 2),
                (-length / 2, thickness / 2),
            ]
        )
    )
    return bar.rotateZ(angle).translateX(place[0]).translateY(place[1]).translateZ(z_bottom)


def add_magnet_hole_on_side(
    obj: OpenSCADObject,
    line_for_hole: Line2D,
    magnet_center_z: float,
    magnet_radius: float = MAGNET_RADIUS,
    magnet_depth: float = MAGNET_DEPTH,
) -> OpenSCADObject:
    """Port of addMagnetHoleOnSide; center Z is offset from FLOWER_BOTTOM_Z."""
    cylinder_tool = cylinder(h=magnet_depth * 4, center=True, r=magnet_radius)
    center_of_line = (
        line_for_hole[1][0] - 0.5 * (line_for_hole[1][0] - line_for_hole[0][0]),
        line_for_hole[1][1] - 0.5 * (line_for_hole[1][1] - line_for_hole[0][1]),
    )
    _, angle_deg = angle_with_x_axis(line_for_hole)
    cylinder_tool = (
        cylinder_tool.rotateX(90)
        .rotateZ(angle_deg)
        .translateX(center_of_line[0])
        .translateY(center_of_line[1])
        .translateZ(magnet_center_z)
    )
    return obj - cylinder_tool


class EdgeGeometry:
    """Apply exterior bevels and magnets to a flower solid."""

    def __init__(self, layout: FlowerLayout, bevel_size: float = HEXAGON_BEVEL_SIZE) -> None:
        self.layout = layout
        self.bevel_size = bevel_size

    def apply_magnets(
        self,
        flower_solid: OpenSCADObject,
        flower: FlowerDef,
        catalog: EdgeProfileCatalog,
    ) -> OpenSCADObject:
        for edge in self.layout.exterior_edges():
            profile = flower.edges[edge.key]
            if not catalog.is_mating_profile(profile):
                continue
            flower_solid = add_magnet_hole_on_side(
                flower_solid, edge.line_2d, MAGNET_CENTER_Z
            )
        return flower_solid

    def apply_bevels(
        self,
        flower_solid: OpenSCADObject,
        max_z: float,
    ) -> OpenSCADObject:
        tools: list[OpenSCADObject] = []
        tool_settings: list[Line3D] = []
        z_anchor = max(max_z, 1.2)

        for edge in self.layout.exterior_edges():
            bevel_settings: Line3D = (
                (edge.line_2d[0][0], edge.line_2d[0][1], z_anchor),
                (edge.line_2d[1][0], edge.line_2d[1][1], z_anchor),
            )
            if bevel_settings not in tool_settings:
                tool = add_bevel(None, bevel_settings, self.bevel_size, z_anchor)
                tools.append(tool)
                tool_settings.append(bevel_settings)

        if tools:
            flower_solid = difference()(flower_solid, union()(tools))
        return flower_solid

    def apply_edges(
        self,
        flower_solid: OpenSCADObject,
        flower: FlowerDef,
        catalog: EdgeProfileCatalog,
        max_z: float,
    ) -> OpenSCADObject:
        flower_solid = self.apply_bevels(flower_solid, max_z)
        return self.apply_magnets(flower_solid, flower, catalog)
