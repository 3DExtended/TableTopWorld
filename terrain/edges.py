"""Exterior edge bevels and magnet holes."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Tuple

import numpy as np
from solid2 import cylinder, difference, polyhedron, union

from terrain.constants import (
    HEXAGON_BEVEL_SIZE,
    MAGNET_DEPTH,
    MAGNET_HEIGHT_OVER_GROUND,
    MAGNET_RADIUS,
)
from terrain.layout import FlowerLayout, Line3D

if TYPE_CHECKING:
    from solid2 import OpenSCADObject

    from terrain.tileset import FlowerDef, Tileset

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


def add_bevel(
    obj: OpenSCADObject | None,
    line: Line3D,
    depth: float,
    z_anchor: float,
) -> OpenSCADObject:
    """Port of addBevel; subtracts bevel prism at z_anchor."""
    line = (
        (
            line[0][0] + 0.0001 * (line[1][0] - line[0][0]),
            line[0][1] + 0.0001 * (line[1][1] - line[0][1]),
            line[0][2] + 0.0001 * (line[1][2] - line[0][2]),
        ),
        (
            line[1][0] + 0.0001 * (line[0][0] - line[1][0]),
            line[1][1] + 0.0001 * (line[0][1] - line[1][1]),
            line[1][2] + 0.0001 * (line[0][2] - line[1][2]),
        ),
    )
    shifted_line_1 = shift_line_3d(line, depth * 0.75, 0, 0)
    shifted_line_2 = shift_line_3d(line, -depth * 0.75, 0, 0)
    shifted_line_3 = shift_line_3d(line, 0, -depth, 0)
    points = [
        shifted_line_1[0],
        shifted_line_1[1],
        shifted_line_2[0],
        shifted_line_2[1],
        shifted_line_3[0],
        shifted_line_3[1],
    ]
    triangle = polyhedron(
        points=points,
        faces=[
            (0, 2, 4),
            (1, 5, 3),
            (0, 4, 1),
            (4, 5, 1),
            (4, 2, 3),
            (5, 4, 3),
            (2, 0, 1),
            (3, 2, 1),
        ],
    ).translateZ(z_anchor)
    if obj is None:
        return triangle
    return obj - triangle


def angle_with_x_axis(line_for_hole: Line2D) -> tuple[float, float]:
    point1, point2 = line_for_hole
    dx = point2[0] - point1[0]
    dy = point2[1] - point1[1]
    angle_rad = math.atan2(dy, dx)
    return angle_rad, math.degrees(angle_rad)


def add_magnet_hole_on_side(
    obj: OpenSCADObject,
    line_for_hole: Line2D,
    magnet_height: float = MAGNET_HEIGHT_OVER_GROUND,
    magnet_radius: float = MAGNET_RADIUS,
    magnet_depth: float = MAGNET_DEPTH,
) -> OpenSCADObject:
    """Port of addMagnetHoleOnSide."""
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
        .translateZ(magnet_height + magnet_radius)
    )
    return obj - cylinder_tool


class EdgeGeometry:
    """Apply exterior bevels and magnets to a flower solid."""

    def __init__(self, layout: FlowerLayout, bevel_size: float = HEXAGON_BEVEL_SIZE) -> None:
        self.layout = layout
        self.bevel_size = bevel_size

    def apply_edges(
        self,
        flower_solid: OpenSCADObject,
        flower: FlowerDef,
        max_z: float,
    ) -> OpenSCADObject:
        tools: list[OpenSCADObject] = []
        tool_settings: list[Line3D] = []
        z_anchor = max(max_z, 1.2)

        for hex_idx in FlowerLayout.RING_HEX_INDICES:
            verts = self.layout.ring_vertices(hex_idx)
            length = len(verts)
            for i in range(length):
                p1 = verts[i]
                p2 = verts[(i + 1) % length]
                angle_of_line = np.cross(
                    np.array([p1[0], p1[1]]), np.array([p2[0], p2[1]])
                )
                if angle_of_line > 0.1:
                    flower_solid = add_magnet_hole_on_side(flower_solid, (p1, p2))
                bevel_settings: Line3D = (
                    (p1[0], p1[1], z_anchor),
                    (p2[0], p2[1], z_anchor),
                )
                if bevel_settings not in tool_settings:
                    tool = add_bevel(None, bevel_settings, self.bevel_size, z_anchor)
                    tools.append(tool)
                    tool_settings.append(bevel_settings)

        if tools:
            flower_solid = difference()(flower_solid, union()(tools))
        return flower_solid
