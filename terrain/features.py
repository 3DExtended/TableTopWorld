"""Road and water channel cutters."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from solid2 import union
from solid2.extensions.bosl2 import beziers, path_sweep, regular_ngon

from terrain.constants import (
    HEXAGON_BEVEL_SIZE,
    STREET_INDENT_HEIGHT,
    STREET_WIDTH_SCALAR,
    WATER_INDENT_HEIGHT,
    WATER_WIDTH_SCALAR,
)
from terrain.edges import add_bevel
from terrain.layout import FlowerLayout

if TYPE_CHECKING:
    from solid2 import OpenSCADObject

    from terrain.tileset import FlowerDef, Tileset


class FeatureCutters:
    """Bezier path_sweep road/water subtraction at host terrain Z."""

    def __init__(
        self,
        layout: FlowerLayout,
        resolution: int = 100,
    ) -> None:
        self.layout = layout
        self.resolution = resolution

    def _path_pair_tools(
        self,
        entry_junction: int,
        exit_junction: int,
        *,
        n_gon: int,
        spin: float,
        width_scalar: float,
        indent_height: float,
        host_z: float,
        counter: int,
    ) -> list[OpenSCADObject]:
        entry_a, entry_b, entry_c = self.layout.junction_lines(entry_junction)
        center_entry = self.layout.center_of_three_lines(entry_a, entry_b, entry_c)
        exit_a, exit_b, exit_c = self.layout.junction_lines(exit_junction)
        center_exit = self.layout.center_of_three_lines(exit_a, exit_b, exit_c)

        entry_dir = np.array(
            [entry_c[1][0] - entry_a[0][0], entry_c[1][1] - entry_a[0][1], 0.0]
        )
        exit_dir = np.array(
            [exit_c[1][0] - exit_a[0][0], exit_c[1][1] - exit_a[0][1], 0.0]
        )
        entry_mag = np.linalg.norm(entry_dir) or 1.0
        exit_mag = np.linalg.norm(exit_dir) or 1.0
        unit_down = np.array([0.0, 0.0, -1.0])
        perp_entry = np.cross(unit_down, entry_dir / entry_mag)
        perp_exit = np.cross(unit_down, exit_dir / exit_mag)

        midpoint = np.array(
            [center_exit[0] - center_entry[0], center_exit[1] - center_entry[1], 0.0]
        )
        midpoint *= 0.1 / (np.linalg.norm(midpoint) or 1.0)

        path_width = (
            np.linalg.norm(
                np.array([entry_c[1][0], entry_c[1][1]])
                - np.array([entry_a[0][0], entry_a[0][1]])
            )
            * width_scalar
        )
        diagonaled_width = math.sqrt(path_width**2 * 2)

        sbez = [
            [center_entry[0], center_entry[1]],
            [-midpoint[0], -midpoint[1]],
            [0.0, 0.0],
            [center_exit[0], center_exit[1]],
        ]
        tool = path_sweep(
            regular_ngon(n=n_gon, d=diagonaled_width, spin=spin),
            beziers.bezpath_curve(
                sbez, N=len(sbez) - 1, splinesteps=self.resolution
            ),
        )

        entry_points = [
            np.array(sbez[0]) * 0.9,
            np.array(center_entry) * 5,
            np.array(center_entry) * 10,
        ]
        tool += path_sweep(
            regular_ngon(n=n_gon, d=diagonaled_width, spin=spin),
            beziers.bezpath_curve(
                entry_points, N=len(entry_points) - 1, splinesteps=self.resolution
            ),
        )
        exit_points = [
            np.array(sbez[-1]) * 0.9,
            np.array(center_exit) * 5,
            np.array(center_exit) * 10,
        ]
        tool += path_sweep(
            regular_ngon(n=n_gon, d=diagonaled_width, spin=spin),
            beziers.bezpath_curve(
                exit_points, N=len(exit_points) - 1, splinesteps=self.resolution
            ),
        )

        bevels_tool = None
        z_anchor = max(host_z, 1.2)
        for hex_idx in FlowerLayout.RING_HEX_INDICES:
            verts = self.layout.ring_vertices(hex_idx)
            for i in range(len(verts)):
                p1, p2 = verts[i], verts[(i + 1) % len(verts)]
                line = ((p1[0], p1[1], z_anchor), (p2[0], p2[1], z_anchor))
                if bevels_tool is None:
                    bevels_tool = add_bevel(None, line, HEXAGON_BEVEL_SIZE, z_anchor)
                else:
                    bevels_tool += add_bevel(None, line, HEXAGON_BEVEL_SIZE, z_anchor)

        bevels_tool = bevels_tool & tool if bevels_tool is not None else None
        tool = tool.translateZ(
            path_width / 2 + host_z - indent_height + counter * 0.00001
        )
        subtractions: list[OpenSCADObject] = [tool]
        if bevels_tool is not None:
            subtractions.append(bevels_tool.translateZ(-indent_height))
        return subtractions

    def apply_features(
        self,
        flower_solid: OpenSCADObject,
        flower: FlowerDef,
        tileset: Tileset,
    ) -> OpenSCADObject:
        max_z = max(tileset.terrain_z(h.terrain) for h in flower.hexes.values())
        host_z = max(max_z, 2.0)
        subtractions: list[OpenSCADObject] = []
        counter = 0
        for entry_j, exit_j in flower.roads:
            counter += 1
            subtractions.extend(
                self._path_pair_tools(
                    entry_j,
                    exit_j,
                    n_gon=4,
                    spin=45,
                    width_scalar=STREET_WIDTH_SCALAR,
                    indent_height=STREET_INDENT_HEIGHT,
                    host_z=host_z,
                    counter=counter,
                )
            )
        for entry_j, exit_j in flower.water:
            counter += 1
            # Channel carved at water plane (one step below dominant ground terrain).
            water_host_z = tileset.water_z("ground") + tileset.meta.model_step
            subtractions.extend(
                self._path_pair_tools(
                    entry_j,
                    exit_j,
                    n_gon=6,
                    spin=60,
                    width_scalar=WATER_WIDTH_SCALAR,
                    indent_height=WATER_INDENT_HEIGHT,
                    host_z=water_host_z,
                    counter=counter,
                )
            )
        if not subtractions:
            return flower_solid
        return flower_solid - union()(subtractions)
