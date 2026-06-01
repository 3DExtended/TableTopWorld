"""Per-hex terrain mesh CSG."""

from __future__ import annotations

from typing import TYPE_CHECKING

from solid2 import cylinder, polygon, union

from terrain.constants import (
    BASE_PLATE_DEPTH,
    FLOWER_BOTTOM_Z,
    MAGNET_DEPTH,
    MAGNET_RADIUS,
    TOPPING_HEX_INDICES,
)
from terrain.layout import FlowerLayout

if TYPE_CHECKING:
    from solid2 import OpenSCADObject

    from terrain.tileset import FlowerDef, Tileset


class FlowerMeshBuilder:
    """Build 7-hex flower solids at terrain Z with standable tops and slopes."""

    def __init__(self, tileset: Tileset, layout: FlowerLayout | None = None) -> None:
        self.tileset = tileset
        self.layout = layout or tileset.layout()

    def terrain_z(self, level: str) -> float:
        return self.tileset.terrain_z(level)

    def hex_top_z(self, flower: FlowerDef, hex_idx: int) -> float:
        """Standable plateau Z before slope ramps (ground uses BASE_PLATE_DEPTH)."""
        z = self.terrain_z(flower.hexes[str(hex_idx)].terrain)
        if z <= 1e-9:
            return BASE_PLATE_DEPTH
        return z

    def slope_ramp_top_z(
        self, flower: FlowerDef, hex_idx: int, *, z_top: float | None = None
    ) -> float | None:
        """Upper Z of a slope wedge when the hex rises toward a taller neighbor."""
        hdef = flower.hexes[str(hex_idx)]
        if hdef.role != "slope":
            return None
        plateau = z_top if z_top is not None else self.hex_top_z(flower, hex_idx)
        neighbor_z = [
            self.hex_top_z(flower, n) for n in self.neighbor_indices(hex_idx)
        ]
        if not neighbor_z:
            return None
        target = max(neighbor_z)
        if target > plateau:
            return target
        return None

    def hex_mesh_top_z(self, flower: FlowerDef, hex_idx: int) -> float:
        """Highest solid Z on the hex, including slope ramps (used for bevel anchors)."""
        plateau = self.hex_top_z(flower, hex_idx)
        ramp_top = self.slope_ramp_top_z(flower, hex_idx, z_top=plateau)
        return ramp_top if ramp_top is not None else plateau

    def hex_height_at(self, flower: FlowerDef, hex_idx: int) -> float:
        return self.hex_mesh_top_z(flower, hex_idx)

    def build_hex_prism(self, hex_idx: int, height: float) -> OpenSCADObject:
        points = self.layout.cell_polygon(hex_idx)
        return polygon(points=points).linear_extrude(height=height)

    def neighbor_indices(self, hex_idx: int) -> list[int]:
        if hex_idx == 0:
            return list(range(1, 7))
        ring = hex_idx - 1
        return [0, ((ring + 5) % 6) + 1, ((ring + 1) % 6) + 1]

    def _hex_prism_to_top(self, hex_idx: int, z_top: float) -> OpenSCADObject:
        """Extrude from FLOWER_BOTTOM_Z (print bed at z=0)."""
        height = z_top - FLOWER_BOTTOM_Z
        return self.build_hex_prism(hex_idx, height).translateZ(FLOWER_BOTTOM_Z)

    def build_hex_solid(self, flower: FlowerDef, hex_idx: int) -> OpenSCADObject:
        hdef = flower.hexes[str(hex_idx)]
        z_top = self.hex_top_z(flower, hex_idx)
        solid = self._hex_prism_to_top(hex_idx, z_top)

        ramp_top = self.slope_ramp_top_z(flower, hex_idx, z_top=z_top)
        if ramp_top is not None:
            ramp_h = ramp_top - z_top
            ramp = self.build_hex_prism(hex_idx, ramp_h).translateZ(z_top)
            solid = solid + ramp
        return solid

    def subtract_topping_holes(
        self, solids: list[OpenSCADObject], flower: FlowerDef
    ) -> list[OpenSCADObject]:
        tools = []
        for hex_idx in flower.topping_hexes:
            if hex_idx not in range(7):
                continue
            z = self.hex_height_at(flower, hex_idx)
            tool = cylinder(h=MAGNET_DEPTH * 4, center=True, r=MAGNET_RADIUS).translateZ(z)
            center = self.layout.cell_center(hex_idx)
            tool = tool.translateX(center[0]).translateY(center[1])
            tools.append(tool)
        if not tools:
            return solids
        remove_tool = union()(tools)
        return [s - remove_tool for s in solids]

    def build_flower(self, flower: FlowerDef) -> OpenSCADObject:
        parts = [self.build_hex_solid(flower, i) for i in range(7)]
        parts = self.subtract_topping_holes(parts, flower)
        return union()(parts)

    def max_flower_z(self, flower: FlowerDef) -> float:
        return max(self.hex_height_at(flower, i) for i in range(7))
