"""Per-hex terrain mesh CSG."""

from __future__ import annotations

from typing import TYPE_CHECKING

from solid2 import cylinder, polygon, union

from terrain.constants import (
    LEGACY_HEXAGON_HEIGHT,
    MAGNET_DEPTH,
    MAGNET_RADIUS,
    TERRAIN_Z,
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

    def hex_height_at(self, flower: FlowerDef, hex_idx: int) -> float:
        hdef = flower.hexes[str(hex_idx)]
        z_top = self.terrain_z(hdef.terrain)
        if hdef.role == "standable":
            return max(z_top, LEGACY_HEXAGON_HEIGHT)
        if hdef.role == "slope":
            return max(z_top, LEGACY_HEXAGON_HEIGHT)
        if hdef.role in ("road_channel", "water"):
            return max(z_top, LEGACY_HEXAGON_HEIGHT)
        return max(z_top, LEGACY_HEXAGON_HEIGHT)

    def build_hex_prism(self, hex_idx: int, height: float) -> OpenSCADObject:
        points = self.layout.cell_polygon(hex_idx)
        return polygon(points=points).linear_extrude(height=height)

    def _neighbor_indices(self, hex_idx: int) -> list[int]:
        if hex_idx == 0:
            return list(range(1, 7))
        ring = hex_idx - 1
        neighbors = [0, ((ring + 5) % 6) + 1, ((ring + 1) % 6) + 1]
        return neighbors

    def build_hex_solid(self, flower: FlowerDef, hex_idx: int) -> OpenSCADObject:
        hdef = flower.hexes[str(hex_idx)]
        z_top = self.terrain_z(hdef.terrain)
        base_z = 0.0
        height = max(z_top - base_z, LEGACY_HEXAGON_HEIGHT)
        solid = self.build_hex_prism(hex_idx, height)

        if hdef.role == "slope":
            neighbor_z = [
                self.terrain_z(flower.hexes[str(n)].terrain)
                for n in self._neighbor_indices(hex_idx)
            ]
            if neighbor_z:
                target = max(neighbor_z)
                if target > z_top:
                    ramp_h = target - z_top
                    ramp = (
                        self.build_hex_prism(hex_idx, ramp_h)
                        .translateZ(z_top)
                    )
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
