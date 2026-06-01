"""Per-hex terrain mesh CSG."""

from __future__ import annotations

from typing import TYPE_CHECKING

from solid2 import cylinder, intersection, polygon, union

from terrain.constants import (
    FLOWER_BOTTOM_Z,
    MAGNET_DEPTH,
    MAGNET_RADIUS,
    MAGNET_WALL_TOP_Z,
    TOPPING_HEX_INDICES,
)
from terrain.edges import magnet_wall_slab
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
        """Standable / slope reference top in model units (terrain_z for the hex)."""
        return self.terrain_z(flower.hexes[str(hex_idx)].terrain)

    def hex_height_at(self, flower: FlowerDef, hex_idx: int) -> float:
        return self.hex_top_z(flower, hex_idx)

    def build_hex_prism(self, hex_idx: int, height: float) -> OpenSCADObject:
        points = self.layout.cell_polygon(hex_idx)
        return polygon(points=points).linear_extrude(height=height)

    def neighbor_indices(self, hex_idx: int) -> list[int]:
        if hex_idx == 0:
            return list(range(1, 7))
        ring = hex_idx - 1
        return [0, ((ring + 5) % 6) + 1, ((ring + 1) % 6) + 1]

    def _hex_prism_to_top(self, hex_idx: int, z_top: float) -> OpenSCADObject:
        """Extrude from FLOWER_BOTTOM_Z so every hex shares one print-bed plane."""
        height = z_top - FLOWER_BOTTOM_Z
        return self.build_hex_prism(hex_idx, height).translateZ(FLOWER_BOTTOM_Z)

    def _exterior_magnet_collar(self, hex_idx: int, z_top: float) -> OpenSCADObject | None:
        """Extra exterior rim inside the hex footprint up to magnet height (ground tiles)."""
        if hex_idx not in FlowerLayout.RING_HEX_INDICES:
            return None
        if z_top >= MAGNET_WALL_TOP_Z - 1e-9:
            return None
        verts = self.layout.ring_vertices(hex_idx)
        center = self.layout.cell_center(hex_idx)
        extended = self._hex_prism_to_top(hex_idx, MAGNET_WALL_TOP_Z)
        bands: list[OpenSCADObject] = []
        for vi in self.layout.exterior_vertex_indices(hex_idx):
            line = (verts[vi], verts[(vi + 1) % 6])
            bands.append(
                magnet_wall_slab(
                    line,
                    FLOWER_BOTTOM_Z,
                    MAGNET_WALL_TOP_Z,
                    center,
                )
            )
        return intersection()(extended, union()(*bands))

    def build_hex_solid(self, flower: FlowerDef, hex_idx: int) -> OpenSCADObject:
        hdef = flower.hexes[str(hex_idx)]
        z_top = self.hex_top_z(flower, hex_idx)
        solid = self._hex_prism_to_top(hex_idx, z_top)

        if hdef.role == "slope":
            neighbor_z = [
                self.hex_top_z(flower, n) for n in self.neighbor_indices(hex_idx)
            ]
            if neighbor_z:
                target = max(neighbor_z)
                if target > z_top:
                    ramp_h = target - z_top
                    ramp = self.build_hex_prism(hex_idx, ramp_h).translateZ(z_top)
                    solid = solid + ramp
        collar = self._exterior_magnet_collar(hex_idx, z_top)
        if collar is not None:
            solid = solid + collar
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
