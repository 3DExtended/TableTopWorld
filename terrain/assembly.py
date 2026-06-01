"""Assemble multiple flowers from preview_map."""

from __future__ import annotations

from typing import TYPE_CHECKING

from solid2 import union

from terrain.catalog import EdgeProfileCatalog
from terrain.edges import EdgeGeometry
from terrain.features import FeatureCutters
from terrain.layout import axial_to_xy
from terrain.mesh import FlowerMeshBuilder
from terrain.render.planner import plan_flower
from terrain.render.spec import FlowerRenderSpec

if TYPE_CHECKING:
    from solid2 import OpenSCADObject

    from terrain.tileset import Tileset


class AssemblyExporter:
    """Place flowers per preview_map and union for combined export."""

    def __init__(self, tileset: Tileset, resolution: int = 100) -> None:
        self.tileset = tileset
        self.resolution = resolution
        self.layout = tileset.layout()
        self.mesh_builder = FlowerMeshBuilder(tileset, self.layout)
        self.edge_geom = EdgeGeometry(self.layout)
        self.catalog = EdgeProfileCatalog()
        self.features = FeatureCutters(self.layout, resolution)

    def describe_flower(self, flower_id: str) -> FlowerRenderSpec:
        if flower_id not in self.tileset.flowers:
            raise KeyError(f"unknown flower {flower_id!r}")
        return plan_flower(
            self.tileset,
            self.tileset.flowers[flower_id],
            resolution=self.resolution,
            catalog=self.catalog,
        )

    def build_flower(self, flower_id: str) -> OpenSCADObject:
        if flower_id not in self.tileset.flowers:
            raise KeyError(f"unknown flower {flower_id!r}")
        flower = self.tileset.flowers[flower_id]
        solid = self.mesh_builder.build_flower(flower)
        max_z = self.mesh_builder.max_flower_z(flower)
        solid = self.edge_geom.apply_bevels(solid, max_z)
        solid = self.features.apply_features(solid, flower, self.tileset)
        solid = self.edge_geom.apply_magnets(solid, flower, self.catalog)
        return solid

    def _place(self, solid: OpenSCADObject, q: int, r: int, rot: int) -> OpenSCADObject:
        spacing = self.layout.flower_center_spacing
        x, y = axial_to_xy(q, r, spacing)
        return solid.rotateZ(rot * 60).translateX(x).translateY(y)

    def build_preview(self) -> OpenSCADObject:
        parts = []
        for placement in self.tileset.preview_map:
            part = self.build_flower(placement.id)
            part = self._place(part, placement.at[0], placement.at[1], placement.rot)
            parts.append(part)
        if not parts:
            raise ValueError("preview_map is empty")
        return union()(parts)
