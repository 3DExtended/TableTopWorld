"""Output backends for flower rendering."""

from __future__ import annotations

from enum import Enum
from typing import TYPE_CHECKING, Protocol

if TYPE_CHECKING:
    from solid2 import OpenSCADObject

    from terrain.render.spec import FlowerRenderSpec
    from terrain.tileset import Tileset


class FlowerOutputFormat(str, Enum):
    SCAD = "scad"
    SPEC = "spec"


class FlowerOutputBackend(Protocol):
    """Abstraction for flower render output (mesh CSG or declarative spec)."""

    def render_flower(self, tileset: Tileset, flower_id: str) -> object:
        """Return backend-specific result (OpenSCADObject or FlowerRenderSpec)."""
        ...


class SpecOutputBackend:
    def __init__(self, resolution: int = 100) -> None:
        self.resolution = resolution

    def render_flower(self, tileset: Tileset, flower_id: str) -> FlowerRenderSpec:
        from terrain.render.planner import plan_flower

        if flower_id not in tileset.flowers:
            raise KeyError(f"unknown flower {flower_id!r}")
        return plan_flower(
            tileset,
            tileset.flowers[flower_id],
            resolution=self.resolution,
        )


class ScadOutputBackend:
    def __init__(self, tileset: Tileset, resolution: int = 100) -> None:
        from terrain.assembly import AssemblyExporter

        self._exporter = AssemblyExporter(tileset, resolution=resolution)

    def render_flower(self, tileset: Tileset, flower_id: str) -> OpenSCADObject:
        return self._exporter.build_flower(flower_id)
