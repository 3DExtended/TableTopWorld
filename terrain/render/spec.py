"""Structured render specification for a single flower."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(frozen=True)
class HexCellSpec:
    hex_idx: int
    terrain: str
    role: str
    bottom_z: float
    top_z: float
    prism_height: float
    polygon_vertices: tuple[tuple[float, float], ...]
    slope_ramp_height: float | None = None
    slope_ramp_base_z: float | None = None


@dataclass(frozen=True)
class ToppingHoleSpec:
    hex_idx: int
    center_x: float
    center_y: float
    center_z: float
    radius: float
    depth: float


@dataclass(frozen=True)
class ExteriorEdgeSpec:
    edge_key: str
    profile: str
    mating: bool


@dataclass(frozen=True)
class MagnetHoleSpec:
    edge_key: str
    center_x: float
    center_y: float
    center_z: float
    radius: float
    depth: float
    angle_deg: float


@dataclass(frozen=True)
class BevelSpec:
    z_anchor: float
    size: float
    segment_count: int


@dataclass(frozen=True)
class PathCutSpec:
    kind: str
    entry_junction: int
    exit_junction: int
    host_z: float
    indent_height: float
    width_scalar: float
    n_gon: int
    spin: float
    entry_center: tuple[float, float]
    exit_center: tuple[float, float]


@dataclass(frozen=True)
class FlowerRenderSpec:
    """Complete declarative description of one flower render pass."""

    flower_id: str
    flower_bottom_z: float
    max_z: float
    hex_outer_width: float
    resolution: int
    hexes: tuple[HexCellSpec, ...]
    topping_holes: tuple[ToppingHoleSpec, ...]
    exterior_edges: tuple[ExteriorEdgeSpec, ...]
    magnet_holes: tuple[MagnetHoleSpec, ...]
    bevel: BevelSpec
    path_cuts: tuple[PathCutSpec, ...] = field(default_factory=tuple)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
