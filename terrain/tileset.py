"""Load and validate tileset YAML."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from terrain.catalog import EdgeProfileCatalog
from terrain.constants import ROLES, TERRAIN_LEVELS, TERRAIN_Z
from terrain.layout import FlowerLayout, axial_to_xy


class TilesetError(Exception):
    """Invalid tileset with human-readable message."""


@dataclass
class HexDef:
    terrain: str
    role: str


@dataclass
class FlowerDef:
    id: str
    hexes: dict[str, HexDef]
    edges: dict[str, str]  # edge_key -> profile name
    roads: list[tuple[int, int]]  # (entry_junction, exit_junction)
    water: list[tuple[int, int]]
    topping_hexes: list[int] = field(default_factory=lambda: [1, 2])


@dataclass
class PreviewPlacement:
    id: str
    at: tuple[int, int]
    rot: int


@dataclass
class TilesetMeta:
    height_step_mm: float
    scale: float
    hex_outer_width: float

    @property
    def model_step(self) -> float:
        return self.height_step_mm / self.scale


@dataclass
class Tileset:
    meta: TilesetMeta
    flowers: dict[str, FlowerDef]
    preview_map: list[PreviewPlacement]
    path: Path | None = None

    def layout(self) -> FlowerLayout:
        return FlowerLayout(self.meta.hex_outer_width)

    def terrain_z(self, level: str) -> float:
        if level not in TERRAIN_Z:
            raise TilesetError(f"unknown terrain level {level!r}")
        return TERRAIN_Z[level]

    def water_z(self, host_level: str) -> float:
        return self.terrain_z(host_level) - self.meta.model_step


def load_tileset(path: str | Path) -> Tileset:
    p = Path(path)
    raw = yaml.safe_load(p.read_text())
    if not isinstance(raw, dict):
        raise TilesetError(f"{p}: expected mapping at root")
    meta = _parse_meta(raw.get("meta") or {}, p)
    catalog = EdgeProfileCatalog()
    layout = FlowerLayout(meta.hex_outer_width)
    flowers_raw = raw.get("flowers") or {}
    if not flowers_raw:
        raise TilesetError(f"{p}: flowers section is required")
    flowers: dict[str, FlowerDef] = {}
    for fid, fdata in flowers_raw.items():
        flowers[fid] = _parse_flower(fid, fdata, layout, catalog, p)
    preview = _parse_preview_map(raw.get("preview_map") or [], flowers, p)
    tileset = Tileset(meta=meta, flowers=flowers, preview_map=preview, path=p)
    validate_tileset(tileset, catalog, layout)
    return tileset


def _parse_meta(data: dict[str, Any], path: Path) -> TilesetMeta:
    return TilesetMeta(
        height_step_mm=float(data.get("height_step_mm", 20)),
        scale=float(data.get("scale", 5)),
        hex_outer_width=float(data.get("hex_outer_width", 5.1961525)),
    )


def _parse_flower(
    fid: str,
    data: dict[str, Any],
    layout: FlowerLayout,
    catalog: EdgeProfileCatalog,
    path: Path,
) -> FlowerDef:
    if not isinstance(data, dict):
        raise TilesetError(f"{path}: flower {fid!r} must be a mapping")
    hexes_raw = data.get("hexes") or {}
    hexes: dict[str, HexDef] = {}
    for idx in range(7):
        key = str(idx)
        h = hexes_raw.get(key)
        if h is None:
            raise TilesetError(f"{path}: flower {fid!r} missing hex {key}")
        terrain = h.get("terrain")
        role = h.get("role")
        if terrain not in TERRAIN_LEVELS:
            raise TilesetError(
                f"{path}: flower {fid!r} hex {key}: terrain must be one of {TERRAIN_LEVELS}"
            )
        if role not in ROLES:
            raise TilesetError(
                f"{path}: flower {fid!r} hex {key}: role must be one of {ROLES}"
            )
        hexes[key] = HexDef(terrain=terrain, role=role)

    edges_raw = data.get("edges") or {}
    edges: dict[str, str] = {}
    valid_keys = set(layout.exterior_edge_keys())
    for ekey, edata in edges_raw.items():
        if ekey not in valid_keys:
            raise TilesetError(
                f"{path}: flower {fid!r} edge {ekey!r} is not a valid exterior side "
                f"(expected one of {sorted(valid_keys)})"
            )
        profile = edata.get("profile") if isinstance(edata, dict) else edata
        try:
            catalog.validate_profile_name(profile)
        except ValueError as exc:
            raise TilesetError(f"{path}: flower {fid!r} edge {ekey!r}: {exc}") from exc
        edges[ekey] = profile

    missing = valid_keys - set(edges.keys())
    if missing:
        raise TilesetError(
            f"{path}: flower {fid!r} missing edge definitions for {sorted(missing)}"
        )

    roads = _parse_junction_paths(data.get("roads") or [], "roads", fid, path)
    water = _parse_junction_paths(data.get("water") or [], "water", fid, path)
    topping = data.get("topping_hexes", [1, 2])
    return FlowerDef(
        id=fid,
        hexes=hexes,
        edges=edges,
        roads=roads,
        water=water,
        topping_hexes=list(topping),
    )


def _parse_junction_paths(
    items: list[Any], label: str, fid: str, path: Path
) -> list[tuple[int, int]]:
    result: list[tuple[int, int]] = []
    for item in items:
        if isinstance(item, dict):
            entry = item.get("entry_junction", item.get("entry"))
            exit_ = item.get("exit_junction", item.get("exit"))
        elif isinstance(item, (list, tuple)) and len(item) == 2:
            entry, exit_ = item[0], item[1]
        else:
            raise TilesetError(
                f"{path}: flower {fid!r} {label} entry must be "
                "{{entry_junction, exit_junction}} or [entry, exit]"
            )
        for j in (entry, exit_):
            if not isinstance(j, int) or j not in range(6):
                raise TilesetError(
                    f"{path}: flower {fid!r} {label} junction {j!r} must be integer 0..5"
                )
        result.append((int(entry), int(exit_)))
    return result


def _parse_preview_map(
    items: list[Any], flowers: dict[str, FlowerDef], path: Path
) -> list[PreviewPlacement]:
    placements: list[PreviewPlacement] = []
    for item in items:
        if not isinstance(item, dict):
            raise TilesetError(f"{path}: preview_map entries must be mappings")
        fid = item.get("id")
        if fid not in flowers:
            raise TilesetError(f"{path}: preview_map references unknown flower {fid!r}")
        at = item.get("at", [0, 0])
        rot = int(item.get("rot", 0))
        if rot not in range(6):
            raise TilesetError(f"{path}: preview_map rot for {fid!r} must be 0..5")
        placements.append(PreviewPlacement(id=fid, at=(int(at[0]), int(at[1])), rot=rot))
    return placements


def validate_tileset(
    tileset: Tileset,
    catalog: EdgeProfileCatalog | None = None,
    layout: FlowerLayout | None = None,
) -> None:
    catalog = catalog or EdgeProfileCatalog()
    layout = layout or tileset.layout()
    step = tileset.meta.model_step
    for i, level in enumerate(TERRAIN_LEVELS):
        expected = i * step
        if TERRAIN_Z[level] != expected:
            raise TilesetError(
                f"TERRAIN_Z[{level!r}] is {TERRAIN_Z[level]}, "
                f"expected {expected} (level_index * model_step)"
            )
    for fid, flower in tileset.flowers.items():
        standable = sum(1 for h in flower.hexes.values() if h.role == "standable")
        if standable < 1:
            raise TilesetError(
                f"flower {fid!r}: requires at least one standable hex, found {standable}"
            )
    _validate_preview_adjacency(tileset, catalog, layout)


def _validate_preview_adjacency(
    tileset: Tileset, catalog: EdgeProfileCatalog, layout: FlowerLayout
) -> None:
    if len(tileset.preview_map) < 2:
        return
    spacing = layout.flower_center_spacing
    tol = spacing * 0.15
    instances: list[tuple[str, str, tuple[float, float], tuple[float, float], str]] = []
    for placement in tileset.preview_map:
        flower = tileset.flowers[placement.id]
        origin_xy = axial_to_xy(placement.at[0], placement.at[1], spacing)
        for edge in layout.exterior_edges():
            mid, normal, _ = layout.transform_edge_to_world(
                edge, origin_xy, placement.rot
            )
            instances.append(
                (placement.id, edge.key, mid, normal, flower.edges[edge.key])
            )

    matched: set[int] = set()
    for i, a in enumerate(instances):
        for j, b in enumerate(instances):
            if j <= i or i in matched or j in matched:
                continue
            if a[0] == b[0]:
                continue
            dist = math.hypot(a[2][0] - b[2][0], a[2][1] - b[2][1])
            if dist > tol:
                continue
            dot = a[3][0] * b[3][0] + a[3][1] * b[3][1]
            if dot > -0.5:
                continue
            matched.add(i)
            matched.add(j)
            if not catalog.is_compatible(a[4], b[4]):
                raise TilesetError(
                    f"preview_map adjacency incompatible: {a[0]!r} edge {a[1]!r} "
                    f"({a[4]!r}) meets {b[0]!r} edge {b[1]!r} ({b[4]!r})"
                )
