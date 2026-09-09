"""Load and validate tileset YAML - explicit-mesh schema.

Replaces the old named-edge-profile-catalog schema: a flower's boundary
is no longer authored as a profile name per exterior edge, but as a
discrete height level per hex cell plus a 4-corner height sequence per
side (design decisions #2-#6). Two flowers may only share a side if the
numbers match according to the reversed-declaration contract (a side k
of one flower always meets side (k+3)%6 of its neighbor, in reversed
corner order - verified in terrain/layout.py and terrain/heightfield.py).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from terrain.heights import HeightLevels
from terrain.layout import FlowerLayout


class TilesetError(Exception):
    """Invalid tileset with human-readable message."""


@dataclass
class HexDef:
    height_level: int


@dataclass(frozen=True)
class PathDecl:
    """A road or river: enters at side `entry`'s middle edge, leaves at
    side `exit`'s, optionally routed through the hex cells in `via`
    (default: the two ring hexes owning those edges, through the centre
    hex when they are not adjacent). See terrain/field.py."""

    entry: int
    exit: int
    via: tuple[int, ...] = ()


@dataclass
class FlowerDef:
    id: str
    hexes: dict[str, HexDef]
    side_corner_heights: dict[int, tuple[int, int, int, int]]
    seed: int
    roads: list[PathDecl] = field(default_factory=list)
    water: list[PathDecl] = field(default_factory=list)


@dataclass
class PreviewPlacement:
    id: str
    at: tuple[int, int]
    rot: int


@dataclass
class TilesetMeta:
    height_step_mm: float
    level_count: int
    scale: float
    hex_outer_width: float
    min_standable_hexes: int

    @property
    def heights(self) -> HeightLevels:
        return HeightLevels(mm_per_level=self.height_step_mm, level_count=self.level_count)


@dataclass
class Tileset:
    meta: TilesetMeta
    flowers: dict[str, FlowerDef]
    preview_map: list[PreviewPlacement]
    path: Path | None = None

    def layout(self) -> FlowerLayout:
        """hex_outer_width x scale, not hex_outer_width alone - decision
        #12: DEFAULT_HEX_OUTER_WIDTH is a small abstract model unit, only
        physically sized once multiplied by meta.scale. Height (via
        meta.heights.z(), height_step_mm) is already authored directly in
        real mm and needs no such multiplier - only the XY footprint does.
        `scale` was previously parsed and validated by _parse_meta() but
        never actually consumed anywhere in the pipeline (the same class
        of gap as HexDef.height_level, flagged separately) - without this,
        a flower's footprint (~20 model units across) would be smaller
        than its own possible height variation (up to height_step_mm x
        (level_count-1), in real mm), an implausible sliver of a shape.
        """
        return FlowerLayout(self.meta.hex_outer_width * self.meta.scale)


def load_tileset(path: str | Path) -> Tileset:
    p = Path(path)
    raw = yaml.safe_load(p.read_text())
    if not isinstance(raw, dict):
        raise TilesetError(f"{p}: expected mapping at root")
    meta = _parse_meta(raw.get("meta") or {})
    layout = FlowerLayout(meta.hex_outer_width)
    flowers_raw = raw.get("flowers") or {}
    if not flowers_raw:
        raise TilesetError(f"{p}: flowers section is required")
    flowers: dict[str, FlowerDef] = {}
    for fid, fdata in flowers_raw.items():
        flowers[fid] = _parse_flower(fid, fdata, meta.heights, p)
    preview = _parse_preview_map(raw.get("preview_map") or [], flowers, p)
    tileset = Tileset(meta=meta, flowers=flowers, preview_map=preview, path=p)
    validate_tileset(tileset, layout)
    return tileset


def _parse_meta(data: dict[str, Any]) -> TilesetMeta:
    return TilesetMeta(
        height_step_mm=float(data.get("height_step_mm", 15.0)),
        level_count=int(data.get("level_count", 4)),
        scale=float(data.get("scale", 5)),
        hex_outer_width=float(data.get("hex_outer_width", 5.1961525)),
        min_standable_hexes=int(data.get("min_standable_hexes", 1)),
    )


def _parse_flower(
    fid: str, data: dict[str, Any], heights: HeightLevels, path: Path
) -> FlowerDef:
    if not isinstance(data, dict):
        raise TilesetError(f"{path}: flower {fid!r} must be a mapping")

    hexes_raw = data.get("hexes") or {}
    hexes: dict[str, HexDef] = {}
    for idx in range(FlowerLayout.HEX_CELL_COUNT):
        key = str(idx)
        h = hexes_raw.get(key)
        if h is None:
            raise TilesetError(f"{path}: flower {fid!r} missing hex {key}")
        level = h.get("height_level") if isinstance(h, dict) else h
        _validate_level(level, heights, f"flower {fid!r} hex {key}", path)
        hexes[key] = HexDef(height_level=int(level))

    sides_raw = data.get("side_corner_heights") or {}
    side_corner_heights: dict[int, tuple[int, int, int, int]] = {}
    for side_idx in range(FlowerLayout.SIDE_COUNT):
        values = sides_raw.get(side_idx, sides_raw.get(str(side_idx)))
        if values is None:
            raise TilesetError(
                f"{path}: flower {fid!r} missing side_corner_heights for side {side_idx}"
            )
        if not isinstance(values, (list, tuple)) or len(values) != 4:
            raise TilesetError(
                f"{path}: flower {fid!r} side {side_idx} must have exactly 4 "
                f"corner heights, got {values!r}"
            )
        for v in values:
            _validate_level(v, heights, f"flower {fid!r} side {side_idx}", path)
        side_corner_heights[side_idx] = tuple(int(v) for v in values)

    for side_idx in range(FlowerLayout.SIDE_COUNT):
        this_side = side_corner_heights[side_idx]
        next_side = side_corner_heights[(side_idx + 1) % FlowerLayout.SIDE_COUNT]
        if this_side[3] != next_side[0]:
            raise TilesetError(
                f"{path}: flower {fid!r} side {side_idx}'s last corner "
                f"({this_side[3]}) must equal side {(side_idx + 1) % 6}'s "
                f"first corner ({next_side[0]}) - they are the same "
                "physical junction"
            )

    seed = data.get("seed")
    if not isinstance(seed, int):
        raise TilesetError(f"{path}: flower {fid!r} requires an integer seed")

    roads = _parse_junction_paths(data.get("roads") or [], "roads", fid, path)
    water = _parse_junction_paths(data.get("water") or [], "water", fid, path)
    return FlowerDef(
        id=fid,
        hexes=hexes,
        side_corner_heights=side_corner_heights,
        seed=seed,
        roads=roads,
        water=water,
    )


def _validate_level(value: Any, heights: HeightLevels, label: str, path: Path) -> None:
    if not isinstance(value, int) or value not in range(heights.level_count):
        raise TilesetError(
            f"{path}: {label}: height level must be an integer 0..{heights.level_count - 1}, "
            f"got {value!r}"
        )


def _parse_junction_paths(
    items: list[Any], label: str, fid: str, path: Path
) -> list[PathDecl]:
    result: list[PathDecl] = []
    for item in items:
        via: tuple[int, ...] = ()
        if isinstance(item, dict):
            entry = item.get("entry_junction", item.get("entry"))
            exit_ = item.get("exit_junction", item.get("exit"))
            raw_via = item.get("via") or []
            if not isinstance(raw_via, list) or not all(
                isinstance(h, int) and h in range(FlowerLayout.HEX_CELL_COUNT) for h in raw_via
            ):
                raise TilesetError(
                    f"{path}: flower {fid!r}: {label} 'via' must be a list of hex indices 0..6"
                )
            via = tuple(raw_via)
        elif isinstance(item, (list, tuple)) and len(item) == 2:
            entry, exit_ = item[0], item[1]
        else:
            raise TilesetError(
                f"{path}: flower {fid!r}: {label} entries must be [entry, exit] or "
                "{entry, exit, via}"
            )
        for value in (entry, exit_):
            if not isinstance(value, int) or value not in range(FlowerLayout.SIDE_COUNT):
                raise TilesetError(
                    f"{path}: flower {fid!r}: {label} sides must be integers 0..5, got {value!r}"
                )
        result.append(PathDecl(entry=entry, exit=exit_, via=via))
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


# Flower-grid axial direction deltas, in the same order as
# FlowerLayout.neighbor_flower_offset()'s side indices 0..5 - verified
# numerically (see terrain/layout.py::flower_grid_to_xy's docstring) that
# these are exactly the (q, r) unit steps reproducing each direction.
_NEIGHBOR_GRID_DELTAS: tuple[tuple[int, int], ...] = (
    (1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1),
)


MAX_ROAD_STEP_LEVELS = 1


def path_stations(flower: FlowerDef, path: PathDecl) -> list[tuple[str, float]]:
    """(name, level) of everything a road or river passes, in order: the
    entry crossing (the mean of the middle edge's two corner levels), each
    hex of its route, the exit crossing."""

    def crossing(side: int) -> float:
        corners = flower.side_corner_heights[side]
        return 0.5 * (corners[1] + corners[2])

    stations = [(f"side {path.entry}", crossing(path.entry))]
    for h in FlowerLayout.path_route(path.entry, path.exit, path.via):
        stations.append((f"hex {h}", float(flower.hexes[str(h)].height_level)))
    stations.append((f"side {path.exit}", crossing(path.exit)))
    return stations


def _check_road_grade(flower: FlowerDef, road: PathDecl) -> None:
    """A road may change at most MAX_ROAD_STEP_LEVELS per hex: its bed
    follows the terrain, so a two-level step is a 30 mm ramp inside one
    blend band, a wall no miniature walks up. Rivers are not held to this."""
    stations = path_stations(flower, road)
    for (a, la), (b, lb) in zip(stations, stations[1:]):
        if abs(lb - la) > MAX_ROAD_STEP_LEVELS + 1e-9:
            raise TilesetError(
                f"flower {flower.id!r}: road from side {road.entry} to side {road.exit} changes "
                f"{abs(lb - la):g} levels between {a} (level {la:g}) and {b} (level {lb:g}); "
                f"a road may change at most {MAX_ROAD_STEP_LEVELS} level per hex"
            )


def validate_tileset(tileset: Tileset, layout: FlowerLayout | None = None) -> None:
    """Structural validation, plus preview_map cross-flower side matching
    (decision #3/#4's reversed-declaration contract) now that
    flower_grid_to_xy() gives a correct placement convention to check
    adjacency against - the old axial_to_xy()/flower_center_spacing bug
    that used to block this is fixed (see terrain/layout.py).

    Deliberately NOT checked here: standability (decision #11 - checked
    post-generation against the real mesh, see terrain/standability.py).

    The cross-flower check only applies to placements with rot=0 on both
    sides: a nonzero `rot` visually rotates a flower in the preview but
    this validation doesn't (yet) remap side indices under rotation, so a
    rotated flower's side-matching correctness against its neighbors is
    simply not verified rather than checked against the wrong side index.
    """
    layout = layout or tileset.layout()
    if tileset.meta.min_standable_hexes < 0:
        raise TilesetError("meta.min_standable_hexes must be >= 0")
    for flower in tileset.flowers.values():
        for road_or_water in (*flower.roads, *flower.water):
            if road_or_water.entry == road_or_water.exit:
                raise TilesetError(
                    f"flower {flower.id!r}: road/water entry and exit junction "
                    f"must differ, got {road_or_water.entry}"
                )
            if (road_or_water.exit - road_or_water.entry) % FlowerLayout.SIDE_COUNT in (1, 5):
                # A path enters along the middle edge's normal and must leave
                # along the exit side's, which for adjacent sides differ by
                # 120 degrees whatever hexes it visits in between - a bend
                # tighter than a road or river is wide (terrain/field.py
                # rounds bends with a 15 mm fillet; the sharpest it can round
                # inside a 22.5 mm leg is 60 degrees).
                raise TilesetError(
                    f"flower {flower.id!r}: road/water cannot leave through side "
                    f"{road_or_water.exit}, adjacent to its entry side {road_or_water.entry} "
                    "(a 120-degree bend); use the opposite side or the two next to it"
                )
            if any(road_or_water is road for road in flower.roads):
                _check_road_grade(flower, road_or_water)

    by_position = {p.at: p for p in tileset.preview_map}
    for placement in tileset.preview_map:
        if placement.rot != 0:
            continue
        flower = tileset.flowers[placement.id]
        q, r = placement.at
        for k, (dq, dr) in enumerate(_NEIGHBOR_GRID_DELTAS):
            neighbor_placement = by_position.get((q + dq, r + dr))
            if neighbor_placement is None or neighbor_placement.rot != 0:
                continue
            neighbor = tileset.flowers[neighbor_placement.id]
            this_side = flower.side_corner_heights[k]
            neighbor_side = neighbor.side_corner_heights[(k + 3) % FlowerLayout.SIDE_COUNT]
            if this_side != tuple(reversed(neighbor_side)):
                raise TilesetError(
                    f"preview_map: flower {flower.id!r} at {placement.at} side "
                    f"{k} ({this_side}) does not match flower "
                    f"{neighbor.id!r} at {neighbor_placement.at} side "
                    f"{(k + 3) % FlowerLayout.SIDE_COUNT} reversed "
                    f"({tuple(reversed(neighbor_side))}) - these are the same "
                    "physical shared boundary"
                )
            # a road/river leaving through this side must continue next door
            for label, mine, theirs in (
                ("road", flower.roads, neighbor.roads),
                ("river", flower.water, neighbor.water),
            ):
                here = any(k in (p.entry, p.exit) for p in mine)
                there = any((k + 3) % FlowerLayout.SIDE_COUNT in (p.entry, p.exit) for p in theirs)
                if here != there:
                    raise TilesetError(
                        f"preview_map: flower {flower.id!r} at {placement.at} side {k} "
                        f"{'has' if here else 'has no'} {label} but flower "
                        f"{neighbor.id!r} at {neighbor_placement.at} side "
                        f"{(k + 3) % FlowerLayout.SIDE_COUNT} {'has none' if here else 'has one'}"
                    )
