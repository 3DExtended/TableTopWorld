"""Minimal tilesets for atomic and composed visual tests."""

from __future__ import annotations

from terrain.layout import FlowerLayout
from terrain.tileset import FlowerDef, HexDef, Tileset, TilesetMeta


def _flat_edges(layout: FlowerLayout) -> dict[str, str]:
    return {key: "flat_ground" for key in layout.exterior_edge_keys()}


def uniform_hexes(
    terrain: str = "ground",
    role: str = "standable",
) -> dict[str, HexDef]:
    return {str(i): HexDef(terrain=terrain, role=role) for i in range(7)}


def make_flower(
    flower_id: str,
    hexes: dict[str, HexDef],
    *,
    roads: list[tuple[int, int]] | None = None,
    water: list[tuple[int, int]] | None = None,
    topping_hexes: list[int] | None = None,
    edges: dict[str, str] | None = None,
) -> FlowerDef:
    layout = FlowerLayout()
    return FlowerDef(
        id=flower_id,
        hexes=hexes,
        edges=edges if edges is not None else _flat_edges(layout),
        roads=roads or [],
        water=water or [],
        topping_hexes=topping_hexes if topping_hexes is not None else [1, 2],
    )


def make_tileset(flower: FlowerDef) -> Tileset:
    return Tileset(
        meta=TilesetMeta(
            height_step_mm=20.0,
            scale=5.0,
            hex_outer_width=5.1961525,
        ),
        flowers={flower.id: flower},
        preview_map=[],
    )


def atom_all_ground() -> Tileset:
    return make_tileset(make_flower("atom_all_ground", uniform_hexes()))


def atom_height_step() -> Tileset:
    """Center middle, one high peak at hex 3, rest ground."""
    hexes = uniform_hexes()
    hexes["0"] = HexDef(terrain="middle", role="standable")
    hexes["3"] = HexDef(terrain="high", role="standable")
    return make_tileset(make_flower("atom_height_step", hexes))


def atom_slope_pad() -> Tileset:
    """Ground slope hex 5 rising toward middle neighbors on hex 0 and 4."""
    hexes = uniform_hexes()
    hexes["0"] = HexDef(terrain="middle", role="standable")
    hexes["4"] = HexDef(terrain="middle", role="standable")
    hexes["5"] = HexDef(terrain="ground", role="slope")
    return make_tileset(make_flower("atom_slope_pad", hexes))


def atom_road_only() -> Tileset:
    flower = make_flower(
        "atom_road_only",
        uniform_hexes(),
        roads=[(0, 4)],
    )
    return make_tileset(flower)


def atom_water_only() -> Tileset:
    flower = make_flower(
        "atom_water_only",
        uniform_hexes(),
        water=[(3, 5)],
    )
    return make_tileset(flower)
