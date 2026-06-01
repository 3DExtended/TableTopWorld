"""Build a FlowerRenderSpec from tileset data (no CSG)."""

from __future__ import annotations

from terrain.catalog import EdgeProfileCatalog
from terrain.constants import (
    FLOWER_BOTTOM_Z,
    HEXAGON_BEVEL_SIZE,
    MAGNET_CENTER_Z,
    MAGNET_DEPTH,
    MAGNET_RADIUS,
    STREET_INDENT_HEIGHT,
    STREET_WIDTH_SCALAR,
    WATER_INDENT_HEIGHT,
    WATER_WIDTH_SCALAR,
)
from terrain.edges import angle_with_x_axis
from terrain.layout import FlowerLayout
from terrain.mesh import FlowerMeshBuilder
from terrain.render.spec import (
    BevelSpec,
    ExteriorEdgeSpec,
    FlowerRenderSpec,
    HexCellSpec,
    MagnetHoleSpec,
    PathCutSpec,
    ToppingHoleSpec,
)
from terrain.tileset import FlowerDef, Tileset


def _round_xy(point: tuple[float, float], places: int = 6) -> tuple[float, float]:
    return (round(point[0], places), round(point[1], places))


def _round_polygon(
    vertices: list[tuple[float, float]], places: int = 6
) -> tuple[tuple[float, float], ...]:
    return tuple(_round_xy(v, places) for v in vertices)


def _magnet_center(line: tuple[tuple[float, float], tuple[float, float]]) -> tuple[float, float]:
    p0, p1 = line
    return (
        p1[0] - 0.5 * (p1[0] - p0[0]),
        p1[1] - 0.5 * (p1[1] - p0[1]),
    )


def _bevel_segment_count(layout: FlowerLayout) -> int:
    tool_settings: list[tuple[tuple[float, float, float], tuple[float, float, float]]] = []
    count = 0
    for hex_idx in FlowerLayout.RING_HEX_INDICES:
        verts = layout.ring_vertices(hex_idx)
        length = len(verts)
        for i in range(length):
            p1 = verts[i]
            p2 = verts[(i + 1) % length]
            bevel_settings = ((p1[0], p1[1], 0.0), (p2[0], p2[1], 0.0))
            if bevel_settings not in tool_settings:
                tool_settings.append(bevel_settings)
                count += 1
    return count


def _path_cut_centers(
    layout: FlowerLayout, entry_j: int, exit_j: int
) -> tuple[tuple[float, float], tuple[float, float]]:
    entry_a, entry_b, entry_c = layout.junction_lines(entry_j)
    entry_center = layout.center_of_three_lines(entry_a, entry_b, entry_c)
    exit_a, exit_b, exit_c = layout.junction_lines(exit_j)
    exit_center = layout.center_of_three_lines(exit_a, exit_b, exit_c)
    return _round_xy(entry_center), _round_xy(exit_center)


def plan_flower(
    tileset: Tileset,
    flower: FlowerDef,
    *,
    resolution: int = 100,
    catalog: EdgeProfileCatalog | None = None,
) -> FlowerRenderSpec:
    """Describe everything the SCAD pipeline would apply to this flower."""
    catalog = catalog or EdgeProfileCatalog()
    layout = tileset.layout()
    mesh = FlowerMeshBuilder(tileset, layout)

    hex_specs: list[HexCellSpec] = []
    for hex_idx in range(7):
        hdef = flower.hexes[str(hex_idx)]
        z_top = mesh.hex_top_z(flower, hex_idx)
        prism_height = z_top - FLOWER_BOTTOM_Z
        slope_ramp_height: float | None = None
        slope_ramp_base_z: float | None = None

        if hdef.role == "slope":
            neighbor_z = [
                mesh.hex_top_z(flower, n) for n in mesh.neighbor_indices(hex_idx)
            ]
            if neighbor_z:
                target = max(neighbor_z)
                if target > z_top:
                    slope_ramp_height = target - z_top
                    slope_ramp_base_z = z_top

        hex_specs.append(
            HexCellSpec(
                hex_idx=hex_idx,
                terrain=hdef.terrain,
                role=hdef.role,
                bottom_z=FLOWER_BOTTOM_Z,
                top_z=z_top,
                prism_height=prism_height,
                polygon_vertices=_round_polygon(layout.cell_polygon(hex_idx)),
                slope_ramp_height=slope_ramp_height,
                slope_ramp_base_z=slope_ramp_base_z,
            )
        )

    topping: list[ToppingHoleSpec] = []
    for hex_idx in flower.topping_hexes:
        if hex_idx not in range(7):
            continue
        cx, cy = layout.cell_center(hex_idx)
        z = mesh.hex_height_at(flower, hex_idx)
        topping.append(
            ToppingHoleSpec(
                hex_idx=hex_idx,
                center_x=round(cx, 6),
                center_y=round(cy, 6),
                center_z=z,
                radius=MAGNET_RADIUS,
                depth=MAGNET_DEPTH,
            )
        )

    exterior: list[ExteriorEdgeSpec] = []
    magnets: list[MagnetHoleSpec] = []
    for edge in layout.exterior_edges():
        profile = flower.edges[edge.key]
        mating = catalog.is_mating_profile(profile)
        exterior.append(
            ExteriorEdgeSpec(edge_key=edge.key, profile=profile, mating=mating)
        )
        if mating:
            cx, cy = _magnet_center(edge.line_2d)
            _, angle_deg = angle_with_x_axis(edge.line_2d)
            magnets.append(
                MagnetHoleSpec(
                    edge_key=edge.key,
                    center_x=round(cx, 6),
                    center_y=round(cy, 6),
                    center_z=MAGNET_CENTER_Z,
                    radius=MAGNET_RADIUS,
                    depth=MAGNET_DEPTH,
                    angle_deg=round(angle_deg, 6),
                )
            )

    max_z = mesh.max_flower_z(flower)
    bevel = BevelSpec(
        z_anchor=max(max_z, 1.2),
        size=HEXAGON_BEVEL_SIZE,
        segment_count=_bevel_segment_count(layout),
    )

    path_cuts: list[PathCutSpec] = []
    host_z = max_z
    for entry_j, exit_j in flower.roads:
        entry_c, exit_c = _path_cut_centers(layout, entry_j, exit_j)
        path_cuts.append(
            PathCutSpec(
                kind="road",
                entry_junction=entry_j,
                exit_junction=exit_j,
                host_z=host_z,
                indent_height=STREET_INDENT_HEIGHT,
                width_scalar=STREET_WIDTH_SCALAR,
                n_gon=4,
                spin=45.0,
                entry_center=entry_c,
                exit_center=exit_c,
            )
        )

    water_host_z = tileset.water_z("ground") + tileset.meta.model_step
    for entry_j, exit_j in flower.water:
        entry_c, exit_c = _path_cut_centers(layout, entry_j, exit_j)
        path_cuts.append(
            PathCutSpec(
                kind="water",
                entry_junction=entry_j,
                exit_junction=exit_j,
                host_z=water_host_z,
                indent_height=WATER_INDENT_HEIGHT,
                width_scalar=WATER_WIDTH_SCALAR,
                n_gon=6,
                spin=60.0,
                entry_center=entry_c,
                exit_center=exit_c,
            )
        )

    return FlowerRenderSpec(
        flower_id=flower.id,
        flower_bottom_z=FLOWER_BOTTOM_Z,
        max_z=max_z,
        hex_outer_width=layout.hex_outer_width,
        resolution=resolution,
        hexes=tuple(hex_specs),
        topping_holes=tuple(topping),
        exterior_edges=tuple(exterior),
        magnet_holes=tuple(magnets),
        bevel=bevel,
        path_cuts=tuple(path_cuts),
    )
