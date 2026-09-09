"""Steps 4+5 of the 2026-09 goal spec: the terrain field (flat pads,
S-curve bands, organic relief) and roads/rivers with real width that
cross the silhouette at a side's middle-edge midpoint.

Every test measures the real field or the real mesh. The numbers are the
ones the approved mockups (docs/mockups) were drawn with: 16 mm roads
1 mm below the terrain, 22 mm rivers 2.5 mm deep with 5 mm banks, pads
covering 55% of a standable hex.
"""

from __future__ import annotations

import math
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial import cKDTree

from dataclasses import replace

from terrain.assembly import build_flower_mesh, pick_standable_hexes
from terrain.constants import BASE_PLATE_DEPTH_MM
from terrain.field import TerrainField
from terrain.layout import FlowerLayout
from terrain.magnets import DEFAULT_MAGNET_BORES
from terrain.standability import hex_cell_z_range
from terrain.tileset import load_tileset

ROOT = Path(__file__).resolve().parent.parent
DEFAULT = ROOT / "tilesets" / "default.yaml"


@pytest.fixture(scope="module")
def tileset():
    return load_tileset(DEFAULT)


def _field(tileset, flower_id: str, plateau=None, **param_overrides) -> TerrainField:
    fl = tileset.flowers[flower_id]
    if plateau is None:
        plateau = pick_standable_hexes(fl.seed, min_standable=tileset.meta.min_standable_hexes)
    return TerrainField(
        tileset.layout(),
        {h: fl.hexes[str(h)].height_level for h in range(FlowerLayout.HEX_CELL_COUNT)},
        plateau,
        fl.side_corner_heights,
        tileset.meta.heights.z,
        fl.seed,
        params=replace(TerrainField.__init__.__kwdefaults__["params"], **param_overrides),
        roads=[(r.entry, r.exit, r.via) for r in fl.roads],
        rivers=[(w.entry, w.exit, w.via) for w in fl.water],
    )


def _leg_point(layout, from_hex: int, to_hex: int, along_mm: float):
    """A point on the straight leg from one hex centre towards another,
    `along_mm` from the first centre, plus the leg's unit normal. Paths
    round their 60-degree bends with a 15 mm fillet whose tangents end
    8.7 mm from a centre, so anything past ~9 mm is on the straight part."""
    a = np.array(layout.cell_center(from_hex))
    b = np.array(layout.cell_center(to_hex))
    u = (b - a) / np.linalg.norm(b - a)
    return a + along_mm * u, np.array([-u[1], u[0]])


def _overused_edges(mesh) -> int:
    counts = Counter(map(tuple, mesh.edges_sorted))
    return sum(1 for c in counts.values() if c > 2)


@pytest.mark.parametrize("flower_id", ["crossroads", "river_bend"])
@pytest.mark.parametrize("n", [4, 8, 16])
def test_road_and_river_flowers_are_printable_solids(tileset, flower_id, n) -> None:
    mesh = build_flower_mesh(tileset, flower_id, subdivisions_per_edge=n)
    assert mesh.is_watertight
    assert mesh.is_winding_consistent
    assert mesh.volume > 0
    assert _overused_edges(mesh) == 0
    assert mesh.area_faces.min() > 1e-7


def test_shared_hex_edges_agree_from_both_sides(tileset) -> None:
    """Each hex cell is triangulated on its own and evaluates the field in
    its own frame, so the field must give the same height on a shared
    edge whichever of the two hexes asks - otherwise the seam is a
    sawtooth. A plateau/organic pair is the case that used to fail."""
    layout = tileset.layout()
    for flower_id in tileset.flowers:
        field = _field(tileset, flower_id)
        for i in range(FlowerLayout.HEX_CELL_COUNT):
            for k in range(FlowerLayout.EDGES_PER_HEX):
                j = field._edges[i][k].neighbour
                if j is None:
                    continue
                p1, p2 = layout.hex_edge_line(i, k)
                for t in np.linspace(0.0, 1.0, 25):
                    x, y = p1[0] + (p2[0] - p1[0]) * t, p1[1] + (p2[1] - p1[1]) * t
                    assert field.z(i, x, y) == pytest.approx(field.z(j, x, y), abs=1e-9), (
                        flower_id, i, j, t)


@pytest.mark.parametrize("flower_id", ["crossroads", "river_bend", "hill_peak"])
def test_field_has_no_steps(tileset, flower_id) -> None:
    """Sampling 12 radial lines per hex every 0.05 mm: the biggest jump
    between neighbouring samples must be far below any designed feature's
    height. The steepest designed feature is a road's 0.6 mm edge, which
    is at most ~0.5 mm per sample; a real discontinuity (the river bank
    used to be a 5 mm wall on cross-sloped ground) is several mm."""
    layout = tileset.layout()
    field = _field(tileset, flower_id)
    step = 0.05
    radii = np.arange(0.0, 14.0, step)
    worst = 0.0
    for h in range(FlowerLayout.HEX_CELL_COUNT):
        cx, cy = layout.cell_center(h)
        for angle in np.arange(0.0, 2 * math.pi, math.pi / 6):
            zs = np.array([field.z(h, cx + r * math.cos(angle), cy + r * math.sin(angle)) for r in radii])
            worst = max(worst, float(np.abs(np.diff(zs)).max()))
    assert worst < 1.0, worst


def test_plateau_pads_are_exactly_flat(tileset) -> None:
    layout = tileset.layout()
    field = _field(tileset, "hill_peak", plateau=range(FlowerLayout.HEX_CELL_COUNT))
    pad_radius = field.params.pad_min_r * field.apothem
    for h in range(FlowerLayout.HEX_CELL_COUNT):
        cx, cy = layout.cell_center(h)
        level = tileset.meta.heights.z(tileset.flowers["hill_peak"].hexes[str(h)].height_level)
        for r in np.linspace(0.0, 0.98 * pad_radius, 8):
            for angle in np.arange(0.0, 2 * math.pi, math.pi / 5):
                z = field.terrain_z(h, cx + r * math.cos(angle), cy + r * math.sin(angle))
                assert z == level, (h, r, angle, z)


def test_road_is_sixteen_wide_and_one_mm_below_a_flat_hex(tileset) -> None:
    """crossroads' road runs side 3 -> hex 4 -> hex 0 -> hex 1 -> side 0
    (the river crosses it at hex 0, so probe the leg inside hex 4). With
    every hex a plateau and the texture off, the pad of hex 4 (level 1,
    15 mm) is exactly flat, so the road bed must be exactly 14 mm across
    its full 16 mm width and the pad back at 15 mm beyond the 0.6 mm edge."""
    layout = tileset.layout()
    field = _field(tileset, "crossroads", plateau=range(7), noise_mm=0.0, organic_relief_mm=0.0, wobble=0.0)
    p = field.params
    assert p.road_half_width_mm * 2 == 16.0
    origin, normal = _leg_point(layout, 4, 0, 11.0)
    road = field.roads[0]
    assert road.dist_s(*origin)[0] < 1e-9
    # The bed is the terrain averaged over 12 mm along the road; hex 5 at
    # level 0 pulls the band beyond the pad down by a hair, so allow 1e-3.
    for d in np.linspace(-10.0, 10.0, 81):
        z = field.z_at(*(origin + d * normal))
        if abs(d) <= p.road_half_width_mm:
            assert z == pytest.approx(15.0 - p.road_depth_mm, abs=1e-3), d
        elif abs(d) >= p.road_half_width_mm + p.road_edge_mm:
            assert z == pytest.approx(15.0, abs=1e-9), d


def test_river_is_22_wide_and_keeps_its_roof_over_the_bores(tileset) -> None:
    """river_bend's river runs side 4 -> hex 5 -> hex 0 -> hex 2 -> side 1,
    through level-0 ground: the deepest case. On the leg between hexes 5
    and 0 (the road crosses at hex 0 itself, so probe 20 mm short of it)
    the flat floor must stop at river_bed_min_z, leaving at least 1 mm of
    material above the wall bores, and the 5 mm banks must rise
    monotonically to the terrain at 11 mm from the centreline."""
    layout = tileset.layout()
    field = _field(tileset, "river_bend", plateau=range(7), noise_mm=0.0, organic_relief_mm=0.0, wobble=0.0)
    p = field.params
    assert p.river_half_width_mm * 2 == 22.0
    origin, normal = _leg_point(layout, 5, 0, 25.0)
    river = field.rivers[0]
    assert river.dist_s(*origin)[0] < 1e-9
    bore_top = -BASE_PLATE_DEPTH_MM + DEFAULT_MAGNET_BORES.top_above_bed_mm
    profile = []
    for d in np.linspace(0.0, 12.0, 121):
        x, y = origin + d * normal
        assert field.roads[0].dist_s(x, y)[0] > p.road_half_width_mm + p.road_edge_mm
        z = field.z_at(x, y)
        profile.append(z)
        if d <= p.river_half_width_mm - p.river_bank_mm:
            assert z == pytest.approx(p.river_bed_min_z, abs=1e-9), d
            assert z - bore_top >= 1.0
        elif d >= p.river_half_width_mm:
            assert z == pytest.approx(0.0, abs=1e-9), d
    assert all(b >= a - 1e-9 for a, b in zip(profile, profile[1:])), "bank not monotone"


def test_river_banks_blend_from_the_local_terrain(tileset) -> None:
    """crossroads hex 3: the river runs along the foot of hex 2's plateau
    band, so the ground rises across the channel. Carving from the
    centreline bed left a 5 mm wall at the channel edge; the profile must
    now reach the terrain smoothly on both sides."""
    layout = tileset.layout()
    field = _field(tileset, "crossroads")
    p = field.params
    origin, normal = _leg_point(layout, 3, 0, 12.0)
    river = field.rivers[0]
    assert river.dist_s(*origin)[0] < 1e-9
    for sign in (-1.0, 1.0):
        prev = None
        for d in np.arange(0.0, 14.0, 0.05):
            x, y = origin + sign * d * normal
            z = field.z_at(x, y)
            if prev is not None:
                assert abs(z - prev) < 0.3, (sign, d, z, prev)
            prev = z
            if river.dist_s(x, y)[0] >= p.river_half_width_mm:
                assert z == pytest.approx(field.terrain_z_at(x, y), abs=1e-9)


def test_road_continues_across_the_crossroads_river_bend_seam(tileset) -> None:
    """The preview places river_bend one flower to the right of crossroads:
    crossroads' side 0 meets river_bend's side 3, both carrying the road.
    The two meshes must coincide vertex for vertex along that side, and
    the road bed must be at the same height (the declared crossing height
    minus the road depth) on both sides across the road's full width."""
    layout = tileset.layout()
    at = {p.id: p.at for p in tileset.preview_map}
    off_a = np.array(layout.flower_grid_to_xy(*at["crossroads"]))
    off_b = np.array(layout.flower_grid_to_xy(*at["river_bend"]))
    mesh_a = build_flower_mesh(tileset, "crossroads")
    mesh_b = build_flower_mesh(tileset, "river_bend")
    va = mesh_a.vertices + np.array([*off_a, 0.0])
    vb = mesh_b.vertices + np.array([*off_b, 0.0])
    corners_a = np.array(layout.side_corners(0)) + off_a
    corners_b = np.array(layout.side_corners(3)) + off_b
    assert np.allclose(corners_a, corners_b[::-1])

    def rim_on_side(mesh, v, corners):
        """The silhouette rim along the side: vertices shared by a vertical
        (wall) face and an upward (top) face, within 1.5 mm of the side's
        three edges (the walls are flat, the tolerance is slack). Not every top
        vertex near the edge: the hex-line strip's inset row runs 0.5 mm
        inside each flower and cannot coincide by design; and not every
        wall vertex: the bore collars' ray hits land where each bore's
        own mouth polygon sends them."""
        nz = mesh.face_normals[:, 2]
        wall_ids = set(np.unique(mesh.faces[np.abs(nz) < 1e-6]))
        top_ids = set(np.unique(mesh.faces[nz > 1e-6]))
        w = v[sorted(wall_ids & top_ids)]
        keep = np.zeros(len(w), dtype=bool)
        for p, q in zip(corners[:-1], corners[1:]):
            d = q - p
            length = np.linalg.norm(d)
            u = d / length
            rel = w[:, :2] - p
            t = rel @ u
            perp = np.abs(rel[:, 0] * u[1] - rel[:, 1] * u[0])
            keep |= (perp < 1.5) & (t > -1e-6) & (t < length + 1e-6)
        return w[keep]

    pa, pb = rim_on_side(mesh_a, va, corners_a), rim_on_side(mesh_b, vb, corners_b)
    assert len(pa) >= 3 * 8 + 1  # 3 edges x subdivisions, plus the far corner
    assert len(pa) == len(pb)
    dist, _ = cKDTree(pb).query(pa)
    assert dist.max() < 1e-8, dist.max()

    field_a = _field(tileset, "crossroads")
    field_b = _field(tileset, "river_bend")
    crossing = np.array(layout.junction_center(0)) + off_a
    edge_dir = corners_a[2] - corners_a[1]
    edge_dir /= np.linalg.norm(edge_dir)
    p = field_a.params
    for u in np.linspace(-p.road_half_width_mm, p.road_half_width_mm, 17):
        world = crossing + u * edge_dir
        za = field_a.z_at(*(world - off_a))
        zb = field_b.z_at(*(world - off_b))
        assert za == pytest.approx(15.0 - p.road_depth_mm, abs=1e-9), u
        assert zb == pytest.approx(za, abs=1e-9), u


def test_sockets_only_on_plateau_hexes_clear_of_roads_and_rivers(tileset) -> None:
    """crossroads seeds plateaus on hexes 0, 2 and 6; the road and the
    river both pass through 0 and the river through 6, so only hex 2 may
    carry a top socket: its centre is recessed by the socket depth, the
    other centres sit on their road/river bed."""
    layout = tileset.layout()
    field = _field(tileset, "crossroads")
    assert sorted(field.plateau_hexes) == [0, 2, 6]
    assert [h for h in range(7) if field.path_crosses_hex(h)] == [0, 1, 3, 4, 6]
    mesh = build_flower_mesh(tileset, "crossroads")
    bores = DEFAULT_MAGNET_BORES
    for h in (0, 2, 6):
        cx, cy = layout.cell_center(h)
        near = mesh.vertices[np.hypot(mesh.vertices[:, 0] - cx, mesh.vertices[:, 1] - cy) <= bores.radius_mm]
        near = near[near[:, 2] > -BASE_PLATE_DEPTH_MM + 1e-6]
        surface = field.z_at(cx, cy)
        if h == 2:
            assert near[:, 2].min() == pytest.approx(surface - bores.depth_mm, abs=1e-6)
        else:
            assert near[:, 2].min() >= surface - 1e-6


def test_plateau_hexes_without_a_path_are_standable(tileset) -> None:
    layout = tileset.layout()
    for flower_id in tileset.flowers:
        field = _field(tileset, flower_id)
        mesh = build_flower_mesh(tileset, flower_id)
        for h in field.plateau_hexes:
            if field.path_crosses_hex(h):
                continue
            assert hex_cell_z_range(mesh, layout, h) <= 1.0, (flower_id, h)
