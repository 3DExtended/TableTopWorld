"""The generated landscapes (tilesets/landscape.yaml on four levels,
tilesets/hills.yaml on three): 19 flowers on a hexagon of the flower grid,
every number derived from one elevation field by scripts/gen_landscape.py,
with a river across six flowers and a road across five that ford.

Built coarse (4 subdivisions per edge) so all 19 flowers build in about a
second; the checks are about seams and paths, not surface detail.
"""

from __future__ import annotations

import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial import cKDTree

from terrain.assembly import build_flower_mesh, build_preview_mesh, pick_standable_hexes
from terrain.field import TerrainField
from terrain.layout import FlowerLayout
from terrain.tileset import _NEIGHBOR_GRID_DELTAS, TilesetError, load_tileset

ROOT = Path(__file__).resolve().parent.parent
GENERATOR = ROOT / "scripts" / "gen_landscape.py"
PRESETS = ("landscape", "hills")
N = 4


@pytest.fixture(scope="module", params=PRESETS)
def tileset_path(request) -> Path:
    return ROOT / "tilesets" / f"{request.param}.yaml"


@pytest.fixture(scope="module")
def tileset(tileset_path):
    return load_tileset(tileset_path)


@pytest.fixture(scope="module")
def meshes(tileset):
    return {p.id: build_flower_mesh(tileset, p.id, subdivisions_per_edge=N) for p in tileset.preview_map}


def _field(tileset, flower_id: str) -> TerrainField:
    fl = tileset.flowers[flower_id]
    plateau = pick_standable_hexes(fl.seed, min_standable=tileset.meta.min_standable_hexes)
    return TerrainField(
        tileset.layout(),
        {h: fl.hexes[str(h)].height_level for h in range(FlowerLayout.HEX_CELL_COUNT)},
        plateau,
        fl.side_corner_heights,
        tileset.meta.heights.z,
        fl.seed,
        roads=[(r.entry, r.exit, r.via) for r in fl.roads],
        rivers=[(w.entry, w.exit, w.via) for w in fl.water],
    )


def _rim(mesh, layout, side: int, offset) -> np.ndarray:
    """World-space silhouette rim vertices along one side: vertices shared by
    a vertical (wall) face and an upward (top) face within 1.5 mm of the
    side's three edges."""
    v = mesh.vertices + np.array([offset[0], offset[1], 0.0])
    nz = mesh.face_normals[:, 2]
    wall_ids = set(np.unique(mesh.faces[np.abs(nz) < 1e-6]))
    top_ids = set(np.unique(mesh.faces[nz > 1e-6]))
    w = v[sorted(wall_ids & top_ids)]
    corners = np.array(layout.side_corners(side)) + np.array(offset)
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


def test_generator_reproduces_the_committed_tileset(tmp_path, tileset_path) -> None:
    out = tmp_path / tileset_path.name
    subprocess.run([sys.executable, str(GENERATOR), "--preset", tileset_path.stem, "--output", str(out)], check=True, cwd=ROOT)
    assert out.read_text() == tileset_path.read_text()


def test_nineteen_flowers_fill_a_hexagon_of_radius_two(tileset) -> None:
    cells = {p.at for p in tileset.preview_map}
    assert len(tileset.preview_map) == len(cells) == 19
    assert all(max(abs(q), abs(r), abs(q + r)) <= 2 for q, r in cells)
    assert {p.id for p in tileset.preview_map} == set(tileset.flowers)
    levels = {fl.hexes[str(h)].height_level for fl in tileset.flowers.values() for h in range(7)}
    assert levels == set(range(tileset.meta.level_count)), "a landscape should use every level it declares"


def test_every_flower_is_a_printable_solid(meshes) -> None:
    for flower_id, mesh in meshes.items():
        assert mesh.is_watertight, flower_id
        assert mesh.is_winding_consistent, flower_id
        assert mesh.volume > 0, flower_id


def test_every_shared_side_coincides_vertex_for_vertex(tileset, meshes) -> None:
    """42 seams in a radius-2 hexagon of flowers; on each, the two rims
    must be the same points (the boundary contract, decision #4, plus the
    road/river crossing rule and the corner-jitter window)."""
    layout = tileset.layout()
    by_position = {p.at: p.id for p in tileset.preview_map}
    seams = 0
    for (q, r), flower_id in by_position.items():
        for k, (dq, dr) in enumerate(_NEIGHBOR_GRID_DELTAS[:3]):  # each seam once
            neighbour_id = by_position.get((q + dq, r + dr))
            if neighbour_id is None:
                continue
            rim_a = _rim(meshes[flower_id], layout, k, layout.flower_grid_to_xy(q, r))
            rim_b = _rim(meshes[neighbour_id], layout, (k + 3) % 6, layout.flower_grid_to_xy(q + dq, r + dr))
            assert len(rim_a) == len(rim_b) >= 3 * N + 1, (flower_id, neighbour_id)
            dist, _ = cKDTree(rim_b).query(rim_a)
            assert dist.max() < 1e-8, (flower_id, neighbour_id, dist.max())
            seams += 1
    assert seams == 42


def test_river_and_road_continue_across_their_seams(tileset) -> None:
    """Wherever a road or river crosses a seam, both flowers' fields give
    exactly the same height across the flat core (a road's full bed, a
    river's floor) along the shared middle edge, agree within the seeded
    texture's amplitude over the banks and the road edge (the silhouette
    vertices themselves are the shared contour - see the rim test), and
    the river is lower than the ground beside it."""
    layout = tileset.layout()
    by_position = {p.at: p.id for p in tileset.preview_map}
    fields = {fid: _field(tileset, fid) for fid in tileset.flowers}
    crossings = 0
    for (q, r), flower_id in by_position.items():
        flower = tileset.flowers[flower_id]
        origin = np.array(layout.flower_grid_to_xy(q, r))
        for k, (dq, dr) in enumerate(_NEIGHBOR_GRID_DELTAS[:3]):
            neighbour_id = by_position.get((q + dq, r + dr))
            if neighbour_id is None:
                continue
            n_origin = np.array(layout.flower_grid_to_xy(q + dq, r + dr))
            for kind, paths in (("road", flower.roads), ("river", flower.water)):
                if not any(k in (p.entry, p.exit) for p in paths):
                    continue
                fa, fb = fields[flower_id], fields[neighbour_id]
                params = fa.params
                half = params.road_half_width_mm if kind == "road" else params.river_half_width_mm
                core = half if kind == "road" else half - params.river_bank_mm
                crossing = np.array(layout.junction_center(k)) + origin
                corners = np.array(layout.side_corners(k)) + origin
                edge_dir = corners[2] - corners[1]
                edge_len = np.linalg.norm(edge_dir)
                edge_dir /= edge_len
                # stay on the shared middle edge: past its corners the line runs into one flower only
                reach = min(half + 3.0, 0.5 * edge_len - 0.5)
                for u in np.linspace(-reach, reach, 25):
                    world = crossing + u * edge_dir
                    za = fa.z_at(*(world - origin))
                    zb = fb.z_at(*(world - n_origin))
                    if abs(u) <= core:
                        assert za == pytest.approx(zb, abs=1e-6), (flower_id, neighbour_id, kind, u)
                    else:
                        assert abs(za - zb) <= 2 * params.noise_mm + 1e-6, (flower_id, neighbour_id, kind, u)
                if kind == "river":
                    # floor = declared crossing height - depth, clamped to keep the roof over the bores
                    h1, h2 = flower.side_corner_heights[k][1], flower.side_corner_heights[k][2]
                    crossing_z = 0.5 * (tileset.meta.heights.z(h1) + tileset.meta.heights.z(h2))
                    floor = max(crossing_z - params.river_depth_mm, params.river_bed_min_z)
                    centre = fa.z_at(*(crossing - origin))
                    assert centre == pytest.approx(floor, abs=1e-6)
                    side = fa.z_at(*(crossing + reach * edge_dir - origin))
                    assert side - centre >= (crossing_z - floor) - 2 * params.noise_mm - 1e-6
                crossings += 1
    assert crossings == 9  # 5 river seams + 4 road seams


def test_render_tileset_writes_one_printable_stl_per_flower(tmp_path, tileset_path) -> None:
    """`render tileset` is how a landscape gets printed: one STL per placed
    flower plus a README that says where each file goes."""
    import trimesh

    from terrain.cli import main

    out_dir = tmp_path / "flowers"
    assert main(["render", "tileset", "--tileset", str(tileset_path), "--output-dir", str(out_dir), "--subdivisions-per-edge", str(N)]) == 0
    ts = load_tileset(tileset_path)
    files = sorted(p.name for p in out_dir.glob("*.stl"))
    assert files == sorted(f"{p.id}.stl" for p in ts.preview_map)
    for name in files[:3]:
        mesh = trimesh.load(str(out_dir / name))
        assert mesh.is_watertight and mesh.volume > 0, name
    readme = (out_dir / "README.md").read_text()
    for p in ts.preview_map:
        assert f"`{p.id}.stl`" in readme and f"({p.at[0]}, {p.at[1]})" in readme
    assert ts.preview_map[0].id in readme.split("```")[1]  # the placement map names every flower


def test_preview_volume_is_the_sum_of_its_flowers(tileset, meshes) -> None:
    preview = build_preview_mesh(tileset, subdivisions_per_edge=N)
    assert preview.volume == pytest.approx(sum(m.volume for m in meshes.values()), rel=1e-6)


def _single_flower_yaml(levels, roads="[]", water="[]") -> str:
    hexes = "\n".join(f'      "{h}": {{ height_level: {lvl} }}' for h, lvl in enumerate(levels))
    sides = "\n".join(f"      {k}: [0, 0, 0, 0]" for k in range(6))
    return "\n".join(
        [
            "meta:",
            "  height_step_mm: 15",
            "  level_count: 4",
            "  scale: 5",
            "  hex_outer_width: 5.1961525",
            "  min_standable_hexes: 1",
            "preview_map: []",
            "flowers:",
            "  steep:",
            "    hexes:",
            hexes,
            "    side_corner_heights:",
            sides,
            "    seed: 7",
            f"    roads: {roads}",
            f"    water: {water}",
            "",
        ]
    )


def test_a_road_climbing_two_levels_in_one_hex_is_rejected(tmp_path) -> None:
    """Peter's rule: a road may change at most one level per hex. The
    centre hex at level 2 between level-0 ring hexes is a 30 mm step."""
    path = tmp_path / "steep.yaml"
    path.write_text(_single_flower_yaml([2, 0, 0, 0, 0, 0, 0], roads="[[0, 3]]"))
    with pytest.raises(TilesetError, match="at most 1 level per hex"):
        load_tileset(path)
    path.write_text(_single_flower_yaml([1, 0, 0, 0, 0, 0, 0], roads="[[0, 3]]"))  # one level per hex: fine
    load_tileset(path)
    path.write_text(_single_flower_yaml([2, 0, 0, 0, 0, 0, 0], water="[[0, 3]]"))  # rivers may fall
    load_tileset(path)


def test_every_generated_road_climbs_at_most_one_level_per_hex(tileset) -> None:
    from terrain.tileset import path_stations

    for flower in tileset.flowers.values():
        for road in flower.roads:
            levels = [level for _, level in path_stations(flower, road)]
            assert max(abs(b - a) for a, b in zip(levels, levels[1:])) <= 1, (flower.id, levels)


def test_a_path_leaving_through_an_adjacent_side_is_rejected(tmp_path) -> None:
    """Entering along one middle edge's normal and leaving along the next
    side's is a 120-degree bend, tighter than any road or river is wide."""
    text = (ROOT / "tilesets" / "default.yaml").read_text()
    flat = text[text.index("  flat_plains:") : text.index("  hill_peak:")]
    bad = flat.replace("roads: []", "roads: [[0, 1]]")
    yaml = "\n".join(
        [
            "meta:",
            "  height_step_mm: 15",
            "  level_count: 4",
            "  scale: 5",
            "  hex_outer_width: 5.1961525",
            "  min_standable_hexes: 1",
            "preview_map: []",
            "flowers:",
            bad,
        ]
    )
    path = tmp_path / "bad.yaml"
    path.write_text(yaml)
    with pytest.raises(TilesetError, match="adjacent"):
        load_tileset(path)
