"""Phase 7: multi-flower preview placement + cross-flower validation.

build_preview_mesh() places every tileset.preview_map entry via
FlowerLayout.flower_grid_to_xy() - the fix for the axial_to_xy()/
flower_center_spacing placement bug flagged since Phase 1 (see
terrain/layout.py). tileset.py's validate_tileset() now also checks that
adjacent rot=0 preview_map placements declare the reversed-matching
side_corner_heights contract (decision #3/#4) - deferred in earlier
phases until the placement math was known correct.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from terrain.assembly import build_preview_mesh
from terrain.tileset import TilesetError, load_tileset

ROOT = Path(__file__).resolve().parent.parent
DEFAULT = ROOT / "tilesets" / "default.yaml"


@pytest.fixture(scope="module")
def tileset():
    return load_tileset(DEFAULT)


def test_default_tileset_preview_map_has_two_adjacent_flowers(tileset) -> None:
    """tilesets/default.yaml places hill_peak at (-1, 1) relative to
    flat_plains at (0, 0) specifically because that's the one delta where
    the two flowers' side_corner_heights actually match (flat_plains side
    2, hill_peak side 5, both all-zero) - see the YAML's own comment."""
    assert len(tileset.preview_map) == 2
    positions = {p.at for p in tileset.preview_map}
    assert positions == {(0, 0), (-1, 1)}


def test_preview_mesh_is_watertight_and_positive_volume(tileset) -> None:
    """Each flower in the scene is its own independently-watertight
    print (decision #13's magnets, not shared mesh topology, do the
    mating) - concatenating them should still read as watertight,
    winding-consistent, positive-volume overall."""
    mesh = build_preview_mesh(tileset)
    assert mesh.is_watertight
    assert mesh.is_winding_consistent
    assert mesh.volume > 0


def test_preview_mesh_volume_is_sum_of_individual_flowers(tileset) -> None:
    from terrain.assembly import build_flower_mesh

    flat = build_flower_mesh(tileset, "flat_plains")
    peak = build_flower_mesh(tileset, "hill_peak")
    preview = build_preview_mesh(tileset)
    assert preview.volume == pytest.approx(flat.volume + peak.volume, rel=1e-6)


def test_adjacent_flowers_boundary_matches_along_its_full_length(tileset) -> None:
    """Regression test for a real bug that ALL the other preview tests in
    this file passed straight through: flat_plains and hill_peak's
    matching side (side 2 / side 5, both declared (0,0,0,0)) only lined
    up at 4 coarse corners, not along the fine jagged contour between
    them - a palindromic corner-height sequence defeated
    canonicalize_sequence_position's direction-disambiguation (see
    tests/test_boundary_determinism.py). Watertightness, volume, and even
    a visual top-down render all looked plausible in isolation; only
    checking actual (x, y, z) coincidence point-by-point along the shared
    edge - not just its 4 declared corners - catches this. This is the
    check that would have caught the original regression automatically.
    """
    import numpy as np

    from terrain.heightfield import build_side_boundary_vertices

    layout = tileset.layout()
    level_z = tileset.meta.heights.z
    flat = tileset.flowers["flat_plains"]
    peak = tileset.flowers["hill_peak"]

    flat_side = build_side_boundary_vertices(
        flat.side_corner_heights[2], layout.side_corners(2), level_z
    )
    peak_side = build_side_boundary_vertices(
        peak.side_corner_heights[5], layout.side_corners(5), level_z
    )

    x0, y0 = layout.flower_grid_to_xy(0, 0)
    x1, y1 = layout.flower_grid_to_xy(-1, 1)
    dx, dy = x1 - x0, y1 - y0
    peak_side_shifted = [(x + dx, y + dy, z) for x, y, z in peak_side]

    assert len(flat_side) == len(peak_side_shifted)
    for a, b in zip(flat_side, reversed(peak_side_shifted)):
        assert np.allclose(a, b, atol=1e-6), (
            f"boundary point mismatch: flat_plains has {a}, "
            f"hill_peak (shifted) has {b} at the same fine-contour position"
        )


def test_flower_grid_to_xy_places_neighbors_at_the_real_adjacency_offset(tileset) -> None:
    """The old axial_to_xy()/flower_center_spacing convention (removed in
    this phase) did not correspond to true edge-sharing adjacency - this
    checks the replacement, flower_grid_to_xy(), against the same
    neighbor_flower_offset() ground truth Phase 1 verified side_corners()
    against."""
    layout = tileset.layout()
    origin = layout.flower_grid_to_xy(0, 0)
    for k in range(6):
        dq, dr = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)][k]
        neighbor = layout.flower_grid_to_xy(dq, dr)
        expected_offset = layout.neighbor_flower_offset(k)
        got_offset = (neighbor[0] - origin[0], neighbor[1] - origin[1])
        assert got_offset == pytest.approx(expected_offset, abs=1e-6)


def test_cross_flower_side_mismatch_is_rejected(tmp_path: Path) -> None:
    """Two adjacent rot=0 placements whose facing sides don't reverse-
    match must be rejected - the deferred check added to
    validate_tileset() once flower_grid_to_xy() gave it a correct
    adjacency convention to check against."""
    bad = tmp_path / "bad_preview_adjacency.yaml"
    bad.write_text(
        """
meta: { height_step_mm: 15, level_count: 4, scale: 5 }
preview_map:
  - { id: a, at: [0, 0], rot: 0 }
  - { id: b, at: [1, 0], rot: 0 }
flowers:
  a:
    hexes:
      "0": { height_level: 0 }
      "1": { height_level: 0 }
      "2": { height_level: 0 }
      "3": { height_level: 0 }
      "4": { height_level: 0 }
      "5": { height_level: 0 }
      "6": { height_level: 0 }
    side_corner_heights:
      0: [0,0,0,0]
      1: [0,0,0,0]
      2: [0,0,0,0]
      3: [0,0,0,0]
      4: [0,0,0,0]
      5: [0,0,0,0]
    seed: 1
  b:
    hexes:
      "0": { height_level: 0 }
      "1": { height_level: 0 }
      "2": { height_level: 0 }
      "3": { height_level: 0 }
      "4": { height_level: 0 }
      "5": { height_level: 0 }
      "6": { height_level: 0 }
    side_corner_heights:
      0: [0,0,0,0]
      1: [0,0,0,0]
      2: [0,0,0,0]
      3: [0,1,1,0]
      4: [0,0,0,0]
      5: [0,0,0,0]
    seed: 1
"""
    )
    with pytest.raises(TilesetError, match="does not match"):
        load_tileset(bad)


def test_cli_render_preview_writes_a_loadable_stl(tmp_path: Path) -> None:
    """Not checked here: reloaded-from-disk watertightness. STL has no
    concept of separate bodies - it's one flat triangle soup - and
    flat_plains/hill_peak are placed deliberately touching along their
    shared, bit-identical (decision #4) boundary. Reloading with
    trimesh's default vertex-merge welds that seam across the two
    otherwise-independent solids, which correctly fails a manifold check
    (a real edge there would need exactly 2 faces, a merged seam between
    two touching solids gives 4) - not a defect in the generated
    geometry, just not the same "watertight" build_preview_mesh's own
    in-memory result already covers (see the test above, which is the
    real invariant: export_stl's own pre-write gate passed on the
    unmerged-across-parts in-memory mesh)."""
    out = tmp_path / "preview.stl"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "terrain.cli",
            "render",
            "preview",
            "--tileset",
            str(DEFAULT),
            "--output",
            str(out),
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert out.exists()

    import trimesh

    mesh = trimesh.load(str(out))
    assert len(mesh.faces) > 0
    assert mesh.volume > 0


def test_no_references_to_discarded_modules() -> None:
    """The old CSG pipeline (terrain/mesh.py, edges.py, features.py,
    catalog.py, and their CLI/assembly wiring) was deleted in Phase 4;
    nothing in terrain/ should still import from it."""
    discarded = {"terrain.mesh", "terrain.edges", "terrain.features", "terrain.catalog"}
    terrain_dir = ROOT / "terrain"
    for path in terrain_dir.rglob("*.py"):
        text = path.read_text()
        for name in discarded:
            assert name not in text, f"{path} still references discarded module {name!r}"
