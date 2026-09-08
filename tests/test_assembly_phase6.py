"""Phase 6: base plate + magnet bores + full assembly, one printable solid.

Verifies the actual constructed geometry (watertight/winding/volume, the
plate depth, the bore positions and the removed bore volume), not just
that the functions run without raising.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from terrain.assembly import (
    build_flower_mesh,
    build_preview_mesh,
    pick_standable_hexes,
    standability_report,
)
from terrain.constants import BASE_PLATE_DEPTH_MM, MAGNET_CENTER_Z_MM
from terrain.layout import FlowerLayout
from terrain.magnets import DEFAULT_MAGNET_BORES, MagnetBores
from terrain.standability import hex_cell_z_range
from terrain.tileset import load_tileset

ROOT = Path(__file__).resolve().parent.parent
DEFAULT = ROOT / "tilesets" / "default.yaml"


@pytest.fixture(scope="module")
def tileset():
    return load_tileset(DEFAULT)


@pytest.mark.parametrize("flower_id", ["flat_plains", "hill_peak"])
def test_full_flower_is_watertight_and_positive_volume(tileset, flower_id: str) -> None:
    mesh = build_flower_mesh(tileset, flower_id)
    assert mesh.is_watertight
    assert mesh.is_winding_consistent
    assert mesh.volume > 0


def test_center_hex_is_not_a_low_poly_fan(tileset) -> None:
    """Regression test: internal (hex-to-hex) edges used to be a single
    raw corner-to-corner segment with no subdivision at all - only the
    flower's own exterior boundary got fine detail. The center hex cell
    (index 0) has ALL 6 of its edges internal, so it rendered as a
    literal 6-triangle flat fan regardless of subdivisions_per_edge - a
    real "too low poly" look a human immediately noticed. Internal edges
    now get the same subdivisions_per_edge granularity as the exterior
    boundary. This checks the *whole* flower's triangle count scales
    with subdivisions_per_edge accordingly, well above what the old
    unsubdivided-interior code could ever produce (872 faces for
    flat_plains at subdivisions_per_edge=8, before this fix)."""
    mesh = build_flower_mesh(tileset, "flat_plains", subdivisions_per_edge=8)
    assert len(mesh.faces) > 900

    finer = build_flower_mesh(tileset, "flat_plains", subdivisions_per_edge=16)
    assert len(finer.faces) > len(mesh.faces)


def test_the_generated_mesh_matches_the_standable_hexes_picked(tileset) -> None:
    """The number of standable hexes varies per flower (drawn from its own
    seed), and the REAL mesh must agree with that choice - every picked
    cell measures flat, every unpicked one does not. This is the check
    that ties the seeded decision to the actual geometry rather than
    trusting one or the other in isolation."""
    layout = tileset.layout()
    for flower_id in ("flat_plains", "hill_peak"):
        flower = tileset.flowers[flower_id]
        picked = pick_standable_hexes(
            flower.seed, min_standable=tileset.meta.min_standable_hexes
        )
        mesh = build_flower_mesh(tileset, flower_id)
        ok, count, min_required = standability_report(mesh, tileset)
        assert count == len(picked), flower_id
        assert count >= min_required, flower_id
        assert ok is True, flower_id
        for hex_idx in range(FlowerLayout.HEX_CELL_COUNT):
            z_range = hex_cell_z_range(mesh, layout, hex_idx)
            if hex_idx in picked:
                assert z_range == pytest.approx(0.0, abs=1e-6), (flower_id, hex_idx)
            else:
                assert z_range > 1.0, (flower_id, hex_idx)


def test_standable_hex_count_is_seeded_varied_and_reproducible() -> None:
    """It must genuinely vary between flowers (not quietly always 7), stay
    within [min_standable, 7], and be a pure function of the seed - the
    same reproducibility contract as every other seeded choice here."""
    counts = {
        len(pick_standable_hexes(seed, min_standable=1)) for seed in range(1, 40)
    }
    assert len(counts) > 1, "count never varies - randomization is not doing anything"
    assert min(counts) >= 1 and max(counts) <= FlowerLayout.HEX_CELL_COUNT

    for seed in range(1, 10):
        assert pick_standable_hexes(seed, min_standable=1) == pick_standable_hexes(
            seed, min_standable=1
        )

    # the declared minimum is always honoured, even when it is the whole flower
    for seed in range(1, 10):
        assert len(pick_standable_hexes(seed, min_standable=7)) == 7


def test_zero_standable_hexes_is_permitted(tileset) -> None:
    """decision #11: a flower may validly have 0 standable hexes, and
    standability_report must report that rather than raise. Exercised by
    pinning standable_hexes=0, which is the purely organic terrain the
    generator produced before height_level was wired up - noise and relief
    alone leave no hex flat within tolerance."""
    mesh = build_flower_mesh(tileset, "hill_peak", standable_hexes=0)
    ok, count, min_required = standability_report(mesh, tileset)
    assert count == 0
    assert min_required == 1
    assert ok is False


def test_plateau_sits_exactly_at_the_declared_height_level(tileset) -> None:
    """height_level must actually drive geometry (it was parsed and
    validated but ignored by the mesh pipeline until now), and a plateau
    must be genuinely flat - not merely flat within tolerance."""
    flower = tileset.flowers["hill_peak"]
    layout = tileset.layout()
    mesh = build_flower_mesh(tileset, "hill_peak")
    picked = pick_standable_hexes(
        flower.seed, min_standable=tileset.meta.min_standable_hexes
    )
    assert picked, "fixture must have at least one plateau for this to test anything"
    for hex_idx in picked:
        declared = tileset.meta.heights.z(flower.hexes[str(hex_idx)].height_level)
        assert hex_cell_z_range(mesh, layout, hex_idx) == pytest.approx(0.0, abs=1e-6)
        # and it is flat AT the declared level, not just flat somewhere
        cx, cy = layout.cell_center(hex_idx)
        centre_z = max(
            float(z)
            for x, y, z in mesh.vertices
            if abs(float(x) - cx) < 1e-6 and abs(float(y) - cy) < 1e-6
        )
        assert centre_z == pytest.approx(declared, abs=1e-6)


def test_flat_plains_center_hex_is_standable(tileset) -> None:
    """decision #11's positive case: with relief/jitter turned off, a
    flower with uniform declared heights is genuinely flat and standable.

    Not built with default jitter_amplitude/interior_relief_mm/groove_depth_mm:
    those (0.3 * 15mm one_level_z = 4.5mm jitter swing, 6mm interior relief
    amplitude, and a real engraved groove reaching its full nominal depth
    now that groove width is decoupled from mesh resolution) apply
    unconditionally to every hex's boundary/interior/edge-adjacent points
    regardless of declared height (decision #5 - interior relief is
    deliberately independent of the boundary contract), and each exceeds
    flatness_tolerance_mm's default of 1.0mm on its own. Under those
    defaults even flat_plains legitimately has 0 standable hexes - that's
    not a bug, just a stronger statement than this test needs to make.
    This test isolates the specific claim decision #11 requires: a
    uniform-height flower CAN be standable, once nothing is sculpting it
    away from flat.
    """
    mesh = build_flower_mesh(
        tileset,
        "flat_plains",
        jitter_amplitude=0.0,
        interior_relief_mm=0.0,
        groove_depth_mm=0.0,
    )
    ok, count, _ = standability_report(mesh, tileset)
    assert count >= 1
    assert ok is True


def _expected_bore_centres(layout, bores: MagnetBores) -> list[tuple[float, float, float]]:
    """The blind end's centre of each of the 18 bores: the silhouette
    edge's midpoint pushed `depth_mm` INTO the flower, at the fixed
    magnet height above the print bed."""
    z = -BASE_PLATE_DEPTH_MM + bores.center_above_bed_mm
    centres = []
    for side in range(6):
        corners = layout.side_corners(side)
        for e in range(3):
            (x0, y0), (x1, y1) = corners[e], corners[e + 1]
            mx, my = (x0 + x1) / 2, (y0 + y1) / 2
            tx, ty = x1 - x0, y1 - y0
            ox, oy = ty, -tx
            if ox * mx + oy * my < 0:
                ox, oy = -ox, -oy
            norm = math.hypot(ox, oy)
            ox, oy = ox / norm, oy / norm
            centres.append((mx - ox * bores.depth_mm, my - oy * bores.depth_mm, z))
    return centres


def test_plate_is_ten_millimetres_deep_in_real_millimetres(tileset) -> None:
    """The print bed sits BASE_PLATE_DEPTH_MM below the level-0 surface -
    a physical distance that must not scale with meta.scale (until 2026-09
    the plate constants were in model units and never scaled, silently
    giving a 2 mm plate)."""
    mesh = build_flower_mesh(tileset, "flat_plains")
    assert BASE_PLATE_DEPTH_MM == 10.0
    assert float(mesh.vertices[:, 2].min()) == pytest.approx(-BASE_PLATE_DEPTH_MM)


@pytest.mark.parametrize("flower_id", ["flat_plains", "hill_peak"])
def test_every_silhouette_edge_carries_one_blind_magnet_bore(tileset, flower_id: str) -> None:
    """18 bores (one per exterior edge), each a blind recess of the right
    depth at the fixed magnet height: the blind end's centre vertex must
    exist in the welded mesh at exactly the expected position, and the
    solid stays watertight around every one of them."""
    layout = tileset.layout()
    mesh = build_flower_mesh(tileset, flower_id)
    assert mesh.is_watertight and mesh.is_winding_consistent
    verts = np.asarray(mesh.vertices)
    expected = _expected_bore_centres(layout, DEFAULT_MAGNET_BORES)
    assert len(expected) == 18
    for cx, cy, cz in expected:
        d = np.linalg.norm(verts - np.array([cx, cy, cz]), axis=1)
        assert d.min() < 1e-6, (flower_id, (cx, cy, cz), d.min())
    assert MAGNET_CENTER_Z_MM == 3.9


def test_magnet_bores_remove_their_cylinder_volume(tileset) -> None:
    """The bores are real cavities, not decoration: the solid built with
    bores is lighter than the same solid without them by 18 blind
    cylinders (a 24-gon prism, hence the tolerance)."""
    with_bores = build_flower_mesh(tileset, "flat_plains")
    without = build_flower_mesh(tileset, "flat_plains", magnet_bores=None)
    assert without.is_watertight and with_bores.is_watertight
    b = DEFAULT_MAGNET_BORES
    expected = 18 * math.pi * b.radius_mm**2 * b.depth_mm
    removed = without.volume - with_bores.volume
    assert removed == pytest.approx(expected, rel=0.05)


def test_neighbouring_flowers_bores_face_each_other(tileset) -> None:
    """Two flowers placed as preview_map neighbours must present their
    bores at the same physical points along the shared side, so the
    magnets actually meet: each of this side's 3 bore mouths (edge
    midpoint at magnet height) coincides for both flowers, and the
    two blind ends sit depth_mm apart on opposite sides of the seam."""
    layout = tileset.layout()
    b = DEFAULT_MAGNET_BORES
    z = -BASE_PLATE_DEPTH_MM + b.center_above_bed_mm
    placements = {p.id: p for p in tileset.preview_map}
    assert len(placements) == 2
    ends: dict[str, list[np.ndarray]] = {}
    for pid, placement in placements.items():
        mesh = build_flower_mesh(tileset, pid)
        x, y = layout.flower_grid_to_xy(*placement.at)
        centres = np.array(_expected_bore_centres(layout, b)) + np.array([x, y, 0.0])
        ends[pid] = list(centres)
    a, c = ends.values()
    pairs = 0
    for pa in a:
        for pc in c:
            d = np.linalg.norm(pa - pc)
            if abs(d - 2 * b.depth_mm) < 1e-6:
                mouth = (pa + pc) / 2
                assert mouth[2] == pytest.approx(z)
                pairs += 1
    assert pairs == 3


def test_too_shallow_plate_is_rejected(tileset) -> None:
    with pytest.raises(ValueError, match="too shallow"):
        build_flower_mesh(tileset, "flat_plains", plate_depth_mm=MAGNET_CENTER_Z_MM)


def test_default_plate_depth_leaves_a_roof_over_the_bore() -> None:
    b = DEFAULT_MAGNET_BORES
    assert b.top_above_bed_mm == pytest.approx(3.9 + 2.65)
    assert BASE_PLATE_DEPTH_MM >= b.min_plate_depth_mm


def test_preview_scene_still_assembles_with_bores(tileset) -> None:
    scene = build_preview_mesh(tileset)
    assert scene.volume > 0
