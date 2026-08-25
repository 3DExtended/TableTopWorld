"""Phase 6: base plate + full assembly, welded into one printable solid.

Verifies the actual constructed geometry (watertight/winding/volume, plus
the base-plate weld itself), not just that the functions run without
raising.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from terrain.assembly import build_flower_mesh, standability_report
from terrain.base_plate import build_base_plate_parts
from terrain.constants import BASE_PLATE_DEPTH, MAGNET_CENTER_Z, MAGNET_RADIUS
from terrain.layout import FlowerLayout
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


def test_hex_plateaus_make_every_hex_standable(tileset) -> None:
    """Each hex is built with a flat, noise-free plateau at its declared
    height_level covering 2/3 of its area, so every hex is standable even
    on hill_peak (a deliberate cliff on every side) - the cliffs live in
    the outer band BETWEEN plateaus, not across the part you stand on.

    Before height_level drove geometry, every hex of both fixtures
    measured as unstandable at production defaults, so min_standable_hexes
    could never be satisfied by any flower at all."""
    for flower_id in ("flat_plains", "hill_peak"):
        mesh = build_flower_mesh(tileset, flower_id)
        ok, count, min_required = standability_report(mesh, tileset)
        assert count == FlowerLayout.HEX_CELL_COUNT, flower_id
        assert min_required == 1
        assert ok is True, flower_id


def test_zero_standable_hexes_is_permitted(tileset) -> None:
    """decision #11: a flower may validly have 0 standable hexes, and
    standability_report must report that rather than raise. Exercised
    with plateaus off (use_hex_plateaus=False), which is the purely
    organic terrain the generator produced before height_level was wired
    up - noise and relief alone leave no hex flat within tolerance."""
    mesh = build_flower_mesh(tileset, "hill_peak", use_hex_plateaus=False)
    ok, count, min_required = standability_report(mesh, tileset)
    assert count == 0
    assert min_required == 1
    assert ok is False


def test_plateau_sits_exactly_at_the_declared_height_level(tileset) -> None:
    """height_level must actually drive geometry (it was parsed and
    validated but ignored by the mesh pipeline until now), and the
    plateau must be genuinely flat - not merely flat within tolerance."""
    flower = tileset.flowers["hill_peak"]
    layout = tileset.layout()
    mesh = build_flower_mesh(tileset, "hill_peak")
    for hex_idx in range(FlowerLayout.HEX_CELL_COUNT):
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


def test_base_plate_rejects_too_shallow_plate_for_the_magnet() -> None:
    with pytest.raises(ValueError, match="plate_depth"):
        build_base_plate_parts(
            [],
            [],
            None,  # type: ignore[arg-type]
            subdivisions_per_edge=8,
            bottom_z=0.0,
            plate_depth=MAGNET_CENTER_Z,  # too shallow: no room for the radius
            magnet_center_z=MAGNET_CENTER_Z,
            magnet_radius=MAGNET_RADIUS,
        )


def test_base_plate_default_depth_fits_the_magnet() -> None:
    assert BASE_PLATE_DEPTH > MAGNET_CENTER_Z + MAGNET_RADIUS
