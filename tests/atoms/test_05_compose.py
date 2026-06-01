"""Layer 5 — composed flowers: mesh + edges + features from atoms upward."""

from __future__ import annotations

from pathlib import Path

from terrain.assembly import AssemblyExporter
from terrain.render.format import format_render_spec
from terrain.render.planner import plan_flower
from tests.helpers.atom_builders import atom_height_step, atom_road_only, atom_slope_pad

ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT = ROOT / "tilesets" / "default.yaml"


def test_22_height_step_mesh_plus_bevels(visual) -> None:
    """Height-step atom: mesh + bevels (still no magnets or features)."""
    tileset = atom_height_step()
    exporter = AssemblyExporter(tileset, resolution=32)
    flower = tileset.flowers["atom_height_step"]
    solid = exporter.mesh_builder.build_flower(flower)
    max_z = exporter.mesh_builder.max_flower_z(flower)
    solid = exporter.edge_geom.apply_bevels(solid, max_z)
    visual(
        "22_height_step_bevels",
        solid,
        "Composed: 11_seven_hex_union + 16_bevels_on_height_step.",
        subdir="05_compose",
    )


def test_23_height_step_full_edges(visual) -> None:
    """Height-step with bevels and all mating magnets."""
    tileset = atom_height_step()
    exporter = AssemblyExporter(tileset, resolution=32)
    solid = exporter.build_flower("atom_height_step")
    spec = exporter.describe_flower("atom_height_step")
    assert len(spec.magnet_holes) == 18
    visual(
        "23_height_step_full_edges",
        solid,
        "Mesh + bevels + 18 magnets on height-step atom.",
        subdir="05_compose",
    )


def test_24_slope_atom_full_export(visual) -> None:
    """Slope-focused atom exported through full assembly pipeline."""
    tileset = atom_slope_pad()
    exporter = AssemblyExporter(tileset, resolution=32)
    solid = exporter.build_flower("atom_slope_pad")
    spec = plan_flower(tileset, tileset.flowers["atom_slope_pad"])
    assert spec.hexes[5].slope_ramp_height is not None
    visual(
        "24_slope_atom_full",
        solid,
        "Full pipeline on slope atom — ramp hex 5 between middle neighbors.",
        subdir="05_compose",
    )


def test_25_road_atom_full_export(visual) -> None:
    """Road atom: mesh + bevels + magnets + road cut."""
    tileset = atom_road_only()
    exporter = AssemblyExporter(tileset, resolution=32)
    solid = exporter.build_flower("atom_road_only")
    assert len(plan_flower(tileset, tileset.flowers["atom_road_only"]).path_cuts) == 1
    visual(
        "25_road_atom_full",
        solid,
        "Full flat flower with road — compare groove to 19_road_on_plateau.",
        subdir="05_compose",
    )


def test_26_flat_plains_production(visual, default_tileset) -> None:
    """Production flat_plains: ground + road + magnets + bevels."""
    exporter = AssemblyExporter(default_tileset, resolution=32)
    solid = exporter.build_flower("flat_plains")
    spec = exporter.describe_flower("flat_plains")
    assert spec.max_z == 0.0
    assert any(p.kind == "road" for p in spec.path_cuts)
    visual(
        "26_flat_plains",
        solid,
        "Reference tile flat_plains — all ground with one road.",
        subdir="05_compose",
    )


def test_27_hill_north_production(visual, default_tileset) -> None:
    """Production hill_north: mixed heights, slopes, no water."""
    exporter = AssemblyExporter(default_tileset, resolution=32)
    solid = exporter.build_flower("hill_north")
    spec = exporter.describe_flower("hill_north")
    assert spec.max_z == 8.0
    assert spec.path_cuts == ()
    visual(
        "27_hill_north",
        solid,
        "Reference hill_north — high peak, slopes, no road.",
        subdir="05_compose",
    )


def test_28_river_grove_production(visual, default_tileset) -> None:
    """Production river_grove: water channel on varied terrain."""
    exporter = AssemblyExporter(default_tileset, resolution=32)
    solid = exporter.build_flower("river_grove")
    spec = exporter.describe_flower("river_grove")
    assert any(p.kind == "water" for p in spec.path_cuts)
    visual(
        "28_river_grove",
        solid,
        "Reference river_grove — includes water cut.",
        subdir="05_compose",
    )


def test_29_preview_map_two_flowers(visual, default_tileset) -> None:
    """Two flowers from preview_map placed side by side."""
    exporter = AssemblyExporter(default_tileset, resolution=32)
    solid = exporter.build_preview()
    visual(
        "29_preview_two_flowers",
        solid,
        "All flowers in default preview_map placed in axial layout.",
        subdir="05_compose",
    )


def test_30_render_spec_snapshot_written(default_tileset, visual) -> None:
    """Text render spec for hill_north saved next to SCAD for cross-check."""
    exporter = AssemblyExporter(default_tileset, resolution=32)
    spec = exporter.describe_flower("hill_north")
    text = format_render_spec(spec)
    out = ROOT / "output" / "visual_atoms" / "05_compose" / "30_hill_north.render.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    assert '"flower_id": "hill_north"' in text
    assert spec.hexes[3].top_z == 8.0
