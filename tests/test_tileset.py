from pathlib import Path

import pytest

from terrain.tileset import TilesetError, load_tileset

ROOT = Path(__file__).resolve().parent.parent
DEFAULT = ROOT / "tilesets" / "default.yaml"


def test_load_default_tileset() -> None:
    tileset = load_tileset(DEFAULT)
    assert "flat_plains" in tileset.flowers
    assert "hill_peak" in tileset.flowers
    assert len(tileset.preview_map) >= 1


def test_heights_from_meta() -> None:
    tileset = load_tileset(DEFAULT)
    heights = tileset.meta.heights
    assert heights.mm_per_level == 15.0
    assert heights.z(0) == 0.0
    assert heights.z(3) == 45.0


def test_hill_peak_side_corner_heights_carry_a_cliff() -> None:
    tileset = load_tileset(DEFAULT)
    side1 = tileset.flowers["hill_peak"].side_corner_heights[1]
    assert side1 == (1, 1, 3, 3)  # no intermediate slope step - decision #7


def test_load_minimal_fixture(fixtures_dir: Path) -> None:
    tileset = load_tileset(fixtures_dir / "minimal_tileset.yaml")
    assert tileset.flowers["alpha"].hexes["0"].height_level == 0


def test_rejects_missing_hex(tmp_path: Path) -> None:
    bad = tmp_path / "missing_hex.yaml"
    bad.write_text(
        """
meta: { height_step_mm: 15, level_count: 4 }
preview_map: []
flowers:
  x:
    hexes:
      "0": { height_level: 0 }
    side_corner_heights:
      0: [0,0,0,0]
      1: [0,0,0,0]
      2: [0,0,0,0]
      3: [0,0,0,0]
      4: [0,0,0,0]
      5: [0,0,0,0]
    seed: 1
"""
    )
    with pytest.raises(TilesetError, match="missing hex"):
        load_tileset(bad)


def test_rejects_out_of_range_height_level(tmp_path: Path) -> None:
    bad = tmp_path / "bad_level.yaml"
    bad.write_text(
        """
meta: { height_step_mm: 15, level_count: 4 }
preview_map: []
flowers:
  x:
    hexes:
      "0": { height_level: 99 }
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
"""
    )
    with pytest.raises(TilesetError, match="height level must be"):
        load_tileset(bad)


def test_rejects_inconsistent_junction(tmp_path: Path) -> None:
    """side k's last corner and side (k+1)'s first corner are the same
    physical point - see the Phase 4 discovery in
    tests/test_flower_mesh_phase4.py's docstring: a mismatched pair here
    would silently produce a real but unintended cliff at that junction.
    """
    bad = tmp_path / "bad_junction.yaml"
    bad.write_text(
        """
meta: { height_step_mm: 15, level_count: 4 }
preview_map: []
flowers:
  x:
    hexes:
      "0": { height_level: 0 }
      "1": { height_level: 0 }
      "2": { height_level: 0 }
      "3": { height_level: 0 }
      "4": { height_level: 0 }
      "5": { height_level: 0 }
      "6": { height_level: 0 }
    side_corner_heights:
      0: [0,0,0,1]
      1: [2,0,0,0]
      2: [0,0,0,0]
      3: [0,0,0,0]
      4: [0,0,0,0]
      5: [0,0,0,0]
    seed: 1
"""
    )
    with pytest.raises(TilesetError, match="same physical junction"):
        load_tileset(bad)


def test_rejects_bad_road_junction(tmp_path: Path) -> None:
    bad = tmp_path / "bad_road.yaml"
    bad.write_text(
        """
meta: { height_step_mm: 15, level_count: 4 }
preview_map: []
flowers:
  x:
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
    roads:
      - { entry_junction: 0, exit_junction: 9 }
"""
    )
    with pytest.raises(TilesetError, match="junction"):
        load_tileset(bad)
