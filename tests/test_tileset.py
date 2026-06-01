from pathlib import Path

import pytest

from terrain.constants import TERRAIN_Z
from terrain.tileset import TilesetError, load_tileset

ROOT = Path(__file__).resolve().parent.parent
DEFAULT = ROOT / "tilesets" / "default.yaml"


def test_load_default_tileset() -> None:
    tileset = load_tileset(DEFAULT)
    assert "flat_plains" in tileset.flowers
    assert len(tileset.preview_map) >= 3


def test_model_step_from_meta() -> None:
    tileset = load_tileset(DEFAULT)
    assert tileset.meta.model_step == 4.0
    assert tileset.terrain_z("middle") == TERRAIN_Z["middle"]
    assert tileset.water_z("ground") == -4.0


def test_load_minimal_fixture(fixtures_dir: Path) -> None:
    tileset = load_tileset(fixtures_dir / "minimal_tileset.yaml")
    assert tileset.flowers["alpha"].hexes["0"].role == "standable"


def test_rejects_no_standable(fixtures_dir: Path) -> None:
    with pytest.raises(TilesetError, match="at least one standable"):
        load_tileset(fixtures_dir / "invalid_no_standable.yaml")


def test_rejects_bad_edge_key(fixtures_dir: Path, tmp_path: Path) -> None:
    bad = tmp_path / "bad_edge.yaml"
    bad.write_text(
        """
meta: { height_step_mm: 20, scale: 5 }
preview_map: []
flowers:
  x:
    hexes:
      "0": { terrain: ground, role: standable }
      "1": { terrain: ground, role: standable }
      "2": { terrain: ground, role: standable }
      "3": { terrain: ground, role: standable }
      "4": { terrain: ground, role: standable }
      "5": { terrain: ground, role: standable }
      "6": { terrain: ground, role: standable }
    edges:
      "9-9": { profile: flat_ground }
    roads: []
    water: []
"""
    )
    with pytest.raises(TilesetError, match="not a valid exterior side"):
        load_tileset(bad)


def test_rejects_bad_junction(tmp_path: Path) -> None:
    bad = tmp_path / "bad_junction.yaml"
    bad.write_text(
        """
meta: { height_step_mm: 20, scale: 5 }
preview_map: []
flowers:
  x:
    hexes:
      "0": { terrain: ground, role: standable }
      "1": { terrain: ground, role: standable }
      "2": { terrain: ground, role: standable }
      "3": { terrain: ground, role: standable }
      "4": { terrain: ground, role: standable }
      "5": { terrain: ground, role: standable }
      "6": { terrain: ground, role: standable }
    edges:
      "1-0": { profile: flat_ground }
      "1-1": { profile: flat_ground }
      "1-2": { profile: flat_ground }
      "2-0": { profile: flat_ground }
      "2-1": { profile: flat_ground }
      "2-2": { profile: flat_ground }
      "3-0": { profile: flat_ground }
      "3-1": { profile: flat_ground }
      "3-2": { profile: flat_ground }
      "4-0": { profile: flat_ground }
      "4-1": { profile: flat_ground }
      "4-2": { profile: flat_ground }
      "5-0": { profile: flat_ground }
      "5-1": { profile: flat_ground }
      "5-2": { profile: flat_ground }
      "6-0": { profile: flat_ground }
      "6-1": { profile: flat_ground }
      "6-2": { profile: flat_ground }
    roads:
      - { entry_junction: 0, exit_junction: 9 }
    water: []
"""
    )
    with pytest.raises(TilesetError, match="junction"):
        load_tileset(bad)
