"""Shared fixtures for atomic visual test suite."""

from __future__ import annotations

from pathlib import Path

import pytest

from terrain.layout import FlowerLayout
from terrain.tileset import load_tileset
from tests.helpers.visual_scad import VISUAL_ATOMS_DIR, write_visual_scad

ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_TILESET = ROOT / "tilesets" / "default.yaml"


@pytest.fixture(scope="session")
def default_tileset():
    return load_tileset(DEFAULT_TILESET)


@pytest.fixture
def layout() -> FlowerLayout:
    return FlowerLayout()


@pytest.fixture
def visual():
    """Export a solid to output/visual_atoms/ and return the path."""

    def _export(
        name: str,
        solid,
        description: str,
        *,
        scale: float = 5.0,
        subdir: str = "",
    ) -> Path:
        return write_visual_scad(name, solid, description, scale=scale, subdir=subdir)

    return _export


@pytest.fixture(scope="session", autouse=True)
def _ensure_visual_output_dir() -> None:
    VISUAL_ATOMS_DIR.mkdir(parents=True, exist_ok=True)
