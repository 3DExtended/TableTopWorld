"""Write OpenSCAD files for manual visual verification alongside unit tests."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from solid2 import OpenSCADObject

ROOT = Path(__file__).resolve().parent.parent.parent
VISUAL_ATOMS_DIR = ROOT / "output" / "visual_atoms"


def write_visual_scad(
    name: str,
    solid: OpenSCADObject,
    description: str,
    *,
    scale: float = 5.0,
    subdir: str = "",
) -> Path:
    """Export one solid to output/visual_atoms/<subdir>/<name>.scad with a header comment."""
    out_dir = VISUAL_ATOMS_DIR / subdir if subdir else VISUAL_ATOMS_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{name}.scad"
    solid.scale(scale).save_as_scad(str(path))
    body = path.read_text()
    header = (
        f"// {description}\n"
        f"// Atom visual: {name}.scad — compare in OpenSCAD with the matching pytest.\n\n"
    )
    path.write_text(header + body)
    return path
