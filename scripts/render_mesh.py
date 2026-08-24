"""Render a mesh (STL or trimesh.Trimesh) to PNG snapshots, headlessly.

No GPU/OpenGL/pyglet dependency (trimesh's default scene viewer needs
those and isn't available in this environment) - uses matplotlib's
mplot3d Poly3DCollection instead, with simple Z-based flat shading so
terrain relief and grooves are visible without real lighting.

Primarily a self-verification tool: render a just-generated STL and look
at the PNG before telling a human the geometry is correct, rather than
relying on watertight/volume checks alone to catch something a human
would immediately see (e.g. a tiling seam mismatch, unexpectedly flat
terrain, a missing groove).

Usage:
    python scripts/render_mesh.py output/flower.stl
    python scripts/render_mesh.py output/flower.stl --out /tmp/flower.png --views top,iso,front
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import trimesh
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

VIEWS = {
    "top": (90, -90),
    "iso": (30, -60),
    "front": (0, -90),
    "side": (0, 0),
}


def render_mesh(
    mesh: trimesh.Trimesh, out_path: str | Path, *, views: list[str] | None = None
) -> Path:
    views = views or ["iso", "top"]
    out_path = Path(out_path)

    tris = mesh.vertices[mesh.faces]
    z = tris[:, :, 2].mean(axis=1)
    z_norm = (z - z.min()) / max(z.max() - z.min(), 1e-9)
    colors = plt.cm.terrain(0.25 + 0.5 * z_norm)

    bounds = mesh.bounds
    center = bounds.mean(axis=0)
    radius = np.linalg.norm(bounds[1] - bounds[0]) / 2

    fig = plt.figure(figsize=(6 * len(views), 6))
    for i, view in enumerate(views):
        elev, azim = VIEWS[view]
        ax = fig.add_subplot(1, len(views), i + 1, projection="3d")
        coll = Poly3DCollection(tris, facecolor=colors, edgecolor="none", linewidths=0)
        ax.add_collection3d(coll)
        ax.set_xlim(center[0] - radius, center[0] + radius)
        ax.set_ylim(center[1] - radius, center[1] + radius)
        ax.set_zlim(center[2] - radius, center[2] + radius)
        ax.view_init(elev=elev, azim=azim)
        ax.set_box_aspect((1, 1, 1))
        ax.set_axis_off()
        ax.set_title(view)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("stl_path", type=Path)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--views", type=str, default="iso,top")
    args = parser.parse_args(argv)

    mesh = trimesh.load(str(args.stl_path))
    out = args.out or args.stl_path.with_suffix(".png")
    views = [v.strip() for v in args.views.split(",")]
    result = render_mesh(mesh, out, views=views)
    print(f"wrote {result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
