"""Annotate a top-down render of a tileset's preview with its flower sides,
coloured by mating class - a picture of how modular the set is.

Two sides can be placed against each other exactly when their four corner
heights match once one of them is read backwards (the boundary contract)
and the same feature (road, river or nothing) crosses their middle edge.
So every side gets a class = (canonical corner profile, crossing feature);
equal letters in the picture can mate, any flower rotated as needed.

Usage:
    .venv/bin/python scripts/annotate_sides.py --tileset tilesets/landscape.yaml --out output/landscape_sides.png
    (add --stl output/landscape.stl to draw over an existing render's mesh)
"""

from __future__ import annotations

import argparse
import string
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import trimesh
from matplotlib.collections import PolyCollection

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from terrain.assembly import build_preview_mesh  # noqa: E402
from terrain.tileset import _NEIGHBOR_GRID_DELTAS, Tileset, load_tileset  # noqa: E402


def side_feature(flower, side: int) -> str:
    if any(side in (p.entry, p.exit) for p in flower.roads):
        return "road"
    if any(side in (p.entry, p.exit) for p in flower.water):
        return "river"
    return ""


def side_classes(tileset: Tileset) -> tuple[dict[tuple, str], list[dict]]:
    """Every placed side with its class key; keys lettered by frequency."""
    placed = {p.at: p.id for p in tileset.preview_map}
    sides = []
    for (q, r), fid in placed.items():
        flower = tileset.flowers[fid]
        for k in range(6):
            profile = tuple(flower.side_corner_heights[k])
            key = (min(profile, profile[::-1]), side_feature(flower, k))
            dq, dr = _NEIGHBOR_GRID_DELTAS[k]
            joined = (q + dq, r + dr) in placed
            sides.append({"id": fid, "at": (q, r), "side": k, "key": key, "joined": joined, "profile": profile})
    counts = Counter(s["key"] for s in sides)
    letters = {}
    names = list(string.ascii_uppercase) + [a + b for a in string.ascii_uppercase for b in string.ascii_uppercase]
    for i, (key, _) in enumerate(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0][1], kv[0][0]))):
        letters[key] = names[i]
    return letters, sides


def draw_top_down(ax, mesh: trimesh.Trimesh) -> None:
    tris = mesh.vertices[mesh.faces]
    z = tris[:, :, 2].mean(axis=1)
    z_norm = (z - z.min()) / max(z.max() - z.min(), 1e-9)
    colors = plt.cm.terrain(0.25 + 0.5 * z_norm)
    light = np.array([0.3, 0.3, 1.0])
    light /= np.linalg.norm(light)
    brightness = np.clip(mesh.face_normals @ light, 0.35, 1.0)
    colors = np.clip(colors[:, :3] * brightness[:, None], 0, 1)
    order = np.argsort(z)  # painter's algorithm: higher faces drawn last
    ax.add_collection(PolyCollection(tris[order][:, :, :2], facecolors=colors[order], edgecolors="none", antialiased=False))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--tileset", type=Path, default=ROOT / "tilesets" / "landscape.yaml")
    parser.add_argument("--stl", type=Path, default=None, help="render to draw over (default: build the preview at n=8)")
    parser.add_argument("--out", type=Path, default=ROOT / "output" / "landscape_sides.png")
    args = parser.parse_args(argv)

    tileset = load_tileset(args.tileset)
    layout = tileset.layout()
    mesh = trimesh.load(str(args.stl)) if args.stl else build_preview_mesh(tileset, subdivisions_per_edge=8)
    letters, sides = side_classes(tileset)
    counts = Counter(s["key"] for s in sides)
    open_counts = Counter(s["key"] for s in sides if not s["joined"])
    palette = plt.cm.tab20(np.linspace(0, 1, 20))

    def style(key):  # 20 distinct colours; classes beyond that are dashed
        i = list(letters).index(key)
        return palette[i % 20], "-" if i < 20 else (0, (2, 1))

    fig, ax = plt.subplots(figsize=(15, 13))
    draw_top_down(ax, mesh)
    labelled_seams = set()
    for s in sides:
        origin = np.array(layout.flower_grid_to_xy(*s["at"]))
        corners = np.array(layout.side_corners(s["side"])) + origin
        colour, dash = style(s["key"])
        letter = letters[s["key"]]
        if s["joined"]:
            dq, dr = _NEIGHBOR_GRID_DELTAS[s["side"]]
            seam = frozenset({s["at"], (s["at"][0] + dq, s["at"][1] + dr)})
            ax.plot(corners[:, 0], corners[:, 1], color=colour, linestyle=dash, linewidth=2.2, solid_capstyle="round", zorder=3)
            if seam in labelled_seams:
                continue
            labelled_seams.add(seam)
            size, weight, edge = 7, "normal", "white"
        else:
            ax.plot(corners[:, 0], corners[:, 1], color=colour, linestyle=dash, linewidth=5, solid_capstyle="round", zorder=4)
            size, weight, edge = 9, "bold", "black"
        mid = np.array(layout.junction_center(s["side"])) + origin
        ax.text(
            mid[0], mid[1], letter, ha="center", va="center", fontsize=size, fontweight=weight, color="black", zorder=6,
            bbox={"boxstyle": "circle,pad=0.25", "facecolor": colour, "edgecolor": edge, "linewidth": 1.2, "alpha": 0.95},
        )
    for p in tileset.preview_map:
        origin = layout.flower_grid_to_xy(*p.at)
        ax.text(origin[0], origin[1], p.id.replace("land_", ""), ha="center", va="center", fontsize=10, color="white",
                zorder=5, bbox={"boxstyle": "round,pad=0.2", "facecolor": "black", "alpha": 0.45, "edgecolor": "none"})

    bounds = mesh.bounds
    ax.set_xlim(bounds[0][0] - 8, bounds[1][0] + 8)
    ax.set_ylim(bounds[0][1] - 8, bounds[1][1] + 8)
    ax.set_aspect("equal")
    ax.set_axis_off()

    n_open = sum(1 for s in sides if not s["joined"])
    partnered = sum(1 for s in sides if counts[s["key"]] >= 2)
    open_partnered = sum(1 for s in sides if not s["joined"] and open_counts[s["key"]] >= 2)
    handles = []
    for key, letter in letters.items():
        profile, feature = key
        label = f"{letter}: {'-'.join(map(str, profile))}" + (f" + {feature}" if feature else "")
        label += f"   {counts[key]} sides ({open_counts[key]} open)"
        colour, dash = style(key)
        handles.append(plt.Line2D([], [], color=colour, linestyle=dash, linewidth=5, label=label))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8, frameon=False,
              title="Side classes: corner heights (levels, 15 mm each) + crossing feature", title_fontsize=9)
    ax.set_title(
        f"{args.tileset.name}: {len(sides)} flower sides, {len(letters)} mating classes. Equal letters fit together "
        f"(any flower may be rotated).\n{n_open} sides are open around the edge; {open_partnered} of them could also go "
        f"against another open side, and {partnered} of all {len(sides)} sides have a partner somewhere in the set. "
        "Thick = open side, thin = joined seam.",
        fontsize=10,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=130, bbox_inches="tight")
    print(f"wrote {args.out}")
    print(f"{len(letters)} classes over {len(sides)} sides; {n_open} open, {open_partnered} open sides with an open partner")
    for key, letter in letters.items():
        print(f"  {letter}: profile {key[0]} {key[1] or 'plain':6s} {counts[key]:3d} sides, {open_counts[key]} open")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
