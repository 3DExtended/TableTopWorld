"""Command-line entry point - written fresh for the explicit-mesh redesign
(the old terrain/cli.py, CSG-based, was deleted in Phase 4 along with the
rest of the discarded pipeline; see the project plan's Phase 4 note).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from terrain.assembly import build_flower_mesh, build_preview_mesh, standability_report
from terrain.export import MeshValidationError, export_stl
from terrain.tileset import Tileset, TilesetError, load_tileset


def _cmd_render_flower(args: argparse.Namespace) -> int:
    try:
        tileset = load_tileset(args.tileset)
    except TilesetError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if args.flower_id not in tileset.flowers:
        print(
            f"error: flower {args.flower_id!r} not found in {args.tileset} "
            f"(available: {', '.join(sorted(tileset.flowers))})",
            file=sys.stderr,
        )
        return 1

    mesh = build_flower_mesh(
        tileset, args.flower_id, subdivisions_per_edge=args.subdivisions_per_edge
    )
    ok, count, min_required = standability_report(mesh, tileset)
    print(
        f"flower {args.flower_id!r}: {count}/7 hex cells standable "
        f"(minimum {min_required}, {'OK' if ok else 'BELOW MINIMUM'})"
    )

    try:
        out_path = export_stl(mesh, args.output)
    except MeshValidationError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"wrote {out_path}")
    return 0


def _cmd_render_preview(args: argparse.Namespace) -> int:
    try:
        tileset = load_tileset(args.tileset)
    except TilesetError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if not tileset.preview_map:
        print(f"error: {args.tileset} has an empty preview_map", file=sys.stderr)
        return 1

    mesh = build_preview_mesh(tileset, subdivisions_per_edge=args.subdivisions_per_edge)
    print(f"preview: {len(tileset.preview_map)} flower(s) placed")

    try:
        out_path = export_stl(mesh, args.output)
    except MeshValidationError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"wrote {out_path}")
    return 0


def _placement_readme(tileset: Tileset, tileset_path: Path, rows: list[dict]) -> str:
    """A README for a folder of per-flower STLs: which file goes where."""
    lines = [
        f"# Flowers of `{tileset_path.name}`",
        "",
        "One STL per flower, each a complete printable solid (base plate, 18 wall",
        "magnet bores, top sockets). Print every file once and place them on the",
        "flower grid below; the two files sharing a seam mate edge for edge.",
        "",
        "| file | grid (q, r) | hex levels 0..6 | roads | water | standable |",
        "|---|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| `{row['file']}` | {row['at']} | {row['levels']} | {row['roads']} | {row['water']} | {row['standable']}/7 |"
        )
    placed = {p.at: p.id for p in tileset.preview_map}
    if placed:
        qs = sorted({q for q, _ in placed})
        rs = sorted({r for _, r in placed})
        lines += [
            "",
            "Placement, seen from above with north up: q counts east, r counts north, and every",
            "row is a third of a flower further east than the row south of it. Side 0 of a flower",
            "faces east and sides count anti-clockwise; side k touches the flower at grid delta",
            "(1,0), (0,1), (-1,1), (-1,0), (0,-1), (1,-1) for k = 0..5.",
            "",
            "```",
        ]
        slot = max(len(v) for v in placed.values()) + 1
        for r in reversed(rs):
            indent = " " * round((r - min(rs)) * slot / 3)
            cells = [placed.get((q, r), "").ljust(slot) for q in qs]
            lines.append((indent + "".join(cells)).rstrip())
        lines.append("```")
    return "\n".join(lines) + "\n"


def _cmd_render_tileset(args: argparse.Namespace) -> int:
    try:
        tileset = load_tileset(args.tileset)
    except TilesetError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if tileset.preview_map and not args.all:
        placements = [(p.id, p.at) for p in tileset.preview_map]
    else:
        placed = {p.id: p.at for p in tileset.preview_map}
        placements = [(fid, placed.get(fid)) for fid in tileset.flowers]

    out_dir: Path = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for flower_id, at in placements:
        mesh = build_flower_mesh(tileset, flower_id, subdivisions_per_edge=args.subdivisions_per_edge)
        ok, count, min_required = standability_report(mesh, tileset)
        try:
            out_path = export_stl(mesh, out_dir / f"{flower_id}.stl")
        except MeshValidationError as exc:
            print(f"error: {flower_id}: {exc}", file=sys.stderr)
            return 1
        flower = tileset.flowers[flower_id]
        rows.append(
            {
                "file": out_path.name,
                "at": f"({at[0]}, {at[1]})" if at is not None else "-",
                "levels": " ".join(str(flower.hexes[str(h)].height_level) for h in range(7)),
                "roads": ", ".join(f"{r.entry}->{r.exit}" for r in flower.roads) or "-",
                "water": ", ".join(f"{w.entry}->{w.exit}" for w in flower.water) or "-",
                "standable": count,
            }
        )
        print(
            f"wrote {out_path} ({count}/7 standable, minimum {min_required}, {'OK' if ok else 'BELOW MINIMUM'})"
        )
    readme = out_dir / "README.md"
    readme.write_text(_placement_readme(tileset, args.tileset, rows))
    print(f"wrote {readme} ({len(rows)} flowers)")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="terrain")
    subparsers = parser.add_subparsers(dest="command", required=True)

    render = subparsers.add_parser("render", help="render flowers to STL")
    render_sub = render.add_subparsers(dest="target", required=True)

    flower = render_sub.add_parser("flower", help="render a single flower")
    flower.add_argument("flower_id", help="flower id as declared in the tileset")
    flower.add_argument(
        "--tileset",
        type=Path,
        default=Path("tilesets/default.yaml"),
        help="path to the tileset YAML (default: tilesets/default.yaml)",
    )
    flower.add_argument(
        "--output",
        type=Path,
        default=Path("output/flower.stl"),
        help="output STL path (default: output/flower.stl)",
    )
    flower.add_argument(
        "--subdivisions-per-edge",
        type=int,
        default=16,
        help="mesh detail along each hex edge, higher = smoother/more triangles (default: 16)",
    )
    flower.set_defaults(func=_cmd_render_flower)

    preview = render_sub.add_parser(
        "preview", help="render the tileset's whole preview_map as one assembled scene"
    )
    preview.add_argument(
        "--tileset",
        type=Path,
        default=Path("tilesets/default.yaml"),
        help="path to the tileset YAML (default: tilesets/default.yaml)",
    )
    preview.add_argument(
        "--output",
        type=Path,
        default=Path("output/preview.stl"),
        help="output STL path (default: output/preview.stl)",
    )
    preview.add_argument(
        "--subdivisions-per-edge",
        type=int,
        default=16,
        help="mesh detail along each hex edge, higher = smoother/more triangles (default: 16)",
    )
    preview.set_defaults(func=_cmd_render_preview)

    every = render_sub.add_parser(
        "tileset", help="render every flower of a tileset as its own STL into a folder, plus a placement README"
    )
    every.add_argument(
        "--tileset",
        type=Path,
        default=Path("tilesets/default.yaml"),
        help="path to the tileset YAML (default: tilesets/default.yaml)",
    )
    every.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/flowers"),
        help="folder for <flower id>.stl files and README.md (default: output/flowers)",
    )
    every.add_argument(
        "--subdivisions-per-edge",
        type=int,
        default=16,
        help="mesh detail along each hex edge, higher = smoother/more triangles (default: 16)",
    )
    every.add_argument(
        "--all",
        action="store_true",
        help="also render flowers that are declared but not placed in preview_map",
    )
    every.set_defaults(func=_cmd_render_tileset)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
