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
from terrain.tileset import TilesetError, load_tileset


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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="terrain")
    subparsers = parser.add_subparsers(dest="command", required=True)

    render = subparsers.add_parser("render", help="render one flower to STL")
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
        default=8,
        help="mesh detail along each hex edge, higher = smoother/more triangles (default: 8)",
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
        default=8,
        help="mesh detail along each hex edge, higher = smoother/more triangles (default: 8)",
    )
    preview.set_defaults(func=_cmd_render_preview)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
