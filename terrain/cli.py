"""CLI: render flower or preview assembly to SCAD."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from solid2 import set_global_fa, set_global_fn, set_global_fs

from terrain.assembly import AssemblyExporter
from terrain.render.format import format_render_spec
from terrain.tileset import TilesetError, load_tileset

DEFAULT_TILESET = Path(__file__).resolve().parent.parent / "tilesets" / "default.yaml"
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "output"


def _configure_resolution(mode: str) -> int:
    if mode == "print":
        resolution = 100
    else:
        resolution = 32
    set_global_fn(resolution)
    set_global_fa(resolution)
    set_global_fs(resolution)
    return resolution


def _export_scad(solid, path: Path, scale: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    solid.scale(scale).save_as_scad(str(path))
    print(f"wrote {path}")


def cmd_render_flower(args: argparse.Namespace) -> int:
    tileset_path = Path(args.tileset)
    resolution = _configure_resolution(args.resolution)
    try:
        tileset = load_tileset(tileset_path)
    except TilesetError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    exporter = AssemblyExporter(tileset, resolution=resolution)
    solid = exporter.build_flower(args.flower_id)
    out = Path(args.output) if args.output else OUTPUT_DIR / f"{args.flower_id}.scad"
    _export_scad(solid, out, tileset.meta.scale)
    return 0


def cmd_render_spec_flower(args: argparse.Namespace) -> int:
    tileset_path = Path(args.tileset)
    try:
        tileset = load_tileset(tileset_path)
    except TilesetError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    exporter = AssemblyExporter(tileset, resolution=32)
    spec = exporter.describe_flower(args.flower_id)
    text = format_render_spec(spec)
    out = Path(args.output) if args.output else OUTPUT_DIR / f"{args.flower_id}.render.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    print(f"wrote {out}")
    return 0


def cmd_render_preview(args: argparse.Namespace) -> int:
    tileset_path = Path(args.tileset)
    resolution = _configure_resolution(args.resolution)
    try:
        tileset = load_tileset(tileset_path)
    except TilesetError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    exporter = AssemblyExporter(tileset, resolution=resolution)
    solid = exporter.build_preview()
    out = Path(args.output) if args.output else OUTPUT_DIR / "preview.scad"
    _export_scad(solid, out, tileset.meta.scale)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="terrain", description="Hex-flower terrain SCAD export")
    sub = parser.add_subparsers(dest="command", required=True)

    render = sub.add_parser("render", help="Export SCAD")
    render_sub = render.add_subparsers(dest="render_target", required=True)

    flower = render_sub.add_parser("flower", help="Single flower by id")
    flower.add_argument("flower_id", help="Flower id from tileset")
    flower.add_argument(
        "--tileset",
        default=str(DEFAULT_TILESET),
        help="Path to tileset YAML",
    )
    flower.add_argument(
        "--resolution",
        choices=("preview", "print"),
        default="preview",
        help="Mesh resolution (preview=32, print=100)",
    )
    flower.add_argument("--output", "-o", help="Output .scad path")
    flower.set_defaults(func=cmd_render_flower)

    preview = render_sub.add_parser("preview", help="Combined preview_map assembly")
    preview.add_argument("--tileset", default=str(DEFAULT_TILESET))
    preview.add_argument(
        "--resolution",
        choices=("preview", "print"),
        default="preview",
    )
    preview.add_argument("--output", "-o")
    preview.set_defaults(func=cmd_render_preview)

    spec = render_sub.add_parser("spec", help="Text render spec for a flower (for tests)")
    spec_flower = spec.add_subparsers(dest="spec_target", required=True)
    spec_one = spec_flower.add_parser("flower", help="Single flower render spec")
    spec_one.add_argument("flower_id", help="Flower id from tileset")
    spec_one.add_argument("--tileset", default=str(DEFAULT_TILESET))
    spec_one.add_argument("--output", "-o", help="Output .json path")
    spec_one.set_defaults(func=cmd_render_spec_flower)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
