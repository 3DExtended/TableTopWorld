"""Slice a folder of flower STLs into a Bambu Studio project (.gcode.3mf)
with Bambu Studio's own command line: the flowers are arranged onto as
many plates as they need, every plate is sliced, and the project opens in
Bambu Studio ready to print plate by plate.

Bambu Studio's bundled profiles inherit from each other (`inherits`), which
its CLI does not resolve - a machine profile loaded raw comes out with a
200 x 200 bed. So the chosen machine, process and filament profiles are
flattened along their inheritance chain into self-contained JSON first.

Usage (after `render tileset` wrote the STLs):

    .venv/bin/python scripts/bambu_project.py --stl-dir output/hills --output output/hills/hills_P1S.gcode.3mf

Defaults are a P1S with the 0.4 mm nozzle, 0.20 mm Standard, Generic PLA
on the textured PEI plate; see --help for the knobs. Plate previews land in
<output dir>/plates/.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

DEFAULT_STUDIO = Path("/Applications/BambuStudio.app")


def flatten_profile(profiles: Path, kind: str, name: str) -> dict:
    path = profiles / kind / f"{name}.json"
    if not path.exists():
        raise SystemExit(f"no {kind} profile {name!r} in {profiles / kind}")
    data = json.loads(path.read_text())
    parent = data.pop("inherits", None)
    merged = flatten_profile(profiles, kind, parent) if parent else {}
    merged.update(data)
    return merged


def run_cli(studio: Path, args: list[str], log: Path) -> None:
    binary = studio / "Contents" / "MacOS" / "BambuStudio"
    with log.open("w") as handle:
        result = subprocess.run([str(binary), "--debug", "2", *args], stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        tail = "".join(log.read_text().splitlines(keepends=True)[-15:])
        raise SystemExit(f"Bambu Studio exited with {result.returncode}; see {log}\n{tail}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--stl-dir", type=Path, required=True, help="folder of <flower>.stl files")
    parser.add_argument("--output", type=Path, required=True, help="the .gcode.3mf to write")
    parser.add_argument("--printer", default="Bambu Lab P1S 0.4 nozzle")
    parser.add_argument("--process", default="0.20mm Standard @BBL X1C")
    parser.add_argument("--filament", default="Generic PLA")
    parser.add_argument("--bed-type", default="Textured PEI Plate")
    parser.add_argument("--studio", type=Path, default=DEFAULT_STUDIO, help="Bambu Studio app bundle")
    args = parser.parse_args(argv)

    stls = sorted(args.stl_dir.glob("*.stl"))
    if not stls:
        raise SystemExit(f"no STL files in {args.stl_dir}")
    profiles = args.studio / "Contents" / "Resources" / "profiles" / "BBL"
    if not profiles.is_dir():
        raise SystemExit(f"Bambu Studio profiles not found under {args.studio}")

    with tempfile.TemporaryDirectory(prefix="bambu_project_") as tmp_name:
        tmp = Path(tmp_name)
        machine = tmp / "machine.json"
        process = tmp / "process.json"
        filament = tmp / "filament.json"
        machine.write_text(json.dumps(flatten_profile(profiles, "machine", args.printer)))
        process_data = flatten_profile(profiles, "process", args.process)
        process_data["curr_bed_type"] = args.bed_type
        process.write_text(json.dumps(process_data))
        filament.write_text(json.dumps(flatten_profile(profiles, "filament", args.filament)))

        # --export-3mf is joined onto --outputdir, so give a bare name.
        run_cli(
            args.studio,
            [
                "--load-settings", f"{machine};{process}",
                "--load-filaments", str(filament),
                "--arrange", "1", "--orient", "0", "--ensure-on-bed",
                "--slice", "0",
                "--export-3mf", "project.gcode.3mf",
                "--outputdir", str(tmp),
                *map(str, stls),
            ],
            tmp / "slice.log",
        )
        png_dir = tmp / "png"
        png_dir.mkdir()
        run_cli(args.studio, ["--export-png", "0", "--outputdir", str(png_dir), str(tmp / "project.gcode.3mf")], tmp / "png.log")

        args.output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(tmp / "project.gcode.3mf", args.output)
        plates_dir = args.output.parent / "plates"
        plates_dir.mkdir(exist_ok=True)
        for png in sorted(png_dir.glob("plate_*_0.png")):
            shutil.copyfile(png, plates_dir / png.name.replace("_0.png", ".png"))
        result = json.loads((tmp / "result.json").read_text())

    print(f"wrote {args.output} for {args.printer}: {args.process}, {args.filament}, {args.bed_type}")
    total_h = total_g = 0.0
    for plate in result.get("sliced_plates", []):
        names = ", ".join(o["name"].removesuffix(".stl") for o in plate["objects"])
        grams = sum(f["total_used_g"] for f in plate["filaments"])
        hours = plate["total_predication"] / 3600
        total_h += hours
        total_g += grams
        warn = f"  WARNING: {plate['warning_message']}" if plate.get("warning_message") else ""
        print(f"  plate {plate['id']:2d}: {hours:4.1f} h {grams:5.0f} g  {names}{warn}")
    print(f"  total: {len(result.get('sliced_plates', []))} plates, {total_h:.1f} h, {total_g:.0f} g")
    return 0


if __name__ == "__main__":
    sys.exit(main())
