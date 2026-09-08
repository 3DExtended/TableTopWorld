# TableTop World — hex flower terrain

Printable **7-hex flower** terrain tiles for D&D / TTRPGs. Each flower is
declared in YAML (per-hex height levels, the height of every silhouette
corner, roads, rivers) and built as one watertight explicit mesh - no CSG,
no OpenSCAD - then written as STL. Neighbouring flowers mate along a
deterministic jagged silhouette and hold together with 5 × 2 mm disc
magnets.

The look the generator targets is in [`docs/mockups/`](docs/mockups/)
(`goal_*.png`, approved 2026-09): flat standable pads with S-curve steps
between height levels, rolling organic hexes, 16 mm roads sunk 1 mm into the
ground, 22 mm rivers 2.5 mm deep, and a thin V-line along every hex edge.

## Physical spec (millimetres)

| | |
|---|---|
| Hex | 45 mm across flats (circumradius 25.98 mm); flower about 130 × 135 mm |
| Height levels | 4 × 15 mm; level 0 at `z = 0`, print bed at `z = -10` (`BASE_PLATE_DEPTH_MM`) |
| Wall magnets | 18 blind bores per flower, one per silhouette edge at its midpoint, 5.3 × 2.2 mm, centred 3.9 mm above the bed |
| Top sockets | same 5.3 × 2.2 mm bore in the centre of every standable hex not crossed by a road or river |
| Hex lines | V-line 1.0 mm wide × 0.6 mm deep; each hex carries half of it, so the V forms only where two hexes (or two flowers) meet |
| Standable hexes | a seeded count between `min_standable_hexes` and 7; a flat pad covers 55% of the hex, a 6 mm S-curve band steps to the neighbours |
| Roads | 16 mm wide, 1 mm below the smoothed terrain, 0.6 mm edge |
| Rivers | 22 mm wide, 2.5 mm deep, 5 mm smoothstep banks; never lower than 1 mm above the wall bores |

Roads and rivers cross a flower side at the midpoint of the side's middle
edge, perpendicular to it, so the neighbour's copy continues exactly.

## Requirements

- **Python 3.11+** (3.14 works with the bundled venv)
- numpy, scipy, trimesh, mapbox_earcut, PyYAML (`requirements.txt`)
- A slicer or mesh viewer to look at the STLs

## Setup

```bash
cd tableTopWorld
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Always run commands with the venv active, or prefix with `.venv/bin/python`.

## Run — export STL

```bash
# One flower from tilesets/default.yaml → output/flower_<id>.stl
.venv/bin/python -m terrain.cli render flower crossroads --output output/flower_crossroads.stl

# The 19-flower landscape (tilesets/landscape.yaml) as one scene, to look at
.venv/bin/python -m terrain.cli render preview --tileset tilesets/landscape.yaml --output output/landscape.stl

# The same landscape as one STL per flower, to print: output/landscape/<id>.stl + README.md saying where each goes
.venv/bin/python -m terrain.cli render tileset --tileset tilesets/landscape.yaml --output-dir output/landscape

# The sample-flower preview scene (every entry in default.yaml's preview_map)
.venv/bin/python -m terrain.cli render preview --output output/preview.stl

# Coarser or finer mesh (default 16 subdivisions per hex edge; a road edge is one lattice step wide)
.venv/bin/python -m terrain.cli render flower hill_peak --subdivisions-per-edge 8 --output output/flower_hill_peak.stl
```

## Print — Bambu Lab project file

With Bambu Studio installed, `scripts/bambu_project.py` turns a folder of
flower STLs into a sliced project (`.gcode.3mf`) using Bambu Studio's own
command line: the flowers are arranged two per 256 mm plate, every plate is
sliced, and plate previews are written next to it. Open the file in Bambu
Studio and print plate by plate (or change the process or filament there
and re-slice).

```bash
# P1S, 0.4 nozzle, 0.20 mm Standard, Generic PLA, textured PEI plate (the defaults)
.venv/bin/python scripts/bambu_project.py --stl-dir output/hills --output output/hills/hills_P1S.gcode.3mf

# Another printer or filament: any profile name from Bambu Studio's bundled BBL profiles
.venv/bin/python scripts/bambu_project.py --stl-dir output/hills --output output/hills/hills_A1.gcode.3mf --printer "Bambu Lab A1 0.4 nozzle" --process "0.20mm Standard @BBL A1" --filament "Bambu PLA Basic @BBL A1"
```

The script flattens the chosen profiles along their `inherits` chain first;
Bambu Studio's CLI does not do that itself and would otherwise slice for a
200 mm bed.

Every export is gated by `terrain/export.py`: the mesh must be watertight,
winding-consistent and of positive volume, or nothing is written. The CLI
also prints how many of the 7 hexes are standable in the real mesh.

`scripts/render_mesh.py` renders a PNG (top-down and isometric, light-shaded)
of an STL for a quick visual check without a GPU. `scripts/annotate_sides.py`
draws a tileset's preview from above with every flower side coloured by its
mating class (corner-height profile plus road/river crossing): equal letters
are sides that could be placed against each other.

## Tilesets

Flower definitions live in **`tilesets/default.yaml`**:

- **`meta`** — `hex_outer_width`, `scale`, `heights` (level count and mm per level), `min_standable_hexes`
- **`flowers.<id>`** — `seed`, seven `hexes` with a `height_level`, six `side_corner_heights` (4 levels per side), optional `roads` / `water` as `[entry_side, exit_side]` or `{entry, exit, via: [hex, ...]}`
- **`preview_map`** — which flowers to place where in `render preview`

Two flowers may share a side only if their corner-height sequences match
once reversed; `load_tileset()` validates that for every adjacent pair in
`preview_map`, that a road or river leaving a flower continues in its
neighbour, and that adjacent sides of one flower agree on their shared
corner. Sample flowers: `flat_plains`, `hill_peak` (a cliff on one side),
`crossroads` (a road and a river that ford at the centre hex),
`river_bend` (a river through level-0 ground plus a road).

**`tilesets/landscape.yaml`** is a 19-flower landscape (a hexagon of
flowers, radius 2 on the flower grid) written by `scripts/gen_landscape.py`:
every hex level and corner height is sampled from one continuous elevation
field, so all 42 seams match by construction. A river meanders from the west
edge to the east edge through a level-0 valley, hills rise to level 3 in the
north-east and west, and a road runs south to north and fords the river in
the centre column. Edit the script (hills, valley, the river and road
routes) and rerun it; a test checks the committed YAML is what the script
produces. `render tileset` writes the 19 printable STLs into a folder.
**`tilesets/hills.yaml`** (`--preset hills`) is the same hexagon on three
levels only: valley, plains and level-2 hills that straddle seams and
three-flower corners, with the road entering from the south-east.

## How a flower is built

1. `terrain/heightfield.py` — the deterministic jagged contour of each of the
   6 sides (`terrain/boundary_noise.py`, bit-identical between the two
   flowers that share it), then a triangular lattice per hex cell with a
   sunken half-V strip along every edge and a top socket where allowed.
2. `terrain/field.py` — the height of every interior vertex: level blend with
   flat pads and S-curve bands, wobble, noise and organic relief, then roads
   and rivers carved along polylines with rounded bends.
3. `terrain/surface_mesh.py` + `terrain/magnets.py` — walls straight down to
   the bed with the 18 blind bores cut in as structured collars.
4. `terrain/base_plate.py` — the floor, a mirror of the top footprint.
5. `terrain/assembly.py` — one `trimesh.Trimesh`; `terrain/standability.py`
   measures which hexes are flat within 1 mm over their central 40%.

## Test

```bash
.venv/bin/python -m pytest -q
```

About 100 tests, all against real constructed geometry: watertightness and
winding across mesh resolutions, the cross-flower boundary contract
(point-by-point along a shared side), bore and socket positions and
volumes, pad flatness, road and river profiles, road continuity across
the crossroads/river_bend seam, and the landscape's 42 seams and 9
road/river crossings.

## Project layout

```
tableTopWorld/
├── terrain/           # Generator (layout, boundary noise, heightfield, field, magnets, assembly, cli)
├── tilesets/          # YAML flower definitions
├── tests/             # pytest
├── scripts/           # render_mesh.py and prototype scripts
├── docs/              # PRD, plan, approved mockups
└── output/            # Generated STLs (local)
```

## Docs

- [`docs/prd/3d-hex-flower-terrain.md`](docs/prd/3d-hex-flower-terrain.md) — product requirements
- [`docs/plans/3d-hex-flower-terrain.md`](docs/plans/3d-hex-flower-terrain.md) — implementation summary
- [`docs/mockups/`](docs/mockups/) — the approved target images and the script that drew them
