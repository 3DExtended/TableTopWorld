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

# The preview scene (every entry in preview_map, placed on the flower grid)
.venv/bin/python -m terrain.cli render preview --output output/preview.stl

# Finer mesh (default 8 subdivisions per hex edge)
.venv/bin/python -m terrain.cli render flower hill_peak --subdivisions-per-edge 16 --output output/flower_hill_peak.stl
```

Every export is gated by `terrain/export.py`: the mesh must be watertight,
winding-consistent and of positive volume, or nothing is written. The CLI
also prints how many of the 7 hexes are standable in the real mesh.

`scripts/render_mesh.py` renders a PNG (top-down and isometric, light-shaded)
of an STL for a quick visual check without a GPU.

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

About 90 tests, all against real constructed geometry: watertightness and
winding across mesh resolutions, the cross-flower boundary contract
(point-by-point along a shared side), bore and socket positions and
volumes, pad flatness, road and river profiles, and road continuity across
the crossroads/river_bend seam.

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
