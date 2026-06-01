# Visual atom tests

Atomic unit tests build the terrain pipeline from small pieces to full flowers. Each test writes an OpenSCAD file under `output/visual_atoms/` so you can compare the mesh to the test description.

## Where to read each test’s description

| Source | What you get |
|--------|----------------|
| **This file — [Per-test catalog](#per-test-catalog)** | ID, SCAD path, one-line expectation (best overview) |
| **SCAD file header** | First `//` lines in `output/visual_atoms/<layer>/NN_*.scad` after running tests |
| **Pytest** | `pytest tests/atoms/ -v` shows each test’s docstring; or open `tests/atoms/test_*.py` |

Numbering note: **test `03`** (junction markers) is in **layer folder `01_layout/`**, not `03_edges/`. Folder `03_edges/` is **layer 3** (magnets/bevels), tests **13–16**.

## Run

From the repo root (see [main README](../../README.md) for venv setup):

```bash
.venv/bin/python -m pytest tests/atoms/ -v
```

OpenSCAD files appear in `output/visual_atoms/<layer>/NN_name.scad`.

All tests (unit + visual): `.venv/bin/python -m pytest tests/ -q`

## Layers (build order)

| Layer | Directory | What it proves |
|-------|-----------|----------------|
| 1 | `01_layout/` | Hex footprints, junctions, exterior edges (geometry only) |
| 2 | `02_mesh/` | Single hex prisms, slopes, topping holes, 7-hex union |
| 3 | `03_edges/` | One magnet, one bevel, all magnets, bevels on mesh |
| 4 | `04_features/` | Road/water cutter tools and cuts on slab/plateau |
| 5 | `05_compose/` | Full pipeline on atoms + production tileset flowers |

## Per-test catalog

### Layer 1 — `01_layout/` (no terrain mesh)

| Test | SCAD | Expect in OpenSCAD |
|------|------|-------------------|
| 01 | `01_center_hex_footprint.scad` | One thin center hex plate (subdivided 18-vertex outline). |
| 02 | `02_seven_footprints_stacked.scad` | Seven hex outlines stacked in Z; ring around center, no gaps. |
| 03 | `03_junction_centers.scad` | **Six vertical pegs only** — see [Why test 03 looks weird](#why-test-03-looks-weird). |
| 04 | `04_exterior_edge_wires.scad` | 18 red bars on outward ring-hex sides only (midpoint farther from center than hex center). |

### Layer 2 — `02_mesh/`

| Test | SCAD | Expect in OpenSCAD |
|------|------|-------------------|
| 05 | `05_hex_prism_ground.scad` | One ground hex slab from z=0 to z=2 (no magnet collar). |
| 06 | `06_hex_prism_middle.scad` | One taller center hex, top at Z=4. |
| 07 | `07_hex_prism_high.scad` | One tall hex, top at Z=8. |
| 08 | `08_slope_ramp.scad` | Ground hex plus wedge ramp up to neighbor height. |
| 09 | `09_topping_hole_tools.scad` | Two vertical cylinders (hole tools, not a printable part). |
| 10 | `10_topping_subtracted.scad` | One hex with circular dimple on top. |
| 11 | `11_seven_hex_union.scad` | Seven hexes fused; center middle, hex 3 high, rest ground. |
| 12 | `12_hill_north_slope_hex5.scad` | Production hill_north slope hex (compare to 08). |

### Layer 3 — `03_edges/` (folder name, not test number)

| Test | SCAD | Expect in OpenSCAD |
|------|------|-------------------|
| 13 | `13_single_magnet.scad` | Ground hex slab z=0–2; fixed-height magnet bore on edge 1-0. |
| 14 | `14_single_bevel.scad` | High hex 3; chamfer on edge 3-0 — side wedge, sloped top wedge, and flat top box (three cutters). |
| 15 | `15_all_magnets_ground_plateau.scad` | Flat 7-hex plate, 18 magnet holes on outer rim. |
| 16 | `16_bevels_on_height_step.scad` | Height-step union with chamfer on all seven hex cells (42 sides). |

### Layer 4 — `04_features/`

| Test | SCAD | Expect in OpenSCAD |
|------|------|-------------------|
| 17 | `17_road_cut_tools.scad` | Curved road cutter solid only (junction 0 → 4). |
| 18 | `18_road_on_block.scad` | Big slab with groove matching 17. |
| 19 | `19_road_on_plateau.scad` | Flat 7-hex with road channel. |
| 20 | `20_water_cut_tools.scad` | Water cutter (hex profile), tools only. |
| 21 | `21_water_on_plateau.scad` | Flat 7-hex with water channel (deeper/wider than road). |

### Layer 5 — `05_compose/`

| Test | SCAD | Expect in OpenSCAD |
|------|------|-------------------|
| 22 | `22_height_step_bevels.scad` | 11 + bevels, no magnets. |
| 23 | `23_height_step_full_edges.scad` | 22 + 18 magnets. |
| 24 | `24_slope_atom_full.scad` | Full pipeline on slope atom. |
| 25 | `25_road_atom_full.scad` | Full flat flower with road. |
| 26 | `26_flat_plains.scad` | Production flat_plains tile. |
| 27 | `27_hill_north.scad` | Production hill_north tile. |
| 28 | `28_river_grove.scad` | Production river_grove + water. |
| 29 | `29_preview_two_flowers.scad` | Entire default preview_map assembly. |
| 30 | *(no SCAD)* | Writes `30_hill_north.render.json` text spec only. |

## Why test 03 looks weird

`03_junction_centers.scad` is **intentionally not a flower**. It only draws six small cylinders at the **road/water junction points** — the spots where three hex edges meet on the ring (where bezier roads start/end).

You see six pegs floating in space because there is **no hex geometry** in that file. That is a layout debug view, not a printable piece.

To make sense of it:

1. Open `02_seven_footprints_stacked.scad` first (the seven hex outlines).
2. Then `03_junction_centers.scad` — each peg should sit on a **outer** corner where three cells meet (roughly on the ring between center and rim).
3. Roads in later tests (17+) run between pairs of these junctions.

If pegs look wrong, the bug is in `terrain/layout.py` junction math, not in mesh or features.

## How to verify

1. Run the tests (generates SCAD).
2. Open the matching `NN_*.scad` in OpenSCAD.
3. Read the `//` comment at the top of the file — it matches the pytest docstring.
4. If the shape looks wrong, find the first layer in the table above where expectation fails.

## Composition chain

```
layout footprints → single hex mesh → 7-hex union
  → + bevels → + magnets → + road/water → production flowers → preview map
```

Examples:

- `08_slope_ramp` + `16_bevels` → `24_slope_atom_full`
- `17_road_cut_tools` → `19_road_on_plateau` → `25_road_atom_full` → `26_flat_plains`
