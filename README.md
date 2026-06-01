# TableTop World — hex flower terrain

Parametric **7-hex flower** terrain tiles for 3D printing. Each flower is defined in YAML (heights, slopes, roads, water, edge profiles), built as CSG in Python ([solidpython2](https://github.com/WolfgangFahl/solidpython2)), and exported to OpenSCAD (`.scad`).

## Requirements

- **Python 3.11+** (3.14 works with the bundled venv)
- **[OpenSCAD](https://openscad.org/)** — to preview `.scad` files (F6 render / F5 preview)
- Optional: `entr` for auto-reexport while editing legacy scripts (see below)

## Setup

```bash
cd tableTopWorld
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Always run commands with the venv active, or prefix with `.venv/bin/python`.

## Run — export SCAD

### CLI (recommended)

From the repo root:

```bash
# One flower from tilesets/default.yaml → output/<id>.scad
.venv/bin/python -m terrain.cli render flower flat_plains

# Preview assembly (all entries in preview_map)
.venv/bin/python -m terrain.cli render preview

# Higher mesh density for slicing (fn/fa/fs = 100)
.venv/bin/python -m terrain.cli render flower hill_north --resolution print

# Custom tileset or output path
.venv/bin/python -m terrain.cli render flower river_grove \
  --tileset tilesets/default.yaml \
  -o output/river_grove.scad

# Text render spec (JSON) for unit tests — no mesh CSG
.venv/bin/python -m terrain.cli render spec flower hill_north
```

| Command | Output (default) |
|---------|------------------|
| `render flower <id>` | `output/<id>.scad` |
| `render preview` | `output/preview.scad` |
| `render spec flower <id>` | `output/<id>.render.json` |

Models are written in **model units**, then scaled by `meta.scale` from the tileset (default **5×** for mm). OpenSCAD files include BOSL2 paths from the solidpython2 install.

### Legacy script

```bash
.venv/bin/python printableFiles/hexagon.py
```

Writes `printableFiles/hexagon.scad` (production resolution, `flat_plains` only).

Watch mode (if you use `entr`):

```bash
ls printableFiles/*.py | entr .venv/bin/python printableFiles/hexagon.py
```

## Tilesets

Flower definitions live in **`tilesets/default.yaml`**:

- **`meta`** — `height_step_mm`, `scale`, `hex_outer_width`
- **`flowers.<id>`** — seven `hexes` (terrain level + role), eighteen `edges` (profile names), optional `roads` / `water`
- **`preview_map`** — which flowers to place in `render preview`

Terrain Z heights: `ground` = 0, `middle` = 4, `high` = 8 (before `meta.scale`). All hexes share a print bed at `FLOWER_BOTTOM_Z` = −2 mm in model units.

## Test

### Quick check (all automated tests)

```bash
.venv/bin/python -m pytest tests/ -q
```

About **85 tests**: layout math, mesh, edges (magnets/bevels), features, tileset loading, render specs, and visual atoms.

Verbose:

```bash
.venv/bin/python -m pytest tests/ -v
```

Run a subset:

```bash
.venv/bin/python -m pytest tests/test_layout.py tests/test_bevel.py -v
.venv/bin/python -m pytest tests/atoms/ -v          # visual atom layer only
```

### Two kinds of tests

| Kind | Location | What it checks | Output |
|------|----------|----------------|--------|
| **Unit** | `tests/test_*.py` (except `tests/atoms/`) | Numbers, specs, SCAD structure strings | None |
| **Visual atoms** | `tests/atoms/test_*.py` | Pipeline slices + writes OpenSCAD | `output/visual_atoms/` |

**Unit tests** avoid parsing full SCAD where possible. **`tests/test_render_spec.py`** compares a structured **`FlowerRenderSpec`** (JSON) from `terrain render spec flower …` — hex tops, magnet count, road junctions, etc.

**Visual atoms** build geometry in **five layers** (layout → mesh → edges → features → full compose). Each test exports one `.scad` file so you can eyeball regressions in OpenSCAD.

### Visual atoms — interpret results

1. Run: `.venv/bin/python -m pytest tests/atoms/ -v`
2. Open the matching file under **`output/visual_atoms/<layer>/NN_*.scad`**
3. Read the **`//` comment** at the top of the file (matches the pytest docstring)
4. Use the **per-test catalog** in [`tests/visual_atoms/README.md`](tests/visual_atoms/README.md) for expected shapes

**Layer folders** (build order):

| Folder | Proves |
|--------|--------|
| `01_layout/` | Footprints, junctions, exterior edge lines (no solid terrain) |
| `02_mesh/` | Hex prisms, slopes, topping holes, 7-hex union |
| `03_edges/` | Magnets and bevels on simple meshes |
| `04_features/` | Road/water cutters and cuts |
| `05_compose/` | Full pipeline on atoms and production flowers |

**Numbering quirk:** pytest **test 03** (junction pegs) lives in **`01_layout/`**, not `03_edges/`. Folder `03_edges/` is layer 3 (tests 13–16).

**If a visual test fails:** find the **first layer** where the shape wrong; fix upstream (e.g. layout before bevels). Composition chain:

```
layout → mesh → + bevels → + magnets → + road/water → production flowers → preview
```

**Test 03** intentionally shows only six junction cylinders — not a printable tile. Compare with `02_seven_footprints_stacked.scad` first.

### Regenerate production outputs

Tests do not overwrite `output/flat_plains.scad` unless you run compose tests or the CLI. To refresh production SCAD:

```bash
.venv/bin/python -m terrain.cli render flower flat_plains
.venv/bin/python -m terrain.cli render flower hill_north
.venv/bin/python -m terrain.cli render preview
```

## Project layout

```
tableTopWorld/
├── terrain/           # Core library (layout, mesh, edges, features, assembly, render, cli)
├── tilesets/          # YAML flower definitions
├── tests/             # pytest (unit + tests/atoms/ visual suite)
├── output/            # Generated .scad and .render.json (gitignored or local)
├── printableFiles/    # Legacy hexagon.py export
└── docs/              # PRD and design notes
```

## OpenSCAD tips

- Open generated `.scad` from the repo root so **include paths** resolve.
- Use **Preview (F5)** for speed; **Render (F6)** before export to STL.
- If the mesh looks faceted, re-export with `--resolution print`.
- Z fighting in preview often means coplanar faces; the atom tests for bevels (`14_single_bevel.scad`) are the reference for chamfer + top rim cuts.

## Docs

- [`docs/prd/3d-hex-flower-terrain.md`](docs/prd/3d-hex-flower-terrain.md) — product requirements
- [`docs/plans/3d-hex-flower-terrain.md`](docs/plans/3d-hex-flower-terrain.md) — implementation plan
- [`tests/visual_atoms/README.md`](tests/visual_atoms/README.md) — full visual test catalog
