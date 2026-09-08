# PRD: 3D modular hex-flower terrain system

**Issue:** https://github.com/3DExtended/TableTopWorld/issues/1
**Labels:** `enhancement`, `ready-for-agent`
**Plan:** [docs/plans/3d-hex-flower-terrain.md](../plans/3d-hex-flower-terrain.md)

> **Superseded design note:** this PRD originally described a CSG/solid2/BOSL2
> pipeline with a named edge-profile catalog (`flat_ground`, `slope_up_*`,
> `cliff_*`, ...). That pipeline broke once variable terrain height was
> introduced (slopes were faked as stacked prisms, road/water cuts used one
> flat Z per flower instead of per-hex Z) and was fully replaced by an
> **explicit-mesh** generator (heightfield-driven vertices/faces via
> `trimesh`, no CSG booleans). This document has been rewritten to describe
> what was actually built; the full decision log and phase-by-phase history
> live in [`.claude/plans/elegant-twirling-scott.md`](../../.claude/plans/elegant-twirling-scott.md).

---

## Problem Statement

TableTopWorld generates flat 7-hex **flower** tiles with carved roads, water
channels, and magnets — but all terrain lived on a single height plane, and
the CSG-based generator could not cleanly support variable per-hex height.
The project needed **3D modular terrain** where hills and cliffs connect
across flowers while preserving the hex grid for miniature placement, driven
by a declarative tileset rather than hand-edited Python constants.

## Solution

A **declarative tileset** (YAML) describes each flower's seven hex cells (a
discrete height level each) and, for each of the flower's 6 physical sides, a
4-corner height sequence. Two flowers may only share a side if their 4
corner-height sequences match numerically (the second flower's sequence must
be the **reverse** of the first's, since the two flowers walk the shared edge
in opposite directions) — there is no named profile catalog; matching is
just equal numbers. The flower's fine, jagged boundary contour between those
corners is a deterministic hash-based noise function of the corner heights
and position only, so two independently-built flowers with matching corner
heights produce bit-identical shared geometry with nothing else exchanged.

The generator (`terrain/`) builds each flower as an explicit heightfield
mesh (`trimesh.Trimesh`) — no boolean CSG anywhere in the pipeline — welds a
flat base plate with magnet bores' worth of clearance underneath, and
exports STL through a pre-export gate that refuses to write a mesh that
isn't watertight, winding-consistent, and positive-volume. Every test
asserts against real constructed geometry (a `Trimesh` or a written/reloaded
`.stl`), not generated source text, closing the gap that let the old
pipeline's breakage go undetected for a long time.

## User Stories

1. As a terrain designer, I want to define flowers in a tileset YAML file, so that I can add new tiles without editing generator code.
2. As a terrain designer, I want each hex cell to carry a discrete declared height level, so that per-cell terrain intent is explicit and validated.
3. As a terrain designer, I want "≥1 standable hex" checked against the real generated mesh after building, so that a flower's authored intent and its actual printed geometry can't silently diverge — and so a flower may validly have 0 standable hexes when that's the deliberate design (e.g. a cliff on every side).
4. As a terrain designer, I want a configurable number of discrete height levels (default 4, 15mm each), so that hills and cliffs read clearly at print scale.
5. As a terrain designer, I want no cap on the height delta between adjacent boundary corners, so that intentional cliffs are supported, not just gentle slopes.
6. As a terrain designer, I want roads carved as channels that follow each hex's real terrain Z, so that roads follow hills instead of only existing on flat ground.
7. As a terrain designer, I want each flower's 6 physical sides described by a 4-corner height sequence, so that two flowers may mate only when the numbers genuinely match.
8. As a terrain designer, I want the shared-edge jaggedness between two matching flowers to be bit-identical, computed from nothing but the declared corner heights and position, so that no extra data has to be exchanged for tiles to fit physically.
9. As a terrain designer, I want the generator to reject tilesets where preview-map neighbors have mismatched corner-height sequences, so that layout mistakes fail at build time, not at the printer.
10. As a terrain designer, I want road paths defined by junction pairs (same convention as before), so that existing road-authoring intuition carries over.
11. As a terrain designer, I want every one of a flower's 7 hex cells to get an engraved groove outlining it, following the real terrain surface however steep, so that hex cells read clearly for gameplay regardless of slope.
12. As a terrain designer, I want a manual preview_map listing flower IDs, grid positions, and rotations, so that I can see a deliberate test layout for the whole tileset.
13. As a terrain designer, I want a combined STL export of the preview_map, so that I can visually judge how tileset edits affect the full assembled layout.
14. As a maker, I want the printable unit to remain a single 7-hex flower plate with edge magnets, so that my existing magnet workflow still applies.
15. As a maker, I want a flat "basement" base plate beneath all terrain variation, with magnets at one fixed Z, so that mating and printing stay consistent regardless of terrain height.
16. As a maker, I want print scale to remain consistent (hex circumradius × scale, same physical size as before), so that existing sizing assumptions for minis and magnets hold.
17. As a developer, I want every exported mesh gated on `is_watertight` / `is_winding_consistent` / positive volume before it's ever written to disk, so that a broken mesh is caught immediately, not discovered later in a slicer.
18. As a player, I want standable hex pads to be genuinely flat at their terrain level in the real generated mesh, so that miniatures sit stable during play.
19. As a player, I want hex tiling preserved for character placement on standable pads, so that movement rules tied to the hex grid still work.
20. As a developer, I want flower layout math isolated in one module, so that coordinate and adjacency logic is testable independent of meshing.
21. As a developer, I want tileset loading and validation isolated, so that schema errors (including cross-flower side mismatches) are caught with clear messages.
22. As a developer, I want mesh building separated into layered modules (heightfield → triangulation → surface mesh → base plate → assembly), so that each stage is independently testable.
23. As a developer, I want a CLI to render a single flower by ID, so that I can iterate on one tile quickly.
24. As a developer, I want a CLI to render the full preview assembly, so that tileset QA is one command.
25. As a terrain designer, I want a flat_plains sample flower (uniform height, all-zero boundary), so that the simplest possible case is available as a smoke test and starting point.
26. As a terrain designer, I want a hill_peak sample flower with a deliberate cliff on one side and no intermediate slope step, so that decision #7 (no cap on cliff delta) and decision #11 (0 standable hexes is valid) are both demonstrated.
27. As a developer, I want the flower-to-flower grid placement math (`flower_grid_to_xy`) verified against the same axial directions the side-adjacency math uses, so that a tileset's `preview_map` positions genuinely correspond to real physical neighbors.
28. As a developer, I want `meta.scale` to actually be applied to the flower's physical footprint, so that a flower's XY size and its own possible height variation stay proportionate.

## Implementation Decisions

### Modules (as built)

| Module | Responsibility |
|--------|----------------|
| `terrain/layout.py` (`FlowerLayout`) | 7-hex flower graph: cell centers/polygons, the 18 true exterior edges grouped into 6 sides of 4 corners, hex-to-hex edges, road/water junction points, `flower_grid_to_xy()` for flower-to-flower placement |
| `terrain/heights.py` (`HeightLevels`) | Configurable discrete height-level system (level count, mm per level) |
| `terrain/boundary_noise.py` | Deterministic hash-based value noise (splitmix64-style bit-mixing, no trig/libm) for the jagged boundary contour and freeform interior relief |
| `terrain/heightfield.py` | Builds the per-side deterministic boundary contour, the whole-flower boundary loop, the freeform interior PSLG, and per-hex-cell groove triangulation |
| `terrain/triangulate.py` | Thin wrapper around `mapbox_earcut` (ear-clipping) and a point-in-polygon test, isolating the third-party triangulation dependency |
| `terrain/roads.py` | Road/river centerlines starting/ending exactly at side-corner positions, split via ear-clipping into two boundary-preserving regions |
| `terrain/surface_mesh.py` | Triangulates the PSLG (or road/groove variants), lifts to Z, returns a `trimesh.Trimesh` — or, for assembly, an open solid with no bottom cap |
| `terrain/base_plate.py` | Floor on the print bed, mirrored from the top surface onto its own bottom vertex copies (no boolean union, no ear-clipping) |
| `terrain/magnets.py` | 18 blind magnet bores (5.3 × 2.2 mm, centred 3.9 mm above the bed) cut into the walls at every silhouette edge's midpoint, as structured radial-sector collars |
| `terrain/standability.py` | Post-generation check of which hex cells are flat within tolerance in the real mesh, vs. the tileset's configured minimum |
| `terrain/assembly.py` | `build_flower_mesh()` (terrain + base plate, one welded solid), `build_preview_mesh()` (every `preview_map` flower placed and concatenated into one scene) |
| `terrain/export.py` | STL writer gated on `is_watertight` / `is_winding_consistent` / positive volume |
| `terrain/tileset.py` | Loads and validates the YAML schema, including cross-flower `preview_map` side matching |
| `terrain/cli.py` | `render flower <id>`, `render preview` |

### Tileset schema (as built)

```yaml
meta:
  height_step_mm: 15
  level_count: 4
  scale: 5
  hex_outer_width: 5.1961525
  min_standable_hexes: 1

preview_map:
  - { id: flat_plains, at: [0, 0], rot: 0 }
  - { id: hill_peak, at: [-1, 1], rot: 0 }

flowers:
  <flower_id>:
    hexes:
      "<0-6>": { height_level: <int> }
    side_corner_heights:
      "<0-5>": [c0, c1, c2, c3]
    seed: <int>
    roads:
      - { entry_junction: 0, exit_junction: 4 }
    water: []
```

- **Hex index 0** = center; **1–6** = ring.
- **`side_corner_heights`**: one 4-corner height-level sequence per physical
  side. Side `k`'s last corner must equal side `(k+1)%6`'s first corner
  (same-flower junction consistency, validated at load time). A neighbor
  sharing side `k` must declare the **reversed** sequence for its own side
  `(k+3)%6` — validated for any two `rot=0` `preview_map` placements that are
  grid-adjacent.
- **`seed`**: per-flower integer seeding the freeform interior relief noise —
  independent of the boundary contract (decision #5).
- **No named edge-profile catalog** — matching is exact numeric equality
  (once correctly oriented), not compatibility-table lookup.

## Testing Decisions

- Every test asserts against **real constructed geometry**: a `trimesh.Trimesh` or an actually-written, reloaded `.stl` — never string-matched generated text.
- Determinism itself is a first-class tested property: independently-built matching boundaries must be bit-identical (`tests/test_boundary_determinism.py`).
- Watertightness/winding/volume are checked at every meshing layer, not only at final export.
- Framework: pytest + YAML fixtures.

## Out of Scope

- Auto-layout / CSP map generator
- Top-face magnet sockets on standable hexes (the wall bores are done; see `terrain/magnets.py`)
- `HexDef.height_level` currently has no effect on generated geometry beyond validation — only `side_corner_heights` (boundary) and `seed` (interior) drive the mesh; flagged, not silently assumed
- Bridges (road + water on the same corner)
- Single-hex prints, web UI, game rules engine

## Further Notes

See [plan](../plans/3d-hex-flower-terrain.md) for the phase-by-phase build log, and [`.claude/plans/elegant-twirling-scott.md`](../../.claude/plans/elegant-twirling-scott.md) for the full decision-by-decision history (bugs found and fixed, library choices, and why).
