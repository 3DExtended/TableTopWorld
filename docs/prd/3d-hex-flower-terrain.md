# PRD: 3D modular hex-flower terrain system (MVP)

**Issue:** https://github.com/3DExtended/TableTopWorld/issues/1  
**Labels:** `enhancement`, `ready-for-agent`  
**Plan:** [docs/plans/3d-hex-flower-terrain.md](../plans/3d-hex-flower-terrain.md)

---

## Problem Statement

TableTopWorld currently generates flat 7-hex **flower** tiles with carved roads, water channels, magnets, and topping holes—but all terrain lives on a single height plane. The project needs **3D modular terrain** where hills and mountains connect across flowers while preserving the hex grid for miniature placement. Today, every new tile variant requires editing Python constants (junction indices, heights, flags), there is no formal model for which flowers may connect to which, and there is no way to preview an entire **tileset** as one assembled map.

## Solution

Introduce a **declarative tileset** (YAML) that describes each flower’s seven hex cells (terrain level + standable vs non-standable role), eighteen exterior **edge profiles** for modular mating, and road/water features at the correct height. A refactored **terrain generator** (solid2/BOSL2, same stack as today) reads tilesets, validates adjacency rules, and exports per-flower and **combined preview** meshes. Terrain uses three dramatic height steps (~2 cm per step at print scale); water sits one step below local terrain. MVP ships three sample flowers and a manual preview map; auto-layout and bridge toppings follow later.

## User Stories

1. As a terrain designer, I want to define flowers in a tileset YAML file, so that I can add new tiles without editing generator code.
2. As a terrain designer, I want each hex cell to be either mini-standable or non-standable, so that gameplay surfaces are unambiguous on the tabletop.
3. As a terrain designer, I want every flower to require at least one standable hex, so that players always have a legal placement option on each tile.
4. As a terrain designer, I want three terrain height levels (ground, middle, high), so that hills and mountains read clearly at ~2 cm per step when printed.
5. As a terrain designer, I want water to sit one terrain step below the hosting hex’s level, so that rivers and lakes feel recessed relative to surrounding land.
6. As a terrain designer, I want roads carved as channels at each hex’s terrain height, so that roads follow hills instead of only existing on flat ground.
7. As a terrain designer, I want eighteen exterior edges per flower each tagged with a profile, so that only compatible tiles physically and visually mate.
8. As a terrain designer, I want a catalog of edge profiles (flat per level, slopes, cliffs), so that transitions between flowers are standardized and repeatable.
9. As a terrain designer, I want the generator to reject tilesets where preview-map neighbors have incompatible edge profiles, so that layout mistakes fail at build time not at the printer.
10. As a terrain designer, I want road paths defined by junction pairs (same convention as today), so that existing road logic remains intuitive.
11. As a terrain designer, I want water paths defined similarly to roads, so that rivers use the same authoring model.
12. As a terrain designer, I want a manual preview_map listing flower IDs, positions, and rotations, so that I can see a deliberate test layout for the whole tileset.
13. As a terrain designer, I want a combined SCAD/STL export of the preview_map, so that I can visually judge how tileset edits affect the full terrain.
14. As a maker, I want the printable unit to remain a single 7-hex flower plate with edge magnets, so that my existing magnet workflow still applies.
15. As a maker, I want magnet holes and bevels preserved on exterior edges, so that tiles still align and interlock as they do today.
16. As a maker, I want print scale to remain consistent (5× export), so that existing sizing assumptions for minis and magnets hold.
17. As a maker, I want configurable mesh resolution for preview vs final print, so that iteration stays fast while final exports stay smooth.
18. As a player, I want standable hex pads to be flat at their terrain level, so that miniatures sit stable during play.
19. As a player, I want non-standable hexes (slopes, cliffs, deep channels) to be obviously not for standing, so that table rules are clear.
20. As a player, I want hex tiling preserved for character placement on standable pads, so that movement rules tied to the hex grid still work.
21. As a developer, I want flower layout math isolated in one module, so that coordinate and edge-ID logic is testable without CSG.
22. As a developer, I want tileset loading and validation isolated, so that schema errors are caught with clear messages.
23. As a developer, I want edge-profile compatibility isolated, so that adjacency rules are data-driven and unit-testable.
24. As a developer, I want mesh building separated from feature cutting (roads/water), so that height geometry can evolve without rewriting bezier cutters.
25. As a developer, I want a CLI to render a single flower by ID, so that I can iterate on one tile quickly.
26. As a developer, I want a CLI to render the full preview assembly, so that tileset QA is one command.
27. As a terrain designer, I want a sample flat_plains flower matching current flat+road behavior, so that the migration path from today’s script is obvious.
28. As a terrain designer, I want a sample hill_north flower with mixed heights and directional rise, so that slope and standable combinations are demonstrated.
29. As a terrain designer, I want a sample river_grove flower with water channels, so that sub-terrain water is demonstrated without bridges yet.
30. As a terrain designer, I want interior hex-to-hex transitions (slope role) within a flower, so that height changes on a single plate look natural.
31. As a terrain designer, I want exterior edge geometry to reflect profile type (flat face, slope, cliff), so that seams between flowers look intentional.
32. As a developer, I want validation to ensure every edge key refers to a real exterior side, so that typos in YAML cannot silently generate wrong geometry.
33. As a developer, I want validation to ensure junction indices are in range 0–5, so that road/water definitions cannot reference invalid paths.
34. As a terrain designer, I want tileset meta to include height_step_mm and scale, so that physical height intent is documented alongside the file.
35. As a maker, I want combined flowers in preview_map positioned in hex-flower coordinates, so that multi-tile maps match how tiles will be placed on the table.
36. As a terrain designer, I want incompatible profile pairs documented in the edge catalog, so that I know which edges can mate before authoring a new flower.
37. As a developer, I want the legacy monolithic generator preserved or callable during refactor, so that regression against today’s output is possible for flat_plains.
38. As a terrain designer, I want non-standable hexes to include road channels and water without being standable, so that features and play surfaces stay distinct.
39. As a player, I want height differences to be dramatic (~2 cm per step), so that terrain reads as mountains and hills not subtle texture.
40. As a terrain designer, I want to express per-hex terrain level independently across the seven cells, so that asymmetric hills are supported on one flower.
41. As a developer, I want assembly export to apply rotation and translation from preview_map, so that combined previews match authored orientation.
42. As a terrain designer, I want clear error output when preview_map adjacency fails validation, so that I can fix edge profiles on the offending pair.
43. As a maker, I want topping magnet holes to remain available on designated hexes (MVP unchanged), so that future bridge and building pieces still attach.
44. As a terrain designer, I want road and water features to respect bevel intersections along their paths, so that carved channels match today’s visual quality where possible.
45. As a developer, I want export to produce OpenSCAD via solid2 save_as_scad, so that the existing OpenSCAD/BOSL2 toolchain remains the source of truth.

## Implementation Decisions

### Deep modules (build or extract)

| Module | Responsibility | Interface (stable) |
|--------|----------------|-------------------|
| **FlowerLayout** | 7-hex flower graph: cell centers, subdivided hex polygons, 18 exterior edge identifiers, 6 three-hex junctions for roads/water | `edge_id(hex_idx, side_idx)`, `junction_lines(junction_idx)`, `cell_polygon(hex_idx)` |
| **EdgeProfileCatalog** | Named profiles, Z semantics, compatibility pairs | `is_compatible(profile_a, profile_b)`, `elevation_at_edge(profile)` |
| **TilesetSchema** | Load YAML tileset; validate flowers, edges, preview_map | `load(path) -> Tileset`, `validate(tileset) -> list[Error]` |
| **FlowerMeshBuilder** | CSG for per-hex slabs at terrain Z; interior slope blending; standable flat tops | `build_flower(flower_def) -> Solid` |
| **EdgeGeometry** | Bevels, magnet pockets, cliff/slope faces along exterior edges from profile + height delta | `apply_edges(flower_solid, flower_def, layout) -> Solid` |
| **FeatureCutters** | Road/water bezier channel subtraction at host hex Z (port of existing path_sweep logic) | `apply_features(flower_solid, flower_def, layout) -> Solid` |
| **AssemblyExporter** | Place multiple flowers per preview_map; union or scene export | `build_preview(tileset) -> Solid` |
| **CLI** | `render flower <id>`, `render preview <tileset>` | argparse entry points |

### Tileset schema (decision-rich shape)

```yaml
meta:
  height_step_mm: 20
  scale: 5
  hex_outer_width: 5.1961525

preview_map:
  - { id: flat_plains, at: [q, r], rot: 0 }

flowers:
  <flower_id>:
    hexes:
      "<0-6>": { terrain: ground|middle|high, role: standable|slope|cliff|road_channel|water|... }
    edges:
      "<hex>-<side>": { profile: <catalog_name> }
    roads:
      - { entry_junction: 0, exit_junction: 4 }
    water: []
```

- **Hex index 0** = center; **1–6** = ring CCW.
- **Terrain Z** in model units: `level_index * (height_step_mm / scale)` → **4.0 units per step** with 20 mm step and scale 5.
- **Water Z** = host hex terrain Z minus one step.
- **Edge profiles (initial)**: `flat_ground`, `flat_middle`, `flat_high`; paired slopes; `cliff_*`; deferred ports for bridges.
- **MVP samples**: `flat_plains`, `hill_north`, `river_grove` + 3–5 flower preview_map.

## Testing Decisions

- Test **external behavior** of pure logic (validation, compatibility, layout IDs, Z math)—not CSG graphs.
- **Test:** TilesetSchema, EdgeProfileCatalog, FlowerLayout; partial AssemblyExporter adjacency checks.
- **Skip MVP:** mesh/STL snapshot tests.
- **Framework:** pytest + YAML fixtures (no prior tests in repo).

## Out of Scope

- Auto-map CSP generator (Milestone 2)
- Bridge/deck magnet toppings and dual-layer standable hexes
- Road + river combo / bridge YAML
- Per-table physical edge matching (build-time preview only in MVP)
- Single-hex prints, web UI, game rules engine

## Further Notes

See [plan](../plans/3d-hex-flower-terrain.md) for architecture diagram and success criteria.
