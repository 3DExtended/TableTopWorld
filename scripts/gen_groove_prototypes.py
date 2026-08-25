"""Generate a handful of genuinely different groove-style prototypes on the
same (now fixed) exact-inset-scale ring engine in terrain/heightfield.py.

Each style is just a different (groove_width_mm, groove_profile) pair - the
profile is a list of (width_fraction, z_offset_mm) control points tracing
the groove's cross-section from the true edge inward. Run directly:

    python scripts/gen_groove_prototypes.py
"""

from __future__ import annotations

from pathlib import Path

import trimesh

from terrain.assembly import build_flower_mesh
from terrain.tileset import load_tileset

OUT_DIR = Path(__file__).resolve().parent.parent / "output" / "groove_prototypes"

STYLES = {
    # No separate entry-wall segment ("outer rim") at all - a single row
    # covering the full (very narrow) width, straight from the true edge
    # down to the floor. Chosen after comparing against flat_channel,
    # narrow_notch, wide_shallow, and terraced (all dropped).
    "thin_line": dict(groove_width_mm=0.5, groove_profile=[(1.0, -0.6)]),
}

# xy_jitter_mm disabled - on hill_peak's steep cliffs, nudging a boundary
# point sideways by up to 1mm combined with a big per-step height change
# creates a real (non-defective, but undesirable) narrow steep sliver in
# the terrain. Not a mesh bug (confirmed: watertight, no inverted normals,
# zero self-intersections found via brute-force triangle-pair testing),
# just an unwanted side effect of XY jitter on very cliffy flowers.
NOISE_OVERRIDES = dict(xy_jitter_mm=0.0)


def main() -> None:
    tileset = load_tileset("tilesets/default.yaml")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for flower_id in ("flat_plains", "hill_peak"):
        for name, kwargs in STYLES.items():
            mesh = build_flower_mesh(
                tileset, flower_id, subdivisions_per_edge=12, **kwargs, **NOISE_OVERRIDES
            )
            assert mesh.is_watertight, f"{flower_id}/{name} not watertight"
            assert mesh.is_winding_consistent, f"{flower_id}/{name} bad winding"
            out_path = OUT_DIR / f"{flower_id}_{name}.stl"
            mesh.export(out_path)
            print(f"wrote {out_path} ({len(mesh.faces)} faces)")


if __name__ == "__main__":
    main()
