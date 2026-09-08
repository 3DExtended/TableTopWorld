"""The flat bottom of a flower: a floor on the print bed.

The "basement" of decision #13 (a flat band below all terrain variation,
where every tile's magnets sit at one fixed height) is simply the part of
the solid between the level-0 surface and the print bed - the terrain's
own walls run all the way down to bottom_z = -BASE_PLATE_DEPTH_MM, and
this module closes them with a floor. The magnet bores are cut into
those same walls by terrain/magnets.py; there is no separate wall band
any more (an earlier version stacked a second, bore-less wall band under
the terrain walls, doubling the base and leaving nowhere for a blind
recess to live).

The floor is the top surface's own triangulation, mirrored onto the flat
bottom copies of its vertices (terrain/surface_mesh.py keeps one bottom
copy per top vertex, at index top_count + i). That is a few hundred more
triangles than the outline strictly needs, but it can never be
degenerate: the top is a lifted planar triangulation, so every mirrored
face has the same non-zero footprint. Ear-clipping the bottom rim on its
own was tried first and emits zero-area slivers wherever the rim runs
exactly straight - which it now does around every magnet bore.

No boolean union anywhere: the floor's rim edges are literally the wall
quads' own bottom vertices, so the seam is shared indices, not two
coincidentally-equal coordinate sets (the failure class debugged at
length in Phase 5).
"""

from __future__ import annotations

from typing import Sequence

Triangle = tuple[int, int, int]


def build_floor_cap(top_triangles: Sequence[Triangle], top_count: int) -> list[Triangle]:
    """Mirror the top surface onto its bottom copies, wound to face -Z."""
    return [(top_count + k, top_count + j, top_count + i) for i, j, k in top_triangles]
