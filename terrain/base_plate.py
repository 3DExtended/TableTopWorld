"""Flat basement layer beneath all terrain variation - design decision #13.

Magnets are meant to sit at one fixed Z (MAGNET_CENTER_Z, measured up
from the print bed) regardless of what height level the terrain above
happens to be - this module is that flat "basement": a floor cap plus 6
solid wall panels (one per flower side, not per raw exterior edge - a
side's 2 end corners are enough for a *structural* panel even though the
terrain surface above traces a finer jagged silhouette between them).

Magnet bore holes are NOT cut into these walls yet - a real bore is a
BLIND recess (MAGNET_DEPTH deep, not a hole clean through), which needs
the wall to have actual thickness (an inner face offset inward from the
outer one, with the recess stopping partway between them) so the solid
stays watertight. This module's walls are currently a single zero-
thickness sheet (extruded only via the floor cap and the terrain's own
walls above, same pattern as every other surface in this codebase), so
there is no "inward" direction to stop a blind recess partway through -
modeling that properly means giving these 6 panels real thickness
(inner + outer faces, connected at all 4 edges, with the recess cut only
partway from the outer face) which is a genuinely separate, larger task
from "weld a flat plate onto the terrain's flat rim". Tracked as a
follow-up; `magnet_center_z`/`magnet_radius` are accepted here (and
plate_depth is still validated against them) so the call signature
doesn't need to change again once that follow-up lands.

This does NOT build a second, separately-capped solid and boolean-union
it onto the terrain surface - this whole redesign deliberately avoids
CSG (see the project plan's library-choice rationale). Instead,
terrain/surface_mesh.py::build_flower_open_solid() leaves the terrain
mesh open at its flat bottom rim (constant Z = bottom_z all the way
around, unlike the jagged top - the rim is genuinely flat, so it's a
legitimate weld seam) and returns that rim's own vertex indices; this
module's build_base_plate_parts() takes those SAME vertices (by index,
not by recomputing matching coordinates and hoping trimesh's merge
tolerance welds them - the exact failure mode debugged at length in
Phase 5) and extends the vertex/face lists with the floor and walls
below. terrain/assembly.py wraps the combined result in one
trimesh.Trimesh.
"""

from __future__ import annotations

from typing import Sequence

from terrain.layout import FlowerLayout
from terrain.triangulate import earcut_triangulate_polygon

Point2D = tuple[float, float]
Vertex3D = tuple[float, float, float]
Triangle = tuple[int, int, int]


def build_base_plate_parts(
    vertices: list[Vertex3D],
    rim_indices: Sequence[int],
    layout: FlowerLayout,
    *,
    subdivisions_per_edge: int,
    bottom_z: float,
    plate_depth: float,
    magnet_center_z: float,
    magnet_radius: float,
) -> tuple[list[Vertex3D], list[Triangle]]:
    """Extend `vertices` with the base plate's own new vertices (floor cap
    + wall panels) and return (extra_vertices, new_faces) - new_faces
    index into vertices + extra_vertices combined (i.e. some face indices
    are < len(vertices), reusing the terrain mesh's own rim vertices
    directly; others are >= len(vertices), into the newly-appended ones).
    The caller concatenates and builds one Trimesh; nothing here calls
    trimesh itself.

    `rim_indices` must be the flat bottom rim of a flower built with this
    same `subdivisions_per_edge`, in boundary order, all at z=bottom_z
    (terrain.surface_mesh.build_flower_open_solid's own return value) -
    this is what lets the plate's top opening weld to the terrain
    without a boolean union.

    Requires plate_depth > magnet_center_z + magnet_radius (a future
    magnet bore must fit entirely within the plate's own vertical band,
    below bottom_z) - raises ValueError otherwise, since a magnet poking
    out of the floor or into the terrain above would be a real, not
    cosmetic, print defect. Checked now even though the bore itself isn't
    cut yet, so a tileset/constants change can't silently violate it
    later.
    """
    if plate_depth <= magnet_center_z + magnet_radius:
        raise ValueError(
            f"plate_depth ({plate_depth}) must exceed magnet_center_z + "
            f"magnet_radius ({magnet_center_z + magnet_radius}) for a "
            "magnet bore to fit inside the base plate"
        )

    floor_z = bottom_z - plate_depth
    edges_per_side = FlowerLayout.EDGES_PER_SIDE * subdivisions_per_edge
    total_rim = len(rim_indices)
    if total_rim != FlowerLayout.SIDE_COUNT * edges_per_side:
        raise ValueError(
            f"rim_indices length {total_rim} doesn't match "
            f"SIDE_COUNT * EDGES_PER_SIDE * subdivisions_per_edge "
            f"({FlowerLayout.SIDE_COUNT * edges_per_side})"
        )

    extra_vertices: list[Vertex3D] = []
    faces: list[Triangle] = []

    def add_vertex(v: Vertex3D) -> int:
        extra_vertices.append(v)
        return len(vertices) + len(extra_vertices) - 1

    # Floor cap: the flower's true footprint (the rim, projected flat),
    # earcut-triangulated (the 18-corner silhouette is not convex, and its
    # fine subdivision points are collinear along each side - both reasons
    # scipy Delaunay is the wrong tool here, per the Phase 5 findings).
    floor_ring_2d = [(vertices[i][0], vertices[i][1]) for i in rim_indices]
    floor_triangles_local = earcut_triangulate_polygon(floor_ring_2d)
    floor_indices = [add_vertex((x, y, floor_z)) for x, y in floor_ring_2d]
    for i, j, k in floor_triangles_local:
        # Reversed winding: floor cap faces -Z (downward, outward from the
        # solid), while the rim's own polygon is wound CCW as seen from +Z.
        faces.append((floor_indices[k], floor_indices[j], floor_indices[i]))

    # 6 wall panels, one per side: a simple quad strip from the floor rim
    # up to the terrain's own bottom rim, both at the SAME fine per-side
    # granularity (not simplified to the 2 coarse end corners - an earlier
    # version of this code did that and left the floor cap's fine edges
    # with only one adjacent face each, a real caught hole: the floor
    # cap's own boundary already uses the full fine rim, so anything
    # welding to it must match that granularity, not a coarser one).
    for side_idx in range(FlowerLayout.SIDE_COUNT):
        start = side_idx * edges_per_side
        end = (start + edges_per_side) % total_rim
        side_rim = list(rim_indices[start : start + edges_per_side])
        side_rim.append(rim_indices[end])
        floor_rim = list(floor_indices[start : start + edges_per_side])
        floor_rim.append(floor_indices[end])

        def resolve(idx: int) -> Vertex3D:
            return vertices[idx] if idx < len(vertices) else extra_vertices[idx - len(vertices)]

        outward = _outward_flip(resolve, floor_rim[0], floor_rim[-1], side_rim[0])

        for e in range(edges_per_side):
            top_a, top_b = side_rim[e], side_rim[e + 1]
            bot_a, bot_b = floor_rim[e], floor_rim[e + 1]
            if outward:
                faces.append((top_a, bot_a, top_b))
                faces.append((top_b, bot_a, bot_b))
            else:
                faces.append((top_b, bot_a, top_a))
                faces.append((bot_b, bot_a, top_b))

    return extra_vertices, faces


def _outward_flip(resolve, floor_a: int, floor_b: int, top_a: int) -> bool:
    """True if the quad-strip winding (top_a, bot_a, top_b) already faces
    outward (away from the flower's center at the origin) for this side,
    checked once via the wall's own triangle normal against the outward
    direction (same convention as the rest of this codebase's wall
    construction: top_a -> bot_a -> top_b, bot_a -> top_b -> bot_b for a
    CCW-as-seen-from-outside boundary)."""
    a = resolve(top_a)
    b = resolve(floor_a)
    c = resolve(floor_b)
    ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    ac = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    normal = (
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    )
    mid_xy = ((b[0] + c[0]) / 2, (b[1] + c[1]) / 2)
    dot = normal[0] * mid_xy[0] + normal[1] * mid_xy[1]
    return dot > 0
