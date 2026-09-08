"""Explicit mesh construction: triangulated top surface + walls + bottom cap.

build_flat_prism_mesh is the first real geometry in the redesigned
pipeline: a flat-topped prism, proving the triangulate -> lift-to-Z ->
wall-stitch -> trimesh.Trimesh round trip before later phases replace the
flat top with a jagged heightfield-driven surface and freeform interior
on top of this same structure.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Collection, Sequence

import numpy as np
import trimesh

from terrain.heightfield import (
    build_flower_boundary_loop,
    build_flower_cells,
    build_flower_pslg,
)
from terrain.layout import FlowerLayout
from terrain.magnets import MagnetBores, cut_magnet_bores
from terrain.roads import build_flower_mesh_with_road
from terrain.triangulate import triangulate_polygon

Point2D = tuple[float, float]


def build_flat_prism_mesh(
    boundary: Sequence[Point2D], top_z: float, bottom_z: float = 0.0
) -> trimesh.Trimesh:
    """A single closed solid: a triangulated flat cap at top_z, a mirrored
    flat cap at bottom_z, and one outward-facing wall quad (2 triangles)
    per boundary edge. `boundary` must be CCW as seen from +Z looking down.
    """
    points_2d, top_triangles = triangulate_polygon(boundary)
    n = len(points_2d)

    top_verts = [(x, y, top_z) for x, y in points_2d]
    bottom_verts = [(x, y, bottom_z) for x, y in points_2d]
    vertices = top_verts + bottom_verts  # bottom vertex i lives at index n + i

    faces: list[tuple[int, int, int]] = []
    faces.extend(top_triangles)  # already CCW -> normal +Z
    for i, j, k in top_triangles:
        faces.append((n + k, n + j, n + i))  # reversed winding -> normal -Z

    boundary_n = len(boundary)
    for e in range(boundary_n):
        i0, i1 = e, (e + 1) % boundary_n
        top_a, top_b = i0, i1
        bot_a, bot_b = n + i0, n + i1
        faces.append((top_a, bot_a, top_b))
        faces.append((top_b, bot_a, bot_b))

    return trimesh.Trimesh(
        vertices=np.array(vertices), faces=np.array(faces), process=True
    )


def _build_flower_top_and_walls(
    side_corner_heights: dict[int, tuple[int, int, int, int]],
    layout: FlowerLayout,
    level_z: Callable[[int], float],
    seed: int,
    *,
    bottom_z: float,
    subdivisions_per_edge: int,
    jitter_amplitude: float,
    xy_jitter_mm: float,
    interior_grid_step: float | None,
    interior_relief_mm: float,
    include_hex_grooves: bool,
    groove_depth_mm: float,
    groove_width_mm: float,
    groove_profile: list[tuple[float, float]] | None = None,
    hex_height_levels: dict[int, int] | None = None,
    road_water_side_pairs: Sequence[tuple[int, int]],
    road_subdivisions: int,
    magnet_bores: MagnetBores | None = None,
    magnet_sockets: MagnetBores | None = None,
    socket_hexes: Collection[int] = (),
) -> tuple[
    list[tuple[float, float, float]],
    list[tuple[int, int, int]],
    list[tuple[int, int, int]],
    list[int],
    int,
    list[tuple[int, int, int]],
]:
    """Shared body behind both build_flower_surface_mesh (closed solid) and
    build_flower_open_solid (open at the bottom, for terrain/assembly.py to
    weld a base plate onto): the top surface (one of 3 branches - road /
    grooves / Phase 4 default) plus the wall quads connecting it down to a
    flat bottom rim at bottom_z, but NOT the bottom cap itself.

    `magnet_bores` cuts one blind bore per silhouette edge into those walls
    (terrain/magnets.py); the boundary's XY jitter is switched off around
    each bore so the wall there is planar. bottom_z is then the print bed.
    `magnet_sockets` (same disc dimensions) sinks a socket into the middle
    of each hex in `socket_hexes` - grooves branch only.

    Returns (vertices, faces, top_triangles, bottom_rim_indices, n,
    footprint_triangles) - n is the vertex count of the top surface alone
    (bottom vertex i lives at index n+i); top_triangles is also returned
    separately (not just folded into faces); bottom_rim_indices are
    (already offset by n) the flat-bottom counterparts of the boundary
    loop, in boundary order; footprint_triangles is the planar XY
    triangulation a bottom cap must mirror (the top surface with any
    magnet sockets capped over - see heightfield.build_flower_cells).
    """
    xy_flat_window_mm = 0.0
    if magnet_bores is not None:
        # every silhouette edge is one hex edge long (= circumradius)
        xy_flat_window_mm = magnet_bores.flat_window_mm(
            layout.hex_outer_width, subdivisions_per_edge
        )

    if magnet_sockets is not None and socket_hexes and not (
        include_hex_grooves and not road_water_side_pairs
    ):
        raise NotImplementedError("magnet sockets are only built by the hex-cell branch")

    if road_water_side_pairs:
        if len(road_water_side_pairs) > 1:
            raise NotImplementedError(
                "only one road/river per flower is currently supported"
            )
        entry_side, exit_side = road_water_side_pairs[0]
        boundary_loop_3d = build_flower_boundary_loop(
            side_corner_heights,
            layout,
            level_z,
            subdivisions_per_edge=subdivisions_per_edge,
            jitter_amplitude=jitter_amplitude,
            xy_jitter_mm=xy_jitter_mm,
            xy_flat_window_mm=xy_flat_window_mm,
        )
        final_vertices_3d, top_triangles = build_flower_mesh_with_road(
            boundary_loop_3d,
            entry_side,
            exit_side,
            side_corner_heights[entry_side][0],
            side_corner_heights[exit_side][0],
            level_z,
            subdivisions_per_edge=subdivisions_per_edge,
            road_subdivisions=road_subdivisions,
        )
        boundary_count = len(boundary_loop_3d)
        final_points_2d = [(x, y) for x, y, _ in final_vertices_3d]
        boundary_indices = list(range(boundary_count))
        footprint = top_triangles
    elif include_hex_grooves:
        final_vertices_3d, boundary_indices, top_triangles, footprint = build_flower_cells(
            side_corner_heights,
            layout,
            level_z,
            seed,
            subdivisions_per_edge=subdivisions_per_edge,
            jitter_amplitude=jitter_amplitude,
            xy_jitter_mm=xy_jitter_mm,
            xy_flat_window_mm=xy_flat_window_mm,
            interior_relief_mm=interior_relief_mm,
            groove_depth_mm=groove_depth_mm,
            groove_width_mm=groove_width_mm,
            groove_profile=groove_profile,
            hex_height_levels=hex_height_levels,
            socket_hexes=socket_hexes if magnet_sockets is not None else (),
            socket_radius_mm=magnet_sockets.radius_mm if magnet_sockets else 0.0,
            socket_depth_mm=magnet_sockets.depth_mm if magnet_sockets else 0.0,
        )
        boundary_count = len(boundary_indices)
        final_points_2d = [(x, y) for x, y, _ in final_vertices_3d]
    else:
        all_points_2d, all_vertices_3d, boundary_count = build_flower_pslg(
            side_corner_heights,
            layout,
            level_z,
            seed,
            subdivisions_per_edge=subdivisions_per_edge,
            jitter_amplitude=jitter_amplitude,
            xy_jitter_mm=xy_jitter_mm,
            xy_flat_window_mm=xy_flat_window_mm,
            interior_grid_step=interior_grid_step,
            interior_relief_mm=interior_relief_mm,
        )
        boundary_2d = all_points_2d[:boundary_count]
        interior_2d = all_points_2d[boundary_count:]
        final_points_2d = all_points_2d
        final_vertices_3d = all_vertices_3d
        _, top_triangles = triangulate_polygon(boundary_2d, interior_2d)
        boundary_indices = list(range(boundary_count))
        footprint = top_triangles

    n = len(final_points_2d)
    top_verts = final_vertices_3d
    bottom_verts = [(x, y, bottom_z) for x, y in final_points_2d]
    vertices = top_verts + bottom_verts

    faces: list[tuple[int, int, int]] = list(top_triangles)
    skipped: set[int] = set()
    if magnet_bores is not None:
        bore_vertices, bore_faces, skipped = cut_magnet_bores(
            vertices,
            n,
            boundary_indices,
            layout,
            subdivisions_per_edge=subdivisions_per_edge,
            bottom_z=bottom_z,
            bores=magnet_bores,
        )
        vertices = vertices + bore_vertices
        faces.extend(bore_faces)
    bottom_rim_indices: list[int] = []
    for e in range(boundary_count):
        i0 = boundary_indices[e]
        i1 = boundary_indices[(e + 1) % boundary_count]
        bottom_rim_indices.append(n + i0)
        if e in skipped:
            continue  # replaced by a bore patch
        faces.append((i0, n + i0, i1))
        faces.append((i1, n + i0, n + i1))

    return vertices, faces, top_triangles, bottom_rim_indices, n, footprint


def build_flower_surface_mesh(
    side_corner_heights: dict[int, tuple[int, int, int, int]],
    layout: FlowerLayout,
    level_z: Callable[[int], float],
    seed: int,
    *,
    bottom_z: float = 0.0,
    subdivisions_per_edge: int = 8,
    jitter_amplitude: float = 0.3,
    xy_jitter_mm: float = 0.0,
    interior_grid_step: float | None = None,
    interior_relief_mm: float = 6.0,
    include_hex_grooves: bool = False,
    groove_depth_mm: float = 0.0,
    groove_width_mm: float = 2.0,
    groove_profile: list[tuple[float, float]] | None = None,
    hex_height_levels: dict[int, int] | None = None,
    road_water_side_pairs: Sequence[tuple[int, int]] = (),
    road_subdivisions: int = 10,
    magnet_bores: MagnetBores | None = None,
    magnet_sockets: MagnetBores | None = None,
    socket_hexes: Collection[int] = (),
) -> trimesh.Trimesh:
    """Full 7-hex flower as one closed solid: the deterministic jagged
    boundary (decision #2-#4) triangulated together with a freeform seeded
    interior (decision #5), extruded straight down to bottom_z with no cap
    on the top surface's own relief (decision #7 - cliffs, not slopes,
    are the expected shape).

    `include_hex_grooves` (decision #10) switches the top surface from one
    whole-flower Delaunay triangulation to 7 independently-triangulated
    hex cells (terrain.heightfield.build_flower_cells) - a convex
    polygon's own hull edges are always part of any Delaunay triangulation
    of its own points, so each hex's 6 edges are guaranteed present with
    no ambiguity, unlike the whole-flower approach at the hex grid's
    highly symmetric shared vertices (see build_flower_cells' docstring).

    `road_water_side_pairs` (decisions #8-#9) embeds a single road/river
    centerline as a real edge chain via
    terrain.roads.build_flower_mesh_with_road, which splits the flower
    into two regions at the road and triangulates each with ear-clipping
    (mapbox_earcut) rather than Delaunay - a guarantee, not a probabilistic
    hope; see terrain/roads.py's module docstring for why Delaunay wasn't
    reliable here. Only one road per flower is supported (splitting into
    more than two regions isn't implemented); a road replaces the
    freeform interior noise with a plain interpolated surface, and combining
    it with `include_hex_grooves` on the same flower isn't wired up yet.

    A flower whose corner heights all sit at the same level as bottom_z
    has zero material thickness there and will NOT be watertight (a real
    hole, not a bug) - terrain/assembly.py always calls this with a
    bottom_z below every configured height level for exactly this reason.
    """
    vertices, faces, _, _, n, footprint = _build_flower_top_and_walls(
        side_corner_heights,
        layout,
        level_z,
        seed,
        bottom_z=bottom_z,
        subdivisions_per_edge=subdivisions_per_edge,
        jitter_amplitude=jitter_amplitude,
        xy_jitter_mm=xy_jitter_mm,
        interior_grid_step=interior_grid_step,
        interior_relief_mm=interior_relief_mm,
        include_hex_grooves=include_hex_grooves,
        groove_depth_mm=groove_depth_mm,
        groove_width_mm=groove_width_mm,
        groove_profile=groove_profile,
        hex_height_levels=hex_height_levels,
        road_water_side_pairs=road_water_side_pairs,
        road_subdivisions=road_subdivisions,
        magnet_bores=magnet_bores,
        magnet_sockets=magnet_sockets,
        socket_hexes=socket_hexes,
    )
    for i, j, k in footprint:
        faces.append((n + k, n + j, n + i))

    return trimesh.Trimesh(
        vertices=np.array(vertices), faces=np.array(faces), process=True
    )


@dataclass(frozen=True)
class OpenSolid:
    """A flower's top surface + walls, open at the bottom (see
    build_flower_open_solid). `vertices` is the top surface's vertices,
    then the flat-bottom copy of every one of them (bottom copy of top
    index i is top_count + i), then any bore-patch vertices;
    bottom_rim_indices are the boundary loop's bottom copies in boundary
    order; top_triangles are the top surface's faces (also the first
    entries of `faces`), indexing the top vertices."""

    vertices: list[tuple[float, float, float]]
    faces: list[tuple[int, int, int]]
    bottom_rim_indices: list[int]
    top_triangles: list[tuple[int, int, int]]
    top_count: int
    # the planar XY triangulation a floor must mirror (sockets capped over)
    footprint_triangles: list[tuple[int, int, int]]


def build_flower_open_solid(
    side_corner_heights: dict[int, tuple[int, int, int, int]],
    layout: FlowerLayout,
    level_z: Callable[[int], float],
    seed: int,
    *,
    bottom_z: float,
    subdivisions_per_edge: int = 8,
    jitter_amplitude: float = 0.3,
    xy_jitter_mm: float = 0.0,
    interior_grid_step: float | None = None,
    interior_relief_mm: float = 6.0,
    include_hex_grooves: bool = False,
    groove_depth_mm: float = 0.0,
    groove_width_mm: float = 2.0,
    groove_profile: list[tuple[float, float]] | None = None,
    hex_height_levels: dict[int, int] | None = None,
    road_water_side_pairs: Sequence[tuple[int, int]] = (),
    road_subdivisions: int = 10,
    magnet_bores: MagnetBores | None = None,
    magnet_sockets: MagnetBores | None = None,
    socket_hexes: Collection[int] = (),
) -> OpenSolid:
    """Same top surface + walls as build_flower_surface_mesh, but with NO
    bottom cap - left open at bottom_z for terrain/assembly.py to close
    with a floor (terrain/base_plate.py) that reuses this solid's own
    vertex indices directly, rather than building two separately-capped
    solids and unioning them with a boolean (which this whole redesign
    deliberately avoids - see the project plan's library-choice rationale).
    """
    vertices, faces, top_triangles, bottom_rim_indices, top_count, footprint = _build_flower_top_and_walls(
        side_corner_heights,
        layout,
        level_z,
        seed,
        bottom_z=bottom_z,
        subdivisions_per_edge=subdivisions_per_edge,
        jitter_amplitude=jitter_amplitude,
        xy_jitter_mm=xy_jitter_mm,
        interior_grid_step=interior_grid_step,
        interior_relief_mm=interior_relief_mm,
        include_hex_grooves=include_hex_grooves,
        groove_depth_mm=groove_depth_mm,
        groove_width_mm=groove_width_mm,
        groove_profile=groove_profile,
        hex_height_levels=hex_height_levels,
        road_water_side_pairs=road_water_side_pairs,
        road_subdivisions=road_subdivisions,
        magnet_bores=magnet_bores,
        magnet_sockets=magnet_sockets,
        socket_hexes=socket_hexes,
    )
    return OpenSolid(vertices, faces, bottom_rim_indices, top_triangles, top_count, footprint)
