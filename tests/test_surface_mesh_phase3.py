"""Phase 3: first real mesh - single hex, flat top, watertight prism.

Proves the triangulate -> lift-to-Z -> wall-stitch -> trimesh round trip
before later phases replace the flat top with a jagged heightfield-driven
surface. Asserts against real constructed geometry (a trimesh.Trimesh, and
an actually-exported/reloaded STL) rather than string-matched text.
"""

from __future__ import annotations

import pytest
import trimesh

from terrain.export import MeshValidationError, export_stl
from terrain.layout import FlowerLayout
from terrain.surface_mesh import build_flat_prism_mesh

ONE_LEVEL_MM = 15.0


def test_single_hex_prism_is_watertight_and_correctly_sized(tmp_path) -> None:
    layout = FlowerLayout()
    boundary = layout.ring_vertices(0)
    mesh = build_flat_prism_mesh(boundary, top_z=ONE_LEVEL_MM, bottom_z=0.0)

    assert mesh.is_watertight
    assert mesh.is_winding_consistent
    assert mesh.volume > 0

    bbox_min, bbox_max = mesh.bounds
    assert bbox_max[2] - bbox_min[2] == pytest.approx(ONE_LEVEL_MM)
    # Vertices sit at angles 0/60/120/...: vertex-to-vertex along X (2r),
    # flat-to-flat along Y (r*sqrt(3)).
    assert bbox_max[0] - bbox_min[0] == pytest.approx(
        2 * layout.hex_outer_width, rel=1e-6
    )
    assert bbox_max[1] - bbox_min[1] == pytest.approx(
        layout.hex_outer_width * 3**0.5, rel=1e-6
    )

    out_path = export_stl(mesh, tmp_path / "phase3_single_hex.stl")
    reloaded = trimesh.load(out_path)
    assert reloaded.is_watertight
    assert reloaded.is_winding_consistent
    assert reloaded.volume == pytest.approx(mesh.volume, rel=1e-6)


def test_export_rejects_a_non_watertight_mesh(tmp_path) -> None:
    layout = FlowerLayout()
    boundary = layout.ring_vertices(0)
    mesh = build_flat_prism_mesh(boundary, top_z=ONE_LEVEL_MM, bottom_z=0.0)
    # Drop one wall triangle to break watertightness.
    broken = mesh.copy()
    broken.update_faces(list(range(len(broken.faces) - 1)))
    broken.remove_unreferenced_vertices()

    with pytest.raises(MeshValidationError):
        export_stl(broken, tmp_path / "should_not_be_written.stl")
    assert not (tmp_path / "should_not_be_written.stl").exists()
