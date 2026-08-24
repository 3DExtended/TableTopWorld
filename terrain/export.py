"""STL export with pre-export sanity checks on real constructed geometry.

This is the check that would have caught the old system's break: the old
suite only string-matched generated SCAD text, so a mesh that wasn't
actually manifold could still pass every test. Every export here refuses
to write a file unless the mesh is watertight, has consistent winding,
and encloses positive volume.
"""

from __future__ import annotations

import pathlib

import trimesh


class MeshValidationError(RuntimeError):
    pass


def export_stl(mesh: trimesh.Trimesh, path: str | pathlib.Path) -> pathlib.Path:
    if not mesh.is_watertight:
        raise MeshValidationError("mesh is not watertight - cannot export")
    if not mesh.is_winding_consistent:
        raise MeshValidationError("mesh has inconsistent face winding")
    if mesh.volume <= 0:
        raise MeshValidationError(f"mesh has non-positive volume: {mesh.volume}")

    out_path = pathlib.Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mesh.export(str(out_path))
    return out_path
