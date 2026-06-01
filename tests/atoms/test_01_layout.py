"""Layer 1 — layout atoms: footprints, junctions, exterior edges (no terrain CSG)."""

from __future__ import annotations

import math

from solid2 import color, cylinder, linear_extrude, polygon, translate, union

from terrain.layout import FlowerLayout


def test_01_center_hex_footprint(layout: FlowerLayout, visual) -> None:
    """Center cell (hex 0) footprint: subdivided hex polygon at Z=0."""
    solid = linear_extrude(height=0.3)(polygon(points=layout.cell_polygon(0)))
    path = visual(
        "01_center_hex_footprint",
        solid,
        "Thin extrusion of center hex 0 polygon — should match inner flower shape.",
        subdir="01_layout",
    )
    assert path.name == "01_center_hex_footprint.scad"
    assert len(layout.cell_polygon(0)) == 18


def test_02_all_seven_footprints_stacked(layout: FlowerLayout, visual) -> None:
    """All seven cell footprints offset in Z so each hex is visible in preview."""
    parts = []
    for hex_idx in range(7):
        foot = linear_extrude(height=0.25)(
            polygon(points=layout.cell_polygon(hex_idx))
        )
        parts.append(translate([0, 0, hex_idx * 0.4])(foot))
    solid = union()(*parts)
    path = visual(
        "02_seven_footprints_stacked",
        solid,
        "Seven hex footprints stacked — verify ring surrounds center with no gaps.",
        subdir="01_layout",
    )
    assert len(parts) == 7
    assert path.exists()


def test_03_junction_centers_marked(layout: FlowerLayout, visual) -> None:
    """Small cylinders at each of the six road/water junction centers."""
    markers = []
    for j in range(6):
        cx, cy = layout.junction_center(j)
        markers.append(
            translate([cx, cy, 0])(cylinder(h=2.0, r=0.35, center=False))
        )
    solid = union()(*markers)
    visual(
        "03_junction_centers",
        solid,
        "Six junction markers — each should sit where three hex sides meet.",
        subdir="01_layout",
    )
    for j in range(6):
        c = layout.junction_center(j)
        assert math.hypot(c[0], c[1]) < layout.hex_outer_width * 3


def test_04_exterior_edges_as_wires(layout: FlowerLayout, visual) -> None:
    """Thin prisms on all 18 outward ring-hex edges (midpoint farther from center than cell)."""
    origin = (0.0, 0.0)
    for edge in layout.exterior_edges():
        center = layout.cell_center(edge.hex_idx)
        p1, p2 = edge.line_2d
        mid = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
        assert (mid[0] - origin[0]) ** 2 + (mid[1] - origin[1]) ** 2 > (
            center[0] - origin[0]
        ) ** 2 + (center[1] - origin[1]) ** 2
    bars = []
    for edge in layout.exterior_edges():
        p1, p2 = edge.line_2d
        dx, dy = p2[0] - p1[0], p2[1] - p1[1]
        length = math.hypot(dx, dy) or 1.0
        angle = math.degrees(math.atan2(dy, dx))
        mid = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
        bar = linear_extrude(height=1.5, center=True)(
            polygon(points=[(-length / 2, -0.08), (length / 2, -0.08), (length / 2, 0.08), (-length / 2, 0.08)])
        )
        bars.append(
            translate([mid[0], mid[1], 0.75])(
                color("red")(bar.rotateZ(angle))
            )
        )
    solid = union()(*bars)
    visual(
        "04_exterior_edge_wires",
        solid,
        "18 red edge bars on ring hex outward sides — magnet/bevel targets.",
        subdir="01_layout",
    )
    assert len(layout.exterior_edges()) == 18
