import pytest

from terrain.layout import FlowerLayout


@pytest.fixture
def layout() -> FlowerLayout:
    return FlowerLayout()


def test_exterior_edge_count(layout: FlowerLayout) -> None:
    assert len(layout.exterior_edge_keys()) == 18


def test_edge_id_roundtrip(layout: FlowerLayout) -> None:
    assert FlowerLayout.edge_id(3, 1) == "3-1"
    assert FlowerLayout.parse_edge_id("3-1") == (3, 1)


def test_junction_lines_count(layout: FlowerLayout) -> None:
    for j in range(6):
        lines = layout.junction_lines(j)
        assert len(lines) == 3
        for line in lines:
            assert len(line[0]) == 2
            assert len(line[1]) == 2


def test_junction_center_inside_flower(layout: FlowerLayout) -> None:
    c = layout.junction_center(0)
    assert abs(c[0]) < 20 and abs(c[1]) < 20


def test_invalid_junction_raises(layout: FlowerLayout) -> None:
    with pytest.raises(ValueError):
        layout.junction_lines(6)


def test_exterior_edges_face_away_from_flower_center(layout: FlowerLayout) -> None:
    """Each exterior edge midpoint is farther from origin than its hex center."""
    origin = (0.0, 0.0)
    for edge in layout.exterior_edges():
        center = layout.cell_center(edge.hex_idx)
        p1, p2 = edge.line_2d
        mid = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
        dist_mid = (mid[0] - origin[0]) ** 2 + (mid[1] - origin[1]) ** 2
        dist_center = (center[0] - origin[0]) ** 2 + (center[1] - origin[1]) ** 2
        assert dist_mid > dist_center + 1e-12, edge.key


def test_exterior_edge_outward_normal(layout: FlowerLayout) -> None:
    """transform_edge_to_world normal points away from the flower origin."""
    for edge in layout.exterior_edges():
        mid, normal, _ = layout.transform_edge_to_world(edge, (0.0, 0.0), 0)
        assert mid[0] * normal[0] + mid[1] * normal[1] > 0, edge.key
