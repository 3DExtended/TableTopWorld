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
