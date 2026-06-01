import pytest

from terrain.catalog import EdgeProfileCatalog


@pytest.fixture
def catalog() -> EdgeProfileCatalog:
    return EdgeProfileCatalog()


def test_flat_profiles_self_compatible(catalog: EdgeProfileCatalog) -> None:
    assert catalog.is_compatible("flat_ground", "flat_ground")
    assert catalog.is_compatible("flat_middle", "flat_middle")
    assert not catalog.is_compatible("flat_ground", "flat_middle")


def test_slope_pairs_compatible(catalog: EdgeProfileCatalog) -> None:
    assert catalog.is_compatible(
        "slope_up_ground_to_middle", "slope_down_middle_to_ground"
    )
    assert not catalog.is_compatible(
        "slope_up_ground_to_middle", "slope_up_ground_to_middle"
    )


def test_elevation_at_edge(catalog: EdgeProfileCatalog) -> None:
    assert catalog.elevation_at_edge("flat_high") == 8.0
    assert catalog.elevation_at_edge("flat_ground") == 0.0


def test_unknown_profile_raises(catalog: EdgeProfileCatalog) -> None:
    with pytest.raises(ValueError, match="unknown edge profile"):
        catalog.validate_profile_name("not_a_profile")
