"""Edge profile catalog and compatibility rules."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from terrain.constants import TERRAIN_Z


@dataclass(frozen=True)
class ProfileInfo:
    name: str
    level: str  # ground | middle | high
    mates_with: frozenset[str]


class EdgeProfileCatalog:
    """Named exterior profiles and symmetric compatibility."""

    FLAT_GROUND = "flat_ground"
    FLAT_MIDDLE = "flat_middle"
    FLAT_HIGH = "flat_high"

    SLOPE_UP_GROUND_TO_MIDDLE = "slope_up_ground_to_middle"
    SLOPE_DOWN_MIDDLE_TO_GROUND = "slope_down_middle_to_ground"
    SLOPE_UP_MIDDLE_TO_HIGH = "slope_up_middle_to_high"
    SLOPE_DOWN_HIGH_TO_MIDDLE = "slope_down_high_to_middle"
    SLOPE_UP_GROUND_TO_HIGH = "slope_up_ground_to_high"
    SLOPE_DOWN_HIGH_TO_GROUND = "slope_down_high_to_ground"

    CLIFF_GROUND = "cliff_ground"
    CLIFF_MIDDLE = "cliff_middle"
    CLIFF_HIGH = "cliff_high"

    ALL_PROFILES: tuple[str, ...] = (
        FLAT_GROUND,
        FLAT_MIDDLE,
        FLAT_HIGH,
        SLOPE_UP_GROUND_TO_MIDDLE,
        SLOPE_DOWN_MIDDLE_TO_GROUND,
        SLOPE_UP_MIDDLE_TO_HIGH,
        SLOPE_DOWN_HIGH_TO_MIDDLE,
        SLOPE_UP_GROUND_TO_HIGH,
        SLOPE_DOWN_HIGH_TO_GROUND,
        CLIFF_GROUND,
        CLIFF_MIDDLE,
        CLIFF_HIGH,
    )

    def __init__(self) -> None:
        self._profiles: dict[str, ProfileInfo] = {}
        self._build()

    def _pair(self, a: str, b: str) -> frozenset[str]:
        return frozenset({a, b})

    def _build(self) -> None:
        defs: list[tuple[str, str, frozenset[str]]] = [
            (self.FLAT_GROUND, "ground", self._pair(self.FLAT_GROUND, self.FLAT_GROUND)),
            (self.FLAT_MIDDLE, "middle", self._pair(self.FLAT_MIDDLE, self.FLAT_MIDDLE)),
            (self.FLAT_HIGH, "high", self._pair(self.FLAT_HIGH, self.FLAT_HIGH)),
            (
                self.SLOPE_UP_GROUND_TO_MIDDLE,
                "ground",
                self._pair(
                    self.SLOPE_UP_GROUND_TO_MIDDLE,
                    self.SLOPE_DOWN_MIDDLE_TO_GROUND,
                ),
            ),
            (
                self.SLOPE_DOWN_MIDDLE_TO_GROUND,
                "middle",
                self._pair(
                    self.SLOPE_DOWN_MIDDLE_TO_GROUND,
                    self.SLOPE_UP_GROUND_TO_MIDDLE,
                ),
            ),
            (
                self.SLOPE_UP_MIDDLE_TO_HIGH,
                "middle",
                self._pair(
                    self.SLOPE_UP_MIDDLE_TO_HIGH,
                    self.SLOPE_DOWN_HIGH_TO_MIDDLE,
                ),
            ),
            (
                self.SLOPE_DOWN_HIGH_TO_MIDDLE,
                "high",
                self._pair(
                    self.SLOPE_DOWN_HIGH_TO_MIDDLE,
                    self.SLOPE_UP_MIDDLE_TO_HIGH,
                ),
            ),
            (
                self.SLOPE_UP_GROUND_TO_HIGH,
                "ground",
                self._pair(
                    self.SLOPE_UP_GROUND_TO_HIGH,
                    self.SLOPE_DOWN_HIGH_TO_GROUND,
                ),
            ),
            (
                self.SLOPE_DOWN_HIGH_TO_GROUND,
                "high",
                self._pair(
                    self.SLOPE_DOWN_HIGH_TO_GROUND,
                    self.SLOPE_UP_GROUND_TO_HIGH,
                ),
            ),
            (self.CLIFF_GROUND, "ground", frozenset({self.CLIFF_GROUND})),
            (self.CLIFF_MIDDLE, "middle", frozenset({self.CLIFF_MIDDLE})),
            (self.CLIFF_HIGH, "high", frozenset({self.CLIFF_HIGH})),
        ]
        for name, level, mates in defs:
            self._profiles[name] = ProfileInfo(name=name, level=level, mates_with=mates)

    def known_profiles(self) -> frozenset[str]:
        return frozenset(self._profiles.keys())

    def is_compatible(self, profile_a: str, profile_b: str) -> bool:
        if profile_a not in self._profiles or profile_b not in self._profiles:
            return False
        info_a = self._profiles[profile_a]
        if profile_b not in info_a.mates_with:
            return False
        # Flat/cliff profiles only mate with themselves; slope pairs do not self-mate.
        if profile_a == profile_b:
            return len(info_a.mates_with) == 1
        return True

    FLAT_MATING_PROFILES: frozenset[str] = frozenset(
        {FLAT_GROUND, FLAT_MIDDLE, FLAT_HIGH}
    )

    def is_flat_mating_profile(self, profile: str) -> bool:
        return profile in self.FLAT_MATING_PROFILES

    def is_mating_profile(self, profile: str) -> bool:
        """Exterior edges that connect to another flower (flat, slope pairs; not cliffs)."""
        if profile not in self._profiles:
            return False
        return not profile.startswith("cliff_")

    def elevation_at_edge(self, profile: str) -> Optional[float]:
        """Nominal terrain Z at the exterior edge for flat profiles."""
        if profile not in self._profiles:
            return None
        if not self.is_flat_mating_profile(profile):
            return None
        level = self._profiles[profile].level
        return TERRAIN_Z.get(level)

    def validate_profile_name(self, profile: str) -> None:
        if profile not in self._profiles:
            known = ", ".join(sorted(self._profiles))
            raise ValueError(f"unknown edge profile {profile!r}; known: {known}")
