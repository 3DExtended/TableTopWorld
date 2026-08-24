"""Configurable discrete terrain height levels - design decision #6.

Replaces the old hardcoded TERRAIN_LEVELS/TERRAIN_Z 3-level dict with a
plain arithmetic level -> Z mapping whose level count is configurable per
tileset (decision #6: "more A but I would like to add more later on").
"""

from __future__ import annotations

from dataclasses import dataclass

DEFAULT_MM_PER_LEVEL = 15.0
DEFAULT_LEVEL_COUNT = 4


@dataclass(frozen=True)
class HeightLevels:
    mm_per_level: float = DEFAULT_MM_PER_LEVEL
    level_count: int = DEFAULT_LEVEL_COUNT

    def __post_init__(self) -> None:
        if self.level_count < 1:
            raise ValueError("level_count must be >= 1")
        if self.mm_per_level <= 0:
            raise ValueError("mm_per_level must be > 0")

    def z(self, level: int) -> float:
        if level not in range(self.level_count):
            raise ValueError(
                f"level {level} out of range 0..{self.level_count - 1}"
            )
        return level * self.mm_per_level
