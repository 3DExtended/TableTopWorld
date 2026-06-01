"""Shared magnet, bevel, and feature constants (from legacy hexagon.py)."""

# Plate / extrusion
DEFAULT_HEX_OUTER_WIDTH = 5.1961525
LEGACY_HEXAGON_HEIGHT = 2.0  # flat_plains parity at ground level

# Magnets
MAGNET_DEPTH = 0.15
MAGNET_RADIUS = 0.53
MAGNET_HEIGHT_OVER_GROUND = 0.25

# Bevels
HEXAGON_BEVEL_SIZE = 1.26

# Roads / water channel carving
STREET_INDENT_HEIGHT = 0.1
STREET_WIDTH_SCALAR = 1.0
WATER_INDENT_HEIGHT = 0.58  # below terrain surface (legacy used -0.58 offset)
WATER_WIDTH_SCALAR = 0.75

# Roles
TERRAIN_LEVELS = ("ground", "middle", "high")
TERRAIN_Z = {"ground": 0.0, "middle": 4.0, "high": 8.0}

ROLES = (
    "standable",
    "slope",
    "cliff",
    "road_channel",
    "water",
)

TOPPING_HEX_INDICES = (1, 2)
