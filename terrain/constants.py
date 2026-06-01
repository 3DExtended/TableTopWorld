"""Shared magnet, bevel, and feature constants (from legacy hexagon.py)."""

# Plate / extrusion
DEFAULT_HEX_OUTER_WIDTH = 5.1961525
# Minimum ground-level hex height from the print bed (z=0) so exterior walls reach magnet height.
BASE_PLATE_DEPTH = 2.0
FLOWER_BOTTOM_Z = 0.0
LEGACY_HEXAGON_HEIGHT = BASE_PLATE_DEPTH  # alias for callers expecting legacy name

# Magnets (horizontal cylinder through outward vertical wall)
MAGNET_DEPTH = 0.15
MAGNET_RADIUS = 0.53
# Fixed center Z for every exterior mating magnet (above print bed z=0).
MAGNET_HEIGHT_OVER_GROUND = 0.25
MAGNET_CENTER_Z = MAGNET_HEIGHT_OVER_GROUND + MAGNET_RADIUS

# Bevels (model units; × meta.scale for mm — 0.7 × 5 ≈ 3.5 mm chamfer depth)
HEXAGON_BEVEL_SIZE = 0.7
# Side wedge sits this far below z_anchor to avoid coplanar fights on the vertical face.
BEVEL_Z_INSET = 0.02
# Top rim box extends ± this depth around z_anchor so the flat face fully clears in CSG.
BEVEL_TOP_SLAB_DEPTH = 0.05

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
