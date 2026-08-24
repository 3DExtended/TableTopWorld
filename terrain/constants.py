"""Shared magnet and plate constants (from legacy hexagon.py).

Bevel/chamfer constants, terrain-level/role dicts, and magnet-"topping"
indices from the old CSG pipeline were removed here along with that
pipeline itself (bevel/chamfer finishing and magnet topping holes are
out of scope for this redesign; terrain levels are now
terrain.heights.HeightLevels, configurable per tileset).
"""

# Plate / extrusion
DEFAULT_HEX_OUTER_WIDTH = 5.1961525
# Minimum ground-level hex height from the print bed (z=0) so exterior walls reach magnet height.
BASE_PLATE_DEPTH = 2.0
FLOWER_BOTTOM_Z = 0.0

# Magnets (horizontal cylinder through outward vertical wall)
MAGNET_DEPTH = 0.15
MAGNET_RADIUS = 0.53
# Fixed center Z for every exterior mating magnet (above print bed z=0).
MAGNET_HEIGHT_OVER_GROUND = 0.25
MAGNET_CENTER_Z = MAGNET_HEIGHT_OVER_GROUND + MAGNET_RADIUS
