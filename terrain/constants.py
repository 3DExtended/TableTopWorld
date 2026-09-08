"""Physical constants, in millimetres.

Peter's original June 2026 design (printableFiles/hexagonWithRoad.stl,
measured directly: 10 mm plate, 18 wall bores of 5.3 mm at 3.9 mm above
the bed) was authored in "model units" x5. The magnet and plate numbers
here are the printed millimetre values - they describe physical parts
(5 x 2 mm disc magnets, a print bed) and must NOT scale with
tileset.meta.scale the way the hex layout does. Until 2026-09 they were
kept in model units and never multiplied by scale, which silently
produced a 2 mm plate with 1 mm "magnets" - the reason this file spells
the unit out in every name.
"""

# Hex layout (model units - multiplied by tileset.meta.scale in Tileset.layout()).
DEFAULT_HEX_OUTER_WIDTH = 5.1961525

# Distance from the level-0 terrain surface down to the print bed: the
# flat "basement" every tile shares, where the mating magnets live at one
# fixed height regardless of the terrain above (decision #13).
BASE_PLATE_DEPTH_MM = 10.0

# The physical magnets Peter uses: 5 mm diameter x 2 mm thick discs.
MAGNET_DIAMETER_MM = 5.0
MAGNET_THICKNESS_MM = 2.0

# Blind bore that receives one disc: 0.3 mm diametral clearance (legacy
# 0.53 x 5 = 5.3 mm, which fit his prints), 0.2 mm extra depth so the
# disc sits flush or just below the wall face.
MAGNET_BORE_RADIUS_MM = (MAGNET_DIAMETER_MM + 0.3) / 2
MAGNET_BORE_DEPTH_MM = MAGNET_THICKNESS_MM + 0.2

# Bore centre height above the print bed (legacy 0.78 x 5 = 3.9 mm) - the
# same on every tile so any two tiles mate at any height difference.
MAGNET_CENTER_Z_MM = 3.9
