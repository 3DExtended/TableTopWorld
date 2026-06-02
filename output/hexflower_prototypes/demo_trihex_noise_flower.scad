// demo_trihex_noise_flower.scad
// Demo: 7-hex flower paddle using fine triangle grid + world-aligned noise.

include <trihex_noise.scad>;

SCALE = 5.0;
HEX_R = 5.1961525;

scale(SCALE) {
  // Flower demo: macro is a world-aligned step band, noise only on that band.
  // (All 7 hexes use the same world-aligned step definition, so edges match.)
  ang = -10;
  off = 0; // step band centered on world origin line

  hf_trihex_noise_flower(
    hex_r = HEX_R,
    // Coarse mesh (≈4× base step); flower enforces a minimum of hex_r/10 for clean junctions.
    tri_step = 0.26 * 4,
    base_z = 1.2,
    thickness = 3.5,
    noise = [14, 2.2, 5, 0.55, 2.05, 1337],
    z_quant = 0,
    macro = ["step", 0.6, 2.6, 2.0, ang, off],
    // Detail only on the slope band; low amp avoids faceted "zipper" on the step face.
    detail = [6, 0.04, 3, 0.55, 2.05, 4242],
    band_fade = 0.28,
    clip_height = 50
  );
}

