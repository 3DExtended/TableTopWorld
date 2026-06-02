// demo_trihex_noise_flower_scurve.scad
// Demo: 7-hex flower with an S-curved height step (sine transition centerline).

include <trihex_noise.scad>;

SCALE = 5.0;
HEX_R = 5.1961525;

scale(SCALE) {
  // Low plateau a, high plateau b; band_w is transition thickness.
  // amp / period / phase_deg bend the step front: wx - amp*sin(wy*360/period + phase) = 0.
  macro_scurve = ["scurve", 0.6, 2.6, 2.2, 5.5, 16, 0];

  hf_trihex_noise_flower(
    hex_r = HEX_R,
    tri_step = 0.26 * 4,
    thickness = 3.5,
    noise = [14, 2.2, 5, 0.55, 2.05, 1337],
    z_quant = 0,
    macro = macro_scurve,
    detail = [6, 0.05, 3, 0.55, 2.05, 5150],
    band_fade = 0.28,
    clip_height = 50,
    bevel_size = HF_BEVEL_SIZE,
    magnets = true
  );
}
