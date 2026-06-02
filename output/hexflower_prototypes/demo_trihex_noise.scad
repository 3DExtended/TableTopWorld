// demo_trihex_noise.scad
// Demo: fine triangle grid + world-aligned noise in a single hex.
//
// Open this file in OpenSCAD and render.

include <trihex_noise.scad>;

SCALE = 5.0;
HEX_R = 5.1961525;

scale(SCALE) {
  // Two hexes: left is standable-flat (no noise), right is stepped (noise only on slope band).
  spacing = HEX_R * 2.15;

  // Left: standable hex (perfectly flat, zero noise anywhere).
  translate([-spacing/2, 0, 0]) {
    color([0.75, 0.75, 0.75])
      hf_trihex_noise_solid(
        hex_r = HEX_R,
        world_xy = [-spacing/2, 0],
        tri_step = 0.28,
        thickness = 3.5,
        noise = [14, 2.2, 5, 0.55, 2.05, 1337],
        z_quant = 0,
        macro = ["standable", 1.0],
        clip_height = 50
      );
  }

  // Right: stepped hex (a->b). Only the transition band gets detail noise.
  translate([spacing/2, 0, 0]) {
    // Step boundary is world-aligned so seams match across tiles.
    // Here we center the band on the local hex center: offset_world = dot(center,[cos,sin]).
    ang = 25;
    c = [spacing/2, 0];
    off = c[0] * cos(ang) + c[1] * sin(ang);

    color([0.75, 0.75, 0.75])
      hf_trihex_noise_solid(
        hex_r = HEX_R,
        world_xy = c,
        tri_step = 0.28,
        thickness = 3.5,
        noise = [14, 2.2, 5, 0.55, 2.05, 1337],
        z_quant = 0,
        macro = ["step", 0.4, 2.4, 1.2, ang, off], // a, b, band_w, angle_deg, offset_world
        detail = [6, 0.16, 3, 0.55, 2.05, 4242],
        band_fade = 0.18,
        clip_height = 50
      );
  }
}

