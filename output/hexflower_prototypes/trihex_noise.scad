// trihex_noise.scad
// Finer low-poly triangle grid inside a hex, displaced by Perlin-like value noise.
//
// Key idea:
// - Generate a fine triangular heightfield over a bounding rectangle.
// - Make it a CLOSED solid (top + bottom + sides).
// - Clip it to a hex prism via intersection() => clean hex boundary.
// - Use world-aligned noise coordinates so adjacent hexes can match along edges.
//
// Usage (single hex):
//   include <trihex_noise.scad>;
//   hf_trihex_noise_solid(
//     hex_r = 5.1961525,
//     world_xy = [0,0],        // hex center in world coords
//     tri_step = 0.35,         // smaller => finer triangles
//     thickness = 2,
//     noise = [18, 1.0, 4, 0.55, 2.0, 1337], // [freq, amp, octaves, persistence, lacunarity, seed]
//     z_quant = 0.25           // optional low-poly banding
//   );
//
// Usage (7-hex flower paddle):
//   include <trihex_noise.scad>;
//   hf_trihex_noise_flower(
//     hex_r = 5.1961525,
//     tri_step = 0.28,
//     thickness = 3,
//     noise = [14, 2.2, 5, 0.55, 2.05, 1337],
//     z_quant = 0.22
//   );

$fn = 18;

// All terrain sits on the print bed; bottom Z is not configurable.
HF_FLOWER_BOTTOM_Z = 0;

include <hf_bevel.scad>
include <hf_magnets.scad>

// --- math helpers ---
function hf_clamp(x, a, b) = x < a ? a : (x > b ? b : x);
function hf_lerp(a, b, t) = a + (b - a) * t;
function hf_fract(x) = x - floor(x);
function hf_smoothstep(t) = t * t * (3 - 2 * t);

// Deterministic hash -> [0..1).
function hf_hash2(ix, iy, seed) =
  hf_fract(sin((ix * 127.1 + iy * 311.7 + seed * 74.7)) * 43758.5453123);

// 2D value noise, smooth interpolation, range roughly [-1..1].
function hf_value_noise2(x, y, seed) =
  let(
    ix = floor(x), iy = floor(y),
    fx = x - ix,   fy = y - iy,
    u = hf_smoothstep(fx),
    v = hf_smoothstep(fy),
    a = hf_hash2(ix,     iy,     seed),
    b = hf_hash2(ix + 1, iy,     seed),
    c = hf_hash2(ix,     iy + 1, seed),
    d = hf_hash2(ix + 1, iy + 1, seed),
    ab = hf_lerp(a, b, u),
    cd = hf_lerp(c, d, u)
  )
  (hf_lerp(ab, cd, v) * 2 - 1);

function hf_fbm2_r(x, y, seed, octaves, persistence, lacunarity, i, amp, freq, acc) =
  (i >= octaves) ? acc :
  hf_fbm2_r(
    x, y, seed, octaves, persistence, lacunarity,
    i + 1,
    amp * persistence,
    freq * lacunarity,
    acc + amp * hf_value_noise2(x * freq, y * freq, seed + i * 101)
  );

function hf_fbm2(x, y, seed, octaves, persistence, lacunarity) =
  let(
    // `sum()` isn't available in some OpenSCAD builds; use recursion instead.
    _o = max(1, floor(octaves))
  )
  hf_fbm2_r(x, y, seed, _o, persistence, lacunarity, 0, 1, 1, 0);

function hf_get(arr, idx, fallback) = (is_undef(arr) || idx >= len(arr)) ? fallback : arr[idx];

function hf_quantize(z, step) =
  is_undef(z) ? 0 :
  (is_undef(step) || step <= 0) ? z :
  round(z / step) * step;

// Hex polygon points (flat-top orientation).
function hf_hex_pts(r) = [
  [ r, 0 ],
  [ r/2,  r*sqrt(3)/2 ],
  [ -r/2, r*sqrt(3)/2 ],
  [ -r, 0 ],
  [ -r/2, -r*sqrt(3)/2 ],
  [ r/2,  -r*sqrt(3)/2 ]
];

module hf_hex_prism(r, h) {
  linear_extrude(height = h, convexity = 10)
    polygon(points = hf_hex_pts(r));
}

// 7-hex flower layout (matches `terrain/layout.py` ring center math).
function hf_flower_inner_r(hex_r) = hex_r * sqrt(3) / 2;
function hf_flower_ring_spacing(hex_r) = hf_flower_inner_r(hex_r) * 2;

function hf_flower_pos_xy(hex_idx, hex_r) =
  (hex_idx == 0) ? [0, 0] :
  let(
    ring_index = hex_idx - 1,
    ang = 60 * ring_index + 30,
    d = hf_flower_ring_spacing(hex_r)
  )
  [d * cos(ang), d * sin(ang)];

function hf_flower_extent(hex_r) =
  hf_flower_ring_spacing(hex_r) + hex_r;

// Point inside flat-top hex centered at (cx, cy) with vertex radius r.
function hf_in_hex_at(px, py, cx, cy, r) =
  let(
    pts = hf_hex_pts(r),
    xs = [for (p = pts) p[0] + cx],
    ys = [for (p = pts) p[1] + cy]
  )
  len([
    for (i = [0 : 5])
      if ((xs[(i + 1) % 6] - xs[i]) * (py - ys[i]) - (ys[(i + 1) % 6] - ys[i]) * (px - xs[i]) < -1e-6)
        i
  ]) == 0;

function hf_in_flower(px, py, hex_r) =
  len([
    for (hi = [0 : 6])
      let(c = hf_flower_pos_xy(hi, hex_r))
      if (hf_in_hex_at(px, py, c[0], c[1], hex_r))
        hi
  ]) > 0;

function hf_macro_mode(macro) = hf_get(macro, 0, "standable");

function hf_macro_is_ruled(macro) =
  let(m = hf_macro_mode(macro))
  m == "standable" || m == "step" || m == "scurve";

// Signed phase across a height transition (0 = center of slope band).
// step:  straight world line; scurve: wx - amp*sin(wy/period + phase).
function hf_macro_phase_u(macro, wx, wy) =
  (hf_macro_mode(macro) == "step")
    ? let(
        ang = hf_get(macro, 4, 0),
        off = hf_get(macro, 5, 0),
        nx = cos(ang),
        ny = sin(ang)
      )
      (wx * nx + wy * ny) - off
    : (hf_macro_mode(macro) == "scurve")
      ? let(
          amp = hf_get(macro, 4, 0),
          period = max(1e-9, hf_get(macro, 5, 1)),
          phase = hf_get(macro, 6, 0)
        )
        wx - amp * sin(wy * 360 / period + phase)
      : 0;

function hf_macro_band_w(macro) =
  (hf_macro_mode(macro) == "step" || hf_macro_mode(macro) == "scurve")
    ? max(1e-9, hf_get(macro, 3, 1))
    : 0;

// Macro height at world XY (no detail noise) — used for terrain and bevel anchors.
function hf_macro_z_world(macro, wx, wy) =
  (hf_macro_mode(macro) == "standable")
    ? hf_get(macro, 1, 0)
    : (hf_macro_mode(macro) == "step" || hf_macro_mode(macro) == "scurve")
      ? let(
          a = hf_get(macro, 1, 0),
          b = hf_get(macro, 2, 0),
          band_w = hf_macro_band_w(macro),
          u = hf_macro_phase_u(macro, wx, wy),
          t = hf_clamp((u + band_w / 2) / band_w, 0, 1),
          s = hf_smoothstep(t)
        )
        hf_lerp(a, b, s)
      : 0;

// Detail noise mask: 1 inside slope band, fades to 0 at band edges.
function hf_macro_band_mask(macro, wx, wy, band_fade = 0.15) =
  (hf_macro_mode(macro) == "step" || hf_macro_mode(macro) == "scurve")
    ? let(
        band_w = hf_macro_band_w(macro),
        u = abs(hf_macro_phase_u(macro, wx, wy)),
        inner = band_w / 2,
        fade = max(1e-9, band_w * hf_clamp(band_fade, 0, 0.49)),
        t = hf_clamp((u - inner) / fade, 0, 1)
      )
      (1 - hf_smoothstep(t))
    : 0;

function hf_hex_edge_line_2d(hex_idx, edge_idx, hex_r) =
  let(
    c = hf_flower_pos_xy(hex_idx, hex_r),
    pts = [for (p = hf_hex_pts(hex_r)) [p[0] + c[0], p[1] + c[1]]]
  )
  [pts[edge_idx], pts[(edge_idx + 1) % 6]];

// --- rectangular triangular heightfield solid ---
// Builds a closed solid from (x0..x1, y0..y1) with triangular tessellation.
// The top surface is displaced by noise; bottom is flat at HF_FLOWER_BOTTOM_Z (z=0).
//
// Parameters:
// - tri_step: approximate triangle edge length on the grid
// - noise: [freq, amp, octaves, persistence, lacunarity, seed]
module hf_tri_heightfield_rect_solid(
  x0, x1, y0, y1,
  tri_step,
  world_xy = [0, 0],
  thickness = 2,
  noise = [12, 1.0, 4, 0.5, 2.0, 0],
  // Legacy: quantize final height directly (kept for compatibility)
  z_quant = 0,
  // --- NEW TERRAIN RULE API ---
  // A hex is either:
  // - standable (perfectly flat, no noise), or
  // - stepped from height a->b where ONLY the slope band gets noise.
  //
  // `macro` encodes the macro height profile (relative to 0):
  // - ["standable", z]                         => z_macro is constant z everywhere, no noise.
  // - ["step", a, b, band_w, angle_deg, offset_world]
  //     transition centered on dot([wx,wy],[cos,sin]) - offset_world = 0, width band_w.
  // - ["scurve", a, b, band_w, amp, period, phase_deg]
  //     transition centered on wx - amp*sin(wy*360/period + phase_deg) = 0 (S-shaped front).
  //
  // Detail noise is applied ONLY within the transition band.
  macro = ["standable", 0],
  detail = [5, 0.22, 3, 0.55, 2.05, 4242], // [freq_mult, amp_mult, oct, pers, lac, seed_offset]
  band_fade = 0.15, // fraction of band_w used to fade noise mask to zero at band edges

  // --- Legacy terrace API (deprecated) ---
  // Kept so existing prototypes still compile, but it violates the new modeling rule because
  // it derives macro levels from noise. Prefer `macro` + `detail`.
  terrace_step = 0,
  detail_noise = [4, 0.35, 3, 0.55, 2.05, 9001],
  mask_power = 1.35,
  // When false, omit the rectangular perimeter walls. Use with a footprint
  // intersection() so CGAL caps the sides cleanly (avoids corner notches).
  include_rect_walls = true,
  // Optional: cull top/bottom triangles outside this footprint (no boolean needed).
  footprint = "none", // "none" | "hex" | "flower"
  footprint_hex_r = 5.1961525
) {
  freq = hf_get(noise, 0, 12);
  amp  = hf_get(noise, 1, 1.0);
  oct  = hf_get(noise, 2, 4);
  pers = hf_get(noise, 3, 0.5);
  lac  = hf_get(noise, 4, 2.0);
  seed = hf_get(noise, 5, 0);

  // New detail noise settings
  d_freq_mult = hf_get(detail, 0, 5);
  d_amp_mult  = hf_get(detail, 1, 0.22);
  d_oct       = hf_get(detail, 2, 3);
  d_pers      = hf_get(detail, 3, 0.55);
  d_lac       = hf_get(detail, 4, 2.05);
  d_seed_off  = hf_get(detail, 5, 4242);

  // Legacy detail noise settings
  legacy_d_freq_mult = hf_get(detail_noise, 0, 4);
  legacy_d_amp_mult  = hf_get(detail_noise, 1, 0.35);
  legacy_d_oct       = hf_get(detail_noise, 2, 3);
  legacy_d_pers      = hf_get(detail_noise, 3, 0.55);
  legacy_d_lac       = hf_get(detail_noise, 4, 2.05);
  legacy_d_seed_off  = hf_get(detail_noise, 5, 9001);

  // Snap horizontal spacing to the bounding rectangle so edges/corners hit grid points.
  span_x = x1 - x0;
  span_y = y1 - y0;
  nx = max(2, ceil(span_x / tri_step) + 1);
  dx = span_x / (nx - 1);
  dy = dx * sqrt(3) / 2;
  ny = max(2, ceil(span_y / dy) + 1);

  function vx(i) = x0 + i * dx;
  function vy(j) = y0 + j * dy;

  // World position (unscaled): used for both macro and detail so tiles align across seams.
  function wpos(x, y) = [x + world_xy[0], y + world_xy[1]];

  // Stagger rows in WORLD space so adjacent flower hexes share identical edge vertices.
  function row_stagger(x, y) =
    let(wy = y + world_xy[1])
    (floor(wy / max(1e-9, dy)) % 2 == 1) ? dx/2 : 0;

  function macro_z(wx, wy) = hf_macro_z_world(macro, wx, wy);

  function band_mask(wx, wy) = hf_macro_band_mask(macro, wx, wy, band_fade);

  function detail_z(wx, wy, mask) =
    (mask <= 0) ? 0 :
    let(
      // Base noise gives a consistent "micro texture" across the world.
      // Note: freq is in world units; we normalize by freq so adjacent hexes match.
      fx = max(1e-9, freq / max(1e-9, d_freq_mult)),
      n = hf_fbm2(wx / fx, wy / fx, seed + d_seed_off, d_oct, d_pers, d_lac)
    )
    (n * (amp * d_amp_mult) * mask);

  function in_footprint(x, y) =
    (footprint == "flower")
      ? hf_in_flower(x, y, footprint_hex_r)
      : (footprint == "hex")
        ? hf_in_hex_at(x, y, world_xy[0], world_xy[1], footprint_hex_r)
        : true;

  function top_z(x, y) =
    let(
      wp = wpos(x, y),
      wx = wp[0],
      wy = wp[1]
    )
    // NEW RULE PATH: explicit macro + masked detail.
    hf_macro_is_ruled(macro)
      ? let(
          mz = macro_z(wx, wy),
          m = band_mask(wx, wy),
          dz = detail_z(wx, wy, m)
        )
        HF_FLOWER_BOTTOM_Z + thickness + mz + dz
      // LEGACY PATH (deprecated): macro derived from noise.
      : let(
          wxn = wx / max(1e-9, freq),
          wyn = wy / max(1e-9, freq),
          n = hf_fbm2(wxn, wyn, seed, oct, pers, lac),
          z0 = n * amp
        )
        (terrace_step > 0)
          ? let(
              t = hf_quantize(z0, terrace_step),
              residual = abs(z0 - t),
              eps = terrace_step * 0.02,
              m0 = (residual <= eps)
                ? 0
                : hf_clamp((residual - eps) / max(1e-9, (terrace_step/2) - eps), 0, 1),
              mask = pow(m0, mask_power),
              fx2 = max(1e-9, freq / max(1e-9, legacy_d_freq_mult)),
              nd = hf_fbm2(wx / fx2, wy / fx2, seed + legacy_d_seed_off, legacy_d_oct, legacy_d_pers, legacy_d_lac),
              detail = nd * (amp * legacy_d_amp_mult) * mask
            )
            HF_FLOWER_BOTTOM_Z + thickness + t + detail
          : HF_FLOWER_BOTTOM_Z + thickness + hf_quantize(z0, z_quant);

  // Vertex indexing: top first, then bottom.
  function idx_top(i, j) = j * nx + i;
  function idx_bot(i, j) = nx * ny + j * nx + i;

  top_verts = [
    for (j = [0 : ny - 1])
      for (i = [0 : nx - 1])
        let(
          y = vy(j),
          x = vx(i) + row_stagger(vx(i), y),
          z = top_z(x, y)
        )
        [x, y, z]
  ];

  bot_verts = [
    for (j = [0 : ny - 1])
      for (i = [0 : nx - 1])
        let(
          y = vy(j),
          x = vx(i) + row_stagger(vx(i), y)
        )
        [x, y, HF_FLOWER_BOTTOM_Z]
  ];

  function tri_keep_verts(idxs) =
    (footprint == "none")
      ? true
      : let(
          v0 = top_verts[idxs[0]],
          v1 = top_verts[idxs[1]],
          v2 = top_verts[idxs[2]]
        )
        in_footprint(v0[0], v0[1])
        || in_footprint(v1[0], v1[1])
        || in_footprint(v2[0], v2[1]);

  // Fixed checkerboard diagonal — height-based flipping leaves a visible seam
  // across sloped macro terrain (looks like holes / missing corners in the flower).
  top_faces = concat(
    [
      for (j = [0 : ny - 2])
        for (i = [0 : nx - 2])
          let(
            a = idx_top(i, j),
            b = idx_top(i + 1, j),
            c = idx_top(i, j + 1),
            d = idx_top(i + 1, j + 1),
            flip = (i + j) % 2 == 1,
            tris = flip
              ? [[a, b, d], [a, d, c]]
              : [[a, b, c], [b, d, c]]
          )
          each [
            for (t = tris)
              if (tri_keep_verts(t))
                t
          ]
    ]
  );

  bottom_faces = [
    for (f = top_faces)
      [ for (k = [len(f) - 1 : -1 : 0]) f[k] + nx * ny ]
  ];

  // Side walls around the rectangular perimeter (as quads split into tris).
  // Left/right
  side_lr = (include_rect_walls) ? concat(
    [
      for (j = [0 : ny - 2])
        let(
          // left edge i=0
          tl0 = idx_top(0, j),
          tl1 = idx_top(0, j + 1),
          bl0 = idx_bot(0, j),
          bl1 = idx_bot(0, j + 1),
          // right edge i=nx-1
          tr0 = idx_top(nx - 1, j),
          tr1 = idx_top(nx - 1, j + 1),
          br0 = idx_bot(nx - 1, j),
          br1 = idx_bot(nx - 1, j + 1)
        )
        each [
          [tl0, tl1, bl1], [tl0, bl1, bl0],
          [tr1, tr0, br0], [tr1, br0, br1]
        ]
    ]
  ) : [];

  // Bottom/top edges (y0/y1)
  side_bt = (include_rect_walls) ? concat(
    [
      for (i = [0 : nx - 2])
        let(
          // bottom row j=0
          tb0 = idx_top(i, 0),
          tb1 = idx_top(i + 1, 0),
          bb0 = idx_bot(i, 0),
          bb1 = idx_bot(i + 1, 0),
          // top row j=ny-1
          tt0 = idx_top(i, ny - 1),
          tt1 = idx_top(i + 1, ny - 1),
          bt0 = idx_bot(i, ny - 1),
          bt1 = idx_bot(i + 1, ny - 1)
        )
        each [
          [tb1, tb0, bb0], [tb1, bb0, bb1],
          [tt0, tt1, bt1], [tt0, bt1, bt0]
        ]
    ]
  ) : [];

  polyhedron(
    points = concat(top_verts, bot_verts),
    faces = concat(top_faces, bottom_faces, side_lr, side_bt),
    convexity = 10
  );
}

// Public module: hex-clipped solid.
module hf_trihex_noise_solid(
  hex_r = 5.1961525,
  world_xy = [0, 0],
  tri_step = 0.35,
  thickness = 2,
  noise = [12, 1.0, 4, 0.5, 2.0, 0],
  z_quant = 0.25,
  macro = ["standable", 0],
  detail = [5, 0.22, 3, 0.55, 2.05, 4242],
  band_fade = 0.15,
  terrace_step = 0,
  detail_noise = [4, 0.35, 3, 0.55, 2.05, 9001],
  mask_power = 1.35,
  clip_height = 100
) {
  // Bounding box generous enough for flat-top hex.
  x0 = -hex_r;
  x1 =  hex_r;
  y0 = -hex_r * sqrt(3)/2;
  y1 =  hex_r * sqrt(3)/2;

  intersection() {
      hf_tri_heightfield_rect_solid(
        x0, x1, y0, y1,
        tri_step = tri_step,
        world_xy = world_xy,
        thickness = thickness,
        noise = noise,
        z_quant = z_quant,
        macro = macro,
        detail = detail,
        band_fade = band_fade,
        terrace_step = terrace_step,
        detail_noise = detail_noise,
        mask_power = mask_power,
        footprint = "none"
      );
      hf_hex_prism(hex_r, clip_height);
  }
}

module hf_flower_clip_prism(hex_r, h) {
  linear_extrude(height = h, convexity = 10)
    union() for (hi = [0 : 6])
      translate(hf_flower_pos_xy(hi, hex_r))
        polygon(points = hf_hex_pts(hex_r));
}

module hf_flower_bevel_cutters(
  hex_r,
  thickness,
  macro,
  bevel_size = HF_BEVEL_SIZE
) {
  union() {
    for (hi = [0 : 6])
      let(center = hf_flower_pos_xy(hi, hex_r))
      for (ei = [0 : 5])
        let(
          edge = hf_hex_edge_line_2d(hi, ei, hex_r),
          mid = [(edge[0][0] + edge[1][0]) / 2, (edge[0][1] + edge[1][1]) / 2],
          mz = hf_macro_z_world(macro, mid[0], mid[1]),
          z_anchor = HF_FLOWER_BOTTOM_Z + thickness + mz,
          line_3d = [
            [edge[0][0], edge[0][1], z_anchor],
            [edge[1][0], edge[1][1], z_anchor]
          ]
        )
        hf_bevel_cutter(line_3d, bevel_size, z_anchor, center);
  }
}

module hf_trihex_flower_body(
  hex_r,
  tri_step,
  thickness,
  noise,
  z_quant,
  macro,
  detail,
  band_fade,
  terrace_step,
  detail_noise,
  mask_power,
  clip_height
) {
  ext = hf_flower_extent(hex_r);
  _tri_step = min(tri_step, hex_r / 12);
  intersection() {
    hf_tri_heightfield_rect_solid(
      -ext, ext, -ext, ext,
      tri_step = _tri_step,
      world_xy = [0, 0],
      thickness = thickness,
      noise = noise,
      z_quant = z_quant,
      macro = macro,
      detail = detail,
      band_fade = band_fade,
      terrace_step = terrace_step,
      detail_noise = detail_noise,
      mask_power = mask_power,
      footprint = "none"
    );
    hf_flower_clip_prism(hex_r, clip_height);
  }
}

module hf_trihex_noise_flower(
  hex_r = 5.1961525,
  tri_step = 0.28,
  thickness = 2,
  noise = [12, 1.0, 4, 0.5, 2.0, 0],
  z_quant = 0.25,
  macro = ["standable", 0],
  detail = [5, 0.22, 3, 0.55, 2.05, 4242],
  band_fade = 0.15,
  terrace_step = 0,
  detail_noise = [4, 0.35, 3, 0.55, 2.05, 9001],
  mask_power = 1.35,
  clip_height = 100,
  bevel_size = 0,
  magnets = false,
  magnet_center_z = HF_MAGNET_CENTER_Z
) {
  render(convexity = 10)
    if (bevel_size > 0 || magnets) {
      difference() {
        hf_trihex_flower_body(
          hex_r, tri_step, thickness,
          noise, z_quant, macro, detail, band_fade,
          terrace_step, detail_noise, mask_power, clip_height
        );
        if (bevel_size > 0)
          hf_flower_bevel_cutters(hex_r, thickness, macro, bevel_size);
        if (magnets)
          hf_flower_magnet_cutters(hex_r, magnet_center_z);
      }
    } else {
      hf_trihex_flower_body(
        hex_r, tri_step, thickness,
        noise, z_quant, macro, detail, band_fade,
        terrace_step, detail_noise, mask_power, clip_height
      );
    }
}

