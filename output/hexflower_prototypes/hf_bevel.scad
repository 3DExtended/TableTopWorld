// hf_bevel.scad — per-hex top-rim chamfers (matches terrain/edges.py cutters).
//
// Three-part cutter per edge: side wedge, 45° top shelf, thin top slab.

HF_BEVEL_SIZE = 0.175;
HF_BEVEL_Z_INSET = 0.02;
HF_BEVEL_TOP_SLAB_DEPTH = 0.05;

_BEVEL_POLY_FACES = [
  [0, 2, 4], [1, 5, 3],
  [0, 4, 1], [4, 5, 1],
  [4, 2, 3], [5, 4, 3],
  [2, 0, 1], [3, 2, 1]
];

function hf_v3_add(a, b) = [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
function hf_v3_sub(a, b) = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
function hf_v3_scale(v, s) = [v[0] * s, v[1] * s, v[2] * s];
function hf_v3_len(v) = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
function hf_v3_norm(v) =
  let(l = hf_v3_len(v))
  (l < 1e-9) ? [0, 0, 1] : hf_v3_scale(v, 1 / l);
function hf_v3_cross(a, b) =
  [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0]
  ];

// Shift segment endpoints: x along first XY perpendicular, y along second, z along edge.
function hf_shift_line_3d(line, x, y, z) =
  let(
    p1 = line[0],
    p2 = line[1],
    u = hf_v3_norm(hf_v3_sub(p2, p1)),
    perp1 = (abs(u[0]) < 1e-9 && abs(u[1]) < 1e-9)
      ? [1, 0, 0]
      : hf_v3_norm([-u[1], u[0], 0]),
    perp2 = hf_v3_norm(hf_v3_cross(u, perp1)),
    off = hf_v3_add(hf_v3_add(hf_v3_scale(perp1, x), hf_v3_scale(perp2, y)), hf_v3_scale(u, z))
  )
  [hf_v3_add(p1, off), hf_v3_add(p2, off)];

function hf_bevel_line_xy(line) =
  let(
    p0 = line[0],
    p1 = line[1]
  )
  [
    [
      p0[0] + 0.0001 * (p1[0] - p0[0]),
      p0[1] + 0.0001 * (p1[1] - p0[1]),
      0
    ],
    [
      p1[0] + 0.0001 * (p0[0] - p1[0]),
      p1[1] + 0.0001 * (p0[1] - p1[1]),
      0
    ]
  ];

function hf_bevel_needs_flip(line_xy, toward_xy) =
  let(
    probe = hf_shift_line_3d(line_xy, 0.075, 0, 0),
    p0 = line_xy[0],
    p1 = line_xy[1],
    mid = [(p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2],
    to_target = [toward_xy[0] - mid[0], toward_xy[1] - mid[1]],
    to_in = [probe[0][0] - mid[0], probe[0][1] - mid[1]]
  )
  (to_target[0] * to_in[0] + to_target[1] * to_in[1]) < 0;

module hf_bevel_polyhedron(line_in, line_out, line_apex) {
  polyhedron(
    points = concat(line_in, line_out, line_apex),
    faces = _BEVEL_POLY_FACES,
    convexity = 4
  );
}

module hf_bevel_side_wedge(line, depth, flip = false) {
  sign = flip ? -1 : 1;
  line_in = hf_shift_line_3d(line, depth * 0.75, 0, 0);
  line_out = hf_shift_line_3d(line, -depth * 0.75, 0, 0);
  line_apex = hf_shift_line_3d(line, 0, sign * -depth, 0);
  hf_bevel_polyhedron(line_in, line_out, line_apex);
}

module hf_bevel_top_shelf_wedge(line, depth, flip = false) {
  sign = flip ? -1 : 1;
  line_in = hf_shift_line_3d(line, depth * 0.75, 0, 0);
  line_apex = hf_shift_line_3d(line, depth * 0.75, sign * depth, 0);
  hf_bevel_polyhedron(line_in, line, line_apex);
}

module hf_bevel_top_flat_box(line, depth, flip = false) {
  inward = flip ? -depth * 0.75 : depth * 0.75;
  half = HF_BEVEL_TOP_SLAB_DEPTH / 2;
  edge_lo = hf_shift_line_3d(line, 0, -half, 0);
  edge_hi = hf_shift_line_3d(line, 0, half, 0);
  in_lo = hf_shift_line_3d(line, inward, -half, 0);
  in_hi = hf_shift_line_3d(line, inward, half, 0);
  polyhedron(
    points = concat(edge_lo, in_lo, edge_hi, in_hi),
    faces = [
      [0, 1, 3], [0, 3, 2],
      [4, 6, 7], [4, 7, 5],
      [0, 4, 5], [0, 5, 1],
      [2, 3, 7], [2, 7, 6],
      [0, 2, 6], [0, 6, 4],
      [1, 5, 7], [1, 7, 3]
    ],
    convexity = 4
  );
}

// Single edge chamfer cutter (union of side + top tools).
module hf_bevel_cutter(line_3d, depth, z_anchor, toward_xy) {
  line_xy = hf_bevel_line_xy(line_3d);
  flip = hf_bevel_needs_flip(line_xy, toward_xy);
  z_side = z_anchor - HF_BEVEL_Z_INSET;
  translate([0, 0, z_side])
    hf_bevel_side_wedge(line_xy, depth, flip = flip);
  translate([0, 0, z_anchor])
    union() {
      hf_bevel_top_shelf_wedge(line_xy, depth, flip = flip);
      hf_bevel_top_flat_box(line_xy, depth, flip = flip);
    }
}
