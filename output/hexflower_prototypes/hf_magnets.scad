// hf_magnets.scad — exterior mating magnet bores (matches terrain/constants.py).

HF_MAGNET_RADIUS = 0.53;
HF_MAGNET_DEPTH = 0.15;
HF_MAGNET_HEIGHT_OVER_GROUND = 0.25;
HF_MAGNET_CENTER_Z = HF_MAGNET_HEIGHT_OVER_GROUND + HF_MAGNET_RADIUS;

function hf_edge_angle_deg(edge_2d) =
  atan2(
    edge_2d[1][1] - edge_2d[0][1],
    edge_2d[1][0] - edge_2d[0][0]
  );

// Outward edge on ring hexes: edge midpoint farther from origin than hex center.
function hf_is_exterior_edge(hex_idx, edge_idx, hex_r) =
  (hex_idx >= 1)
    ? let(
        c = hf_flower_pos_xy(hex_idx, hex_r),
        dist_c = sqrt(c[0] * c[0] + c[1] * c[1]),
        edge = hf_hex_edge_line_2d(hex_idx, edge_idx, hex_r),
        mid = [
          (edge[0][0] + edge[1][0]) / 2,
          (edge[0][1] + edge[1][1]) / 2
        ],
        dist_mid = sqrt(mid[0] * mid[0] + mid[1] * mid[1])
      )
      dist_mid > dist_c + 1e-6
    : false;

module hf_magnet_hole_cutter(
  edge_2d,
  magnet_center_z = HF_MAGNET_CENTER_Z,
  magnet_radius = HF_MAGNET_RADIUS,
  magnet_depth = HF_MAGNET_DEPTH
) {
  center = [
    (edge_2d[0][0] + edge_2d[1][0]) / 2,
    (edge_2d[0][1] + edge_2d[1][1]) / 2
  ];
  // Match terrain/edges.py: rotateX(90) then rotateZ(edge_angle) — NOT rotate([0,90,ang])
  // (that form applies Y then Z and leaves the bore skewed / non-circular on the wall).
  ang = hf_edge_angle_deg(edge_2d);
  translate([center[0], center[1], magnet_center_z])
    rotate([0, 0, ang])
      rotate([90, 0, 0])
        cylinder(h = magnet_depth * 4, r = magnet_radius, center = true, $fn = 32);
}

// 18 horizontal bores on ring hex exterior sides (hex 1–6, outward edges only).
module hf_flower_magnet_cutters(
  hex_r,
  magnet_center_z = HF_MAGNET_CENTER_Z,
  magnet_radius = HF_MAGNET_RADIUS,
  magnet_depth = HF_MAGNET_DEPTH
) {
  union() {
    for (hi = [1 : 6])
      for (ei = [0 : 5])
        if (hf_is_exterior_edge(hi, ei, hex_r))
          hf_magnet_hole_cutter(
            hf_hex_edge_line_2d(hi, ei, hex_r),
            magnet_center_z,
            magnet_radius,
            magnet_depth
          );
  }
}
