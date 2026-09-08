"""Idealised GOAL geometry for mockups - deliberately NOT the generator.

Numbers come from the code/legacy STL where they exist (hex size, plate,
magnet position) and from this session's agreed choices (plateau 55%,
5x2 magnets, thin V-lines); the rest are proposals to be confirmed.
"""
import math, sys
import numpy as np
from scipy.ndimage import map_coordinates
sys.path.insert(0, "/Users/peteresser/Developer/projects/archive/tableTopWorld")
from terrain.layout import FlowerLayout

R = 5.1961525 * 5              # hex circumradius 25.98 mm (code: hex_outer_width * scale)
A = R * math.sqrt(3) / 2       # apothem 22.5 mm -> 45 mm across flats
PLATE = 10.0                   # legacy hexagon_height 2.0 * 5
STEP = 15.0                    # tileset height_step_mm
LEVELS = 4
LINE_W, LINE_D = 1.0, 0.6      # engraved V-line: total width / depth (each hex carries half)
BORE_R, BORE_DEPTH, BORE_Z = 2.65, 2.2, 3.9   # 5x2 disc -> 5.3 x 2.2 bore, centre 3.9 above bed (legacy)
SOCKET_R, SOCKET_DEPTH = 2.65, 2.2
ROAD_HALF_W, ROAD_DEPTH = 8.0, 1.0        # 16 mm wide
RIVER_HALF_W, RIVER_DEPTH = 11.0, 2.5     # 22 mm wide, flat bed with 5 mm banks
RIVER_BANK = 5.0
PLATEAU_R = math.sqrt(0.55)    # flat inner hexagon: 55% of area -> 0.742 of the size
ORGANIC_R = 0.25               # organic hexes: broad rolling blend

LAYOUT = FlowerLayout(R)
EDGE_NORMALS = [(math.cos(math.radians(60 * k + 30)), math.sin(math.radians(60 * k + 30))) for k in range(6)]
HEX_VERTS = [(R * math.cos(math.radians(60 * i)), R * math.sin(math.radians(60 * i))) for i in range(6)]


def level_z(level):
    return PLATE + STEP * level


class Noise:
    def __init__(self, seed, cell=9.0, octaves=2):
        rng = np.random.default_rng(seed)
        self.grids = [rng.random((64, 64)) for _ in range(octaves)]
        self.cell = cell

    def __call__(self, P):
        out = np.zeros(len(P)); amp = 1.0; tot = 0.0
        for i, g in enumerate(self.grids):
            c = self.cell / (2 ** i)
            coords = np.vstack([(P[:, 1] / c) % 64, (P[:, 0] / c) % 64])
            out += amp * map_coordinates(g, coords, order=3, mode="wrap")
            tot += amp; amp *= 0.5
        return (out / tot - 0.5) * 2.0


def catmull_rom(points, samples_per_seg=40):
    pts = np.asarray(points, float)
    ext = np.vstack([pts[0] - (pts[1] - pts[0]), pts, pts[-1] + (pts[-1] - pts[-2])])
    out = []
    for i in range(1, len(ext) - 2):
        p0, p1, p2, p3 = ext[i - 1], ext[i], ext[i + 1], ext[i + 2]
        for t in np.linspace(0, 1, samples_per_seg, endpoint=False):
            out.append(0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t * t
                              + (-p0 + 3 * p1 - 3 * p2 + p3) * t ** 3))
    out.append(ext[-2])
    return np.array(out)


class Path:
    """A road or river centreline: dense polyline + arc length + smoothed bed."""

    def __init__(self, control_points):
        self.poly = catmull_rom(control_points)
        seg = np.linalg.norm(np.diff(self.poly, axis=0), axis=1)
        self.s = np.concatenate([[0.0], np.cumsum(seg)])
        self.bed_z = None

    def set_bed(self, terrain_z_fn, window_mm=12.0):
        z = terrain_z_fn(self.poly)
        ds = np.mean(np.diff(self.s))
        k = max(1, int(window_mm / ds))
        kernel = np.ones(2 * k + 1) / (2 * k + 1)
        zp = np.pad(z, k, mode="edge")
        self.bed_z = np.convolve(zp, kernel, mode="valid")

    def dist_s(self, P):
        a = self.poly[:-1]; b = self.poly[1:]
        best = np.full(len(P), np.inf); best_s = np.zeros(len(P))
        for i in range(len(a)):
            ab = b[i] - a[i]; L2 = ab @ ab
            t = np.clip(((P - a[i]) @ ab) / max(L2, 1e-9), 0, 1)
            d = np.linalg.norm(P - (a[i] + t[:, None] * ab), axis=1)
            better = d < best
            best[better] = d[better]; best_s[better] = self.s[i] + t[better] * math.sqrt(L2)
        return best, best_s

    def bed(self, s):
        return np.interp(s, self.s, self.bed_z)


class Scene:
    def __init__(self, flowers, roads=(), rivers=(), seed=3):
        """flowers: list of dict(origin=(x,y), levels={0..6:int}, standable=set)."""
        centers, levels, standable, flower_id, local = [], [], [], [], []
        for fi, f in enumerate(flowers):
            ox, oy = f["origin"]
            for h in range(7):
                cx, cy = LAYOUT.cell_center(h)
                centers.append((cx + ox, cy + oy)); levels.append(f["levels"][h])
                standable.append(h in f["standable"]); flower_id.append(fi); local.append(h)
        self.flowers = flowers
        self.centers = np.array(centers); self.levels = np.array(levels, float)
        self.standable = np.array(standable); self.flower_id = np.array(flower_id); self.local = np.array(local)
        H = len(centers)
        self.nbr_level = np.tile(self.levels[:, None], (1, 6))
        for i in range(H):
            for k, (nx, ny) in enumerate(EDGE_NORMALS):
                target = self.centers[i] + 2 * A * np.array([nx, ny])
                d = np.linalg.norm(self.centers - target, axis=1)
                j = int(d.argmin())
                if d[j] < 1.0:
                    self.nbr_level[i, k] = self.levels[j]
        self.noise_a = Noise(seed, cell=8.0, octaves=3)
        self.noise_b = Noise(seed + 11, cell=30.0, octaves=2)
        self.roads = [Path(p) for p in roads]
        self.rivers = [Path(p) for p in rivers]
        for p in self.roads + self.rivers:
            p.set_bed(lambda Q: self.terrain_z(Q)[0])
        # socket rule: flat hex whose centre is NOT crossed by a road/river
        self.socket_hexes = []
        for i in np.where(self.standable)[0]:
            c = self.centers[i][None, :]
            clear = all(p.dist_s(c)[0][0] > 8.0 for p in self.roads + self.rivers)
            if clear:
                self.socket_hexes.append(int(i))

    def hex_of(self, P):
        d2 = ((P[:, None, :] - self.centers[None]) ** 2).sum(-1)
        return d2.argmin(1)

    def terrain_z(self, P):
        P = np.asarray(P, float)
        hi = self.hex_of(P)
        c = self.centers[hi]; lvl = self.levels[hi]; st = self.standable[hi]
        r0 = np.where(st, PLATEAU_R, ORGANIC_R)
        wob = 0.12 * self.noise_b(P)
        num = lvl.copy(); den = np.ones(len(P)); wmax = np.zeros(len(P))
        for k, (nx, ny) in enumerate(EDGE_NORMALS):
            n = ((P[:, 0] - c[:, 0]) * nx + (P[:, 1] - c[:, 1]) * ny) / A
            n = n + wob * (1.0 - np.clip(n, 0, 1) ** 2)
            t = np.clip((n - r0) / (1.0 - r0), 0, 1)
            w = t * t * t * (t * (t * 6 - 15) + 10)
            num += w * self.nbr_level[hi, k]; den += w; wmax = np.maximum(wmax, w)
        z = level_z(num / den)
        amp = np.where(st, 0.5 * wmax, 0.5)
        z = z + amp * self.noise_a(P)
        return z, hi, wmax

    def surface(self, P):
        """z, material (0 terrain, 1 plateau, 2 road, 3 river, 4 socket)."""
        P = np.asarray(P, float)
        z, hi, wmax = self.terrain_z(P)
        mat = np.where(self.standable[hi] & (wmax <= 1e-6), 1, 0)
        for road in self.roads:
            d, s = road.dist_s(P)
            bed = road.bed(s) - ROAD_DEPTH
            inside = d <= ROAD_HALF_W
            z[inside] = bed[inside]; mat[inside] = 2
            edge = (d > ROAD_HALF_W) & (d < ROAD_HALF_W + 0.6)
            t = (d[edge] - ROAD_HALF_W) / 0.6
            z[edge] = bed[edge] * (1 - t) + z[edge] * t
        for river in self.rivers:
            d, s = river.dist_s(P)
            bank = np.clip((RIVER_HALF_W - d) / RIVER_BANK, 0, 1)
            bank = bank * bank * (3 - 2 * bank)
            prof = river.bed(s) - RIVER_DEPTH * bank
            inside = d <= RIVER_HALF_W
            z[inside] = np.minimum(z[inside], prof[inside]); mat[inside] = 3
        for i in self.socket_hexes:
            d = np.linalg.norm(P - self.centers[i], axis=1)
            sock = d <= SOCKET_R
            z[sock] -= SOCKET_DEPTH; mat[sock] = 4
        return z, mat


def side_mid(flower_origin, side_idx, outward_mm=0.0):
    """Midpoint of the middle edge of side k (== legacy junction_center): the
    road/river crossing point shared by exactly two flowers."""
    m = LAYOUT.junction_center(side_idx)
    ang = math.atan2(*reversed(LAYOUT.neighbor_flower_offset(side_idx)))
    return (m[0] + flower_origin[0] + outward_mm * math.cos(ang), m[1] + flower_origin[1] + outward_mm * math.sin(ang))


def hex_center(flower_origin, h):
    c = LAYOUT.cell_center(h)
    return (c[0] + flower_origin[0], c[1] + flower_origin[1])


# ---------- meshing ----------

def wedge_lattice(center, n):
    """All 6 wedges of one hex as (points, triangles)."""
    cx, cy = center
    pts = []; tris = []
    for k in range(6):
        v0 = np.array(HEX_VERTS[k]); v1 = np.array(HEX_VERTS[(k + 1) % 6])
        base = len(pts)
        index = {}
        for I in range(n + 1):
            for J in range(n + 1 - I):
                index[(I, J)] = len(pts)
                pts.append((cx + v0[0] * I / n + v1[0] * J / n, cy + v0[1] * I / n + v1[1] * J / n))
        for I in range(n):
            for J in range(n - I):
                tris.append((index[(I, J)], index[(I + 1, J)], index[(I, J + 1)]))
                if I + J < n - 1:
                    tris.append((index[(I + 1, J)], index[(I + 1, J + 1)], index[(I, J + 1)]))
    return np.array(pts), np.array(tris)


def build_scene_mesh(scene, n=32):
    """Returns list of (polys Nx3x3, colors Nx3) groups + line segments for V-lines."""
    groups = []
    lines = []
    for fi, f in enumerate(scene.flowers):
        for h in range(7):
            c = hex_center(f["origin"], h)
            P, T = wedge_lattice(c, n)
            z, mat = scene.surface(P)
            V = np.c_[P, z]
            tri = V[T]
            m = mat[T].max(axis=1)
            groups.append((tri, m))
            # V-lines: every hex edge, lifted to the surface
            for k in range(6):
                a = np.array(HEX_VERTS[k]) + c; b = np.array(HEX_VERTS[(k + 1) % 6]) + c
                t = np.linspace(0, 1, 60)[:, None]
                Q = a + (b - a) * t
                zq, _ = scene.surface(Q)
                lines.append(np.c_[Q, zq + 0.12])
    return groups, lines


def exterior_edges_world(flower_origin):
    out = []
    for e in LAYOUT.exterior_edges():
        (x1, y1), (x2, y2) = e.line_2d
        out.append(((x1 + flower_origin[0], y1 + flower_origin[1]), (x2 + flower_origin[0], y2 + flower_origin[1])))
    return out


def build_walls(scene, m=24):
    """Wall quads (plate band + terrain band) and magnet bore discs."""
    plate_quads, terrain_quads, bores = [], [], []
    for f in scene.flowers:
        ox, oy = f["origin"]
        for (x1, y1), (x2, y2) in exterior_edges_world((ox, oy)):
            t = np.linspace(0, 1, m + 1)
            Q = np.c_[x1 + (x2 - x1) * t, y1 + (y2 - y1) * t]
            z, _ = scene.surface(Q)
            for i in range(m):
                zt0, zt1 = z[i], z[i + 1]
                p0, p1 = Q[i], Q[i + 1]
                zp0, zp1 = min(PLATE, zt0), min(PLATE, zt1)
                plate_quads.append([(p0[0], p0[1], 0), (p1[0], p1[1], 0), (p1[0], p1[1], zp1), (p0[0], p0[1], zp0)])
                if zt0 > PLATE or zt1 > PLATE:
                    terrain_quads.append([(p0[0], p0[1], zp0), (p1[0], p1[1], zp1), (p1[0], p1[1], zt1), (p0[0], p0[1], zt0)])
            mid = np.array([(x1 + x2) / 2, (y1 + y2) / 2])
            tang = np.array([x2 - x1, y2 - y1]); tang /= np.linalg.norm(tang)
            outward = np.array([tang[1], -tang[0]])
            if outward @ (mid - np.array([ox, oy])) < 0: outward = -outward
            ang = np.linspace(0, 2 * np.pi, 28)
            circ = [(mid[0] + 0.25 * outward[0] + BORE_R * math.cos(a) * tang[0],
                     mid[1] + 0.25 * outward[1] + BORE_R * math.cos(a) * tang[1],
                     BORE_Z + BORE_R * math.sin(a)) for a in ang]
            bores.append(circ)
    return plate_quads, terrain_quads, bores
