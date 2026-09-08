import sys, math, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
sys.path.insert(0, "/private/tmp/claude-501/-Users-peteresser-Developer-projects-archive-tableTopWorld/987a149e-0f53-4c5f-952a-baba80f83122/scratchpad/mock")
from mock_common import *

TERRAIN_CMAP = LinearSegmentedColormap.from_list("terr", [
    (0.00, (0.33, 0.50, 0.24)),   # level 0: moss green
    (0.35, (0.58, 0.66, 0.34)),   # level 1: grass/olive
    (0.70, (0.78, 0.70, 0.50)),   # level 2: dry / rocky
    (1.00, (0.88, 0.86, 0.80)),   # level 3: rock
])
COL_ROAD = np.array((0.74, 0.62, 0.44))
COL_RIVER = np.array((0.36, 0.56, 0.80))
COL_SOCKET = np.array((0.12, 0.12, 0.12))
COL_PLATE = np.array((0.84, 0.84, 0.82))
COL_TWALL = np.array((0.56, 0.46, 0.36))
COL_BORE = np.array((0.06, 0.06, 0.06))
LIGHT = np.array([-0.45, 0.35, 0.82]); LIGHT /= np.linalg.norm(LIGHT)


def shade(polys, base_cols):
    polys = np.asarray(polys, float)
    n = np.cross(polys[:, 1] - polys[:, 0], polys[:, 2] - polys[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1)[:, None], 1e-12)
    lam = np.clip(n @ LIGHT, 0, 1)
    b = 0.45 + 0.6 * lam
    return np.clip(base_cols * b[:, None], 0, 1)


def face_colors(tri, mat):
    zc = tri[:, :, 2].mean(axis=1)
    lv = np.clip((zc - PLATE) / (STEP * (LEVELS - 1)), 0, 1)
    cols = TERRAIN_CMAP(lv)[:, :3]
    cols[mat == 1] = cols[mat == 1] * 0.85 + 0.15 * np.array((0.75, 0.85, 0.55))  # plateau: slightly lighter meadow
    cols[mat == 2] = COL_ROAD
    cols[mat == 3] = COL_RIVER
    cols[mat == 4] = COL_SOCKET
    return cols


def render(scene, out, n=32, views=(("iso", 32, -55), ("top", 90, -90)), title="", figsize=(11, 11), legend=True, cutaway=None):
    t0 = time.time()
    groups, lines = build_scene_mesh(scene, n=n)
    plate_quads, terrain_quads, bores = build_walls(scene)
    all_polys, all_cols = [], []
    for tri, mat in groups:
        all_polys.append(tri); all_cols.append(face_colors(tri, mat))
    tri = np.concatenate(all_polys); cols = np.concatenate(all_cols)
    tri_cols = shade(tri, cols)

    def quads_to_tris(quads, col):
        q = np.asarray(quads, float)
        t = np.concatenate([q[:, [0, 1, 2]], q[:, [0, 2, 3]]])
        return t, shade(t, np.tile(col, (len(t), 1)))
    pt, pc = quads_to_tris(plate_quads, COL_PLATE)
    tt, tc = quads_to_tris(terrain_quads, COL_TWALL)
    print(f"mesh: {len(tri)} surface tris, {len(pt)+len(tt)} wall tris, built in {time.time()-t0:.1f}s")

    allx = tri[:, :, 0]; ally = tri[:, :, 1]
    cx, cy = (allx.min() + allx.max()) / 2, (ally.min() + ally.max()) / 2
    rad = max(allx.max() - allx.min(), ally.max() - ally.min()) / 2 * 1.02

    fig = plt.figure(figsize=(figsize[0] * len(views), figsize[1]))
    for vi, (name, elev, azim) in enumerate(views):
        ax = fig.add_subplot(1, len(views), vi + 1, projection="3d")
        ax.set_proj_type("ortho")
        polys = np.concatenate([tri, pt, tt]); pcols = np.concatenate([tri_cols, pc, tc])
        coll = Poly3DCollection(polys, facecolors=pcols, edgecolors="none", linewidths=0, shade=False)
        ax.add_collection3d(coll)
        bc = Poly3DCollection(bores, facecolors=[COL_BORE] * len(bores), edgecolors="none")
        ax.add_collection3d(bc)
        ax.add_collection3d(Line3DCollection(lines, colors=(0.12, 0.10, 0.08, 0.9), linewidths=0.9))
        zmax = 70.0
        ax.set_xlim(cx - rad, cx + rad); ax.set_ylim(cy - rad, cy + rad); ax.set_zlim(0, zmax)
        ax.set_box_aspect((1, 1, zmax / (2 * rad)))
        ax.view_init(elev=elev, azim=azim)
        ax.set_axis_off()
    if title:
        fig.suptitle(title, fontsize=15, y=0.97)
    if legend:
        fig.text(0.5, 0.02,
                 "green = terrain (colour = height level, 15 mm per level) | lighter flat pads = standable plateaus (55% of hex, dead flat)\n"
                 "tan = road (smooth bed, 1 mm below terrain, 16 mm wide) | blue = river channel (2.5 mm deep, 22 mm wide, flat bed) | black dots = 5.3 mm magnet sockets (flat hexes not crossed by a road/river) / wall bores\n"
                 "thin dark lines = engraved V-lines on every hex edge (1.0 x 0.6 mm) | grey band = 10 mm base plate with 18 wall bores at 3.9 mm above the bed",
                 ha="center", fontsize=10.5, family="monospace")
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out, f"({time.time()-t0:.1f}s)")


# ---------------- scenes ----------------
def flower_a(origin=(0.0, 0.0)):
    return dict(origin=origin,
                levels={0: 1, 1: 1, 2: 2, 3: 3, 4: 2, 5: 1, 6: 0},
                standable={0, 2, 3, 5})


def flower_b(origin):
    return dict(origin=origin,
                levels={0: 1, 1: 2, 2: 1, 3: 1, 4: 0, 5: 0, 6: 1},
                standable={0, 1, 4, 5})


def single_flower_scene():
    O = (0.0, 0.0)
    road = [side_mid(O, 3, 14), hex_center(O, 4), hex_center(O, 0), hex_center(O, 1), side_mid(O, 0, 14)]
    river = [side_mid(O, 4, 14), hex_center(O, 5), hex_center(O, 6), side_mid(O, 5, 14)]
    return Scene([flower_a(O)], roads=[road], rivers=[river], seed=3)


def two_flower_scene():
    O = (0.0, 0.0)
    B = LAYOUT.neighbor_flower_offset(0)
    road = [side_mid(O, 3, 14), hex_center(O, 4), hex_center(O, 0), hex_center(O, 1), side_mid(O, 0, 0),
            hex_center(B, 4), hex_center(B, 0), hex_center(B, 2), side_mid(B, 1, 14)]
    river = [side_mid(O, 4, 14), hex_center(O, 5), hex_center(O, 6), side_mid(O, 5, 14)]
    return Scene([flower_a(O), flower_b(B)], roads=[road], rivers=[river], seed=3)


if __name__ == "__main__":
    which = sys.argv[1]; n = int(sys.argv[2]); out = sys.argv[3]
    if which == "single":
        render(single_flower_scene(), out, n=n,
               title="GOAL: one printable hex flower (130 x 135 mm, 45 mm hexes, 4 height levels x 15 mm)")
    elif which == "single_iso":
        render(single_flower_scene(), out, n=n, views=(("iso", 32, -55),), figsize=(13, 11),
               title="GOAL: one printable hex flower (130 x 135 mm, 45 mm hexes, 4 height levels x 15 mm)")
    elif which == "two":
        render(two_flower_scene(), out, n=n, views=(("iso", 40, -50),), figsize=(16, 11),
               title="GOAL: two flowers joined - matching edge contour, V-line closes across the seam, road continues")
