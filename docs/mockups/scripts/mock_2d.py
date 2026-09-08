import sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle, Rectangle
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch
sys.path.insert(0, "/private/tmp/claude-501/-Users-peteresser-Developer-projects-archive-tableTopWorld/987a149e-0f53-4c5f-952a-baba80f83122/scratchpad/mock")
from mock_common import *
from mock_3d import single_flower_scene, flower_a, TERRAIN_CMAP, COL_ROAD, COL_RIVER

LEVEL_COL = {l: TERRAIN_CMAP(l / (LEVELS - 1))[:3] for l in range(LEVELS)}
INK = (0.12, 0.10, 0.08)


def hex_poly(center, scale=1.0):
    return [(center[0] + v[0] * scale, center[1] + v[1] * scale) for v in HEX_VERTS]


def silhouette():
    pts = []
    for k in range(6):
        pts.extend(LAYOUT.side_corners(k)[:3])
    return pts


def band_polygon(poly, half_w):
    P = np.asarray(poly)
    d = np.gradient(P, axis=0); d /= np.maximum(np.linalg.norm(d, axis=1)[:, None], 1e-9)
    nrm = np.c_[-d[:, 1], d[:, 0]]
    return np.vstack([P + nrm * half_w, (P - nrm * half_w)[::-1]])


def dim_h(ax, x0, x1, y, text, above=True, fs=9):
    ax.annotate("", xy=(x0, y), xytext=(x1, y), arrowprops=dict(arrowstyle="<->", color=INK, lw=0.9))
    ax.text((x0 + x1) / 2, y + (1.0 if above else -1.0), text, ha="center", va="bottom" if above else "top", fontsize=fs, color=INK)


def dim_v(ax, x, y0, y1, text, fs=9, side=1, dx=1.0):
    ax.annotate("", xy=(x, y0), xytext=(x, y1), arrowprops=dict(arrowstyle="<->", color=INK, lw=0.9))
    ax.text(x + dx * side, (y0 + y1) / 2, text, ha="left" if side > 0 else "right", va="center", fontsize=fs, color=INK)


# ------------------------------------------------------------------ top view
def top_view(out):
    scene = single_flower_scene()
    f = flower_a()
    fig, ax = plt.subplots(figsize=(15, 14))
    ax.set_aspect("equal")
    clip = PathPatch(MPath(silhouette() + [silhouette()[0]]), transform=ax.transData, visible=False)
    ax.add_patch(clip)
    for h in range(7):
        c = hex_center((0, 0), h); lv = f["levels"][h]
        ax.add_patch(Polygon(hex_poly(c), closed=True, facecolor=LEVEL_COL[lv], edgecolor="none", alpha=0.75))
        if h in f["standable"]:
            ax.add_patch(Polygon(hex_poly(c, PLATEAU_R), closed=True, facecolor=(1, 1, 1, 0.35), edgecolor=INK, ls="--", lw=0.9))
    for path, col, hw in ((scene.roads[0], COL_ROAD, ROAD_HALF_W), (scene.rivers[0], COL_RIVER, RIVER_HALF_W)):
        p = Polygon(band_polygon(path.poly, hw), closed=True, facecolor=col, edgecolor=INK, lw=0.6, alpha=0.9)
        ax.add_patch(p); p.set_clip_path(clip)
    for h in range(7):
        c = hex_center((0, 0), h); lv = f["levels"][h]
        if h in scene.socket_hexes:
            ax.add_patch(Circle(c, SOCKET_R, facecolor=(0.1, 0.1, 0.1), edgecolor="none", zorder=4))
        kind = "FLAT pad" if h in f["standable"] else "organic ~"
        extra = "" if h in scene.socket_hexes or h not in f["standable"] else "\n(no socket: road/river)"
        ax.text(c[0], c[1] - (5 if h in scene.socket_hexes else 2), f"hex {h} - level {lv}\n{kind}{extra}", ha="center", va="top", fontsize=8.5, color=INK, zorder=5,
                bbox=dict(boxstyle="round,pad=0.15", fc=(1, 1, 1, 0.55), ec="none"))
    for h in range(7):
        p = hex_poly(hex_center((0, 0), h)); p.append(p[0])
        ax.plot(*zip(*p), color=INK, lw=1.0)
    for (x1, y1), (x2, y2) in exterior_edges_world((0, 0)):
        ax.plot([x1, x2], [y1, y2], color=INK, lw=2.2)
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        ang = math.atan2(y2 - y1, x2 - x1)
        ax.add_patch(Rectangle((0, 0), 2 * BORE_R, BORE_DEPTH, facecolor=(0.1, 0.1, 0.1),
                               transform=matplotlib.transforms.Affine2D().translate(-BORE_R, -BORE_DEPTH / 2).rotate(ang).translate(mx, my) + ax.transData, zorder=6))
    for k in range(6):
        for p in LAYOUT.side_corners(k):
            ax.plot(p[0], p[1], "o", ms=5, color="white", mec=INK, mew=1.0, zorder=7)
        m = side_mid((0, 0), k, 0)
        ax.plot(m[0], m[1], marker="x", ms=9, color="crimson", mew=1.8, zorder=8)
        lab = side_mid((0, 0), k, 16)
        ax.text(lab[0], lab[1], f"side {k}", ha="center", va="center", fontsize=9, color=INK,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=INK, lw=0.5), zorder=8)
    # dimensions - all outside the flower
    c2 = hex_center((0, 0), 2)
    dim_h(ax, c2[0] - A, c2[0] + A, 74, "45.0 mm across flats (hex 2)", above=True)
    ax.plot([c2[0] - A, c2[0] - A], [c2[1] + R * 0.5, 74], color=INK, lw=0.5, ls=":"); ax.plot([c2[0] + A, c2[0] + A], [c2[1] + R * 0.5, 74], color=INK, lw=0.5, ls=":")
    ax.annotate("", xy=(0, 0), xytext=(HEX_VERTS[0][0], HEX_VERTS[0][1]), arrowprops=dict(arrowstyle="<->", color=INK, lw=0.9))
    ax.text(13, 2.0, "R 26.0", fontsize=8.5, color=INK)
    c5 = hex_center((0, 0), 5)
    dim_h(ax, c5[0] - A * PLATEAU_R, c5[0] + A * PLATEAU_R, -78, "flat pad: 33.4 mm across flats = 55% of the hex area", above=False)
    ax.plot([c5[0] - A * PLATEAU_R] * 2, [c5[1] - R * 0.55, -78], color=INK, lw=0.5, ls=":"); ax.plot([c5[0] + A * PLATEAU_R] * 2, [c5[1] - R * 0.55, -78], color=INK, lw=0.5, ls=":")
    dim_h(ax, -64.95, 64.95, -90, "flower footprint: 129.9 mm (x) by 135.0 mm (y)", above=False)
    dim_v(ax, 92, -67.5, 67.5, "135.0", side=1)
    notes = (
        "GOAL - top view of one flower, true layout coordinates (mm)\n"
        "\n"
        "thin black lines   engraved V-line on EVERY hex edge (1.0 wide x 0.6 deep); a lone edge is '/', two hexes make the 'V'\n"
        "thick outline      flower silhouette: 18 edges = 6 sides x 3 edges; the outer edge carries half a V too, so the seam reads as a line\n"
        "black bars         18 magnet bores, one per silhouette edge, at the edge midpoint (5.3 dia x 2.2 deep, centre 3.9 mm above the bed)\n"
        "black discs        top magnet sockets (5.3 dia x 2.2 deep) at the centre of every FLAT hex not crossed by a road/river\n"
        "dashed hexagon     the flat, noise-free pad of a standable hex (0.742 x hex size = 55% of the area); the 6 mm band outside it\n"
        "                   is the S-curve transition to the neighbour's height; organic hexes have no pad, just rolling terrain + noise\n"
        "red x              road / river crossing point = midpoint of the side's MIDDLE edge (shared by exactly two flowers)\n"
        "tan band           road 16 mm wide, smooth bed 1 mm below the terrain, runs hex-centre to hex-centre\n"
        "blue band          river 22 mm wide, 2.5 mm deep channel, flat bed with 5 mm sloped banks\n"
        "white dots         the 4 corners of each side: declared height levels live here (matching contract with the neighbour)\n"
    )
    ax.text(-105, -105, notes, fontsize=9.2, family="monospace", va="top", ha="left", color=INK,
            bbox=dict(boxstyle="round,pad=0.6", fc=(1, 1, 0.96), ec=INK, lw=0.6))
    ax.set_xlim(-108, 108); ax.set_ylim(-200, 92)
    ax.set_axis_off()
    fig.savefig(out, dpi=120, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)


# ------------------------------------------------------------------ cross-section
def cross_section(out):
    scene = single_flower_scene()
    d = np.array([math.cos(math.radians(-30)), math.sin(math.radians(-30))])   # hex 3 -> 0 -> 6 axis
    t = np.arange(-67.5, 67.5 + 1e-9, 0.1)
    P = np.c_[t * d[0], t * d[1]]
    z, mat = scene.surface(P)
    # carve the V-lines the section crosses: internal edges at +-22.5, silhouette at +-67.5 (half V)
    for te in (-22.5, 22.5):
        z -= LINE_D * np.clip(1 - np.abs(t - te) / (LINE_W / 2), 0, 1)
    z -= LINE_D * np.clip(1 - (t + 67.5) / (LINE_W / 2), 0, 1)
    z -= LINE_D * np.clip(1 - (67.5 - t) / (LINE_W / 2), 0, 1)
    fig, ax = plt.subplots(figsize=(21, 8.5))
    ax.set_aspect("equal")
    ax.add_patch(Rectangle((-67.5, 0), 135, PLATE, facecolor=(0.84, 0.84, 0.82), edgecolor=INK, lw=1.0))
    ax.fill_between(t, PLATE, z, where=z >= PLATE, color=(0.58, 0.48, 0.38), alpha=0.95, lw=0)
    ax.fill_between(t, z, PLATE, where=z < PLATE, color="white", lw=0)
    ax.plot(t, z, color=INK, lw=1.3)
    for m, col in ((2, COL_ROAD), (3, COL_RIVER)):
        sel = mat == m
        if sel.any():
            ax.fill_between(t[sel], z[sel], z[sel] + 0.5, color=col, lw=0, zorder=4)
    for i in scene.socket_hexes:
        pass
    for x in (-67.5, 67.5):
        ax.plot([x, x], [0, np.interp(x, t, z)], color=INK, lw=1.3)
    for x, sgn in ((-67.5, 1), (67.5, -1)):
        ax.add_patch(Rectangle((x if sgn > 0 else x - BORE_DEPTH, BORE_Z - BORE_R), BORE_DEPTH, 2 * BORE_R,
                               facecolor="white", edgecolor=INK, lw=1.0, hatch="////", zorder=3))
    for l in range(LEVELS):
        ax.axhline(level_z(l), color=INK, lw=0.6, ls=":", alpha=0.7)
        ax.text(69, level_z(l), f"L{l} = {level_z(l):.0f} mm", va="center", fontsize=9, color=INK)
    ax.axhline(0, color=INK, lw=0.6, ls=":", alpha=0.7); ax.text(69, 0, "print bed z = 0", va="center", fontsize=9, color=INK)
    dim_v(ax, -71, 0, PLATE, "10 mm base plate", side=-1)
    ax.text(-72, 2.5, "magnets live here, same\nheight on every tile", ha="right", va="center", fontsize=8.5, color=INK)
    dim_v(ax, -71, PLATE, level_z(3), "3 levels x 15 mm", side=-1)
    ax.annotate("wall bore 5.3 dia x 2.2 deep,\ncentre 3.9 above the bed", xy=(-66.4, BORE_Z), xytext=(-64, -6.5), fontsize=9, color=INK, arrowprops=dict(arrowstyle="->", color=INK))
    ax.annotate("flat pad: dead flat at L3 (55 mm), no noise", xy=(-45, 55), xytext=(-40, 64), fontsize=9, color=INK, arrowprops=dict(arrowstyle="->", color=INK))
    ax.annotate("6 mm S-curve band down to the shared\nedge; the edge sits at the mean of both\nhexes (L3 + L1)/2 = 40 mm, then hex 0's\nown band continues down to its pad", xy=(-26, 46), xytext=(-8, 56), fontsize=9, color=INK, arrowprops=dict(arrowstyle="->", color=INK))
    ax.annotate("V-line on the shared edge\n(1.0 x 0.6 mm, to scale)", xy=(-22.5, np.interp(-22.5, t, z)), xytext=(-40, 30), fontsize=9, color=INK, arrowprops=dict(arrowstyle="->", color=INK))
    ax.annotate("road crossing hex 0 (tan): smooth bed\n1 mm below the pad, vertical edges,\n~18.5 mm here because it crosses at 60 deg", xy=(0, np.interp(0, t, z)), xytext=(-2, 36), fontsize=9, color=INK, arrowprops=dict(arrowstyle="->", color=INK))
    ax.annotate("organic hex 6 (L0): no pad, rolling\nterrain with +/-0.5 mm noise", xy=(28, np.interp(28, t, z)), xytext=(22, 27), fontsize=9, color=INK, arrowprops=dict(arrowstyle="->", color=INK))
    riv = np.where(mat == 3)[0]
    if len(riv):
        xr = t[riv].mean()
        ax.annotate("river channel 22 wide x 2.5 deep, flat bed;\nat L0 the bed (7.5) stays >= 1 mm above\nthe bore top (6.55)", xy=(xr, z[riv].min()), xytext=(30, 17), fontsize=9, color=INK, arrowprops=dict(arrowstyle="->", color=INK))
    ax.annotate("silhouette edge: half a V-line ('/'),\nthe neighbour closes it", xy=(67.3, np.interp(67.3, t, z)), xytext=(40, -6.5), fontsize=9, color=INK, arrowprops=dict(arrowstyle="->", color=INK))
    ax.annotate("V-line", xy=(22.5, np.interp(22.5, t, z)), xytext=(14, 11), fontsize=9, color=INK, arrowprops=dict(arrowstyle="->", color=INK))
    ax.set_xlim(-100, 90); ax.set_ylim(-9, 70)
    ax.set_title("GOAL - cross-section through hex 3 (L3, flat) -> hex 0 (L1, flat, road) -> hex 6 (L0, organic, river); true scale, mm", fontsize=13)
    ax.set_xlabel("mm along the section"); ax.set_yticks([]); ax.spines[["top", "right", "left"]].set_visible(False)
    fig.savefig(out, dpi=120, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)


# ------------------------------------------------------------------ details sheet
def details(out):
    fig, axs = plt.subplots(2, 3, figsize=(21, 12))
    for ax in axs.flat:
        ax.set_aspect("equal"); ax.set_axis_off()
    G1, G2 = (0.58, 0.66, 0.34), (0.50, 0.60, 0.30)

    # 1. V-line: single hex vs joined (drawn at 1:1 in a 16 x 5 mm window)
    ax = axs[0, 0]
    ax.set_title("1. Hex line = half a V per hex (skirt), 1:1 in a 16 x 5 mm window", fontsize=11)
    # single hex: surface, 0.5 mm skirt down to -0.6 at the edge, then the flower wall
    xs = [-7.5, -1.5, -1.0, -1.0]; zs = [0, 0, -LINE_D, -3]
    ax.fill_between([-7.5, -1.5, -1.0], -3, [0, 0, -LINE_D], color=G1); ax.plot(xs, zs, color=INK, lw=1.5)
    ax.text(-4.5, 0.6, "one hex alone: '/'", ha="center", fontsize=10, color=INK)
    ax.annotate("flower wall at the edge", xy=(-1.0, -2.0), xytext=(-6.5, -2.6), fontsize=8, color=INK, arrowprops=dict(arrowstyle="->", color=INK))
    # joined: two half-skirts meeting at x = 4
    xa = [1.0, 3.5, 4.0]; za = [0, 0, -LINE_D]
    xb = [4.0, 4.5, 7.0]; zb = [-LINE_D, 0, 0]
    ax.fill_between(xa, -3, za, color=G1); ax.fill_between(xb, -3, zb, color=G2)
    ax.plot(xa + xb, za + zb, color=INK, lw=1.5); ax.plot([4, 4], [-LINE_D, -3], color=INK, lw=0.8, ls="--")
    ax.text(4, 0.6, "two hexes (or two flowers) joined: 'V'", ha="center", fontsize=10, color=INK)
    dim_h(ax, 3.5, 4.5, -1.4, "1.0 wide", above=False, fs=9); dim_v(ax, 5.6, -LINE_D, 0, "0.6 deep", side=1, dx=0.3)
    ax.text(-7.5, -4.2, "legacy June design: 0.45 x 0.3  |  this session's approved prototype: 1.0 x 0.6 (used here)", fontsize=8.5, color=INK)
    ax.set_xlim(-8, 8); ax.set_ylim(-4.8, 2.8)

    # 2. seam between two flowers with facing bores + magnets
    ax = axs[0, 1]
    ax.set_title("2. Two flowers meeting (section through a wall bore)", fontsize=11)
    for x0, col in ((-30, (0.84, 0.84, 0.82)), (0, (0.80, 0.82, 0.86))):
        ax.add_patch(Rectangle((x0, 0), 30, PLATE, facecolor=col, edgecolor=INK, lw=1.0))
        ax.add_patch(Rectangle((x0, PLATE), 30, 8, facecolor=(0.58, 0.48, 0.38), edgecolor=INK, lw=1.0))
    for x0 in (-BORE_DEPTH, 0):
        ax.add_patch(Rectangle((x0, BORE_Z - BORE_R), BORE_DEPTH, 2 * BORE_R, facecolor="white", edgecolor=INK, lw=1.0, hatch="////", zorder=3))
    for x0 in (-2.05, 0.05):
        ax.add_patch(Rectangle((x0, BORE_Z - 2.5), 2.0, 5.0, facecolor=(0.3, 0.3, 0.35), edgecolor=INK, lw=0.8, zorder=4))
    ax.annotate("two 5 x 2 mm disc magnets, N/S facing\n(bore 5.3 x 2.2 leaves 0.15 clearance)", xy=(0.3, BORE_Z - 2.6), xytext=(2, -7.0), fontsize=9, color=INK, arrowprops=dict(arrowstyle="->", color=INK))
    xa = [-8, -0.5, 0]; za = [18, 18, 18 - LINE_D]; xb = [0, 0.5, 8]; zb = [18 - LINE_D, 18, 18]
    ax.plot(xa + xb, za + zb, color=INK, lw=1.5)
    ax.annotate("V-line closes across the seam", xy=(0, 17.4), xytext=(-28, 22.5), fontsize=9, color=INK, arrowprops=dict(arrowstyle="->", color=INK))
    dim_v(ax, -33, 0, PLATE, "10", side=1, dx=0.4); dim_v(ax, -38, BORE_Z - BORE_R, BORE_Z + BORE_R, "5.3", side=-1)
    ax.text(-20, -2.5, "flower A", ha="center", fontsize=10, color=INK); ax.text(28, -2.5, "flower B", ha="center", fontsize=10, color=INK)
    ax.text(-40, 29.5, "plate walls straight and flush. The terrain silhouette above may keep its +/-1 mm jitter,\nbut it must fade to 0 within +/-4 mm of every bore so each bore sits in a flat wall.", fontsize=8.5, color=INK, va="top")
    ax.set_xlim(-42, 34); ax.set_ylim(-9.5, 31)

    # 3. wall elevation
    ax = axs[0, 2]
    ax.set_title("3. One silhouette edge seen from outside (18 per flower)", fontsize=11)
    ax.add_patch(Rectangle((-R / 2, 0), R, PLATE, facecolor=(0.84, 0.84, 0.82), edgecolor=INK, lw=1.0))
    xx = np.linspace(-R / 2, R / 2, 80); top = PLATE + 6 + 1.2 * np.sin(xx / 3.1) + 0.6 * np.sin(xx / 1.3 + 1)
    ax.fill_between(xx, PLATE, top, color=(0.58, 0.48, 0.38)); ax.plot(xx, top, color=INK, lw=1.2)
    ax.add_patch(Circle((0, BORE_Z), BORE_R, facecolor=(0.1, 0.1, 0.1), edgecolor=INK, zorder=3))
    dim_h(ax, -R / 2, R / 2, -1.5, "edge 26.0 mm", above=False)
    dim_v(ax, R / 2 + 2, 0, BORE_Z, "3.9", side=1)
    dim_h(ax, -BORE_R, BORE_R, PLATE + 0.8, "5.3 dia", above=True)
    ax.text(-R / 2, PLATE + 12.5, "bore at the edge MIDPOINT so the neighbour's bore on the same\nphysical edge lines up 1:1 (3 per side, 18 per flower);\nthe plate is solid, so a 2.2 mm blind bore needs no thin wall", fontsize=8.5, color=INK, va="bottom")
    ax.set_xlim(-18, 20); ax.set_ylim(-5, 29)

    # 4. socket on a flat pad
    ax = axs[1, 0]
    ax.set_title("4. Top socket on flat hexes (trees, ruins, markers snap on)", fontsize=11)
    ax.add_patch(Rectangle((-20, 0), 40, 6, facecolor=G1, edgecolor=INK, lw=1.2))
    ax.add_patch(Rectangle((-SOCKET_R, 6 - SOCKET_DEPTH), 2 * SOCKET_R, SOCKET_DEPTH, facecolor="white", edgecolor=INK, lw=1.0, hatch="////", zorder=3))
    ax.add_patch(Rectangle((-6, 10.5), 12, 3, facecolor=(0.75, 0.65, 0.5), edgecolor=INK, lw=1.0))
    ax.add_patch(Rectangle((-2.5, 8.5), 5.0, 2.0, facecolor=(0.3, 0.3, 0.35), edgecolor=INK, lw=0.8))
    ax.annotate("", xy=(0, 6.4), xytext=(0, 8.2), arrowprops=dict(arrowstyle="->", color=INK, lw=1.2))
    ax.text(0, 14.2, "scatter piece with a glued 5 x 2 magnet", ha="center", fontsize=9, color=INK)
    dim_h(ax, -SOCKET_R, SOCKET_R, -1.2, "5.3 dia", above=False); dim_v(ax, 4.5, 6 - SOCKET_DEPTH, 6, "2.2", side=1)
    ax.text(-20, -4.6, "socket at the pad centre, only on hexes that are flat AND not crossed by a road/river\n(1-6 flat hexes per flower, seeded)", fontsize=8.5, color=INK, va="top")
    ax.set_xlim(-22, 22); ax.set_ylim(-8.5, 16)

    # 5. road profile
    ax = axs[1, 1]
    ax.set_title("5. Road: smooth bed with a distinct 1 mm step", fontsize=11)
    xx = np.linspace(-22, 22, 440)
    terr = 8 + 0.5 * np.sin(xx * 1.7) + 0.3 * np.sin(xx * 4.1 + 2)
    road = np.where(np.abs(xx) <= ROAD_HALF_W, 8 - ROAD_DEPTH, terr)
    ax.fill_between(xx, 0, road, color=G1); ax.plot(xx, road, color=INK, lw=1.4)
    sel = np.abs(xx) <= ROAD_HALF_W
    ax.fill_between(xx[sel], 0, road[sel], color=COL_ROAD)
    dim_h(ax, -ROAD_HALF_W, ROAD_HALF_W, 9.8, "16 mm wide", above=True); dim_v(ax, ROAD_HALF_W + 1.5, 8 - ROAD_DEPTH, 8, "1.0", side=1)
    ax.text(-22, -1.8, "bed = terrain smoothed along the path (no noise), 1 mm below it, vertical edges;\nhex V-lines continue across the road; ramps follow slopes; enters/exits at the side's middle edge", fontsize=8.5, color=INK, va="top")
    ax.set_xlim(-23, 23); ax.set_ylim(-5, 12)

    # 6. river profile
    ax = axs[1, 2]
    ax.set_title("6. River: carved channel, flat bed, 5 mm banks", fontsize=11)
    bank = np.clip((RIVER_HALF_W - np.abs(xx)) / RIVER_BANK, 0, 1); bank = bank * bank * (3 - 2 * bank)
    chan = 8 - RIVER_DEPTH * bank
    riv = np.where(np.abs(xx) <= RIVER_HALF_W, np.minimum(terr, chan), terr)
    ax.fill_between(xx, 0, riv, color=G1); ax.plot(xx, riv, color=INK, lw=1.4)
    sel = np.abs(xx) <= RIVER_HALF_W
    ax.fill_between(xx[sel], riv[sel], riv[sel] + 0.35, color=COL_RIVER)
    dim_h(ax, -RIVER_HALF_W, RIVER_HALF_W, 9.8, "22 mm wide", above=True); dim_v(ax, RIVER_HALF_W + 1.5, 8 - RIVER_DEPTH, 8, "2.5", side=1)
    ax.text(-22, -1.8, "sloped banks, smooth bed; same crossing points as roads; a river across a flat pad\nstill leaves the pad flat outside the channel; bridges (road x river) stay out of scope", fontsize=8.5, color=INK, va="top")
    ax.set_xlim(-23, 23); ax.set_ylim(-5, 12)

    fig.suptitle("GOAL - detail sheet (all mm; panels 1, 2 and 4 are enlarged)", fontsize=14, y=0.99)
    fig.savefig(out, dpi=120, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    S = "/private/tmp/claude-501/-Users-peteresser-Developer-projects-archive-tableTopWorld/987a149e-0f53-4c5f-952a-baba80f83122/scratchpad/mock/"
    top_view(S + "goal_3_top_view_dimensions.png")
    cross_section(S + "goal_4_cross_section.png")
    details(S + "goal_5_details.png")
