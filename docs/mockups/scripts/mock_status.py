import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
S = "/private/tmp/claude-501/-Users-peteresser-Developer-projects-archive-tableTopWorld/987a149e-0f53-4c5f-952a-baba80f83122/scratchpad/mock/"

def crop(img, pad=10, drop_top=0.0):
    a = np.asarray(img)
    a = a[int(a.shape[0] * drop_top):]
    mask = (a[:, :, :3] < 250).any(axis=2)
    ys, xs = np.where(mask)
    return a[max(0, ys.min() - pad):ys.max() + pad, max(0, xs.min() - pad):xs.max() + pad]

cols = [
    ("legacy", "JUNE 2026 - Peter's own flat tile\n(printableFiles/hexagonWithRoad.stl)",
     "has: 10 mm plate, 45 mm hexes, 130 x 135 mm\n"
     "has: 18 wall bores (5.3 dia x 1.5 deep, z 3.9)\n"
     "has: hairline V-lines 0.45 x 0.3 on every edge,\n"
     "     half-V on the silhouette\n"
     "has: 2 top holes (5.3 dia x 1.5)\n"
     "has: road sunk 0.5 mm, flat bed (width bug:\n"
     "     covers half the tile)\n"
     "lacks: any height, water not cut"),
    ("scurve", "JUNE 2026 - Peter's organic prototype\n(demo_trihex_noise_flower_scurve.stl)",
     "has: two dead-flat plateaus (20.5 / 30.5)\n"
     "has: smooth sinuous S-curve step front\n"
     "has: crisp hex lines everywhere, also\n"
     "     across the slope\n"
     "has: magnets in the wall\n"
     "lacks: per-hex levels, cross-tile matching,\n"
     "       roads / rivers"),
    ("cli_today", "TODAY - generator, shipped defaults\n(python -m terrain.cli render flower hill_peak)",
     "has: watertight mesh, per-side height contract\n"
     "     (bit-identical seams), 2-flower preview\n"
     "has: plateaus (55%) on a seeded subset of hexes\n"
     "BUT: 6 mm wide / 1.5 mm deep ditch on every\n"
     "     hex edge (defaults not the 1.0 x 0.6 line)\n"
     "BUT: plate 2 mm not 10 (scale never applied),\n"
     "     no bores, no top sockets\n"
     "BUT: interior = blurry IDW gradient, no S-curve\n"
     "     step fronts; roads = zero-width crease"),
    ("current", "TODAY - prototype script\n(scripts/gen_groove_prototypes.py, thin skirt)",
     "has: 1.0 x 0.6 V-lines that read as '/' alone\n"
     "     and 'V' when joined (this session)\n"
     "has: plateaus dead flat at the declared level\n"
     "same gaps as the CLI: 2 mm plate, no magnets,\n"
     "     no sockets, no real roads/rivers,\n"
     "     gradient instead of S-curve bands,\n"
     "     coarse triangles (~3 mm)"),
]
fig, axs = plt.subplots(3, 4, figsize=(26, 15), gridspec_kw=dict(height_ratios=[1, 0.62, 0.5]))
for j, (key, title, notes) in enumerate(cols):
    hm = plt.imread(S + f"{key}_heightmap.png"); iso = plt.imread(S + f"{key}_render.png")
    axs[0, j].imshow(crop(hm)); axs[0, j].set_title(title, fontsize=12); axs[0, j].set_axis_off()
    axs[1, j].imshow(crop(iso, drop_top=0.12)); axs[1, j].set_axis_off()
    axs[2, j].set_axis_off()
    axs[2, j].text(0.0, 1.0, notes, fontsize=11, family="monospace", va="top", ha="left",
                   bbox=dict(boxstyle="round,pad=0.5", fc=(1, 1, 0.96), ec=(0.2, 0.2, 0.2), lw=0.6))
fig.suptitle("WHERE WE STAND - measured top-surface heightmaps (top row) and shaded renders (middle row) of the real STL files", fontsize=16, y=0.995)
fig.tight_layout()
fig.savefig(S + "status_0_where_we_stand.png", dpi=100, bbox_inches="tight")
print("wrote status sheet")
