import sys, numpy as np, trimesh, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
path, out, step = sys.argv[1], sys.argv[2], float(sys.argv[3])
m = trimesh.load(path)
b = m.bounds
xs = np.arange(b[0][0], b[1][0]+step, step); ys = np.arange(b[0][1], b[1][1]+step, step)
H = np.full((len(ys), len(xs)), np.nan)
tris = m.vertices[m.faces]
for t in tris:
    (x0,y0,z0),(x1,y1,z1),(x2,y2,z2) = t
    det = (y1-y2)*(x0-x2) + (x2-x1)*(y0-y2)
    if abs(det) < 1e-12: continue
    i0 = max(0, int((min(x0,x1,x2)-xs[0])/step)); i1 = min(len(xs)-1, int((max(x0,x1,x2)-xs[0])/step)+1)
    j0 = max(0, int((min(y0,y1,y2)-ys[0])/step)); j1 = min(len(ys)-1, int((max(y0,y1,y2)-ys[0])/step)+1)
    if i1 < i0 or j1 < j0: continue
    gx, gy = np.meshgrid(xs[i0:i1+1], ys[j0:j1+1])
    l0 = ((y1-y2)*(gx-x2) + (x2-x1)*(gy-y2))/det
    l1 = ((y2-y0)*(gx-x2) + (x0-x2)*(gy-y2))/det
    l2 = 1-l0-l1
    inside = (l0>=-1e-9)&(l1>=-1e-9)&(l2>=-1e-9)
    if not inside.any(): continue
    z = l0*z0 + l1*z1 + l2*z2
    sub = H[j0:j1+1, i0:i1+1]
    upd = inside & (np.isnan(sub) | (z > sub))
    sub[upd] = z[upd]
np.save(out.replace('.png','.npy'), H)
fig, ax = plt.subplots(figsize=(9,9))
im = ax.imshow(H, origin='lower', extent=[xs[0],xs[-1],ys[0],ys[-1]], cmap='viridis', interpolation='nearest')
ax.set_title(f"{path.split('/')[-1]}  bounds x{np.round(b[:,0],1)} y{np.round(b[:,1],1)} z{np.round(b[:,2],1)}", fontsize=9)
plt.colorbar(im, ax=ax, label='top-surface z (mm)')
ax.set_xlabel('x mm'); ax.set_ylabel('y mm'); ax.grid(True, alpha=.3)
fig.savefig(out, dpi=110, bbox_inches='tight')
vals = H[~np.isnan(H)]
print(path, "z min/max", vals.min().round(2), vals.max().round(2))
u, c = np.unique(np.round(vals,1), return_counts=True)
top = sorted(zip(c,u), reverse=True)[:12]
print("most common z (count, z):", [(int(a), float(b)) for a,b in top])
