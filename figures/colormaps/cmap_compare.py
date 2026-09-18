"""Side-by-side of the computed power-dissipation spectrum in candidate colormaps."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from cmcrameri import cm as cmc

M = np.load('pds_M.npy'); u = np.load('pds_u.npy'); dE = np.load('pds_d.npy')
j = u <= 1.2
M, u = M[:, j], u[j]
M = np.clip(M, 1e-5, None)
norm = LogNorm(vmin=3e-4, vmax=1.0)

cands = [('jet', plt.get_cmap('jet'), 'current (rainbow) — avoid'),
         ('viridis', plt.get_cmap('viridis'), 'viridis'),
         ('inferno', plt.get_cmap('inferno'), 'inferno  ← recommended'),
         ('magma', plt.get_cmap('magma'), 'magma'),
         ('batlow', cmc.batlow, 'batlow (Crameri)'),
         ('lipari', cmc.lipari, 'lipari (Crameri)  ← recommended')]

fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.4))
for ax, (name, cmap, title) in zip(axes.ravel(), cands):
    im = ax.pcolormesh(u, dE, M, cmap=cmap, norm=norm, shading='nearest', rasterized=True)
    ax.set_title(title, fontsize=11, fontweight='bold' if '←' in title else 'normal',
                 color='#b00' if 'avoid' in title else 'k')
    ax.set_xlim(0, 1.2); ax.set_ylim(10, 500)
    ax.set_xlabel('Normalized in-plane wavevector  u'); ax.set_ylabel('d$_{ETL}$ (nm)')
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label('power dissipation (norm.)', fontsize=8.5)
    ax.annotate('SPP', xy=(1.133, 60), xytext=(0.80, 120), color='w', fontsize=9, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='w', lw=1.2))
fig.suptitle('Same data, six colormaps — log colour scale, Ag / ETL(1.8) / EML / HTL / TCO / substrate, 550 nm',
             fontsize=12.5, y=0.98)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('cmap_compare.png', dpi=140)
print('saved compare')

# --- ramp swatches + luminance + CVD (deuteranopia) ---
def lum(rgb):
    r, g, b = [np.where(c <= 0.04045, c/12.92, ((c+0.055)/1.055)**2.4) for c in rgb.T[:3]]
    return 0.2126*r + 0.7152*g + 0.0722*b
def deuter(rgb):  # Brettel/Vienot-style LMS simulation
    M1 = np.array([[17.8824, 43.5161, 4.11935], [3.45565, 27.1554, 3.86714], [0.0299566, 0.184309, 1.46709]])
    M2 = np.array([[0.080944, -0.13050, 0.116721], [-0.0102485, 0.0540194, -0.113615], [-0.000365294, -0.00412163, 0.693513]])
    D = np.array([[1, 0, 0], [0.49421, 0, 1.24827], [0, 0, 1]])
    lms = rgb[:, :3] @ M1.T
    return np.clip((lms @ D.T) @ M2.T, 0, 1)

fig2, axs = plt.subplots(len(cands), 3, figsize=(13, 1.05*len(cands)+1.4),
                         gridspec_kw={'width_ratios': [3, 3, 2]})
x = np.linspace(0, 1, 256)
for row, (name, cmap, title) in enumerate(cands):
    rgb = cmap(x)
    axs[row, 0].imshow([rgb], aspect='auto', extent=[0, 1, 0, 1]); axs[row, 0].set_yticks([])
    axs[row, 1].imshow([deuter(rgb)], aspect='auto', extent=[0, 1, 0, 1]); axs[row, 1].set_yticks([])
    L = lum(rgb)
    mono = np.all(np.diff(L) > -0.004)
    axs[row, 2].plot(x, L, color='k', lw=1.4); axs[row, 2].set_ylim(0, 1); axs[row, 2].set_yticks([0, 1])
    axs[row, 2].text(0.03, 0.88, 'monotonic ✓' if mono else 'NOT monotonic ✗', ha='left', va='top', fontsize=8.5,
                     color='#0a0' if mono else '#b00', fontweight='bold', transform=axs[row, 2].transAxes,
                     bbox=dict(boxstyle='round,pad=0.2', fc='w', ec='none', alpha=0.85))
    axs[row, 0].set_ylabel(name, rotation=0, ha='right', va='center', fontsize=10.5, labelpad=8)
    for c in range(3):
        axs[row, c].set_xticks([] if row < len(cands)-1 else [0, 0.5, 1])
for c, t in enumerate(['normal vision', 'deuteranopia (simulated)', 'perceived lightness']):
    axs[0, c].set_title(t, fontsize=10)
fig2.tight_layout()
fig2.savefig('cmap_ramps.png', dpi=140)
print('saved ramps')
