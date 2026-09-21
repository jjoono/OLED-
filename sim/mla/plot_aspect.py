# -*- coding: utf-8 -*-
"""Microlens aspect-ratio tolerance, generic stack.

Left: the lens shape.  Up to AR = 1 it is a spherical cap, which is what a
reflowed lens is and which close-packs without overlap; beyond AR = 1 a cap of
the same base would cut into its neighbours, so a taller lens is a half-
ellipsoid.  The two are the same hemisphere at AR = 1.

Right: EQE against AR for the two reflectors, with the flat (lens-free)
interface of the same stack as the dashed reference.
"""
import numpy as np, matplotlib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
C_AL, C_AG = '#1A1A1A', '#C00000'
FS_LAB, FS_TICK, FS_NOTE = 7.4, 6.8, 6.4
plt.rcParams.update({'font.size': FS_LAB, 'axes.linewidth': 0.8,
                     'xtick.major.width': 0.8, 'ytick.major.width': 0.8})
fig = plt.figure(figsize=(6.6, 2.7))
gs = fig.add_gridspec(1, 2, width_ratios=[0.62, 1.0], wspace=0.34,
                      left=0.02, right=0.975, top=0.93, bottom=0.18)

# ---- geometry ---------------------------------------------------------------
ax = fig.add_subplot(gs[0, 0]); ax.axis('off')
ax.add_patch(Rectangle((-0.1, -0.55), 3.2, 0.55, facecolor='#BBD7EA',
                       edgecolor='0.3', lw=0.8))


def profile(AR, r=0.45, npt=200):
    """Half-profile of the lens: spherical cap up to AR = 1, half-ellipsoid above."""
    h = AR * r
    if AR <= 1.0:                       # spherical cap, centre below the base plane
        R = (r * r + h * h) / (2 * h)
        zc = h - R
        t = np.linspace(-1, 1, npt) * np.arcsin(min(1.0, r / R))
        return R * np.sin(t), zc + R * np.cos(t)
    t = np.linspace(0, np.pi, npt)      # half-ellipsoid, widest at the base circle
    return r * np.cos(t), h * np.sin(t)


for k, AR in enumerate((0.3, 0.7, 1.3)):
    cx = 0.5 + k
    x, z = profile(AR)
    ax.plot(cx + x, z, color='0.15', lw=1.2)
    ax.plot([cx - 0.45, cx + 0.45], [0, 0], color='0.3', lw=0.8)
    ax.text(cx, -0.30, 'AR = %.1f' % AR, ha='center', va='center',
            fontsize=FS_NOTE, color='0.25')
    ax.text(cx, -0.44, 'cap' if AR <= 1 else 'ellipsoid', ha='center', va='center',
            fontsize=FS_NOTE - 0.6, color='0.45', style='italic')
ax.annotate('', xy=(0.5, 0.135), xytext=(0.5, 0.0),
            arrowprops=dict(arrowstyle='<->', lw=0.8, color=C_AG, mutation_scale=6))
ax.text(0.42, 0.07, '$h$', ha='right', va='center', fontsize=FS_NOTE, color=C_AG)
ax.annotate('', xy=(0.95, -0.12), xytext=(0.5, -0.12),
            arrowprops=dict(arrowstyle='<->', lw=0.8, color=C_AG, mutation_scale=6))
ax.text(0.72, -0.19, '$r$', ha='center', va='top', fontsize=FS_NOTE, color=C_AG)
ax.set_xlim(-0.15, 3.15); ax.set_ylim(-0.62, 0.95)
ax.text(1.5, 0.84, 'hexagonally close-packed, index-matched\n'
                   'to the substrate,  AR = $h/r$',
        ha='center', va='center', fontsize=FS_NOTE, color='0.3', linespacing=1.4)

# ---- EQE vs aspect ratio ----------------------------------------------------
bx = fig.add_subplot(gs[0, 1])
A = np.genfromtxt(os.path.join(HERE, 'mla_aspect.csv'), delimiter=',', names=True,
                  dtype=None, encoding='utf-8')
F = np.genfromtxt(os.path.join(HERE, 'mla_flat_reference.csv'), delimiter=',',
                  names=True, dtype=None, encoding='utf-8')
FLAT = {str(r['reflector']): float(r['EQE_flat']) for r in F}
m = A[A['aspect_ratio'] >= 0.09]
for k, c, lab in (('EQE_Ag', C_AG, 'Ag reflector'), ('EQE_Al', C_AL, 'Al reflector')):
    bx.plot(m['aspect_ratio'], 100 * m[k], color=c, lw=1.7, label=lab)
    flat = FLAT[k.split('_')[1]]
    bx.axhline(100 * flat, color=c, lw=0.7, dashes=(1.6, 2.0), zorder=1)
    i = int(np.argmax(m[k]))
    thr = m['aspect_ratio'][m[k] >= 0.95 * m[k][i]][0]
    bx.plot([thr], [100 * m[k][m['aspect_ratio'] == thr][0]], marker='o', ms=4,
            color=c, mew=0, zorder=5)
bx.text(0.03, 100 * FLAT['Ag'] + 1.0, 'flat interface (Ag / Al)', color='0.35',
        fontsize=FS_NOTE - 0.4, ha='left', va='bottom')
bx.axvline(1.0, color='0.75', lw=0.7, dashes=(2.5, 2.5), zorder=0)
bx.text(1.015, 59.0, 'hemisphere', fontsize=FS_NOTE - 0.4, color='0.5',
        ha='left', va='bottom', rotation=90)
bx.set_xlim(0.0, 1.55); bx.set_ylim(20, 95)
bx.set_xlabel('aspect ratio   $h_{\\rm lens}/r_{\\rm lens}$')
bx.set_ylabel('EQE  (%)')
bx.tick_params(labelsize=FS_TICK)
bx.legend(loc='center right', bbox_to_anchor=(1.0, 0.60), frameon=False,
          fontsize=FS_NOTE, handlelength=1.8, borderaxespad=0.6)
bx.text(0.02, 0.965, '$n_{\\rm sub}$ = $n_{\\rm MLA}$ = 1.80', transform=bx.transAxes,
        ha='left', va='top', fontsize=FS_NOTE, color='0.3')
bx.text(0.40, 36.0, 'dots: from here on, within 5 % of the best', fontsize=FS_NOTE - 0.4,
        color='0.45', ha='left', va='center')
for s in ('top', 'right'):
    bx.spines[s].set_visible(False)
for e in ('png', 'pdf'):
    fig.savefig(os.path.join(HERE, 'mla_aspect.' + e), dpi=400)
print('mla_aspect.png written')
