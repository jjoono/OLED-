# -*- coding: utf-8 -*-
"""Microlens aspect-ratio tolerance, generic stack."""
import numpy as np, matplotlib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Wedge

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
for k, AR in enumerate((0.3, 0.7, 1.2)):
    cx = 0.5 + k
    r, h = 0.45, 0.45 * AR
    R = (r * r + h * h) / (2 * h)
    zc = h - R
    th0 = np.degrees(np.arcsin(min(1.0, r / R)))
    a = np.linspace(-th0, th0, 120) + 90
    ax.plot(cx + R * np.cos(np.radians(a)), zc + R * np.sin(np.radians(a)),
            color='0.15', lw=1.2)
    ax.plot([cx - r, cx + r], [0, 0], color='0.3', lw=0.8)
    ax.text(cx, -0.30, 'AR = %.1f' % AR, ha='center', va='center',
            fontsize=FS_NOTE, color='0.25')
ax.annotate('', xy=(0.5, 0.135), xytext=(0.5, 0.0),
            arrowprops=dict(arrowstyle='<->', lw=0.8, color='#C00000',
                            mutation_scale=6))
ax.text(0.42, 0.07, '$h$', ha='right', va='center', fontsize=FS_NOTE, color='#C00000')
ax.annotate('', xy=(0.95, -0.12), xytext=(0.5, -0.12),
            arrowprops=dict(arrowstyle='<->', lw=0.8, color='#C00000',
                            mutation_scale=6))
ax.text(0.72, -0.19, '$r$', ha='center', va='top', fontsize=FS_NOTE, color='#C00000')
ax.set_xlim(-0.15, 3.15); ax.set_ylim(-0.62, 0.95)
ax.text(1.5, 0.82, 'hexagonally close-packed spherical caps,\n'
                   'index-matched to the substrate,  AR = $h/r$',
        ha='center', va='center', fontsize=FS_NOTE, color='0.3', linespacing=1.4)

# ---- EQE vs aspect ratio ----------------------------------------------------
bx = fig.add_subplot(gs[0, 1])
A = np.genfromtxt(os.path.join(HERE, 'mla_aspect.csv'), delimiter=',', names=True)
m = A[A['aspect_ratio'] >= 0.09]
flat = A[A['aspect_ratio'] < 0.09]
for k, c, lab in (('EQE_Ag', C_AG, 'Ag reflector'), ('EQE_Al', C_AL, 'Al reflector')):
    bx.plot(m['aspect_ratio'], 100 * m[k], color=c, lw=1.7, label=lab)
    bx.axhline(100 * flat[k][0], color=c, lw=0.7, dashes=(1.6, 2.0), zorder=1)
    bx.text(0.03, 100 * flat[k][0] + 1.2, 'flat interface', color=c,
            fontsize=FS_NOTE - 0.4, ha='left', va='bottom')
    i = int(np.argmax(m[k]))
    thr = m['aspect_ratio'][m[k] >= 0.95 * m[k][i]][0]
    bx.plot([thr], [100 * m[k][m['aspect_ratio'] == thr][0]], marker='o', ms=4,
            color=c, mew=0, zorder=5)
bx.set_xlim(0.0, 1.55); bx.set_ylim(20, 95)
bx.set_xlabel('aspect ratio   $h_{\\rm lens}/r_{\\rm lens}$')
bx.set_ylabel('EQE  (%)')
bx.tick_params(labelsize=FS_TICK)
bx.legend(loc='center right', frameon=False, fontsize=FS_NOTE, handlelength=1.8,
          borderaxespad=0.6)
bx.text(0.02, 0.965, '$n_{\\rm sub}$ = $n_{\\rm MLA}$ = 1.80', transform=bx.transAxes,
        ha='left', va='top', fontsize=FS_NOTE, color='0.3')
bx.text(0.80, 21.5, 'dots: within 5 % of the best', fontsize=FS_NOTE - 0.4, color='0.45')
for s in ('top', 'right'):
    bx.spines[s].set_visible(False)
for e in ('png', 'pdf'):
    fig.savefig(os.path.join(HERE, 'mla_aspect.' + e), dpi=400)
print('mla_aspect.png written')
