# -*- coding: utf-8 -*-
"""Generic stack + EQE against the horizontal dipole ratio."""
import numpy as np, matplotlib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
C_PRO, C_REF = '#C00000', '#1F3B73'
FS_LAB, FS_TICK, FS_NOTE, FS_LET = 7.4, 6.8, 6.4, 9.0
plt.rcParams.update({'font.size': FS_LAB, 'axes.linewidth': 0.8,
                     'xtick.major.width': 0.8, 'ytick.major.width': 0.8})
fig = plt.figure(figsize=(6.9, 2.9))
gs = fig.add_gridspec(1, 2, width_ratios=[0.72, 1.0], wspace=0.42,
                      left=0.015, right=0.975, top=0.94, bottom=0.17)

# ---- generic stack ---------------------------------------------------------
ax = fig.add_subplot(gs[0, 0]); ax.axis('off')
LAYERS = [('reflector  (100 nm)', '#C9C9C9'),
          ('electron transport  $n$ = 1.8', '#D6E6C8'),
          ('emissive layer  20 nm', '#F6C9A8'),
          ('hole transport  $n$ = 1.8', '#FBE3A0'),
          ('transparent electrode  50 nm', '#BBD7EA'),
          ('substrate + microlens array', '#C5CCE6')]
H = 0.86 / len(LAYERS)
for i, (t, c) in enumerate(LAYERS):
    y = 0.86 - (i + 1) * H
    ax.add_patch(Rectangle((0.04, y), 0.92, H, facecolor=c, edgecolor='0.25', lw=0.7))
    ax.text(0.5, y + H / 2, t, ha='center', va='center', fontsize=FS_NOTE)
ax.set_xlim(0, 1); ax.set_ylim(-0.30, 1.10)
ax.annotate('', xy=(0.5, -0.18), xytext=(0.5, -0.01),
            arrowprops=dict(arrowstyle='-|>', lw=1.0, color='0.3', mutation_scale=8))
ax.text(0.5, -0.22, 'light out', ha='center', va='top', fontsize=FS_NOTE, color='0.3')
ax.text(0.5, 1.02, 'reference:  Al,  $n_{\\rm sub}$ = 1.50,  cavity-optimised\n'
                   'proposed:  Ag,  $n_{\\rm sub}$ = 1.80,  un-tuned',
        ha='center', va='top', fontsize=FS_NOTE, linespacing=1.5)

# ---- EQE vs theta ----------------------------------------------------------
bx = fig.add_subplot(gs[0, 1])
A = np.genfromtxt(os.path.join(HERE, 'fig3_theta.csv'), delimiter=',',
                  names=True, dtype=None, encoding='utf-8')
for name, c, lab in (('proposed', C_PRO, r'proposed  ($n_{\rm sub}$ = 1.80, Ag)'),
                     ('reference', C_REF, r'reference  ($n_{\rm sub}$ = 1.50, Al)')):
    m = A[A['device'] == name]
    o = np.argsort(m['theta_horizontal'])
    bx.plot(m['theta_horizontal'][o], 100 * m['EQE'][o], color=c, lw=1.8, label=lab)
for t, lab in ((2 / 3, 'isotropic\n$\\Theta$ = 0.67'), (0.83, 'oriented\n$\\Theta$ = 0.83')):
    bx.axvline(t, color='0.55', lw=0.7, dashes=(1.6, 1.8), zorder=1)
    bx.text(t, 27, lab, ha='center', va='top', fontsize=FS_NOTE, color='0.45',
            linespacing=1.2)
bx.set_xlim(0.5, 1.0); bx.set_ylim(20, 100)
bx.set_xlabel(r'horizontal dipole ratio  $\Theta$')
bx.set_ylabel('EQE  (%)')
bx.tick_params(labelsize=FS_TICK)
bx.legend(loc='center left', frameon=False, fontsize=FS_NOTE, handlelength=1.8,
          borderaxespad=0.4)
for name, c, dy in (('proposed', C_PRO, -13), ('reference', C_REF, 7)):
    m = A[A['device'] == name]
    o = np.argsort(m['theta_horizontal'])
    lo, hi = 100 * m['EQE'][o][0], 100 * m['EQE'][o][-1]
    bx.annotate('%.0f → %.0f %%' % (lo, hi), (1.0, hi),
                textcoords='offset points', xytext=(-4, dy), ha='right',
                fontsize=FS_NOTE, color=c, fontweight='bold')
for s in ('top', 'right'):
    bx.spines[s].set_visible(False)
for e in ('png', 'pdf'):
    fig.savefig(os.path.join(HERE, 'fig3_theta.' + e), dpi=400)
print('fig3_theta.png written')
