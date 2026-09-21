# -*- coding: utf-8 -*-
"""The three lens constructions side by side: where the maximum sits depends on how
the lens is defined when its aspect ratio changes."""
import numpy as np, matplotlib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
C_AL, C_AG = '#1A1A1A', '#C00000'
FS = 7.4
plt.rcParams.update({'font.size': FS, 'axes.linewidth': 0.8})
A = np.genfromtxt(os.path.join(HERE, 'mla_aspect.csv'), delimiter=',', names=True, dtype=None, encoding='utf-8')
F = np.genfromtxt(os.path.join(HERE, 'mla_aspect_families.csv'), delimiter=',', names=True, dtype=None, encoding='utf-8')
A = A[A['aspect_ratio'] >= 0.09]
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5), sharey=True)
fig.subplots_adjust(left=0.07, right=0.985, top=0.84, bottom=0.2, wspace=0.12)
titles = {'A': 'A  cap of fixed base radius (close-packed at every AR)\nhalf-ellipsoid above AR = 1',
          'B': 'B  sphere of fixed radius = half pitch, cut at height h\n(base shrinks, gap grows); cylinder + hemisphere above 1',
          'C': 'C  hemisphere scaled in height at every AR\n(half-ellipsoid)'}
for ax, fam in zip(axs, ('A', 'B', 'C')):
    d = A if fam == 'A' else F[F['family'] == fam]
    o = np.argsort(d['aspect_ratio']); d = d[o]
    for k, c, lab in (('EQE_Ag', C_AG, 'Ag'), ('EQE_Al', C_AL, 'Al')):
        ax.plot(d['aspect_ratio'], 100 * d[k], color=c, lw=1.6, label=lab)
        i = int(np.argmax(d[k]))
        ax.plot([d['aspect_ratio'][i]], [100 * d[k][i]], marker='v', ms=5, color=c, mew=0)
        ax.annotate('max %.1f %% at AR %.2f' % (100 * d[k][i], d['aspect_ratio'][i]),
                    (d['aspect_ratio'][i], 100 * d[k][i]), xytext=(0, 6), textcoords='offset points',
                    ha='center', fontsize=FS - 1.2, color=c)
    ax.axvline(1.0, color='0.8', lw=0.7, dashes=(2.5, 2.5), zorder=0)
    ax.set_xlim(0.1, 1.55); ax.set_ylim(30, 100)
    ax.set_title(titles[fam], fontsize=FS - 0.8, loc='left')
    ax.set_xlabel('aspect ratio  $h/r$')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
axs[0].set_ylabel('EQE (%)')
axs[0].legend(loc='lower right', frameon=False, fontsize=FS - 0.6)
fig.savefig(os.path.join(HERE, 'mla_aspect_families.png'), dpi=300)
print('mla_aspect_families.png written')
