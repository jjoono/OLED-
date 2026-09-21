# -*- coding: utf-8 -*-
"""EQE against aspect ratio for three lens spacings: touching (pitch 2r, fill factor 0.907),
pitch 2.4r (0.63) and 2.8r (0.46).  The hemisphere is penalised only when the lenses touch."""
import numpy as np, matplotlib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
C_AL, C_AG = '#1A1A1A', '#C00000'
FS = 7.4
plt.rcParams.update({'font.size': FS, 'axes.linewidth': 0.8})
A = np.genfromtxt(os.path.join(HERE, 'mla_aspect.csv'), delimiter=',', names=True, dtype=None, encoding='utf-8')
A = A[A['aspect_ratio'] >= 0.19]
sets = [('touching, pitch 2.0 r (fill factor 0.91)', A, '-')]
for pf, ff, ls in ((2.4, 0.63, '--'), (2.8, 0.46, ':')):
    fn = os.path.join(HERE, 'mla_aspect_pitch%d.csv' % round(pf * 10))
    if os.path.exists(fn):
        D = np.genfromtxt(fn, delimiter=',', names=True, dtype=None, encoding='utf-8')
        sets.append(('pitch %.1f r (fill factor %.2f)' % (pf, ff), D, ls))
fig, ax = plt.subplots(figsize=(3.6, 2.6))
fig.subplots_adjust(left=0.14, right=0.97, top=0.93, bottom=0.17)
for lab, D, ls in sets:
    o = np.argsort(D['aspect_ratio']); D = D[o]
    for k, c in (('EQE_Ag', C_AG), ('EQE_Al', C_AL)):
        ax.plot(D['aspect_ratio'], 100 * D[k], color=c, lw=1.4, ls=ls, label=lab if k == 'EQE_Ag' else None)
        i = int(np.argmax(D[k])); ax.plot([D['aspect_ratio'][i]], [100 * D[k][i]], marker='v', ms=4, color=c, mew=0)
ax.axvline(1.0, color='0.8', lw=0.7, dashes=(2.5, 2.5), zorder=0)
ax.set_xlim(0.1, 1.55); ax.set_ylim(30, 95)
ax.set_xlabel('aspect ratio  $h/r$'); ax.set_ylabel('EQE (%)')
ax.text(0.02, 0.97, 'red Ag, black Al;  triangles mark the maximum', transform=ax.transAxes, va='top', fontsize=FS - 1)
ax.legend(loc='lower right', frameon=False, fontsize=FS - 1)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
fig.savefig(os.path.join(HERE, 'mla_aspect_pitch.png'), dpi=300)
print('mla_aspect_pitch.png written')
