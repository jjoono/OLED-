import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # tooling.xlsx
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fits, tables, figures
import sys, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(4, 2, figsize=(13, 13), sharex=True)
groups = [('1_2', 'run 1, Si'), ('2_2', 'run 2, Si'), ('1_3', 'run 1, glass'), ('2_3', 'run 2, glass')]
for g, (pre, title) in enumerate(groups):
    for k, suf in enumerate('1234'):
        sh = pre + suf
        wl, P, D, ang = L.se(sh, 245, 1100)
        N, C, Sx = L.ncs(P, D)
        c = plt.cm.tab10(k)
        ax[g, 0].plot(wl, P[:, 2], color=c, lw=1, label=sh)
        ax[g, 1].plot(wl, Sx[:, 2], color=c, lw=1, label=sh)
    ax[g, 0].set_title('%s  Psi @55 deg' % title); ax[g, 1].set_title('%s  S = sin2Psi sinDelta @55 deg' % title)
    for a in ax[g]: a.grid(alpha=.3); a.legend(fontsize=7)
for a in ax[-1]: a.set_xlabel('wavelength (nm)')
fig.tight_layout(); fig.savefig(_os.path.join(ELLIPS_OUT, 'overview.png'), dpi=110)
