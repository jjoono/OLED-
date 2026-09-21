# -*- coding: utf-8 -*-
"""The converged planar map that feeds the MLA step."""
import numpy as np, matplotlib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
A = np.genfromtxt(os.path.join(HERE, 'mla_planar_map.csv'), delimiter=',', names=True)
dE = np.unique(A['d_ETL_nm']); dH = np.unique(A['d_HTL_nm'])
def grid(k): return A[k].reshape(len(dE), len(dH)).T
plt.rcParams.update({'font.size': 7.2, 'axes.linewidth': 0.7,
                     'xtick.major.width': 0.7, 'ytick.major.width': 0.7})
fig, ax = plt.subplots(1, 4, figsize=(11.2, 2.7))
for a, (k, ttl, cm) in zip(ax, [('eta_sub', r'$\eta_{\rm sub}$  (air + substrate)', 'viridis'),
                                ('spp', 'SPP', 'magma'),
                                ('wg', 'waveguided', 'cividis'),
                                ('abs', 'absorbed', 'inferno')]):
    m = a.pcolormesh(dE, dH, 100 * grid(k), cmap=cm, shading='nearest', rasterized=True)
    cb = fig.colorbar(m, ax=a, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=6)
    cb.set_label('%', fontsize=6.4)
    a.set_title(ttl, fontsize=7.2, pad=3)
    a.set_xlabel(r'$d_{\rm ETL}$  (nm)')
    a.tick_params(labelsize=6.6)
ax[0].set_ylabel(r'$d_{\rm HTL}$  (nm)')
fig.tight_layout(pad=0.4, w_pad=1.2)
fig.savefig(os.path.join(HERE, 'mla_planar_map.png'), dpi=300)
print('mla_planar_map.png written')
