# -*- coding: utf-8 -*-
"""Fig 3(b) mock: where the power goes as the organic layer thickens and the
substrate index is raised, for  Ag | organic | ITO 50 nm | glass  at 550 nm.

Left   the two loss channels: SPP (solid) dies with distance from the Ag,
       waveguided power (dashed) dies as n_sub approaches n_organic.
Right  what is left, the power delivered into the substrate mode.
Colour is the substrate index, light to dark with increasing n."""
import numpy as np, matplotlib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
C = {1.50: '#9ECAE1', 1.65: '#4292C6', 1.80: '#08519C'}
DASH = (0, (3.2, 1.8))
FS_LAB, FS_TICK, FS_NOTE, FS_LET = 7.2, 6.6, 6.2, 9.0
plt.rcParams.update({'font.size': FS_LAB, 'axes.linewidth': 0.7,
                     'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
                     'xtick.major.size': 2.4, 'ytick.major.size': 2.4})

A = np.genfromtxt(os.path.join(HERE, 'fig3b_modes.csv'), delimiter=',', names=True)
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.5))

for ns in (1.50, 1.65, 1.80):
    m = A[np.isclose(A['n_substrate'], ns)]
    d = m['d_organic_nm']
    ax[0].plot(d, 100 * m['EQE_spp'], color=C[ns], lw=1.3)
    ax[0].plot(d, 100 * m['EQE_wg'], color=C[ns], lw=1.3, linestyle=DASH)
    ax[1].plot(d, 100 * m['eta_sub_total'], color=C[ns], lw=1.4)

ax[0].axhline(10, color='0.6', lw=0.6, dashes=(1.2, 1.8))
ax[0].text(495, 11.5, '10 %', ha='right', va='bottom', fontsize=FS_NOTE, color='0.45')
ax[0].set_ylabel('power fraction  (%)')
ax[0].set_ylim(0, 100)
ax[0].legend(handles=[Line2D([], [], color='0.35', lw=1.3, label='SPP'),
                      Line2D([], [], color='0.35', lw=1.3, linestyle=DASH,
                             label='waveguided')],
             loc='upper right', frameon=False, fontsize=FS_NOTE,
             handlelength=2.2, borderaxespad=0.3)

ax[1].set_ylabel(r'into the substrate mode  $\eta_{\rm sub}$  (%)')
ax[1].set_ylim(0, 100)
ax[1].legend(handles=[Line2D([], [], color=C[n], lw=1.4,
                             label=r'$n_{\rm sub}$ = %.2f' % n)
                      for n in (1.50, 1.65, 1.80)],
             loc='upper left', frameon=False, fontsize=FS_NOTE,
             handlelength=1.6, borderaxespad=0.3)

for k, a in enumerate(ax):
    a.set_xlim(0, 500)
    a.set_xlabel('organic thickness  (nm)')
    a.tick_params(labelsize=FS_TICK)
    a.text(-0.20, 1.06, 'b' + ('i' if k == 0 else 'ii'), transform=a.transAxes,
           fontsize=FS_LET, fontweight='bold', va='top')

fig.tight_layout(pad=0.4, w_pad=1.6)
for e in ('png', 'pdf'):
    fig.savefig(os.path.join(HERE, 'fig3b_mock.' + e), dpi=400)
print('fig3b_mock.png / .pdf written')

# ---- supplementary: the u-grid artefact -----------------------------------
G = np.genfromtxt(os.path.join(HERE, 'fig3b_ugrid_demo.csv'),
                  delimiter=',', names=True)
fig2, bx = plt.subplots(1, 2, figsize=(6.6, 2.5))
for a, key, ttl in ((bx[0], 'wg', 'waveguided'), (bx[1], 'spp', 'SPP')):
    a.plot(G['d_organic_nm'], 100 * G[key + '_N1000'], color='#D55E00', lw=0.8,
           marker='o', ms=1.6, mew=0, label='uniform $u$ grid, $N$ = 1000')
    a.plot(G['d_organic_nm'], 100 * G[key + '_quad'], color='#08519C', lw=1.5,
           label='substitution + Simpson')
    a.set_xlim(0, 500)
    a.set_xlabel('organic thickness  (nm)')
    a.set_ylabel('%s power  (%%)' % ttl)
    a.tick_params(labelsize=FS_TICK)
    a.legend(loc='best', frameon=False, fontsize=FS_NOTE, handlelength=2.0)
bx[0].set_title(r'$n_{\rm sub}$ = 1.50, $\lambda$ = 550 nm', fontsize=FS_NOTE)
fig2.tight_layout(pad=0.4, w_pad=1.6)
fig2.savefig(os.path.join(HERE, 'fig3b_ugrid_demo.png'), dpi=400)
print('fig3b_ugrid_demo.png written')
