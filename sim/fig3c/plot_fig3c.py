# -*- coding: utf-8 -*-
"""Fig 3(c): a low out-of-plane index in the ETL buys the same plasmon
suppression at a thinner layer.

ci   the plasmon sits at a lower in-plane index and carries less power when the
     ETL's out-of-plane index is low.
cii  so the same SPP loss is reached with a thinner ETL, which is also what the
     drive voltage wants."""
import numpy as np, matplotlib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
D_SPEC, N_SUB, TARGET = 60.0, 1.80, 0.20
SER = [('isotropic ETL', 'n_e=1.80', 1.80, '#4D4D4D', '-'),
       ('B3PyMPM', 'B3PyMPM', 1.609, '#D55E00', '-'),
       ('B4PyMPM', 'B4PyMPM', 1.560, '#009E73', '-'),
       (r'$n_e$ = 1.50 (model)', 'n_e=1.50', 1.50, '#0072B2', (0, (3.0, 1.6)))]
FS_LAB, FS_TICK, FS_NOTE, FS_LET = 7.2, 6.6, 6.2, 9.0
plt.rcParams.update({'font.size': FS_LAB, 'axes.linewidth': 0.7,
                     'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
                     'xtick.major.size': 2.4, 'ytick.major.size': 2.4})

fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.5))

# ---- ci : where the plasmon sits ------------------------------------------
S = np.genfromtxt(os.path.join(HERE, 'fig3c_spectrum.csv'), delimiter=',',
                  names=True, dtype=None, encoding='utf-8')
S = S[(S['d_ETL_nm'] == D_SPEC) & (np.abs(S['n_eff'] - N_SUB) > 0.009)]
a = ax[0]
a.axvspan(N_SUB, 2.6, color='#F0F0F0', zorder=0, lw=0)
a.axvline(N_SUB, color='0.45', lw=0.8, dashes=(2.4, 1.6), zorder=3)
for lab, key, ne, c, ls in SER:
    m = S[S['series'] == key]
    a.semilogy(m['n_eff'], m['TM'], color=c, lw=1.3, linestyle=ls, zorder=4)
a.set_xlim(1.5, 2.35)
a.set_ylim(1e-2, 60)
a.set_xlabel(r'in-plane effective index   $k_x/k_0$')
a.set_ylabel('TM power dissipation density')
a.text(N_SUB - 0.015, 40, r'$n_{\rm sub}$', fontsize=FS_NOTE, color='0.45',
       ha='right', va='top')
a.text(0.97, 0.96, r'$d_{\rm ETL}$ = %.0f nm' % D_SPEC, transform=a.transAxes,
       ha='right', va='top', fontsize=FS_NOTE, color='0.3')
a.legend(handles=[Line2D([], [], color=c, lw=1.3, linestyle=ls, label=lab)
                  for lab, _, _, c, ls in SER],
         loc='lower left', frameon=False, fontsize=FS_NOTE,
         handlelength=1.9, labelspacing=0.3, borderaxespad=0.25)

# ---- cii : the thickness it buys ------------------------------------------
A = np.genfromtxt(os.path.join(HERE, 'fig3c_sweep.csv'), delimiter=',',
                  names=True, dtype=None, encoding='utf-8')
b = ax[1]
b.axhline(100 * TARGET, color='0.65', lw=0.7, dashes=(1.2, 1.8), zorder=2)
cross = {}
for lab, key, ne, c, ls in SER:
    sel = ((A['series'] == key.split('=')[0] if 'Py' in key else
            (A['series'] == 'no1.80') & np.isclose(A['n_e'], ne))
           if 'Py' not in key else (A['series'] == key))
    m = A[sel]
    o = np.argsort(m['d_ETL_nm'])
    d, s = m['d_ETL_nm'][o], m['spp'][o]
    b.plot(d, 100 * s, color=c, lw=1.4, linestyle=ls, zorder=4)
    i = np.where(s <= TARGET)[0][0]
    x = np.interp(TARGET, [s[i], s[i - 1]], [d[i], d[i - 1]])
    cross[key] = x
    b.plot([x], [100 * TARGET], marker='o', ms=3.6, color=c, mew=0, zorder=6)
    b.plot([x, x], [0, 100 * TARGET], color=c, lw=0.7, dashes=(1.2, 1.6),
           zorder=3)
b.set_xlim(0, 180)
b.set_ylim(0, 82)
b.set_xlabel(r'ETL thickness  $d_{\rm ETL}$  (nm)')
b.set_ylabel('power lost to the SPP  (%)')
b.text(178, 100 * TARGET + 1.5, 'same SPP loss', ha='right', va='bottom',
       fontsize=FS_NOTE, color='0.55')
b.text(0.97, 0.955, 'reached at  %s nm'
       % ' / '.join('%.0f' % cross[k] for _, k, _, _, _ in SER),
       transform=b.transAxes, ha='right', va='top', fontsize=FS_NOTE,
       color='0.3')

for k, a_ in enumerate(ax):
    a_.tick_params(labelsize=FS_TICK)
    a_.text(-0.19, 1.07, 'c' + ('i' if k == 0 else 'ii'), transform=a_.transAxes,
            fontsize=FS_LET, fontweight='bold', va='top')

fig.tight_layout(pad=0.4, w_pad=1.8)
for e in ('png', 'pdf'):
    fig.savefig(os.path.join(HERE, 'fig3c_mock.' + e), dpi=400)
print('fig3c_mock.png written;  d at SPP = %.0f%% : ' % (100 * TARGET)
      + ', '.join('%s %.0f nm' % (k, v) for k, v in cross.items()))
