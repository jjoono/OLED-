"""Where HATCN-seeded Ag lands on the T(550)-vs-Rs map of Zhang et al., Adv. Mater. 2017, 29, 1605177 (Fig. 5c).
T(550) is absolute (glass substrate incl. its back face, as in that paper), bare film and ZnS/HATCN/Ag/ZnS DMD optimised
for T(550), from dmd_transmittance.TRA with the joint nk of each Ag thickness.  Rs: measured on HATCN seeds
(3 / 5 nm on HATCN 4 nm, 8 nm on HATCN 5 nm); thicknesses without a measurement are not plotted.
The Zhang star and ITO point are read off their Fig. 5c by eye (approximate)."""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import minimize
import dmd_transmittance as M
FILES = {5: 'Ag_HATCN_1-7_5nm_nk.csv', 8: 'Ag_HATCN_2-6_8nm_nk.csv'}
RS = {5: 25.0, 8: 9.1}
i5 = list(M.L).index(550.)
rows = []
for d, fn in FILES.items():
    a = np.genfromtxt('/home/user/OLED-/sim/teoled_opt/nk_joint/' + fn, delimiter=',', names=True)
    M.AG = np.interp(M.L, a['wavelength_nm'], a['n_joint']) + 1j * np.interp(M.L, a['wavelength_nm'], a['k_joint'])
    bare = M.TRA(0, 0, dag=d)[i5]
    r = min((minimize(lambda x: -M.TRA(*np.abs(x), dag=d)[i5, 0], x0, method='Nelder-Mead', options=dict(xatol=0.5, fatol=1e-5))
             for x0 in ([30, 30], [50, 40], [60, 20])), key=lambda r: r.fun)
    d1, d2 = np.abs(r.x); o = M.TRA(d1, d2, dag=d)[i5]
    rows.append((d, RS[d], bare[0], o[0], o[2], d1, d2))
    print('Ag %d nm  Rs %.1f ohm/sq | bare T550 %.3f | DMD ZnS %.0f/HATCN/Ag/ZnS %.0f: T550 %.3f, A550 %.3f' % (d, RS[d], bare[0], d1, d2, o[0], o[2]))
fig, ax = plt.subplots(figsize=(5.2, 4.2))
for d, rs, tb, td, *_ in rows:
    ax.annotate('', (rs, 100 * td), (rs, 100 * tb), arrowprops=dict(arrowstyle='->', color='0.6'))
    ax.plot(rs, 100 * tb, 'o', mfc='white', mec='C0', ms=8); ax.plot(rs, 100 * td, 'o', color='C0', ms=8)
    ax.text(rs * (1.12 if d != 8 else 0.89), 100 * td + (0 if d != 8 else 1.6), 'Ag %d nm' % d, va='center', ha='left' if d != 8 else 'right', fontsize=9)
ax.plot(20, 92.8, '*', color='navy', ms=13, label='Zhang 2017 Al-doped Ag DMD (approx.)')
ax.plot(13.5, 91.8, 's', color='k', ms=6, label='ITO in their Fig. 5c (approx.)')
ax.plot([], [], 'o', mfc='white', mec='C0', label='HATCN/Ag, bare (calc. T, measured Rs)')
ax.plot([], [], 'o', color='C0', label='ZnS/HATCN/Ag/ZnS DMD (calc.)')
ax.set_xscale('log'); ax.set_xlim(1, 1000); ax.set_ylim(30, 100)
ax.set_xlabel('Sheet resistance (ohm/sq)'); ax.set_ylabel('Transmittance at 550 nm (%, absolute)')
ax.legend(fontsize=7.5, loc='lower right', frameon=False); ax.grid(alpha=0.3, which='both')
fig.tight_layout(); fig.savefig('dmd_rs_plot.png', dpi=160)
