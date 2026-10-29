"""Re-plot of last week's devices (same evaluation as run_reference.py / run_afterAgMLA.py, -90..90 deg angle correction):
planar glass reference d32_p5 (2026-10-22) and MLA + Ag back reflector d44_p3 (2026-10-23, two angular datasets).
-> out/lastweek_devices.png"""
import os, sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import pixel_report as P
H = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(H, '..', 'data')
C1, C2, C3 = '#2a78d6', '#eb6834', '#1baf7a'                 # categorical slots 1-3 (validated default palette)
runs = [('Planar glass (reference)', os.path.join(H, 'out', 'ref_d32_p5', '2026-10-22_glass_d32_p5_jvl_dequoted.csv'),
         os.path.join(D, '2026-10-22_glass_d32_p5_gon-spec-deg.csv'), C1, '-'),
        ('MLA + Ag reflector (angles d44_p1)', os.path.join(D, '2026-10-23_afterAgMLA_d44_p3_jvl.csv'),
         os.path.join(D, '2026-10-23_afterAgMLA_d44_p1_gon-spec-deg.csv'), C2, '-'),
        ('MLA + Ag reflector (angles d45_p3)', os.path.join(D, '2026-10-23_afterAgMLA_d44_p3_jvl.csv'),
         os.path.join(D, '2026-10-23_afterAgMLA_d45_p3_gon-spec-deg.csv'), C2, '--')]
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#888',
                     'axes.labelcolor': '#333', 'xtick.color': '#555', 'ytick.color': '#555'})
fig, ax = plt.subplots(1, 3, figsize=(14, 4.3))
for name, jvl, deg, col, ls in runs:
    r = P.evaluate(jvl, deg); f = r['full']; ok = (f.luminance > 1) & np.isfinite(f.eqe)
    L, E = np.asarray(f.luminance[ok], float), np.asarray(f.eqe[ok], float)
    ax[0].semilogx(L, E, ls, color=col, lw=2, label=name)
    ax[1].semilogy(np.asarray(f.voltage), np.clip(np.asarray(f.current_density, float), 1e-6, None), ls, color=col, lw=2)
    ax[2].semilogx(L, np.asarray(f.current_efficiency[ok], float), ls, color=col, lw=2)
    k = np.nanargmax(E); ax[0].annotate('max %.1f %%' % E[k], (L[k], E[k]), textcoords='offset points',
                                       xytext=(6, -13 if ls == '--' else 5), fontsize=9, color='#333')
ax[0].set(xlabel='Luminance (cd/m$^2$)', ylabel='EQE (%)', ylim=(0, 70), title='EQE')
ax[1].set(xlabel='Voltage (V)', ylabel='Current density (mA/cm$^2$)', title='J–V')
ax[2].set(xlabel='Luminance (cd/m$^2$)', ylabel='Current efficiency (cd/A)', title='Current efficiency')
for a in ax: a.grid(alpha=0.25, lw=0.6)
ax[0].legend(frameon=False, fontsize=8.5, loc='lower left')
fig.suptitle('Green Ir(ppy)$_2$acac OLEDs, EQE from −90…90° angle-resolved spectra', fontsize=11, color='#222')
fig.tight_layout(); fig.savefig(os.path.join(H, 'out', 'lastweek_devices.png'), dpi=160)
print('saved')
