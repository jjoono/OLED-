"""ED Fig. 7d from the MATLAB results (ED7d_etl_absorption.csv): figure (png/pdf/svg) and raw-data xlsx.
Usage: python3 plot_ED7d.py <out_dir>"""
import sys, os, csv, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from openpyxl import Workbook
from openpyxl.styles import Font
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = sys.argv[1] if len(sys.argv) > 1 else HERE
R = list(csv.DictReader(open(os.path.join(HERE, 'ED7d_etl_absorption.csv'))))
KS = (0.0, 3e-4, 1e-3, 3e-3); LAB = ('0', r'3$\times$10$^{-4}$', r'1$\times$10$^{-3}$', r'3$\times$10$^{-3}$')
COL = ('#440154', '#31688e', '#35b779', '#b8a100')
D = np.array(sorted({int(r['d_ETL_nm']) for r in R}))
def col(k, q): return np.array([float(r[q]) for r in sorted((r for r in R if float(r['k_ETL']) == k), key=lambda r: int(r['d_ETL_nm']))])
E = {k: 100 * col(k, 'EQE') for k in KS}

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 6.5, 'axes.linewidth': 0.6, 'mathtext.default': 'regular',
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': False, 'ytick.right': True,
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 3, 'ytick.major.size': 3,
                     'svg.fonttype': 'none', 'pdf.fonttype': 42})
fig, ax = plt.subplots(figsize=(2.283, 1.81))
for k, lab, c in zip(KS, LAB, COL):
    ax.plot(D, E[k], color=c, lw=1.0, label=lab, zorder=2)
    if k > 0:
        i = int(np.argmax(E[k])); ax.plot(D[i], E[k][i], 'o', ms=3.2, mfc='white', mec=c, mew=0.8, zorder=3)
ax.set_xlim(0, 400); ax.set_ylim(50, 92); ax.set_xticks([0, 100, 200, 300, 400]); ax.set_yticks([50, 60, 70, 80, 90])
ax.set_xlabel(r'ETL thickness, $\it{d}$$_{ETL}$ (nm)'); ax.set_ylabel('EQE (%)')
lg = ax.legend(title=r'$\it{k}$$_{ETL}$', ncol=2, frameon=False, loc='lower right', handlelength=1.6, columnspacing=1.0, borderaxespad=0.3, fontsize=5.6)
lg.get_title().set_fontsize(6)
fig.text(0.01, 0.97, 'd', fontweight='bold', fontsize=8, va='top')
fig.subplots_adjust(left=0.135, right=0.95, bottom=0.175, top=0.985)
for ext in ('png', 'pdf', 'svg'): fig.savefig(os.path.join(OUT, 'ED7d.' + ext), dpi=600)

wb = Workbook(); ws = wb.active; ws.title = 'EQE (%)'
ws.append(['Extended Data Fig. 7d: Ag 100 nm / B3PyMPM d_ETL (k_ETL added to n_o and n_e) / TCTA:B3PyMPM 25 nm (emitter at the centre) / TAPC 200 nm / ITO 50 nm / '
           'n_sub = 1.8 + hemispherical MLA (LightTools BSDF, n_MLA = 1.80); Ir(ppy)2acac spectrum 400-700 nm (1-nm steps), Theta = 2/3, PLQY = 1; '
           'angle-resolved model in MATLAB, series summed to convergence (last term < 1e-6 of the sum)'])
ws['A1'].font = Font(bold=True)
ws.append(['d_ETL (nm)'] + ['k_ETL = %g' % k for k in KS])
for j, d in enumerate(D): ws.append([int(d)] + [round(float(E[k][j]), 3) for k in KS])
ws.append([]); ws.append(['Maximum EQE (open circles for k_ETL > 0)'] + ['%.1f%% at %d nm' % (E[k].max(), D[np.argmax(E[k])]) for k in KS])
ws2 = wb.create_sheet('all quantities'); ws2.append(['k_ETL', 'd_ETL (nm)', 'eta_sub', 'eta_sta', 'EQE', 'round trips summed'])
for k in KS:
    for j, d in enumerate(D):
        ws2.append([k, int(d), float(col(k, 'eta_sub')[j]), float(col(k, 'eta_ext')[j]), float(col(k, 'EQE')[j]), int(col(k, 'n_terms')[j])])
wb.save(os.path.join(OUT, 'ED7d_rawdata.xlsx'))
print('\n'.join('k %-6g max %.1f%% @ %d nm, 400 nm %.1f%%' % (k, E[k].max(), D[np.argmax(E[k])], E[k][-1]) for k in KS))
