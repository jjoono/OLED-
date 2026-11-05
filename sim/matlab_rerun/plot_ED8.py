"""ED Fig. 8 from ED8_eqe.csv (MATLAB series + LightTools BSDFs): rawdata xlsx and a preview figure.
Usage: python3 plot_ED8.py <out_dir>"""
import sys, os, csv, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from openpyxl import Workbook
from openpyxl.styles import Font
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = sys.argv[1]
R = list(csv.DictReader(open(os.path.join(HERE, 'ED8_eqe.csv'))))
pl = [r for r in R if r['id'] == 'planar'][0]
A = sorted([r for r in R if r['family'] == 'aspect'], key=lambda r: float(r['aspect_ratio']))
def ffset(ar): return sorted([r for r in R if r['family'] == 'fill_factor' and abs(float(r['aspect_ratio']) - ar) < 1e-6], key=lambda r: float(r['fill_factor']))
F5, F10 = ffset(0.5), ffset(1.0)
g = lambda rows, k: np.array([100 * float(r[k]) for r in rows])
ar = np.array([float(r['aspect_ratio']) for r in A])
cond = ('Stack: Al or Ag (100 nm)/ETL (200 nm)/EML (20 nm, emitter at the centre)/HTL (200 nm), all n = 1.8/ITO (50 nm, 1.864+0.0032i)/substrate n = 1.8 '
        '+ hexagonal array of half-spheroidal lenses (n = 1.8); 550 nm, PLQY = 1, Theta = 2/3; BSDFs from LightTools; trans-scale model in MATLAB, series summed to convergence.')
wb = Workbook(); ws = wb.active; ws.title = 'ED8a aspect ratio'
ws.append(['Extended Data Fig. 8a. ' + cond + ' Close-packed (fill factor 0.907).']); ws['A1'].font = Font(bold=True)
ws.append(['Aspect ratio (h/r)', 'EQE Al (%)', 'EQE Ag (%)', 'p (single-pass escape, cos·sin weighted)'])
for r in A: ws.append([float(r['aspect_ratio']), round(100 * float(r['EQE_Al']), 3), round(100 * float(r['EQE_Ag']), 3), round(float(r['p_single_pass']), 4)])
ws.append([]); ws.append(['Planar substrate (no lenses)', round(100 * float(pl['EQE_Al']), 3), round(100 * float(pl['EQE_Ag']), 3)])
for k, lab in (('EQE_Al', 'Al'), ('EQE_Ag', 'Ag')):
    e = g(A, k); m = e.max(); ok = ar[e >= 0.95 * m]
    ws.append(['%s: max %.2f%% at AR %.2f; >=95%% of max for AR %.2f-%.2f' % (lab, m, ar[e.argmax()], ok.min(), ok.max())])
ws2 = wb.create_sheet('ED8b fill factor')
ws2.append(['Extended Data Fig. 8b. ' + cond + ' Base radius fixed, pitch increased; fill factor 0.907 = close-packed.']); ws2['A1'].font = Font(bold=True)
ws2.append(['Fill factor', 'EQE Al, AR 0.5 (%)', 'EQE Ag, AR 0.5 (%)', 'EQE Al, AR 1.0 (%)', 'EQE Ag, AR 1.0 (%)'])
for a, b in zip(F5, F10):
    ws2.append([round(float(a['fill_factor']), 4), round(100 * float(a['EQE_Al']), 3), round(100 * float(a['EQE_Ag']), 3), round(100 * float(b['EQE_Al']), 3), round(100 * float(b['EQE_Ag']), 3)])
wb.save(os.path.join(OUT, 'ED8_rawdata.xlsx'))

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'axes.linewidth': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in',
                     'xtick.top': True, 'ytick.right': True, 'mathtext.default': 'regular'})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(5.6, 2.3))
a1.plot(ar, g(A, 'EQE_Al'), color='#e8731c', lw=1.0, label='Al'); a1.plot(ar, g(A, 'EQE_Ag'), color='#2a55e5', lw=1.0, label='Ag')
a1.set_xlim(0.1, 1.5); a1.set_ylim(30, 90); a1.set_xlabel(r'Aspect ratio, $h_{lens}/r_{lens}$'); a1.set_ylabel('EQE (%)'); a1.legend(frameon=False, loc='lower right')
ff = np.array([float(r['fill_factor']) for r in F5])
for rows, ls, lab in ((F5, '-', 'AR 0.5'), (F10, '--', 'AR 1.0')):
    a2.plot(ff, g(rows, 'EQE_Al'), ls, color='#e8731c', lw=1.0); a2.plot(ff, g(rows, 'EQE_Ag'), ls, color='#2a55e5', lw=1.0, label=lab)
a2.set_xlim(0.08, 0.92); a2.set_ylim(30, 90); a2.set_xlabel('Fill factor'); a2.set_ylabel('EQE (%)'); a2.legend(frameon=False, loc='lower right')
for ax, l in ((a1, 'a'), (a2, 'b')): ax.text(-0.18, 1.02, l, transform=ax.transAxes, fontweight='bold', fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'ED8_preview.png'), dpi=300)
