"""Raw data for eta_ext vs n_sub: matrix series vs equation (3), and their difference.

Source: sim/audit/fig2def_fine.csv (same stack/conditions as main-text Fig. 2d-f:
substrate / ITO 50 nm / 420 nm organic (n = 1.8) / 100 nm reflector, 550 nm,
isotropic emitter, PLQY = 1, hexagonal hemispherical MLA BSDF with n_MLA = n_sub).
The eq. (3) column is identical to figures/update/figure_update_rawdata.xlsx, sheet Fig2e.
"""
import csv, os
import openpyxl
from openpyxl.styles import Font
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
rows = list(csv.DictReader(open(os.path.join(ROOT, 'sim/audit/fig2def_fine.csv'))))
D = {R: sorted(((float(r['n_sub']), float(r['eta_ext_series']), float(r['eta_ext_closed']))
                for r in rows if r['reflector'] == R)) for R in ('Al', 'Ag')}
n = [x[0] for x in D['Al']]
CUTOFF = {'Ag': [1.46, 1.59, 1.76, 1.78], 'Al': [1.42, 1.58, 1.75, 1.78]}
BANDS = [(1.42, 1.46), (1.58, 1.59), (1.75, 1.78), (1.80, 1.80)]

wb = openpyxl.Workbook()
ws = wb.active; ws.title = 'README'
for line in [
    'eta_ext vs n_sub (= n_MLA): matrix series vs equation (3)',
    'Stack: substrate / ITO 50 nm / organic 420 nm (n = 1.8, lossless) / Al or Ag 100 nm; 550 nm; isotropic emitter (Theta = 2/3); PLQY = 1.',
    'Outcoupling: hexagonally close-packed hemispherical microlens array, BSDF with n_MLA = n_sub (ray tracing, 0.05 grid, interpolated to 0.01).',
    'Series: eta_ext = Psub (I - B_R R)^-1 B_T^T with the device angular distribution Psub(theta).',
    'Eq. (3): eta_ext = p / (p + (1 - p) A\'), p = cos(theta)sin(theta)-weighted B_T, A\' = cos(theta)sin(theta)-weighted (1 - R_LED).',
    'Eq. (3) values are identical to main-text Fig. 2e (figure_update_rawdata.xlsx, sheet Fig2e).',
    'Sheet c: panel (c), solid = series, dashed = eq. (3).  Sheet d: panel (d), difference series - eq. (3) in percentage points.',
    'Sheet cutoffs: substrate indices at which a waveguided mode reaches cut-off (n_eff = n_sub) at 550 nm; grey bands in panel (d).',
    'Source data: sim/audit/fig2def_fine.csv. Script: figures/fig2e_series_vs_eq3/make_rawdata.py',
]:
    ws.append([line])
ws.column_dimensions['A'].width = 140

c = wb.create_sheet('c_eta_ext')
c.append(['n_sub (= n_MLA)', 'Al series', 'Al eq. (3)', 'Ag series', 'Ag eq. (3)'])
for i, x in enumerate(n):
    c.append([x, D['Al'][i][1], D['Al'][i][2], D['Ag'][i][1], D['Ag'][i][2]])

d = wb.create_sheet('d_difference')
d.append(['n_sub (= n_MLA)', 'Al: series - eq.(3) (%p)', 'Ag: series - eq.(3) (%p)'])
diff = {R: [100 * (s - e) for _, s, e in D[R]] for R in D}
for i, x in enumerate(n):
    d.append([x, round(diff['Al'][i], 3), round(diff['Ag'][i], 3)])

k = wb.create_sheet('cutoffs')
k.append(['mode order', 'Al cathode n_sub cut-off', 'Ag cathode n_sub cut-off'])
for j, lab in enumerate(['TM0', 'TE0', 'TM1', 'TE1']):
    k.append([lab, CUTOFF['Al'][j], CUTOFF['Ag'][j]])
k.append(['n_sub = n_org', 1.80, 1.80])
k.append([]); k.append(['grey bands for panel (d)', 'from', 'to'])
for a, b in BANDS: k.append(['', a, b])

s = wb.create_sheet('summary')
s.append(['quantity', 'Al', 'Ag'])
for R in ('Al', 'Ag'):
    pass
mins = {R: min(zip(diff[R], n)) for R in D}
s.append(['largest drop (%p)', round(mins['Al'][0], 1), round(mins['Ag'][0], 1)])
s.append(['at n_sub', mins['Al'][1], mins['Ag'][1]])
def outside(R):
    v = [abs(diff[R][i]) for i, x in enumerate(n) if not any(a - 0.005 <= x <= b + 0.035 for a, b in BANDS)]
    return round(max(v), 1), round(sorted(v)[len(v)//2], 1)
for lab, j in (('max |diff| outside cut-off regions (%p)', 0), ('median |diff| outside cut-off regions (%p)', 1)):
    s.append([lab, outside('Al')[j], outside('Ag')[j]])
for sh in wb.worksheets[1:]:
    for cell in sh[1]: cell.font = Font(bold=True)
    sh.column_dimensions['A'].width = 22
    for col in 'BCDE': sh.column_dimensions[col].width = 24
wb.save(os.path.join(HERE, 'fig2e_series_vs_eq3_rawdata.xlsx'))

# preview
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
col = {'Al': '#c0392b', 'Ag': '#2c6fbb'}
for R in ('Al', 'Ag'):
    ax[0].plot(n, [x[1] for x in D[R]], color=col[R], lw=1.6, label=f'{R} series')
    ax[0].plot(n, [x[2] for x in D[R]], color=col[R], lw=1.4, ls='--', label=f'{R} eq. (3)')
    ax[1].plot(n, diff[R], color=col[R], lw=1.6, label=R)
for a, b in BANDS:
    ax[1].axvspan(a - 0.004, b + 0.004, color='0.85', zorder=0)
ax[1].axhline(0, color='0.4', lw=0.6)
ax[0].set(xlabel='$n_{sub}$ = $n_{MLA}$', ylabel=r'$\eta_{ext}$', xlim=(1.3, 2.0)); ax[0].legend(fontsize=7, frameon=False)
ax[1].set(xlabel='$n_{sub}$ = $n_{MLA}$', ylabel=r'series $-$ eq. (3) (%p)', xlim=(1.3, 2.0)); ax[1].legend(fontsize=7, frameon=False)
ax[0].text(0.02, 0.04, 'c', transform=ax[0].transAxes, weight='bold'); ax[1].text(0.02, 0.04, 'd', transform=ax[1].transAxes, weight='bold')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'fig2e_series_vs_eq3_preview.png'), dpi=200)
print(mins, outside('Al'), outside('Ag'))
