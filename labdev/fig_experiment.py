"""Main-text Fig. 5: fabricated devices versus the finite-geometry model (no fitted parameter except the common
internal efficiency taken from the bare devices). Measured values: 260328 measurement series (photodiode goniometer,
IOEL efficiency template; EQE = 0-deg EQE x angular factor), relayed from the lab analysis; S4D3 excluded
(intensity maximum at -40 deg, single file, EQE above the optical ceiling)."""
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyBboxPatch
from openpyxl import Workbook
plt.rcParams.update({'font.size': 7, 'axes.titlesize': 7.5, 'axes.spines.top': False, 'axes.spines.right': False})
C = {'bare': '#7d7c76', 'hemisphere': '#2a78d6', 'freeform': '#eb6834'}; INK2 = '#52514e'
# ---- measured
EQE34 = {'bare': [29.7, 27.7, 30.7, 31.4, 30.8], 'freeform': [22.9, 26.7, 23.7, 28.7], 'hemisphere': [22.9, 27.1, 22.5, 25.3, 24.1]}
PEm = {'bare': [1.01, 1.08, 0.95, 1.02, 0.97], 'freeform': [1.13, 0.99, 1.11, 1.03], 'hemisphere': [0.63, 0.69, 0.63, 0.71, 0.66]}
A10 = np.arange(0, 90, 10)
PROF = {'S3D2': [1.00, 0.97, 0.86, 0.71, 0.49, 0.32, 0.28, 0.13, 0.08], 'S2D4': [1.00, 0.96, 0.84, 0.69, 0.52, 0.37, 0.27, 0.20, 0.13],
        'S1D4': [1.00, 0.97, 0.85, 0.65, 0.44, 0.29, 0.21, 0.17, 0.15], 'S4D2': [1.00, 0.95, 0.84, 0.69, 0.53, 0.42, 0.33, 0.22, 0.10],
        'S5D4': [1.00, 0.96, 0.83, 0.64, 0.44, 0.32, 0.26, 0.23, 0.15]}
# ---- simulated
R = json.load(open('layout_results.json')); SW = json.load(open('layout_sweep.json'))
ETA = np.mean(EQE34['bare']) / (R['bare|patch0mm|rectAl']['EQE'] * 100)       # common internal efficiency from bare devices
E = np.arange(91); dOm = 2 * np.pi * (np.cos(np.radians(E[:-1])) - np.cos(np.radians(E[1:])))
def norm_prof(p, win=5):
    p = np.array(p); i0 = (p[:5] * dOm[:5]).sum() / dOm[:5].sum(); th = np.arange(90) + 0.5
    k = np.ones(win); num = np.convolve(p * dOm, k, 'same'); den = np.convolve(dOm, k, 'same')
    return th, num / den / i0
fig, ax = plt.subplots(2, 2, figsize=(7.1, 5.4), gridspec_kw=dict(hspace=0.45, wspace=0.3)); ax = ax.ravel(); sheets = []
# (a) layout
a = ax[0]; a.set_aspect('equal'); a.set_xlim(-0.5, 26); a.set_ylim(-0.5, 26); a.axis('off')
a.add_patch(Rectangle((0, 0), 25.4, 25.4, fill=False, lw=0.8, ec='k'))
a.add_patch(Rectangle((0, 12.7), 25.4, 12.7, color='#d9e6f2', lw=0, label='ITO'))
for xc in 25.4 - np.array([4.7, 8.6898, 12.6796, 16.6694, 20.6592]):
    a.add_patch(FancyBboxPatch((xc - 1.3919, 4.5618 - 1.3919), 2 * 1.3919, 9.9263 + 2 * 1.3919, boxstyle='round,pad=0,rounding_size=1.3919', fc='#c9c9c4', ec='#8a8a85', lw=0.4))
a.add_patch(Rectangle((1.7, 25.4 - 1.7497 - 4.57), 22.0, 4.57, fc='#c9c9c4', ec='#8a8a85', lw=0.4))
px, py = 12.7, 14.4882
a.add_patch(Rectangle((px - 5, py - 5), 10, 10, fill=False, ec='k', ls='--', lw=0.7))
a.add_patch(Rectangle((px - 1.5, py - 1.5), 3, 3, fill=False, ec=C['hemisphere'], lw=1.2))
a.add_patch(Circle((px, py), 1.0, color='#f2b134'))
a.text(px + 5.3, py - 6.6, '10 × 10 mm (prediction)', fontsize=5.5, color='k'); a.text(px + 1.8, py + 1.6, '3 × 3 mm\nlens patch', fontsize=5.5, color=C['hemisphere'])
a.annotate('pixel (r = 1 mm)', (px - 0.7, py - 0.7), (px - 11.5, py - 3.5), fontsize=5.5, arrowprops=dict(arrowstyle='-', lw=0.5)); a.text(0.6, 23.9, 'ITO (upper half)', fontsize=5.5, color='#4a6a8a'); a.text(17.5, 1.3, 'Al fingers', fontsize=5.5, color='#55554f')
a.text(0.0, -0.4, '25.4 mm glass; bare 0.7 mm, MLA 0.7 + 0.7 mm (oil-coupled)', fontsize=5.2, va='top')
a.set_title('(a) device layout', loc='left')
# (b) angular profiles
b = ax[1]; rows_b = []
for k, v in PROF.items(): b.plot(A10, v, 'o', ms=3, mfc='none', color=C['hemisphere'], mew=0.8)
b.plot([], [], 'o', ms=3, mfc='none', color=C['hemisphere'], label='hemispherical MLA, measured (5 devices)')
for key, col, ls, lab in (('hemisphere|patch3mm|rectAl', C['hemisphere'], '-', 'hemispherical, 3 mm patch'), ('hemisphere|patch10mm|rectAl', C['hemisphere'], '--', 'hemispherical, 10 mm patch'),
                          ('freeform|patch3mm|rectAl', C['freeform'], '-', 'freeform, 3 mm patch'), ('bare|patch0mm|rectAl', C['bare'], '-', 'bare')):
    th, y = norm_prof(R[key]['profile']); b.plot(th, y, ls, color=col, lw=1.0, label=lab + ' (simulated)'); rows_b += [[key, float(t), float(v)] for t, v in zip(th, y)]
b.set_xlim(0, 85); b.set_ylim(0, 1.55); b.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); b.set_xlabel('emission angle (°)'); b.set_ylabel('I(θ) / I(0°)'); b.legend(fontsize=5, frameon=False, loc='upper right', ncol=1)
b.set_title('(b) normalized angular intensity', loc='left')
sheets.append(('b_sim_profiles', ['case', 'angle_deg', 'I/I0'], rows_b)); sheets.append(('b_measured_hemi', ['device'] + [f'{a}deg' for a in A10], [[k] + v for k, v in PROF.items()]))
# (c) angular factor
c = ax[2]; names = ['bare', 'freeform', 'hemisphere']; rng = np.random.default_rng(1); rows_c = []
for i, nm in enumerate(names):
    v = PEm[nm]; c.plot(i + rng.uniform(-0.12, 0.12, len(v)), v, 'o', ms=3.5, color=C[nm], alpha=0.85)
    s3 = R['bare|patch0mm|rectAl']['PE'] if nm == 'bare' else R[f'{nm}|patch3mm|rectAl']['PE']
    c.plot([i - 0.3, i + 0.3], [s3, s3], '-', color='k', lw=1.4)
    if nm != 'bare':
        s10 = R[f'{nm}|patch10mm|rectAl']['PE']; c.plot(i + 0.38, s10, 'D', ms=4, mfc='white', mec='k', mew=0.8)
    rows_c.append([nm, s3, (R[f'{nm}|patch10mm|rectAl']['PE'] if nm != 'bare' else None)] + v)
c.plot([], [], '-', color='k', lw=1.4, label='simulated, 3 mm patch'); c.plot([], [], 'D', ms=4, mfc='white', mec='k', label='simulated, 10 mm patch'); c.plot([], [], 'o', color='#999', ms=3.5, label='measured devices')
c.axhline(1, color=INK2, lw=0.4, ls=':'); c.set_xticks(range(3)); c.set_xticklabels(['bare', 'freeform', 'hemispherical']); c.set_xlim(-0.6, 2.7)
c.set_ylim(0.5, 1.25); c.set_ylabel('angular factor  EQE / EQE$_{0°}$'); c.legend(fontsize=5.5, frameon=False, loc='lower left')
c.set_title('(c) angular factor', loc='left'); sheets.append(('c_angular_factor', ['structure', 'sim_3mm', 'sim_10mm', 'measured...'], rows_c))
# (d) EQE vs patch width
d = ax[3]; W = [2, 3, 4, 5, 7, 10, 14, 20]; rows_d = []
for nm in ('hemisphere', 'freeform'):
    y = [SW[f'{nm}|{w}']['EQE'] * 100 * ETA for w in W]; d.plot(W, y, '-o', color=C[nm], ms=2.5, lw=1.0, label=f'{nm} MLA, simulated')
    v = EQE34[nm]; d.plot(3 * np.exp(rng.uniform(-0.06, 0.06, len(v))), v, 's', ms=3.5, mfc='none', mec=C[nm], mew=0.9)
    rows_d += [[nm, w, yy] for w, yy in zip(W, y)]
bs = R['bare|patch0mm|rectAl']['EQE'] * 100 * ETA; d.axhline(bs, color=C['bare'], lw=1.0, label='bare, simulated')
d.axhspan(min(EQE34['bare']), max(EQE34['bare']), color=C['bare'], alpha=0.15, lw=0, label='bare, measured range')
d.plot([], [], 's', ms=3.5, mfc='none', mec='#777', label='MLA devices, measured (3 mm patch)')
d.set_xscale('log'); d.set_xticks(W); d.set_xticklabels([str(w) for w in W]); d.minorticks_off(); d.set_xlabel('lens patch width (mm)'); d.set_ylabel('EQE (%)')
d.set_ylim(18, 52); d.legend(fontsize=5.2, frameon=False, loc='upper left'); d.set_title(f'(d) EQE vs lens patch width (simulation × {ETA:.2f})', loc='left')
sheets.append(('d_EQE_vs_patch', ['structure', 'patch_width_mm', 'EQE_sim_scaled_%'], rows_d)); sheets.append(('d_measured_EQE', ['structure', 'EQE_% (>=3.4 V)...'], [[k] + v for k, v in EQE34.items()]))
fig.savefig('../figures_oe/fig5_experiment.png', dpi=300, bbox_inches='tight'); fig.savefig('../figures_oe/fig5_experiment.pdf', bbox_inches='tight')
wb = Workbook(); wb.remove(wb.active)
for nm, hdr, rows in sheets:
    ws = wb.create_sheet(nm); ws.append(hdr)
    for r in rows: ws.append(r)
wb.save('../figures_oe/fig5_experiment.xlsx')
print('ETA', round(ETA, 3), 'bare sim scaled', round(bs, 1))
for nm in ('hemisphere', 'freeform'):
    print(nm, '3mm sim scaled', round(R[f'{nm}|patch3mm|rectAl']['EQE'] * 100 * ETA, 1), 'meas', round(np.mean(EQE34[nm]), 1), '±', round(np.std(EQE34[nm], ddof=1), 1), 'ratio meas/sim', round(np.mean(EQE34[nm]) / (R[f'{nm}|patch3mm|rectAl']['EQE'] * 100 * ETA), 2))
