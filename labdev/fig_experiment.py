"""Main-text Fig. 5: fabricated devices of two emitters versus the finite-geometry model with the goniometer geometry.
Protocol: the cavity and the internal efficiency are fixed by the bare devices of each emitter; the lens devices are
predicted without further parameters. Simulated bands span the lattice orientation (0/30 deg), scan plane
(0/45/90 deg), rotation axis (device plane, mid-substrate, lens surface) and photodiode size (5/10 mm); see gonio_az.py.
Measured values (summary statistics only):
  orange (Ir(dmppy-ph)2tmd, 260328 series): EQE at >= 3.4 V and angular factor per device, hemispherical profiles
         (negative-angle side); S4D3 excluded (intensity maximum at -40 deg, single file, EQE above the optical ceiling).
  green (Ir(ppy)2acac, 251120 series): tilt-corrected profiles (tilt fitted to the +/- asymmetry over 5-30 deg),
         angular factor per trial with the lab definition (70-80 deg taken from the simulation where not measured;
         freeform 70-80 deg excluded because of stray light, Note S13); Lambertian-equivalent on-axis EQE per device."""
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyBboxPatch
from openpyxl import Workbook
plt.rcParams.update({'font.size': 7, 'axes.titlesize': 7.5, 'axes.spines.top': False, 'axes.spines.right': False})
C = {'bare': '#7d7c76', 'hemisphere': '#2a78d6', 'freeform': '#eb6834'}; INK2 = '#52514e'
MK = {'orange': 'o', 'green': '^'}
A10 = np.arange(0, 90, 10)
# ---------------- measured: orange emitter
EQE_O = {'bare': [29.7, 27.7, 30.7, 31.4, 30.8], 'freeform': [22.9, 26.7, 23.7, 28.7], 'hemisphere': [22.9, 27.1, 22.5, 25.3, 24.1]}
PE_O = {'bare': [1.01, 1.08, 0.95, 1.02, 0.97], 'freeform': [1.13, 0.99, 1.11, 1.03], 'hemisphere': [0.63, 0.69, 0.63, 0.71, 0.66]}
PROF_O = {'S3D2': [1.00, 0.97, 0.86, 0.71, 0.49, 0.32, 0.28, 0.13, 0.08], 'S2D4': [1.00, 0.96, 0.84, 0.69, 0.52, 0.37, 0.27, 0.20, 0.13],
          'S1D4': [1.00, 0.97, 0.85, 0.65, 0.44, 0.29, 0.21, 0.17, 0.15], 'S4D2': [1.00, 0.95, 0.84, 0.69, 0.53, 0.42, 0.33, 0.22, 0.10],
          'S5D4': [1.00, 0.96, 0.83, 0.64, 0.44, 0.32, 0.26, 0.23, 0.15]}
# ---------------- measured: green emitter (tilt-corrected; mean and s.d. over trials)
PROF_G = {'bare': ([1.0, 0.980, 0.931, 0.854, 0.768, 0.679, 0.554, 0.339, 0.162], [0, 0.005, 0.006, 0.006, 0.015, 0.024, 0.026, 0.032, 0.039]),
          'hemisphere': ([1.0, 0.954, 0.866, 0.754, 0.544, 0.428, 0.247], [0, 0.002, 0.005, 0.010, 0.001, 0.015, 0.005]),
          'freeform': ([1.0, 0.943, 0.844, 0.732, 0.609, 0.524, 0.463], [0, 0.010, 0.025, 0.045, 0.052, 0.035, 0.025])}
PE_G = {'bare': [1.028, 1.030, 0.990, 1.005, 1.029, 1.020, 1.009], 'hemisphere': [0.685, 0.690, 0.700],
        'freeform': [0.819, 0.830, 0.827, 0.820, 0.870, 0.875, 0.849, 0.899, 0.899, 0.902, 0.893, 0.888]}
FRONT_G = {'bare': [22.38, 21.90], 'hemisphere': [26.89], 'freeform': [21.09, 22.40, 22.70]}          # per device
PEDEV_G = {'bare': [1.013, 1.019], 'hemisphere': [0.692], 'freeform': [0.824, 0.878, 0.894]}           # device means
EQE_G = {k: [f * p for f, p in zip(FRONT_G[k], PEDEV_G[k])] for k in FRONT_G}
FRONT_O = {k: [e / p for e, p in zip(EQE_O[k], PE_O[k])] for k in EQE_O}
# ---------------- simulated
GZ = {'orange': json.load(open('gonio_az_orange.json')), 'green': json.load(open('gonio_az_e65h80.json'))}
P10 = {'orange': {nm: json.load(open('layout_results.json'))[f'{nm}|patch{"0" if nm == "bare" else "10"}mm|rectAl'] for nm in ('bare', 'hemisphere', 'freeform')},
       'green': {nm: (lambda d: dict(EQE=d['EQE'], front_EQE=d['front'], PE=d['PE']))(json.load(open('green_patch.json'))[f'{nm}|{"0" if nm == "bare" else "10"}mm']) for nm in ('bare', 'hemisphere', 'freeform')}}


def sim(em, nm):
    """EQE and the ranges of the lab-defined profile and angular factor over the measurement geometries."""
    ks = [k for k in GZ[em] if k.startswith(nm + '|')]; S = [sm for k in ks for sm in GZ[em][k]['scans'].values()]
    eqe = float(np.mean([GZ[em][k]['EQE'] for k in ks])); prof = np.array([s['sym'] for s in S]); pe = np.array([s['PE'] for s in S])
    return dict(EQE=eqe * 100, pmin=prof.min(0), pmax=prof.max(0), pmid=0.5 * (prof.min(0) + prof.max(0)), PE=(pe.min(), pe.max()),
                front=(eqe * 100 / pe.max(), eqe * 100 / pe.min()))


SIM = {em: {nm: sim(em, nm) for nm in ('bare', 'hemisphere', 'freeform')} for em in GZ}
ETA = {'orange': np.mean(EQE_O['bare']) / SIM['orange']['bare']['EQE'], 'green': np.mean(EQE_G['bare']) / SIM['green']['bare']['EQE']}
fig = plt.figure(figsize=(7.1, 5.6)); gs = fig.add_gridspec(2, 6, hspace=0.5, wspace=1.1)
ax_a = fig.add_subplot(gs[0, 0:2]); ax_b = fig.add_subplot(gs[0, 2:4]); ax_c = fig.add_subplot(gs[0, 4:6]); ax_d = fig.add_subplot(gs[1, 0:3]); ax_e = fig.add_subplot(gs[1, 3:6])
sheets = []
# (a) layout
a = ax_a; a.set_aspect('equal'); a.set_xlim(-0.5, 26); a.set_ylim(-3.5, 26); a.axis('off')
a.add_patch(Rectangle((0, 0), 25.4, 25.4, fill=False, lw=0.8, ec='k'))
a.add_patch(Rectangle((0, 12.7), 25.4, 12.7, color='#d9e6f2', lw=0))
for xc in 25.4 - np.array([4.7, 8.6898, 12.6796, 16.6694, 20.6592]):
    a.add_patch(FancyBboxPatch((xc - 1.3919, 4.5618 - 1.3919), 2 * 1.3919, 9.9263 + 2 * 1.3919, boxstyle='round,pad=0,rounding_size=1.3919', fc='#c9c9c4', ec='#8a8a85', lw=0.4))
a.add_patch(Rectangle((1.7, 25.4 - 1.7497 - 4.57), 22.0, 4.57, fc='#c9c9c4', ec='#8a8a85', lw=0.4))
px, py = 12.7, 14.4882
a.add_patch(Rectangle((px - 5, py - 5), 10, 10, fill=False, ec='k', ls='--', lw=0.6))
a.add_patch(Rectangle((px - 1.5, py - 1.5), 3, 3, fill=False, ec=C['hemisphere'], lw=1.1))
a.add_patch(Circle((px, py), 1.0, color='#f2b134'))
a.text(px + 5.3, py - 6.4, '10 mm\n(model)', fontsize=5, color='k'); a.text(px + 1.8, py + 1.6, '3 mm lens\npatch', fontsize=5, color=C['hemisphere'])
a.annotate('pixel\nr = 1 mm', (px - 0.7, py - 0.7), (px - 12.3, py - 4.2), fontsize=5, arrowprops=dict(arrowstyle='-', lw=0.5)); a.text(0.6, 23.6, 'ITO', fontsize=5, color='#4a6a8a'); a.text(17.0, 1.3, 'Al fingers', fontsize=5, color='#55554f')
a.text(0.0, -0.6, 'glass 25.4 mm; bare 0.7 mm,\nMLA 0.7 + 0.7 mm (oil-coupled)', fontsize=4.8, va='top')
a.set_title('(a) layout', loc='left')
# (b, c) normalized angular intensity
rows_b = []
for axp, em, title in ((ax_b, 'orange', '(b) orange emitter'), (ax_c, 'green', '(c) green emitter')):
    shown = ('bare', 'hemisphere', 'freeform')
    for nm in shown:
        s = SIM[em][nm]; axp.fill_between(A10, s['pmin'], s['pmax'], color=C[nm], alpha=0.22, lw=0)
        axp.plot(A10, s['pmid'], '-', color=C[nm], lw=0.8)
        rows_b += [[em, nm, 'sim_min', *s['pmin'].tolist()], [em, nm, 'sim_max', *s['pmax'].tolist()]]
    if em == 'orange':
        for k, v in PROF_O.items(): axp.plot(A10, v, MK[em], ms=3, mfc='none', color=C['hemisphere'], mew=0.7)
        rows_b += [[em, 'hemisphere', f'measured {k}', *v] for k, v in PROF_O.items()]
    else:
        for nm, (m, sd) in PROF_G.items():
            n = len(m); axp.errorbar(A10[:n], m, yerr=sd, fmt=MK[em], ms=3, color=C[nm], mfc=C[nm] if nm != 'freeform' else 'white', mew=0.7, elinewidth=0.6, capsize=1.2)
            rows_b += [[em, nm, 'measured mean', *m], [em, nm, 'measured sd', *sd]]
    axp.set_xlim(0, 82); axp.set_ylim(0, 1.08); axp.set_xticks([0, 20, 40, 60, 80]); axp.set_xlabel('angle (°)'); axp.set_ylabel('I(θ) / I(0°)'); axp.set_title(title, loc='left')
for nm, lab in (('bare', 'bare'), ('hemisphere', 'hemispherical'), ('freeform', 'freeform')):
    ax_b.fill_between([], [], [], color=C[nm], alpha=0.35, lw=0, label=lab)
ax_b.legend(fontsize=5, frameon=False, loc='lower left', title='model (3 mm patch)', title_fontsize=5)
ax_c.plot([], [], '^', color='#777', ms=3, label='measured (mean ± s.d.)'); ax_c.plot([], [], 'o', mfc='none', color='#777', ms=3, label='measured, orange (b)')
ax_c.legend(fontsize=5, frameon=False, loc='lower left')
sheets.append(('bc_profiles', ['emitter', 'structure', 'series'] + [f'{x}deg' for x in A10], rows_b))
# (d) angular factor
d = ax_d; names = ['bare', 'freeform', 'hemisphere']; rng = np.random.default_rng(1); rows_d = []
for i, nm in enumerate(names):
    for em, dx in (('orange', -0.2), ('green', 0.2)):
        lo, hi = SIM[em][nm]['PE']; d.add_patch(Rectangle((i + dx - 0.13, lo), 0.26, max(hi - lo, 0.004), fc='#d8d7d2', ec='k', lw=0.7))
        v = (PE_O if em == 'orange' else PE_G)[nm]
        d.plot(i + dx + rng.uniform(-0.07, 0.07, len(v)), v, MK[em], ms=3.2, color=C[nm], mfc=C[nm] if em == 'orange' else 'white', mew=0.8, alpha=0.9)
        if nm != 'bare': d.plot(i + dx, P10[em][nm]['PE'], 'D', ms=3.6, mfc='white', mec='k', mew=0.7)
        rows_d.append([em, nm, lo, hi, P10[em][nm]['PE'] if nm != 'bare' else None] + list(v))
d.add_patch(Rectangle((0, 0), 0, 0, fc='#d8d7d2', ec='k', lw=0.7, label='model, 3 mm patch (geometry range)'))
d.plot([], [], 'D', ms=3.6, mfc='white', mec='k', label='model, 10 mm patch'); d.plot([], [], 'o', color='#777', ms=3, label='measured, orange'); d.plot([], [], '^', mfc='white', color='#777', ms=3, label='measured, green')
d.axhline(1, color=INK2, lw=0.4, ls=':'); d.set_xticks(range(3)); d.set_xticklabels(names[:2] + ['hemispherical']); d.set_xlim(-0.55, 2.55)
d.set_ylim(0.55, 1.25); d.set_ylabel('angular factor  EQE / EQE$_{0°}$'); d.legend(fontsize=5, frameon=False, loc='upper right', ncol=2)
d.set_title('(d) angular factor', loc='left'); sheets.append(('d_angular_factor', ['emitter', 'structure', 'sim3_min', 'sim3_max', 'sim10', 'measured...'], rows_d))
# (e) lens/bare ratios of on-axis (Lambertian-equivalent) and total EQE
e = ax_e; groups = [('orange', 'hemisphere'), ('green', 'hemisphere'), ('orange', 'freeform')]; rows_e = []; xt = []; xl = []
for g, (em, nm) in enumerate(groups):
    mb = np.mean((EQE_O if em == 'orange' else EQE_G)['bare']); fb = np.mean((FRONT_O if em == 'orange' else FRONT_G)['bare'])
    meas_tot = np.array((EQE_O if em == 'orange' else EQE_G)[nm]) / mb; meas_ax = np.array((FRONT_O if em == 'orange' else FRONT_G)[nm]) / fb
    sb = SIM[em]['bare']; sl = SIM[em][nm]
    sim_tot = sl['EQE'] / sb['EQE']; sim_ax = (sl['front'][0] / np.mean(sb['front']), sl['front'][1] / np.mean(sb['front']))
    t10 = P10[em][nm]['EQE'] / P10[em]['bare']['EQE']; a10 = P10[em][nm]['front_EQE'] / P10[em]['bare']['front_EQE']
    for j, (lab, smin, smax, s10, mv) in enumerate((('on-axis', *sim_ax, a10, meas_ax), ('total', sim_tot, sim_tot, t10, meas_tot))):
        x = 3 * g + j; xt.append(x); xl.append(lab)
        e.add_patch(Rectangle((x - 0.3, min(smin, smax)), 0.6, max(abs(smax - smin), 0.006), fc='#d8d7d2', ec='k', lw=0.7))
        e.plot(x, s10, 'D', ms=3.6, mfc='white', mec='k', mew=0.7)
        e.plot(x + rng.uniform(-0.12, 0.12, len(mv)), mv, MK[em], ms=3.2, color=C[nm], mfc=C[nm] if em == 'orange' else 'white', mew=0.8)
        rows_e.append([em, nm, lab, smin, smax, s10] + mv.tolist())
    e.text(3 * g + 0.5, 1.93, f"{nm[:4] + '.' if nm == 'hemisphere' else 'freeform'}\n{em}", ha='center', va='top', fontsize=5.2, color=INK2)
e.axhline(1, color=INK2, lw=0.5, ls=':'); e.set_xticks(xt); e.set_xticklabels(xl, fontsize=5.5); e.set_xlim(-0.7, 7.7); e.set_ylim(0.6, 1.95)
e.set_ylabel('lens device / bare'); e.set_title('(e) on-axis versus total EQE', loc='left')
sheets.append(('e_ratios', ['emitter', 'structure', 'quantity', 'sim3_min', 'sim3_max', 'sim10', 'measured...'], rows_e))
fig.savefig('../figures_oe/fig5_experiment.png', dpi=300, bbox_inches='tight'); fig.savefig('../figures_oe/fig5_experiment.pdf', bbox_inches='tight')
wb = Workbook(); wb.remove(wb.active)
for nm, hdr, rows in sheets:
    ws = wb.create_sheet(nm); ws.append(hdr)
    for r in rows: ws.append([float(x) if isinstance(x, (np.floating,)) else x for x in r])
wb.save('../figures_oe/fig5_experiment.xlsx')
# ---------------- numbers quoted in the text
for em in ('orange', 'green'):
    print(f'== {em}: eta {ETA[em]:.3f}')
    for nm in ('bare', 'hemisphere', 'freeform'):
        s = SIM[em][nm]; M = (EQE_O if em == 'orange' else EQE_G)[nm]; P = (PE_O if em == 'orange' else PE_G)[nm]
        print(f"  {nm:10s} sim EQE {s['EQE']*ETA[em]:.1f} (unscaled {s['EQE']:.1f})  meas {np.mean(M):.1f} ± {np.std(M, ddof=1) if len(M) > 1 else 0:.1f}  ratio {np.mean(M)/(s['EQE']*ETA[em]):.2f} | "
              f"PE sim {s['PE'][0]:.2f}-{s['PE'][1]:.2f} meas {min(P):.2f}-{max(P):.2f} | 10 mm EQE {P10[em][nm]['EQE']*100*ETA[em]:.1f} PE {P10[em][nm]['PE']:.2f}")
for r in rows_e: print('  ', r[:3], 'sim', round(r[3], 2), round(r[4], 2), '10mm', round(r[5], 2), 'meas', np.round(np.mean(r[6:]), 2))
