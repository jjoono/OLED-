"""Main-text Fig. 5: fabricated devices of two emitters versus the finite-geometry model with the goniometer geometry.
Protocol: the cavity and the internal efficiency are fixed by the bare devices of each emitter; the lens devices are
predicted without further parameters. Simulated bands span the lattice orientation (0/30 deg), scan plane
(0/45/90 deg), rotation axis (device plane, mid-substrate, lens surface) and photodiode size (5/10 mm); see gonio_az.py.
Measured values (summary statistics only):
  orange (Ir(dmppy-ph)2tmd, 260328 series): EQE at >= 3.4 V and angular factor per device, hemispherical profiles
         (negative-angle side); S4D3 excluded (intensity maximum at -40 deg, single file, EQE above the optical ceiling).
  green (Ir(ppy)2acac, 251120 series): tilt-corrected profiles (tilt fitted to the +/- asymmetry over 5-30 deg),
         angular factor per trial with the lab definition (70-80 deg taken from the simulation where not measured;
         freeform 70-80 deg excluded because of stray light, Note S13); hemispherical scan 7 excluded (outliers at 40 and
         60 deg); Lambertian-equivalent on-axis EQE per device."""
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyBboxPatch
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from openpyxl import Workbook
import sys; sys.path.insert(0, '..')
import figstyle as FS
FS.apply()
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
          'hemisphere': ([1.0, 0.958, 0.859, 0.768, 0.547, 0.450, 0.266], [0, 0.004, 0.008, 0.016, 0.005, 0.024, 0.021]),
          'freeform': ([1.0, 0.943, 0.844, 0.732, 0.609, 0.524, 0.463], [0, 0.010, 0.025, 0.045, 0.052, 0.035, 0.025])}
PE_G = {'bare': [1.028, 1.030, 0.990, 1.005, 1.029, 1.020, 1.009], 'hemisphere': [0.685, 0.690, 0.700, 0.724, 0.729, 0.733],
        'freeform': [0.819, 0.830, 0.827, 0.820, 0.870, 0.875, 0.849, 0.899, 0.899, 0.902, 0.893, 0.888]}
FRONT_G = {'bare': [22.38, 21.90], 'hemisphere': [26.89], 'freeform': [21.09, 22.40, 22.70]}          # per device
PEDEV_G = {'bare': [1.013, 1.019], 'hemisphere': [0.710], 'freeform': [0.824, 0.878, 0.894]}           # device means
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
fig = plt.figure(figsize=(FS.W, 6.0)); gs = fig.add_gridspec(3, 2, height_ratios=[0.86, 1, 1], hspace=0.5, wspace=0.36)
ax_a = fig.add_subplot(gs[0, 0]); lg = fig.add_subplot(gs[0, 1]); ax_b = fig.add_subplot(gs[1, 0]); ax_c = fig.add_subplot(gs[1, 1]); ax_d = fig.add_subplot(gs[2, 0]); ax_e = fig.add_subplot(gs[2, 1])
sheets = []; LAB = {'bare': 'Bare', 'hemisphere': 'Hemispherical MLA', 'freeform': 'Freeform MLA'}
# (a) layout
a = ax_a; a.set_aspect('equal'); a.set_xlim(-9, 34.5); a.set_ylim(-1.0, 26.4); a.axis('off')
a.add_patch(Rectangle((0, 0), 25.4, 25.4, fill=False, lw=0.9, ec='k'))
a.add_patch(Rectangle((0, 12.7), 25.4, 12.7, color='#d9e6f2', lw=0))
for xc in 25.4 - np.array([4.7, 8.6898, 12.6796, 16.6694, 20.6592]):
    a.add_patch(FancyBboxPatch((xc - 1.3919, 4.5618 - 1.3919), 2 * 1.3919, 9.9263 + 2 * 1.3919, boxstyle='round,pad=0,rounding_size=1.3919', fc='#c9c9c4', ec='#8a8a85', lw=0.4))
a.add_patch(Rectangle((1.7, 25.4 - 1.7497 - 4.57), 22.0, 4.57, fc='#c9c9c4', ec='#8a8a85', lw=0.4))
px, py = 12.7, 14.4882
a.add_patch(Rectangle((px - 5, py - 5), 10, 10, fill=False, ec='k', ls='--', lw=0.7))
a.add_patch(Rectangle((px - 1.5, py - 1.5), 3, 3, fill=False, ec=C['hemisphere'], lw=1.2))
a.add_patch(Circle((px, py), 1.0, color='#f2b134'))
kw = dict(fontsize=6.4, va='center', arrowprops=dict(arrowstyle='-', lw=0.5, shrinkA=1, shrinkB=0))
a.annotate('3 mm lens patch', (px + 1.5, py + 0.9), (26.6, 19.5), color=C['hemisphere'], ha='left', **kw)
a.annotate('Pixel, $r$ = 1 mm', (px + 0.7, py - 0.3), (26.6, 14.4), ha='left', **kw)
a.annotate('10 mm patch\n(model only)', (px + 5, py - 4.0), (26.6, 8.6), ha='left', **kw)
a.annotate('ITO', (1.0, 23.0), (-1.2, 23.0), ha='right', **kw); a.annotate('Al', (4.74, 6.0), (-1.2, 6.0), ha='right', **kw)
a.plot([0, 25.4], [-0.9, -0.9], color='k', lw=0.6); a.text(12.7, -1.5, '25.4 mm', fontsize=6.4, ha='center', va='top')
FS.label(a, 'a', x=0.0, y=1.0)
# (b, c) normalized angular intensity
rows_b = []
for axp, em, letter, title in ((ax_b, 'orange', 'b', 'Orange emitter'), (ax_c, 'green', 'c', 'Green emitter')):
    for nm in ('bare', 'hemisphere', 'freeform'):
        s = SIM[em][nm]; axp.fill_between(A10, s['pmin'], s['pmax'], color=C[nm], alpha=0.25, lw=0)
        axp.plot(A10, s['pmid'], '-', color=C[nm], lw=0.9)
        rows_b += [[em, nm, 'sim_min', *s['pmin'].tolist()], [em, nm, 'sim_max', *s['pmax'].tolist()]]
    if em == 'orange':
        PO = np.array(list(PROF_O.values())); mo, so = PO.mean(0), PO.std(0, ddof=1)
        axp.errorbar(A10, mo, yerr=so, fmt=MK[em], ms=3.4, color=C['hemisphere'], mec='white', mew=0.3, elinewidth=0.7, capsize=1.4, zorder=3)
        rows_b += [[em, 'hemisphere', 'measured mean', *mo.tolist()], [em, 'hemisphere', 'measured sd', *so.tolist()]] + [[em, 'hemisphere', f'measured {k}', *v] for k, v in PROF_O.items()]
    else:
        for nm, (m, sd) in PROF_G.items():
            n = len(m); axp.errorbar(A10[:n], m, yerr=sd, fmt=MK[em], ms=3.4, color=C[nm], mfc='white', mew=0.9, elinewidth=0.7, capsize=1.4, zorder=3)
            rows_b += [[em, nm, 'measured mean', *m], [em, nm, 'measured sd', *sd]]
    axp.set_xlim(-2.5, 83); axp.set_ylim(0, 1.08); axp.set_xticks([0, 20, 40, 60, 80]); axp.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); axp.set_xlabel('Viewing angle (°)'); axp.set_ylabel(r'Normalized intensity')
    axp.text(0.96, 0.94, title, transform=axp.transAxes, ha='right', va='top', fontsize=7); FS.label(axp, letter)
sheets.append(('bc_profiles', ['emitter', 'structure', 'series'] + [f'{x}deg' for x in A10], rows_b))
def mean_sd(ax, x, v, em, col):
    """Mean with the standard deviation (n >= 3) of the measured values; orange filled circles, green open triangles."""
    v = np.asarray(v, float); m = v.mean(); sd = v.std(ddof=1) if len(v) >= 3 else 0.0
    ax.errorbar(x, m, yerr=sd if sd else None, fmt=MK[em], ms=4.2, color=col, mfc=col if em == 'orange' else 'white', mec='white' if em == 'orange' else col,
                mew=0.3 if em == 'orange' else 0.9, ecolor=col, elinewidth=0.9, capsize=2.0, zorder=3)
    return m, sd


# (d) angular factor
d = ax_d; names = ['bare', 'freeform', 'hemisphere']; rng = np.random.default_rng(1); rows_d = []
for i, nm in enumerate(names):
    for em, dx in (('orange', -0.2), ('green', 0.2)):
        lo, hi = SIM[em][nm]['PE']; d.add_patch(Rectangle((i + dx - 0.14, lo), 0.28, max(hi - lo, 0.004), fc='#dcdcdc', ec='k', lw=0.7))
        v = (PE_O if em == 'orange' else PE_G)[nm]                       # per device (orange) or per scan (green)
        vs = v if (em == 'orange' or len(PEDEV_G[nm]) == 1) else PEDEV_G[nm]  # statistics as in Table 2: over devices, or over scans for a single device
        m_, sd_ = mean_sd(d, i + dx, vs, em, C[nm])
        if nm != 'bare': d.plot(i + dx, P10[em][nm]['PE'], 'D', ms=3.8, mfc='white', mec='k', mew=0.8, zorder=2)
        rows_d.append([em, nm, lo, hi, P10[em][nm]['PE'] if nm != 'bare' else None, m_, sd_, len(vs)] + list(v))
d.axhline(1, color='k', lw=0.6, ls=':'); d.set_xticks(range(3)); d.set_xticklabels(['Bare', 'Freeform\nMLA', 'Hemispherical\nMLA'], fontsize=7); d.set_xlim(-0.55, 2.55); d.tick_params(axis='x', which='both', top=False, length=0)
d.set_ylim(0.55, 1.3); d.set_ylabel('Angular factor'); FS.label(d, 'd')
sheets.append(('d_angular_factor', ['emitter', 'structure', 'sim3_min', 'sim3_max', 'sim10', 'measured_mean', 'measured_sd', 'n', 'individual devices or scans...'], rows_d))
# (e) lens/bare ratios of on-axis (Lambertian-equivalent) and total EQE
# Lens-glass loss (lens_loss_test.py): transmission 0.925 per crossing of the oil-coupled interface, fitted to the total EQE of the
# orange hemispherical devices; the other three device groups are predictions.  f = EQE with loss / lossless EQE; pf = lossless / lossy angular factor.
LL = json.load(open('lens_loss_test.json'))
LOSS = {('orange', nm): (LL['interface'][nm]['EQE'] / LL['lossless'][nm]['EQE'], LL['lossless'][nm]['PE'] / LL['interface'][nm]['PE']) for nm in ('hemisphere', 'freeform')}
LOSS.update({('green', nm): (LL['interface_green'][nm]['ratio'], LL['interface_green'][nm]['PE0'] / LL['interface_green'][nm]['PE']) for nm in ('hemisphere', 'freeform')})
e = ax_e; groups = [('orange', 'hemisphere'), ('green', 'hemisphere'), ('orange', 'freeform'), ('green', 'freeform')]; rows_e = []; xt = []; xl = []
for g, (em, nm) in enumerate(groups):
    mb = np.mean((EQE_O if em == 'orange' else EQE_G)['bare']); fb = np.mean((FRONT_O if em == 'orange' else FRONT_G)['bare'])
    meas_tot = np.array((EQE_O if em == 'orange' else EQE_G)[nm]) / mb; meas_ax = np.array((FRONT_O if em == 'orange' else FRONT_G)[nm]) / fb
    sb = SIM[em]['bare']; sl = SIM[em][nm]
    sim_tot = sl['EQE'] / sb['EQE']; sim_ax = (sl['front'][0] / np.mean(sb['front']), sl['front'][1] / np.mean(sb['front']))
    t10 = P10[em][nm]['EQE'] / P10[em]['bare']['EQE']; a10 = P10[em][nm]['front_EQE'] / P10[em]['bare']['front_EQE']
    f, pf = LOSS[(em, nm)]
    for j, (lab, smin, smax, s10, mv, fac) in enumerate((('Axial', *sim_ax, a10, meas_ax, f * pf), ('Total', sim_tot, sim_tot, t10, meas_tot, f))):
        x = 2.6 * g + j; xt.append(x); xl.append(lab if j == 0 else '\n' + lab); lo, hi = min(smin, smax), max(smin, smax)
        e.add_patch(Rectangle((x - 0.3, lo), 0.6, max(hi - lo, 0.008), fc='none', ec='#6a6a6a', lw=0.8, ls=(0, (2.5, 1.6))))          # lossless model
        e.add_patch(Rectangle((x - 0.3, lo * fac), 0.6, max((hi - lo) * fac, 0.008), fc='#dcdcdc', ec='k', lw=0.7))                 # with the lens-glass loss
        m_, sd_ = mean_sd(e, x, mv, em, C[nm])
        rows_e.append([em, nm, lab.lower(), lo, hi, lo * fac, hi * fac, s10, m_, sd_, len(mv)] + mv.tolist())
    if g: e.axvline(2.6 * g - 0.8, color='#b0b0b0', lw=0.6)
    e.text(2.6 * g + 0.5, 1.71, em.capitalize(), ha='center', va='top', fontsize=6.4)
e.axhline(1, color='k', lw=0.6, ls=':'); e.set_xticks(xt); e.set_xticklabels(xl, fontsize=6.2); e.set_xlim(-0.8, 2.6 * 3 + 1.8); e.set_ylim(0.6, 1.75); e.tick_params(axis='x', which='both', top=False, length=0)
e.set_ylabel('Lens device / bare device'); FS.label(e, 'e')
sheets.append(('e_ratios', ['emitter', 'structure', 'quantity', 'lossless_min', 'lossless_max', 'with_loss_min', 'with_loss_max', 'sim10_lossless (not plotted)', 'measured_mean', 'measured_sd', 'n', 'individual devices...'], rows_e))
# legend cell
lg.axis('off')
lg.legend(handles=[Patch(fc=C[k], ec='none', label=LAB[k]) for k in ('bare', 'hemisphere', 'freeform')] +
                  [Patch(fc='#b9b9b9', ec='none', alpha=0.6, label='Model range, 3 mm patch, in (b), (c)'), Patch(fc='#dcdcdc', ec='k', lw=0.7, label='Model range, 3 mm patch, in (d)'),
                   Line2D([], [], ls='none', marker='D', ms=3.8, mfc='white', mec='k', mew=0.8, label='Model, 10 mm patch, in (d)'),
                   Patch(fc='none', ec='#6a6a6a', lw=0.8, ls=(0, (2.5, 1.6)), label='Lossless model, in (e)'), Patch(fc='#dcdcdc', ec='k', lw=0.7, label='Model with lens-glass loss, in (e)'),
                   Line2D([], [], ls='none', marker='o', ms=3.6, color='#4d4d4d', label='Measured, orange emitter (mean ± s.d.)'), Line2D([], [], ls='none', marker='^', ms=3.6, mfc='white', color='#4d4d4d', mew=0.9, label='Measured, green emitter (mean ± s.d.)')],
          loc='center left', bbox_to_anchor=(-0.06, 0.5), fontsize=6.6, labelspacing=0.32, handlelength=1.1)
FS.save(fig, '../figures_oe/fig5_experiment'); plt.close(fig)
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
for r in rows_e: print('  ', r[:3], 'lossless %.2f-%.2f' % (r[3], r[4]), 'with loss %.2f-%.2f' % (r[5], r[6]), 'meas %.2f ± %.2f (n = %d)' % (r[8], r[9], r[10]))
print('loss factors (EQE, axial):', {k: (round(v[0], 3), round(v[0] * v[1], 3)) for k, v in LOSS.items()})
