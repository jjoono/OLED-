# -*- coding: utf-8 -*-
"""EQE against aspect ratio from the same ray-traced BSDFs, evaluated two ways: the matrix
series of eq. (1) (what the device predictions use) and the closed form of eq. (3) with the
single-pass escape probability p.  The closed form assumes the returned light is angularly
randomised at every pass; the hemisphere's returned light is not, so the two part company."""
import numpy as np, matplotlib, os, sys, csv
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/home/user/OLED-/sim/fig3c'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cps2, materials as M, series

HERE = os.path.dirname(os.path.abspath(__file__))
AL = 0.958336812040867 + 6.68678039868499j
Z = np.load(os.path.join(HERE, 'mla_bsdf.npz')); ar = np.asarray(Z['aspect'], float); TH = Z['theta_deg']
W = np.cos(np.radians(TH)) * np.sin(np.radians(TH))
rows = []
curves = {}
for name, refl, c in (('Ag', M.AG, '#C00000'), ('Al', AL, '#1A1A1A')):
    S = cps2.Stack(550.0, (1.8, 1.8), 20.0, 10.0, above=[(1.8, 1.8, 200.0), (refl, refl, 100.0)],
                   below=[(1.8, 1.8, 200.0), (M.ITO, M.ITO, 50.0)], n_sub=1.8)
    r = cps2.solve_pol(S, npts=12000); es = r['air'] + r['sub']
    R = cps2.stack_reflectance(S, TH); P = cps2.sub_angular(S, TH); A = 1 - np.sum(R * W) / np.sum(W)
    xs, ser, cf1, cf2 = [], [], [], []
    for k in np.argsort(ar):
        if ar[k] < 0.09: continue
        BT, BR = Z['BT'][k], Z['BR'][k]
        pc = np.sum(BT * W) / np.sum(W); pp = BT @ (P / P.sum())
        e, _ = series.eta_ext(BT, BR, R, P, n_term=80)
        xs.append(ar[k]); ser.append(es * e); cf1.append(es * pc / (pc + (1 - pc) * A)); cf2.append(es * pp / (pp + (1 - pp) * A))
        rows.append([name, ar[k], pc, pp, e, es * e, es * pc / (pc + (1 - pc) * A), es * pp / (pp + (1 - pp) * A)])
    curves[name] = (np.array(xs), np.array(ser), np.array(cf1), np.array(cf2), c)
with open(os.path.join(HERE, 'mla_aspect_series_vs_closedform.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['reflector', 'aspect_ratio', 'p_cos_sin', 'p_Psub_weighted', 'eta_ext_series', 'EQE_series', 'EQE_closedform_p_cos_sin', 'EQE_closedform_p_Psub'])
    for r_ in rows: w.writerow([r_[0], '%.2f' % r_[1]] + ['%.6f' % v for v in r_[2:]])

FS = 7.4
plt.rcParams.update({'font.size': FS, 'axes.linewidth': 0.8})
fig, ax = plt.subplots(figsize=(3.8, 2.7)); fig.subplots_adjust(left=0.14, right=0.97, top=0.9, bottom=0.17)
for name, (x, s_, c1, c2, c) in curves.items():
    ax.plot(x, 100 * s_, color=c, lw=1.7, label='%s — series, eq. (1)' % name)
    ax.plot(x, 100 * c2, color=c, lw=1.2, ls='--', label='%s — closed form, eq. (3), p weighted by P_sub' % name)
    ax.plot(x, 100 * c1, color=c, lw=1.0, ls=':', label='%s — closed form, p flux-weighted' % name)
ax.axvline(1.0, color='0.8', lw=0.7, dashes=(2.5, 2.5), zorder=0)
ax.set_xlim(0.05, 1.55); ax.set_ylim(30, 95); ax.set_xlabel('aspect ratio  $h/r$'); ax.set_ylabel('EQE (%)')
ax.set_title('same BSDFs, close-packed lenses, $n_{\\rm sub}$ = $n_{\\rm MLA}$ = 1.8', fontsize=FS - 0.6, loc='left')
ax.legend(loc='lower right', frameon=False, fontsize=FS - 1.6)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
fig.savefig(os.path.join(HERE, 'mla_aspect_series_vs_closedform.png'), dpi=300)
print('written')
