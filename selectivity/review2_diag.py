"""Diagnostics requested by the second technical review.  Output: data/review2_diag.json, ../figures_oe/figS2_composition_window.*"""
import json, glob, os, numpy as np, mixing
src = open('oe_figs.py').read(); exec(src[:src.index('# ------------------------------------------------------------------ Fig. 1')])
import sys; sys.path.insert(0, '../multipass_reset'); import hftrace_gen as G
TH = np.arange(90) + .5; w = W1 / W1.sum(); out = {}
# C1: column-sum bookkeeping and detailed-balance residual (theta < 85)
k = TH < 85; rows = []
for nm, (c, BT, BR) in LIB.items():
    pT = BT.sum(1); cs = BR.sum(0); lost = 1 - pT - cs
    F = W1[k][None, :] * BR[np.ix_(k, k)]                   # F[j,i] = w_i B(i->j)
    rel = np.linalg.norm(F - F.T) / np.linalg.norm(F)
    m = (F + F.T) > 1e-4 * F.max(); e = np.abs(F - F.T)[m] / (F + F.T)[m]
    rows.append(dict(structure=nm, cls=c, lost_max=float(lost.max()), lost_mean_w=float(lost @ w), db_rel_frob=float(rel), db_entry_median=float(np.median(e)), db_entry_p90=float(np.percentile(e, 90))))
out['C1'] = rows
nonpl = [r for r in rows if r['cls'] != 'planar']
print('C1 lost max %.2e; detailed-balance rel. Frobenius %.3f-%.3f (median %.3f); entry-wise median %.3f, p90 %.3f' % (
    max(r['lost_max'] for r in rows), min(r['db_rel_frob'] for r in nonpl), max(r['db_rel_frob'] for r in nonpl), np.median([r['db_rel_frob'] for r in nonpl]),
    np.median([r['db_entry_median'] for r in nonpl]), np.median([r['db_entry_p90'] for r in nonpl])))
# C3: hemispherical return matrix: off-diagonal bands, unit cell vs 6x6 supercell vs 40k rays
def Mof(BR): return BR / np.maximum(BR.sum(0, keepdims=True), 1e-12)
Mh = Mof(LIB['MLA hemi_total'][2])
z0 = np.load('../multipass_reset/random_sweep/A_shape0.0.npz'); M6 = Mof(z0['BR'].T)
z4 = np.load('../multipass_reset/review_traces/A_shape0.0_N40000.npz'); M40 = Mof(z4['BR'].T)
def band(M):  # mass returned more than 10 deg from the incident angle, for incidence 45-80 deg
    I = np.arange(45, 80); J = np.arange(90)[:, None]
    return float(np.mean([(M[np.abs(np.arange(90) - i) > 10, i]).sum() for i in I]))
peaks = {}
for i in (50, 60, 70):
    col = M6[:, i].copy(); col[max(0, i - 10):i + 11] = 0; peaks[i] = int(np.argmax(col))
out['C3'] = dict(offdiag_mass_hemi_unitcell=band(Mh), offdiag_mass_supercell6=band(M6), offdiag_mass_supercell6_40k=band(M40),
                 rel_frob_6x6_20k_vs_40k=float(np.linalg.norm((M6 - M40)[:85, :85]) / np.linalg.norm(M40[:85, :85])),
                 off_band_peak_return_angle={str(a): b for a, b in peaks.items()})
print('C3', out['C3'])
# C4: p, sigma, p_mix over 0-80 / 0-85 / 0-90
c4 = []
for nm, (c, BT, BR) in LIB.items():
    pT = BT.sum(1); r = dict(structure=nm)
    for cut in (80, 85, 90):
        kk = TH < cut; ww = W1[kk] / W1[kk].sum(); p = float(pT[kk] @ ww); v = float(((pT[kk] - p) ** 2) @ ww)
        r[f'p_{cut}'] = p; r[f'pmix_{cut}'] = p - v / (1 - p)
    c4.append(r)
d = [max(abs(r['pmix_80'] - r['pmix_90']), abs(r['pmix_85'] - r['pmix_90'])) for r in c4]
dr = [max(abs(r['pmix_80'] / r['p_80'] - r['pmix_90'] / r['p_90']), abs(r['pmix_85'] / r['p_85'] - r['pmix_90'] / r['p_90'])) for r in c4]
out['C4'] = dict(max_abs_change_pmix=float(max(d)), max_change_ratio=float(max(dr)))
print('C4 p_mix: max change between domains %.4f (ratio p_mix/p %.4f)' % (max(d), max(dr)))
# C5: signed errors by class / device / wavelength, worst case
c5 = {}; worst = None
for (nm, dev, lam), e in EV.items():
    err = (e['red'] - e['eta']) * 100; c = LIB[nm][0]
    c5.setdefault(('class', c), []).append(err); c5.setdefault(('device', dev), []).append(err); c5.setdefault(('lambda', lam), []).append(err); c5.setdefault(('all', ''), []).append(err)
    if worst is None or abs(err) > abs(worst[0]):
        worst = (err, nm, dev, lam)
out['C5'] = {f'{a}:{b}': dict(mean_signed=float(np.mean(v)), mae=float(np.mean(np.abs(v))), max_abs=float(np.max(np.abs(v)))) for (a, b), v in c5.items()}
out['C5']['worst'] = dict(err=worst[0], structure=worst[1], device=worst[2], lam=worst[3])
for kk_, v in out['C5'].items(): print('C5', kk_, {a: (round(b, 2) if isinstance(b, float) else b) for a, b in v.items()})
# C6: xi versus cutoff
def xi_cut(BR, cut):
    kk = TH < cut; B = BR[np.ix_(kk, kk)]; rho = float(np.max(np.abs(np.linalg.eigvals(B))))
    pt = 1 - float((B.sum(0) * W1[kk]).sum() / W1[kk].sum()); return (1 - rho) / pt
c6 = {}
for nm in ('planar', 'MLA hemi_total', 'W_slope0.2', 'scatter g0.8_t4.0', 'scatter g0.0_t0.25', 'W_slope1.0', 'A_shape1.0'):
    c6[nm] = {str(cut): xi_cut(LIB[nm][2], cut) for cut in (80, 85, 88, 90)}
    print('C6', nm, {a: round(b, 3) for a, b in c6[nm].items()})
out['C6'] = c6
# C7: wavy geometry
c7 = {}
for sg in (0.2, 0.4, 0.7, 1.0, 1.5):
    amp, gx, gy, ph, L = G.wave_surface(sg, seed=0)
    c7[str(sg)] = dict(rms_height=float(np.sqrt(0.5 * np.sum(amp ** 2))), corr_length=1.0, box=L, n_modes=int(len(amp)))
print('C7', {a: round(b['rms_height'], 3) for a, b in c7.items()}); out['C7'] = c7
# C8: one-step escape autocorrelation under the stationary chain (C < 0: failed light sent to high-escape angles)
c8 = []
for nm, (c, BT, BR) in LIB.items():
    if c == 'planar':
        continue
    pT = BT.sum(1); q = W1 * (1 - pT); q /= q.sum(); pm = float(q @ pT); f = pT - pm; M = Mof(BR)
    C = float((q * f * (M.T @ f)).sum() / (q * f * f).sum())
    e = EV[(nm, 'Al', 550.0)]; c8.append(dict(structure=nm, cls=c, C=C, gap_pp=(e['eta'] - e['m1']) * 100))
cc = np.corrcoef([r['C'] for r in c8], [r['gap_pp'] for r in c8])[0, 1]
for r in sorted(c8, key=lambda r: r['C'])[:4]: print('C8', r['structure'], 'C %.3f  eta-eta_mix %+.2f pp' % (r['C'], r['gap_pp']))
print('C8 correlation(C, eta - eta_mix) = %.2f' % cc); out['C8'] = dict(rows=c8, corr=float(cc))
# C9: S40-60 uncertainty (ray seeds, realizations) and small multiples
c9 = {}
for f in sorted(glob.glob('../multipass_reset/random_sweep_unc/*.npz')):
    nm = os.path.basename(f)[:-4]; base, var = nm.split('__'); z = np.load(f)
    P = DV['Al_P'][75].reshape(90, 2).sum(1); R = DV['Al_R'][75].reshape(90, 2).mean(1)
    air = np.linalg.solve(np.eye(90) - R[:, None] * z['BR'].T, P) @ z['BT']
    c9.setdefault(base, {'lib': EV[(base, 'Al', 550.0)]['S'][2]})[var] = float(air[40:60].sum() / air.sum())
for b, v in c9.items():
    ray = [v['lib'], v.get('ray2'), v.get('ray3')]; real = [v['lib'], v.get('real1'), v.get('real2')]
    v['sd_ray'] = float(np.std([x for x in ray if x is not None], ddof=1)); v['sd_real'] = float(np.std([x for x in real if x is not None], ddof=1))
    print('C9', b, 'S40-60 %.4f  sd(ray) %.4f  sd(realization) %.4f' % (v['lib'], v['sd_ray'], v['sd_real']))
out['C9'] = c9
fig, axs = plt.subplots(3, 3, figsize=(7.2, 6.4), sharex=False, sharey=True, gridspec_kw=dict(hspace=0.55, wspace=0.12)); rows9 = []
for a, dev in enumerate(DEVS):
    for b, (i, lam) in enumerate(WLS):
        ax = axs[a, b]
        for nm, (c, BT, BR) in LIB.items():
            e = EV[(nm, dev, lam)]; col, mk, _ = CLS[c]; ax.plot(e['eta'], e['S'][2], mk, color=col, ms=3, mew=0); rows9.append([dev, lam, nm, c, e['eta'], e['S'][2]])
        ax.axhline(0.337, color=INK2, ls=':', lw=0.6); ax.set_title(f"{ {'Al': 'Al', 'Ag': 'Ag', 'LL15': 'low loss'}[dev]}, {int(lam)} nm", loc='left', fontsize=6.5)
        if b == 0: ax.set_ylabel('S₄₀₋₆₀')
        if a == 2: ax.set_xlabel('η_sta')
for c, (col, mk, lab) in CLS.items():
    axs[0, 0].plot([], [], mk, color=col, label=lab)
axs[0, 0].legend(fontsize=4.6, frameon=False, loc='lower left')
fig.savefig(f'{OUT}/figS2_composition_window.png', dpi=300, bbox_inches='tight'); fig.savefig(f'{OUT}/figS2_composition_window.pdf', bbox_inches='tight')
save_xlsx('figS2_composition_window', [('points', ['device', 'lambda_nm', 'structure', 'class', 'eta_sta', 'S40_60'], rows9)])
json.dump(out, open('data/review2_diag.json', 'w'), indent=1, default=float)
