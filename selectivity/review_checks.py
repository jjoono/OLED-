"""Robustness checks requested by the technical review.  Writes data/review_checks.json and prints a summary."""
import json, numpy as np
src = open('oe_figs.py').read(); exec(src[:src.index('LIB = library()')])
LIB = library()
W90 = W1.copy(); TH1 = np.arange(90) + 0.5
OUTJ = {}


def lam2_signed(BR, cut=85.0):
    k = TH1 < cut; B = BR[np.ix_(k, k)]; M = B / np.maximum(B.sum(0, keepdims=True), 1e-12)
    ev = np.linalg.eigvals(M); ev = ev[np.argsort(-np.abs(ev))]; return ev[1]


def pinf(BR, cut=85.0):
    k = TH1 < cut; B = BR[np.ix_(k, k)]; return 1 - float(np.max(np.abs(np.linalg.eigvals(B))))


def p_on(BT, cut=85.0):
    k = TH1 < cut; pT = BT.sum(1)[k]; w = W90[k]; return float(pT @ w / w.sum())


def dev(d, i):
    return DV[f'{d}_P'][i].reshape(90, 2).sum(1), DV[f'{d}_R'][i].reshape(90, 2).mean(1)


def red(pT, R, P, m):
    return mixing.reduced(pT, W1, R, P, m)


# (a) sign of lambda2 ------------------------------------------------------------
sg = {nm: lam2_signed(v[2]) for nm, v in LIB.items() if nm != 'planar'}
OUTJ['lam2_signed'] = {k: [float(np.real(v)), float(np.imag(v))] for k, v in sg.items()}
print('(a) lambda2 real part min %.3f, max |imag| %.1e, negatives: %s' % (min(np.real(list(sg.values()))), max(abs(np.imag(list(sg.values())))), [k for k, v in sg.items() if np.real(v) < 0]))

# (b) symmetrized (detailed balance enforced) -------------------------------------
rows = []
for nm, (c, BT, BR) in LIB.items():
    if nm == 'planar': continue
    F = BR * W1[None, :]; Fs = 0.5 * (F + F.T); Bs = Fs / W1[None, :]
    # keep each column's return fraction (hence p(theta)) exactly as traced
    Bs *= (BR.sum(0) / np.maximum(Bs.sum(0), 1e-15))[None, :]
    l0, l1 = mixing.lam2_cut(BR), mixing.lam2_cut(Bs); p0, p1 = pinf(BR), pinf(Bs)
    e = []
    for d in DEVS:
        P, R = dev(d, 75); pT = BT.sum(1)
        e.append((mixing.full(pT, Bs, R, P) - mixing.full(pT, BR, R, P)) * 100)
        e.append((abs(red(pT, R, P, 1 - l1) - mixing.full(pT, Bs, R, P)) - abs(red(pT, R, P, 1 - l0) - mixing.full(pT, BR, R, P))) * 100)
    rows.append((nm, 1 / (1 - l0), 1 / (1 - l1), p0 / p_on(BT), p1 / p_on(BT), max(abs(np.array(e[0::2]))), max(abs(np.array(e[1::2])))))
a = np.array([r[1:] for r in rows])
OUTJ['symmetrized'] = dict(max_dNmix=float(np.max(abs(a[:, 1] - a[:, 0]))), max_rel_dNmix=float(np.max(abs(a[:, 1] / a[:, 0] - 1))), max_dxi=float(np.max(abs(a[:, 3] - a[:, 2]))),
                           max_deta_pp=float(a[:, 4].max()), max_derr_pp=float(a[:, 5].max()))
print('(b) symmetrized:', OUTJ['symmetrized'])

# (c) bin width and cutoff ----------------------------------------------------------
def rebin(BT, BR, P, R, g):
    n = 90 // g; w = W1[:n * g].reshape(n, g)
    pT = BT.sum(1)[:n * g].reshape(n, g); pc = (pT * w).sum(1) / w.sum(1)
    Bc = (BR[:n * g, :n * g] * W1[None, :n * g]).reshape(n, g, n, g).sum(axis=(1, 3)) / w.sum(1)[None, :]
    return pc, Bc, P[:n * g].reshape(n, g).sum(1), (R[:n * g] * W1[:n * g]).reshape(n, g).sum(1) / w.sum(1), w.sum(1)


def metrics_grid(pc, Bc, Pc, Rc, wc, cut):
    nb = len(pc); th = (np.arange(nb) + 0.5) * 90 / nb; k = th < cut
    B = Bc[np.ix_(k, k)]; M = B / np.maximum(B.sum(0, keepdims=True), 1e-12)
    ev = np.linalg.eigvals(M); ev = ev[np.argsort(-np.abs(ev))]; l2 = float(abs(ev[1]))
    rho = float(np.max(np.abs(np.linalg.eigvals(B)))); pk = float(pc[k] @ wc[k] / wc[k].sum())
    x = np.linalg.solve(np.eye(nb) - Rc[:, None] * Bc, Pc); eta = float(pc @ x)
    q = wc * (1 - pc); q /= q.sum(); Br = (1 - pc)[None, :] * (l2 * np.eye(nb) + (1 - l2) * q[:, None])
    er = float(pc @ np.linalg.solve(np.eye(nb) - Rc[:, None] * Br, Pc))
    return dict(Nmix=1 / (1 - l2) if l2 < 1 - 1e-9 else float("inf"), xi=(1 - rho) / pk, eta=eta, err_pp=(er - eta) * 100)


REP = ['planar', 'MLA hemi_total', 'scatter g0.8_t1.0', 'A_shape1.0', 'W_slope0.2', 'W_slope1.0']
conv = {}
P, R = dev('Al', 75)
for nm in REP:
    c, BT, BR = LIB[nm]; conv[nm] = {}
    for g in (1, 2, 3):
        for cut in (80.0, 85.0, 88.0, 90.0):
            if g > 1 and cut != 85.0: continue
            conv[nm][f'{g}deg_cut{int(cut)}'] = metrics_grid(*rebin(BT, BR, P, R, g), cut)
# native 0.5 deg for a scattering layer
Ts, Rs = SL['g0.8_t1.0_w1']; I = np.eye(sel.NB); Rf, Tf = sel.tmm_RT([], 550.0, sel.TH, N_SUB)
U = np.linalg.solve(I - Rs * Rf[None, :], Ts); BRh = Rs + (Ts * Rf[None, :]) @ U; Th = Tf[:, None] * U
pTh = Th.sum(0); Ph = DV['Al_P'][75]; Rh = DV['Al_R'][75]; wh = np.diff(sel.S2)
conv['scatter g0.8_t1.0']['0.5deg_cut85'] = metrics_grid(pTh, BRh, Ph, Rh, wh, 85.0)
OUTJ['convergence'] = conv
print('(c) convergence (Al 550):')
for nm, d in conv.items():
    print('  ', nm, ' | '.join(f"{k}: N {v['Nmix']:.2f} xi {v['xi']:.3f} eta {v['eta']:.4f} err {v['err_pp']:+.2f}" for k, v in d.items()))

# (d) bounds -----------------------------------------------------------------------
viol = []; marg = []
for nm, (c, BT, BR) in LIB.items():
    k = TH1 < 85; pk = BT.sum(1)[k]; pi = pinf(BR); pp = p_on(BT)
    lo, hi = float(pk.min()), pp
    marg.append((nm, lo, pi, hi)); 
    if pi < lo - 1e-9 or pi > hi + 1e-9: viol.append(nm)
OUTJ['bounds'] = dict(violations=viol, min_upper_margin=float(min(h - p for _, _, p, h in marg)), min_lower_margin=float(min(p - l for _, l, p, _ in marg)))
print('(d) bounds: violations', viol, 'min margins upper %.4f lower %.4f' % (OUTJ['bounds']['min_upper_margin'], OUTJ['bounds']['min_lower_margin']))

# (e) device-aware late-pass escape -----------------------------------------------
from scipy.stats import spearmanr
xs, xd = [], []
for nm, (c, BT, BR) in LIB.items():
    if nm == 'planar': continue
    k = TH1 < 85; B = (R[:, None] * BR)[np.ix_(k, k)]
    xd.append(1 - float(np.max(np.abs(np.linalg.eigvals(B))))); xs.append(pinf(BR))
rs = spearmanr(xs, xd).correlation
OUTJ['device_aware'] = dict(spearman=float(rs), n=len(xs))
print('(e) rank correlation p_inf(B_R) vs 1-rho(diag R B_R) (Al 550): %.3f over %d structures' % (rs, len(xs)))

# (f) stratified errors ---------------------------------------------------------------
d = json.load(open('data/oe_eval.json'))
strat = {}
for key in ('Al', 'Ag', 'LL15', '500.0', '550.0', '600.0'):
    e = [abs(v['red'] - v['eta']) * 100 for k, v in d.items() if key in k.split('|')[1:]]
    ew = [abs(v['red'] - v['eta']) * 100 for k, v in d.items() if key in k.split('|')[1:] and not k.startswith(('W_slope0.2', 'W_slope0.4'))]
    strat[key] = dict(n=len(e), mean=float(np.mean(e)), max=float(np.max(e)), mean_wo_shallow=float(np.mean(ew)), max_wo_shallow=float(np.max(ew)))
OUTJ['stratified'] = strat
print('(f) stratified:', {k: f"{v['mean']:.2f}/{v['max']:.2f}" for k, v in strat.items()})

# (g) contour templates ---------------------------------------------------------------
tmpl = {}
P0, R1 = dev('Al', 75)
for nm in ('MLA hemi_total', 'scatter g0.8_t1.0', 'W_slope1.0'):
    pT0 = LIB[nm][1].sum(1); pmean = pT0 @ (W1 / W1.sum()); out = []
    for p in (0.30, 0.38, 0.42):
        pT = np.clip(pT0 * p / pmean, 0, 1)
        etas = [red(pT, R1, P0, 1 / Nm) for Nm in (1.0, 2.0, 4.0, 8.0, 16.0)]
        out.append(dict(p=p, eta_N1=etas[0], eta_N16=etas[-1], span_pp=(max(etas) - min(etas)) * 100))
    # slope in p at N=4 for comparison
    e1 = red(np.clip(pT0 * 0.40 / pmean, 0, 1), R1, P0, 0.25); e2 = red(np.clip(pT0 * 0.44 / pmean, 0, 1), R1, P0, 0.25)
    tmpl[nm] = dict(rows=out, deta_dp_per_0p01_pp=(e2 - e1) / 4 * 100)
OUTJ['templates'] = tmpl
for nm, v in tmpl.items():
    print('(g)', nm, ' '.join(f"p={r['p']}: span N1-16 {r['span_pp']:.2f} pp" for r in v['rows']), '| d eta per 0.01 p: %.2f pp' % v['deta_dp_per_0p01_pp'])

# (h) counterfactual decomposition along the shape sweep -----------------------------
base = LIB['A_shape0.0']; pT_a = base[1].sum(1); m_a = 1 - mixing.lam2_cut(base[2])
dec = []
for s in ('0.1', '0.2', '0.4', '0.7', '1.0'):
    c, BT, BR = LIB[f'A_shape{s}']; pT_b = BT.sum(1); m_b = 1 - mixing.lam2_cut(BR)
    for dv in DEVS:
        P, R = dev(dv, 75)
        e_aa = red(pT_a, R, P, m_a); e_ab = red(pT_a, R, P, m_b); e_ba = red(pT_b, R, P, m_a); e_bb = red(pT_b, R, P, m_b)
        dec.append(dict(s=float(s), dev=dv, total=(e_bb - e_aa) * 100, mixing=(e_ab - e_aa) * 100, shape=(e_ba - e_aa) * 100, interaction=(e_bb - e_ab - e_ba + e_aa) * 100,
                        full=(mixing.full(pT_b, BR, R, P) - mixing.full(pT_a, base[2], R, P)) * 100))
OUTJ['decomposition'] = dec
for r in dec:
    if r['dev'] == 'Al' or r['s'] == 1.0:
        print('(h) s=%.1f %-4s total %+.2f = mixing %+.2f + shape %+.2f + int %+.2f  (full series %+.2f)' % (r['s'], r['dev'], r['total'], r['mixing'], r['shape'], r['interaction'], r['full']))

# (i) composition window per device / wavelength ---------------------------------------
win = {}
for key in ('Al|500.0', 'Al|550.0', 'Al|600.0', 'Ag|550.0', 'LL15|550.0'):
    S = [v['S'][2] for k, v in d.items() if k.endswith(key)]
    win[key] = [float(min(S)), float(max(S))]
OUTJ['window'] = win
print('(i) S40-60 window:', {k: f'{a:.3f}-{b:.3f}' for k, (a, b) in win.items()})

# (j) solver diagnostics ---------------------------------------------------------------
rhomax, resmax, trunc = 0, 0, 0
for nm, (c, BT, BR) in LIB.items():
    for dv in DEVS:
        for i, lam in WLS:
            P, R = dev(dv, i); A = np.eye(90) - R[:, None] * BR
            rho = float(np.max(np.abs(np.linalg.eigvals(R[:, None] * BR)))); rhomax = max(rhomax, rho)
            x = np.linalg.solve(A, P); resmax = max(resmax, float(np.abs(A @ x - P).max()))
            v = P.copy(); acc = np.zeros(90)
            for _ in range(5000):
                acc += v; v = R * (BR @ v)
                if v.sum() < 1e-14: break
            trunc = max(trunc, abs(float(BT.sum(1) @ acc) - float(BT.sum(1) @ x)))
OUTJ['solver'] = dict(rho_max=rhomax, residual_max=resmax, sum_vs_solve_max=trunc)
print('(j) solver: max rho(diag R B_R) %.4f, max residual %.1e, |pass sum - solve| max %.1e' % (rhomax, resmax, trunc))

# (l) emission above 85 deg ------------------------------------------------------------------
fr = {dv: float(dev(dv, 75)[0][85:].sum() / dev(dv, 75)[0].sum()) for dv in DEVS}
OUTJ['P_above85'] = fr; print('(l) fraction of v0 above 85 deg:', fr)
json.dump(OUTJ, open('data/review_checks.json', 'w'), indent=1, default=float)
