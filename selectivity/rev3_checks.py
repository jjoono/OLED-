import numpy as np, json
src = open('oe_figs.py').read(); exec(src[:src.index('# ------------------------------------------------------------------ Fig. 1')])
from scipy.stats import pearsonr, spearmanr
# item 9: correlation of Fig. 2(b)
PLC = json.load(open('data/placement.json')); xs, ys, cl = [], [], []
for nm, (c, BT, BR) in LIB.items():
    v = PLC['structures'][f'{nm}|Al|550']; xs.append(v['p'] - v['p_mix']); ys.append((v['dev_mix'] - v['exit_mix']) * 100); cl.append(c)
xs, ys = np.array(xs), np.array(ys); m = np.array([c != 'planar' for c in cl])
print('Fig2b N', len(xs), 'Pearson %.3f, without planar %.3f (N %d), Spearman %.3f, without planar %.3f' % (pearsonr(xs, ys)[0], pearsonr(xs[m], ys[m])[0], m.sum(), spearmanr(xs, ys)[0], spearmanr(xs[m], ys[m])[0]))
# item 8: bin width for p, p_mix and eta (planar at native 0.5 deg; others merged from 1 deg)
def pm(pT, w):
    wn = w / w.sum(); p = pT @ wn; v = ((pT - p) ** 2) @ wn; return p, p - v / (1 - p)
def merge(BT, BR, w, k):
    n = len(w) // k; W = w.reshape(n, k); f = (W / W.sum(1, keepdims=True)).ravel()
    BTm = (BT * f[:, None]).reshape(n, k, -1).sum(1); BRm = (BR * f[None, :]).reshape(BR.shape[0] // k, k, n, k).sum(axis=(1, 3))
    return BTm, BRm, W.sum(1)
res = {}
P05 = DV['Al_P'][75]; R05 = DV['Al_R'][75]
for nm in ('planar', 'MLA hemi_total', 'W_slope0.2', 'A_shape1.0'):
    c, BT, BR = LIB[nm]; rows = []
    for k in (1, 2, 3, 5):
        BTm, BRm, wm = merge(BT, BR, W1, k)
        P = DV['Al_P'][75].reshape(90 // k, 2 * k).sum(1); R = (DV['Al_R'][75].reshape(90 // k, 2 * k) * sel.WL.reshape(90 // k, 2 * k)).sum(1) / sel.WL.reshape(90 // k, 2 * k).sum(1)
        p, pmx = pm(BTm.sum(1), wm); eta = mixing.full(BTm.sum(1), BRm, R, P)
        rows.append((k, round(float(p), 4), round(float(pmx), 4), round(eta, 4)))
    res[nm] = rows; print(nm, rows)
# planar at native 0.5 deg (exact Fresnel)
Rf, Tf = sel.tmm_RT([], 550.0, sel.TH, N_SUB); w05 = sel.WL
p, pmx = pm(Tf, w05); print('planar 0.5 deg native: p %.4f p_mix %.4f' % (p, pmx))
json.dump(dict(fig2b=dict(pearson=pearsonr(xs, ys)[0], pearson_noplanar=pearsonr(xs[m], ys[m])[0], spearman=spearmanr(xs, ys)[0], spearman_noplanar=spearmanr(xs[m], ys[m])[0], N=len(xs)),
               bins=res, planar05=dict(p=float(p), p_mix=float(pmx))), open('data/rev3_checks.json', 'w'), indent=1)
