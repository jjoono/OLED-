"""Cross-seeded refinement of lensfree_limit.csv: every case is re-evaluated (fine spectrum) at the optima of all cases
that share its bounds (same metal / order), the three best seeds are polished with Nelder-Mead, and the best result is kept.
Removes local optima that the per-case differential evolution can fall into. -> lensfree_refined.csv"""
import csv, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
import lensfree_limit as L

FIX = lambda v: v if v in ('Ir', 'bulk', 'joint8', 'free', 'order1') else float(v)
rows = [{k: FIX(v) for k, v in r.items()} for r in csv.DictReader(open('lensfree_limit.csv'))]
for r in rows:
    if r['fwhm'] != 'Ir': r['fwhm'] = int(r['fwhm']) if r['fwhm'] == int(r['fwhm']) else r['fwhm']
    r['lam0'] = int(r['lam0'])
KEYS = ('theta', 'fwhm', 'lam0', 'metal', 'n_t', 'n_e', 'n_c', 'order')
SEEDS = [np.array([r['dE'], r['dH'], r['dAg'], r['dC']]) for r in rows] + [np.array([60., 70., 25., 70.])]


def refine(r):
    c = {k: r[k] for k in KEYS}
    dmax = 150. if c['order'] == 'order1' else 400.
    agb = (7.87, 7.87) if c['metal'] == 'joint8' else (5., 60.)
    bnds = [(15., dmax), (15., dmax), agb, (0., 250.)]
    lo, hi = np.array([b[0] for b in bnds]), np.array([b[1] for b in bnds])
    l, w = L.spec(c)
    cand = []
    for s in SEEDS:
        s = np.clip(s, lo, hi)
        cand.append((L.eqe(s, c, l, w), tuple(s)))
    cand = sorted(set(cand), reverse=True)[:3]
    x0 = np.array([r['dE'], r['dH'], r['dAg'], r['dC']]); best, xb = r['EQE'], x0
    for _, s in cand:
        rr = minimize(lambda x: -L.eqe(x, c, l, w), np.array(s), method='Nelder-Mead', bounds=bnds,
                      options=dict(maxiter=300, xatol=0.2, fatol=1e-5))
        if -rr.fun > best + 1e-4: best, xb = -rr.fun, rr.x
    out = dict(r, EQE_dE=r['EQE'], EQE=best, dE=xb[0], dH=xb[1], dAg=xb[2], dC=xb[3])
    lam = c['lam0'] if c['fwhm'] != 'Ir' else 615
    i = L.T.idx(lam); ab = L.T.AG[i]; at = L.AGJ[i] if c['metal'] == 'joint8' else ab
    St = L.cps2.Stack(lam, (c['n_e'],) * 2, 2.0, 1.0, above=[(c['n_t'],) * 2 + (xb[0],), (ab, ab, 100.)],
                      below=[(c['n_t'],) * 2 + (xb[1],), (at, at, xb[2])] + ([(c['n_c'],) * 2 + (xb[3],)] if xb[3] > 0.5 else []),
                      n_sub=1.0, n_top=1.0, h=c['theta'])
    b = L.cps2.solve(St, npts=3000)
    out.update(air_lam0=float(np.real(b['air'])), wg=float(np.real(b['wg'])), spp_evan=float(np.real(b['spp'])), abs_cone=float(np.real(b['abs'])))
    print(c, 'DE %.3f -> refined %.3f' % (r['EQE'], best), flush=True)
    return out


if __name__ == '__main__':
    with Pool(4) as p:
        res = p.map(refine, rows, chunksize=1)
    with open('lensfree_refined.csv', 'w', newline='') as f:
        wr = csv.DictWriter(f, list(res[0].keys())); wr.writeheader(); wr.writerows(res)
