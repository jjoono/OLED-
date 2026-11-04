"""Baseline: planar strong-microcavity NIR OLED (no scattering, flat substrate) optimised for forward radiance.
Ag 100 | ETL | EML | HTL | TE | internal DBR (SiO2/TiO2, free) | flat substrate n_sub -> air (Fresnel, specular recycling).
Usage: python3 planar_cavity.py <TE> <n_sub> <npairs> <objective> <lam0> <fwhm> <q> <Theta> <tag>"""
import sys, json, os, time, numpy as np
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize
sys.path.insert(0, '/home/user/OLED-/sim/nir')
import nir_opt as N, cps2, collim as C
te, nsub, npr, obj = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
lam0, fwhm, q, th, tag = float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8]), sys.argv[9]
N.THETA = th; Cn = C.to_air(nsub)
def one(a):
    x, lam, npts = a
    S = cps2.Stack(float(lam), (N.NORG, N.NORG), N.DEML, N.DEML / 2, above=[(N.NORG, N.NORG, x[0]), (N.AG[int(lam - 600)],) * 2 + (100.0,)],
                   below=N.stack(x, te, npr, lam), n_sub=nsub, h=th)
    r = cps2.solve_pol(S, npts=npts); P = np.real(cps2.sub_angular(S, N.TH)); R = np.real(cps2.stack_reflectance(S, N.TH))
    Rf, Tf = C.tmm_RT(nsub, [], float(lam), C.TH)
    Dl = C.air_dist(P, R, Cn * Tf[None, :], np.diag(Rf))
    return float(np.real(r['P_tot'])), float(np.real(r['air'] + r['sub'])), Dl
def evaluate(x, step, npts, mapper):
    lams = np.arange(600, 901, step); w = C.spectrum(lams, lam0, fwhm); keep = w > 1e-4 * w.max(); lams, w = lams[keep], w[keep] / w[keep].sum()
    res = list(mapper(one, [(tuple(x), l, npts) for l in lams])); D = np.zeros(90); Fb = 0.0
    for wi, (F, es, Dl) in zip(w, res): D += wi * F * es * Dl; Fb += wi * F
    return q * D / (1 - q + q * Fb)
class Obj:
    def __init__(self, lo, hi): self.lo, self.hi = lo, hi
    def __call__(self, x): return -C.metrics(evaluate(np.clip(x, self.lo, self.hi), 20, 1000, map))[obj]
if __name__ == '__main__':
    bd = N.bounds(te, npr); lo, hi = np.array([b[0] for b in bd]), np.array([b[1] for b in bd]); t0 = time.time()
    with Pool(2) as pool:
        f = lambda x, st=20, npts=1000: C.metrics(evaluate(np.clip(x, lo, hi), st, npts, map))[obj]
        r = differential_evolution(Obj(lo, hi), bd, seed=3, popsize=8, maxiter=40, tol=1e-4, polish=False, workers=pool.map, updating='deferred')
        m = minimize(lambda x: -C.metrics(evaluate(np.clip(x, lo, hi), 10, 1500, pool.map))[obj], r.x, method='Nelder-Mead',
                     options=dict(maxiter=200, xatol=0.5, fatol=1e-5))
        xb = np.clip(m.x, lo, hi); D = evaluate(xb, 2, 3000, pool.map)
    out = dict(tag=tag, TE=te, n_sub=nsub, npairs=npr, objective=obj, lam0=lam0, fwhm=fwhm, q=q, Theta=th,
               x=[round(float(v), 1) for v in xb], metrics=C.metrics(D), D=[float(v) for v in D], minutes=round((time.time() - t0) / 60, 1))
    os.makedirs('results', exist_ok=True); json.dump(out, open('results/%s.json' % tag, 'w'), indent=1)
    print(tag, {k: round(v, 4) for k, v in out['metrics'].items()}, flush=True)
