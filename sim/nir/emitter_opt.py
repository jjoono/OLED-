"""Optimise fabricable NIR stacks for literature emitters (peak, FWHM, PLQY q, horizontal dipole ratio Theta).
Same stack/BSDF model as nir_opt.py; TE limited to measured n,k (ITO, IZO, Ag on HATCN).
dbr0 cases: DE + Nelder-Mead.  dbr4 cases: every refined q = 0.3 structure (results_refined/*dbr4*) as a seed
(own TE thickness), best two polished with Nelder-Mead.  -> results_emitter/<emitter>_<TE>_n<n>_dbr<nd>_<kind>.json
Usage: python3 emitter_opt.py <emitter>"""
import sys, json, os, glob, numpy as np
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize
import nir_opt as N

EMITTERS = {   # peak nm, FWHM nm, PLQY, Theta (horizontal fraction)
    'Pt-fprpz2': (740, 80, 0.80, 0.78),     # neat film; Ly et al. Nat. Photon. 2017 (PLQY 0.81), Theta 0.78 (Nat. Commun. 2024)
    'PyCNTPA-12': (747, 100, 0.786, 0.80),  # 12 wt% in CBP, TADF; Theta 'high' (not given) -> 0.80 assumed
    'PyCNTPA-24': (761, 100, 0.467, 0.80),  # 24 wt% in CBP
    'Ir-tmd': (760, 80, 0.032, 0.80),       # Ir(mCPDTiq)2tmd, Sci. Rep. 2023 (PLQY 3.2 %, Theta 0.80)
}
CASES = [('ITO', 1.5, 0, 'flat'), ('ITO', 1.8, 0, 'mla'), ('ITO', 1.8, 4, 'mla'), ('IZO', 1.8, 0, 'mla'), ('IZO', 1.8, 4, 'mla'),
         ('Ag', 1.5, 4, 'mla'), ('Ag', 1.8, 4, 'mla')]

class Obj:                                  # picklable DE objective
    def __init__(self, *a): self.a = a
    def __call__(self, x):
        te, nd, nsub, kind, q = self.a
        return -N.evaluate(x, te, nd, nsub, step=10, npts=1500, kind=kind, qs=(q,))['EQE_q%s' % q]

def set_emitter(lam0, fwhm, theta):
    N.THETA = theta
    def spectrum(step):
        l = np.arange(max(600, lam0 - 2 * fwhm), min(900, lam0 + 2 * fwhm) + 1, step)
        w = np.exp(-0.5 * ((l - lam0) / (fwhm / 2.3548)) ** 2); return l, w / w.sum()
    N.spectrum = spectrum

if __name__ == '__main__':
    em = sys.argv[1]; lam0, fwhm, q, theta = EMITTERS[em]; set_emitter(lam0, fwhm, theta)
    os.makedirs('results_emitter', exist_ok=True); key = 'EQE_q%s' % q
    seeds4 = [np.array(list(json.load(open(f))['structure'].values()), float) for f in sorted(glob.glob('results_refined/*_dbr4_*.json'))]
    with Pool(4) as pool:
        for te, nsub, nd, kind in CASES:
            fn = 'results_emitter/%s_%s_n%.1f_dbr%d_%s.json' % (em, te, nsub, nd, kind)
            if os.path.exists(fn): continue
            bd = N.bounds(te, nd); lo, hi = np.array([b[0] for b in bd]), np.array([b[1] for b in bd])
            ev = lambda x, st: N.evaluate(np.clip(x, lo, hi), te, nd, nsub, step=st, npts=1500, kind=kind, mapper=pool.map, qs=(q,))[key]
            if nd == 0:
                r = differential_evolution(Obj(te, nd, nsub, kind, q), bd, seed=1,
                                           popsize=13, maxiter=30, tol=1e-4, polish=False, workers=pool.map, updating='deferred')
                starts = [r.x]
            else:
                ss = []
                for s in seeds4:
                    s = s.copy(); s[2] = (s[2] if te == 'Ag' and 8 <= s[2] <= 12 else lo[2]); ss.append(np.clip(s, lo, hi))
                starts = [np.array(s) for _, s in sorted(((ev(s, 10), tuple(s)) for s in ss), reverse=True)[:2]]
            best = None
            for s in starts:
                m = minimize(lambda x: -ev(x, 5), s, method='Nelder-Mead', options=dict(maxiter=300 if nd else 150, xatol=0.5, fatol=1e-5))
                xr = np.clip(m.x, lo, hi)
                fin = N.evaluate(xr, te, nd, nsub, step=2, npts=3000, kind=kind, mapper=pool.map, qs=(q, 1.0))
                if best is None or fin[key] > best[0][key]: best = (fin, xr)
            fin, xb = best
            names = ['d_ETL', 'd_HTL', 'd_TE'] + [('SiO2_' if k % 2 == 0 else 'TiO2_') + str(k // 2 + 1) for k in range(2 * nd)]
            out = dict(emitter=em, lam0=lam0, fwhm=fwhm, PLQY=q, Theta=theta, TE=te, n_sub=nsub, ndbr=nd, extractor=kind,
                       EQE=fin[key], EQE_q1=fin['EQE_q1.0'], **{k: v for k, v in fin.items() if not k.startswith('EQE')},
                       structure={k: round(float(v), 1) for k, v in zip(names, xb)})
            json.dump(out, open(fn, 'w'), indent=1); print(em, te, nsub, nd, kind, 'EQE %.4f' % fin[key], flush=True)
    print('EMDONE', em, flush=True)
