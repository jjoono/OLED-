"""Cross-seeded refinement of the DBR cases (local optima of the per-case DE): each dbr4 case is re-evaluated with the
DBR + organic thicknesses of every other dbr4 case (own TE thickness kept, clipped to its bounds), the best three seeds
are polished with Nelder-Mead on a 5 nm spectral grid, and the result replaces the case if better.  -> results_refined/"""
import json, glob, os, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
import nir_opt as N

def load(f):
    d = json.load(open(f)); return d, np.array(list(d['structure'].values()), float)

if __name__ == '__main__':
    files = sorted(glob.glob('results/*_dbr4_q0.3_mla.json')); cases = [load(f) for f in files]
    os.makedirs('results_refined', exist_ok=True)
    with Pool(4) as pool:
        for f, (d, x0) in zip(files, cases):
            te, nsub = d['TE'], d['n_sub']; bd = N.bounds(te, 4); lo, hi = np.array([b[0] for b in bd]), np.array([b[1] for b in bd])
            seeds = []
            for _, xs in cases:
                s = xs.copy(); s[2] = x0[2]; seeds.append(np.clip(s, lo, hi))
            f10 = lambda x: N.evaluate(np.clip(x, lo, hi), te, 4, nsub, step=10, npts=1500, mapper=pool.map)['EQE_q0.3']
            sc = sorted(((f10(s), tuple(s)) for s in seeds), reverse=True)[:3]
            best, xb = d['EQE_q0.3'], x0
            for _, s in sc:
                f5 = lambda x: -N.evaluate(np.clip(x, lo, hi), te, 4, nsub, step=5, npts=1500, mapper=pool.map)['EQE_q0.3']
                r = minimize(f5, np.array(s), method='Nelder-Mead', options=dict(maxiter=400, xatol=0.5, fatol=1e-5))
                xr = np.clip(r.x, lo, hi)
                fin = N.evaluate(xr, te, 4, nsub, step=2, npts=3000, mapper=pool.map, qs=(0.1, 0.3, 0.5, 1.0))
                if fin['EQE_q0.3'] > best + 1e-4: best, xb, bf = fin['EQE_q0.3'], xr, fin
            out = dict(d)
            if best > d['EQE_q0.3'] + 1e-4:
                out.update(bf); out['structure'] = {k: round(float(v), 1) for k, v in zip(d['structure'], xb)}; out['refined'] = True
            else: out['refined'] = False
            json.dump(out, open(f.replace('results/', 'results_refined/'), 'w'), indent=1)
            print(f, 'before %.4f after %.4f' % (d['EQE_q0.3'], out['EQE_q0.3']), flush=True)
