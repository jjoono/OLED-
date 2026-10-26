"""Effect of the number of DBR pairs.  (1) quarter-wave TiO2/SiO2 DBR at the emission peak (TiO2 next to the TE),
N = 0, 2, 4, 6, 8 pairs, ETL/HTL/TE optimised (DE + NM);  (2) free chirped DBR with 6 and 8 pairs, Nelder-Mead from the
best 4-pair structure padded with zero-thickness layers and from the quarter-wave optimum.  -> results_pairs/
Usage: python3 dbr_pairs.py <emitter> <TE> <n_sub>"""
import sys, json, os, glob, numpy as np
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize
import nir_opt as N, emitter_opt as E

class Obj:
    def __init__(self, *a): self.a = a
    def __call__(self, y):
        te, nd, nsub, q, dbr = self.a
        return -N.evaluate(np.r_[y, dbr], te, nd, nsub, step=10, npts=1500, qs=(q,))['EQE_q%s' % q]

def qw(nd, lam0):
    j = int(lam0 - 600); dH, dL = lam0 / 4 / N.TIO2[j].real, lam0 / 4 / N.SIO2[j].real
    return np.array([0.0 if k == 0 else (dH if k % 2 else dL) for k in range(2 * nd)])

if __name__ == '__main__':
    em, te, nsub = sys.argv[1], sys.argv[2], float(sys.argv[3]); lam0, fwhm, q, th = E.EMITTERS[em]; E.set_emitter(lam0, fwhm, th)
    key = 'EQE_q%s' % q; os.makedirs('results_pairs', exist_ok=True); out = {}
    with Pool(4) as pool:
        for nd in (0, 2, 4, 6, 8):
            bd = N.bounds(te, nd)[:3]; lo, hi = np.array([b[0] for b in bd]), np.array([b[1] for b in bd]); dbr = qw(nd, lam0)
            r = differential_evolution(Obj(te, nd, nsub, q, dbr), bd, seed=1, popsize=13, maxiter=30, tol=1e-4, polish=False,
                                       workers=pool.map, updating='deferred')
            m = minimize(lambda y: -N.evaluate(np.r_[np.clip(y, lo, hi), dbr], te, nd, nsub, step=5, npts=1500, mapper=pool.map, qs=(q,))[key],
                         r.x, method='Nelder-Mead', options=dict(maxiter=150, xatol=0.5, fatol=1e-5))
            x = np.r_[np.clip(m.x, lo, hi), dbr]
            fin = N.evaluate(x, te, nd, nsub, step=2, npts=3000, mapper=pool.map, qs=(q, 1.0))
            out['qw%d' % nd] = dict(x=[round(float(v), 1) for v in x], **fin); print(em, te, nsub, 'QW', nd, 'EQE %.4f F %.3f eta_sub %.3f eta_ext %.3f A %.3f'
                  % (fin[key], fin['F'], fin['eta_sub'], fin['eta_ext'], fin['A_round']), flush=True)
        f4 = 'results_emitter/%s_%s_n%.1f_dbr4_mla.json' % (em, te, nsub); x4 = np.array(list(json.load(open(f4))['structure'].values()), float)
        for nd in (6, 8):
            bd = N.bounds(te, nd); lo, hi = np.array([b[0] for b in bd]), np.array([b[1] for b in bd])
            ev = lambda x: -N.evaluate(np.clip(x, lo, hi), te, nd, nsub, step=5, npts=1500, mapper=pool.map, qs=(q,))[key]
            seeds = [np.r_[x4, np.zeros(2 * (nd - 4))], np.array(out['qw%d' % nd]['x'])]
            best = None
            for s in seeds:
                m = minimize(ev, s, method='Nelder-Mead', options=dict(maxiter=500, xatol=0.5, fatol=1e-5)); xr = np.clip(m.x, lo, hi)
                fin = N.evaluate(xr, te, nd, nsub, step=2, npts=3000, mapper=pool.map, qs=(q, 1.0))
                if best is None or fin[key] > best[1][key]: best = (xr, fin)
            out['free%d' % nd] = dict(x=[round(float(v), 1) for v in best[0]], **best[1]); fin = best[1]
            print(em, te, nsub, 'FREE', nd, 'EQE %.4f F %.3f eta_sub %.3f eta_ext %.3f A %.3f'
                  % (fin[key], fin['F'], fin['eta_sub'], fin['eta_ext'], fin['A_round']), flush=True)
    json.dump(out, open('results_pairs/%s_%s_n%.1f.json' % (em, te, nsub), 'w'), indent=1); print('PAIRSDONE', flush=True)
