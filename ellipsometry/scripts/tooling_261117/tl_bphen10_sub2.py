"""Round 3.  (1) Does the Bphen 50 nm pair from the same run (2_21/2_22) carry
the same substrate anomaly?  (2) If the rough Si's valleys are filled with
Bphen rather than air, how much Bphen does the data allow?"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, json, zlib, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_ce_genosc as C, tl_glass_global as GG, tl_bphen10 as B
import jnk_model as M
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')

def fit(sheet, model, n_start=20):
    CE = json.load(open(W + 'ce_genosc.json')); p = np.array(CE['Bphen']['with_roughness']['p'])
    wl, P, D, ang = L.se(sheet, 245.0, 1000.0, 2); meas = L.ncs(P, D)
    air = np.ones(len(wl), complex); ox = C.N_ce('NTVE_JAW', wl)**2; si = C.N_ce('SI_JAW', wl)**2
    eb, _ = B.eps_bphen(p, wl)
    rB = lambda e, rr: B.bruggeman_eps(e, air, 0.5)
    if model == 'film + rough top':                    # the clean-substrate reference model of the 50 nm pieces
        names, x0, lo, hi = ['d', 'r'], [45, 2], [0, 0], [100, 20]
        st = lambda v: ([air, rB(eb, 0), eb, ox, si], [v['r'], v['d'], C.D_OX])
    elif model == 'film + rough top + Si damage':
        names, x0, lo, hi = ['d', 'r', 't', 'fv'], [45, 2, 2, 0.3], [0, 0, 0, 0], [100, 20, 30, 0.9]
        st = lambda v: ([air, rB(eb, 0), eb, ox, B.bruggeman_eps(si, air, 1 - v['fv']), si], [v['r'], v['d'], C.D_OX, v['t']])
    elif model == 'rough Si filled with Bphen':        # EMA(Si, Bphen) under the oxide? no: Bphen sits ON the oxide-capped roughness
        names, x0, lo, hi = ['t', 'fsi', 'd'], [13, 0.77, 1.0], [0, 0.1, 0], [40, 1.0, 40]
        st = lambda v: ([air, eb, B.bruggeman_eps(si, eb, v['fsi']), ox, si], [v['d'], v['t'], C.D_OX])
    elif model == 'rough Si (air) + Bphen':
        names, x0, lo, hi = ['t', 'fsi', 'd'], [13, 0.77, 1.0], [0, 0.1, 0], [40, 1.0, 40]
        st = lambda v: ([air, eb, B.bruggeman_eps(si, air, v['fsi']), ox, si], [v['d'], v['t'], C.D_OX])
    lo, hi = np.array(lo, float), np.array(hi, float)
    def res(x):
        v = dict(zip(names, x)); eo, ds = st(v)
        if 'r' in v and v['r'] < 1e-4:                  # drop a vanishing roughness layer
            eo, ds = eo[:1] + eo[2:], ds[1:]
        N, C_, S_ = B.ncs_regions(wl, ang, [(1.0, eo, eo, ds)], GG.CONV['normalised'])
        return np.concatenate([(N - meas[0]).ravel(), (C_ - meas[1]).ravel(), (S_ - meas[2]).ravel()])
    rng = np.random.default_rng(zlib.crc32(repr((sheet, model)).encode())); best = None
    for k in range(n_start):
        s = np.array(x0, float) if k == 0 else lo + rng.random(len(lo)) * (hi - lo)
        try:
            r = least_squares(res, np.clip(s, lo + 1e-6, hi - 1e-6), bounds=(lo, hi), max_nfev=800)
        except Exception:
            continue
        if best is None or r.cost < best.cost: best = r
    rr = best.fun.reshape(3, len(wl), len(ang))
    band = lambda a, b: float(1e3 * np.sqrt(np.mean(rr[:, (wl >= a) & (wl < b), :]**2)))
    v = dict(zip(names, [float(x) for x in best.x]))
    if 'fsi' in v and 'filled' in model:
        v['Bphen equivalent nm'] = v['d'] + v['t'] * (1 - v['fsi'])
    return dict(sheet=sheet, model=model, x=v, mse=float(M.mse_ncs(best.fun, len(names))),
                bands={'UV<350': band(200, 350), '350-500': band(350, 500), '500-1000': band(500, 1001)})

def main():
    from multiprocessing import Pool
    jobs = [(sh, m) for sh in ('2_21', '2_22') for m in ('film + rough top', 'film + rough top + Si damage')]
    jobs += [(sh, m) for sh in ('2_23', '2_24') for m in ('rough Si filled with Bphen', 'rough Si (air) + Bphen')]
    with Pool(4) as pool:
        res = pool.starmap(fit, jobs)
    json.dump(res, open(W + 'bphen10_sub2.json', 'w'), indent=1)
    for r in res:
        print('%s  %-30s MSE %6.2f | UV %5.1f 350-500 %5.1f 500-1000 %5.1f | %s'
              % (r['sheet'], r['model'], r['mse'], r['bands']['UV<350'], r['bands']['350-500'], r['bands']['500-1000'],
                 '  '.join('%s=%.3g' % kv for kv in r['x'].items())))

if __name__ == '__main__':
    main()
