"""Bphen 10 nm on Si, round 2: is the leftover UV residual the SUBSTRATE?

The measured-minus-bare signal of 2_23 has its extrema at Si's E1 (365 nm) and
E2 (~290 nm) critical points, not at Bphen's 277 nm band.  A Si surface that is
not the reference SI_JAW - damaged, roughened, an interfacial mix - broadens and
weakens exactly those features.  Tested here, alone and with the island picture:

  oxide free         NTVE_JAW thickness free
  Si damage layer    Bruggeman(SI_JAW, void) between the oxide and the Si,
                     thickness and void fraction free
  ...+ islands       the same substrate under the area-mixed island film
and the island picture again under the other two readings of a depolarized export."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, json, zlib, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_ce_genosc as C, tl_glass_global as GG, tl_bphen10 as B
import jnk_model as M
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')

MODELS = {
    #                        names                              x0                          lo                           hi
    'bare, oxide free':      (['dox'],                          [3.4],                      [0.0],                       [30.0]),
    'film, oxide free':      (['d', 'dox'],                     [5.0, 3.4],                 [0.0, 0.0],                  [80.0, 30.0]),
    'Si damage only':        (['t', 'fv'],                      [2.0, 0.3],                 [0.0, 0.0],                  [30.0, 0.9]),
    'film + Si damage':      (['d', 't', 'fv'],                 [5.0, 2.0, 0.3],            [0.0, 0.0, 0.0],             [80.0, 30.0, 0.9]),
    'islands + Si damage':   (['a', 'h', 'w', 't', 'fv'],       [0.25, 29.0, 0.5, 2.0, 0.3],[0.0, 1.0, 0.0, 0.0, 0.0],   [1.0, 300.0, 20.0, 30.0, 0.9]),
    'islands + oxide free':  (['a', 'h', 'w', 'dox'],           [0.25, 29.0, 0.5, 3.4],     [0.0, 1.0, 0.0, 0.0],        [1.0, 300.0, 20.0, 30.0]),
    'islands (psi-from-N)':  (['a', 'h', 'w'],                  [0.25, 29.0, 0.5],          [0.0, 1.0, 0.0],             [1.0, 300.0, 20.0]),
    'islands (raw NCS)':     (['a', 'h', 'w'],                  [0.25, 29.0, 0.5],          [0.0, 1.0, 0.0],             [1.0, 300.0, 20.0]),
}


def substrate(wl, v):
    air = np.ones(len(wl), complex)
    ox = C.N_ce('NTVE_JAW', wl)**2; si = C.N_ce('SI_JAW', wl)**2
    lay = [ox]; ds = [v.get('dox', C.D_OX)]
    if v.get('t', 0) > 1e-3:
        lay.append(B.bruggeman_eps(si, air, 1 - v['fv'])); ds.append(v['t'])
    return lay + [si], ds


def fit(sheet, model, n_start=24):
    CE = json.load(open(W + 'ce_genosc.json')); p = np.array(CE['Bphen']['with_roughness']['p'])
    wl, P, D, ang = L.se(sheet, 245.0, 1000.0, 2); meas = L.ncs(P, D)
    names, x0, lo, hi = MODELS[model]; lo, hi = np.array(lo, float), np.array(hi, float)
    conv = GG.CONV['psi-from-N' if 'psi-from-N' in model else 'raw NCS' if 'raw' in model else 'normalised']
    air = np.ones(len(wl), complex)
    def regions(x):
        v = dict(zip(names, x)); sub, dsub = substrate(wl, v)
        eb, _ = B.eps_bphen(p, wl)
        def stack(film):
            return [air] + [f for f, _ in film] + sub, [air] + [f for f, _ in film] + sub, [d for _, d in film] + dsub
        if 'islands' in model:
            isl = stack([(eb, v['h'])]); bare = stack([(eb, v['w'])] if v['w'] > 1e-3 else [])
            return [(v['a'], *isl), (1 - v['a'], *bare)]
        film = [(eb, v['d'])] if v.get('d', 0) > 1e-3 else []
        return [(1.0, *stack(film))]
    def res(x):
        N, C_, S_ = B.ncs_regions(wl, ang, regions(x), conv)
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
    return dict(sheet=sheet, model=model, x=dict(zip(names, [float(v) for v in best.x])),
                mse=float(M.mse_ncs(best.fun, len(names))),
                bands={'UV<350': band(200, 350), '350-500': band(350, 500), '500-1000': band(500, 1001)})


def main():
    from multiprocessing import Pool
    jobs = [(sh, m) for sh in ('2_23', '2_24') for m in MODELS]
    with Pool(4) as pool:
        res = pool.starmap(fit, jobs)
    json.dump(res, open(W + 'bphen10_sub.json', 'w'), indent=1)
    for sh in ('2_23', '2_24'):
        print('== %s' % sh)
        for r in sorted([r for r in res if r['sheet'] == sh], key=lambda r: r['mse']):
            print('   %-22s MSE %7.2f | UV %6.1f  350-500 %6.1f  500-1000 %6.1f | %s'
                  % (r['model'], r['mse'], r['bands']['UV<350'], r['bands']['350-500'], r['bands']['500-1000'],
                     '  '.join('%s=%.3g' % kv for kv in r['x'].items())))


if __name__ == '__main__':
    main()
