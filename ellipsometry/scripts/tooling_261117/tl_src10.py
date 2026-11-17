"""Src10_10nm: a nominal 10 nm organic film on Si, material not given.

SRC10_XLSX=/path/to/Src10_10nm.xlsx (one sheet, same export layout as
tooling.xlsx; default $ELLIPS_DATA/Src10_10nm.xlsx).  CompleteEASE frame as the
tooling set (SI_JAW, NTVE_JAW); the NPB tooling pieces are refitted alongside
as the like-for-like reference.
Every picture is fitted on the same data and footing:
  * bare Si with a free oxide          - is there a film signal at all?
  * Si damage only                     - the reading that won for Bphen 2_23 / 2_24
  * each known dispersion (NPB, TAPC, Bphen) with d free, oxide 33.68 A or 0 A,
    with and without a top roughness layer
  * a material-agnostic transparent film (Cauchy A, B) with d free
The UV residual is reported separately: a 10 nm organic film shows itself
mainly through its absorption bands below 400 nm.

Second pass (CHECKS), on Src10 and on the NPB tooling pieces 1_21 / 1_22 / 1_23:
  * NPB + Si damage   - does the extra interface layer help every NPB piece?
  * oxide free        - d and the oxide together
  * band scale        - Gaussian amplitudes x s: does a 10 nm film have the
                        band strength of the 30 / 50 nm pieces the NPB
                        dispersion was fitted to?"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, json, zlib, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_ce_genosc as C, tl_glass_global as GG, tl_bphen10 as B, tl_genosc as GO
import jnk_model as M
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')
SHEET = 'Sheet1'
SRC10 = _os.environ.get('SRC10_XLSX', _os.path.join(ELLIPS_DATA, 'Src10_10nm.xlsx'))
MATS = ('NPB', 'TAPC', 'Bphen')


def eps_mat(mat, wl, s=1.0):
    """The material's Gen-Osc; s scales the Gaussian amplitudes (UV pole kept)."""
    CE = json.load(open(W + 'ce_genosc.json'))
    p = np.array(CE[mat]['pick']['p'], float); p[2::3] *= s
    return GO.eps_genosc(p, wl)


def model_spec(model, wl):
    """-> names, x0, lo, hi, stack(v) -> (eps layers top-down incl. air and Si, thicknesses)."""
    air = np.ones(len(wl), complex); ox = C.N_ce('NTVE_JAW', wl)**2; si = C.N_ce('SI_JAW', wl)**2
    rB = lambda e: B.bruggeman_eps(e, air, 0.5)
    if model == 'bare Si, oxide free':
        return ['dox'], [3.4], [0.0], [40.0], lambda v: ([air, ox, si], [v['dox']])
    if model == 'Si damage only':
        return (['t', 'fv'], [2.0, 0.3], [0.0, 0.0], [30.0, 0.9],
                lambda v: ([air, ox, B.bruggeman_eps(si, air, 1 - v['fv']), si], [C.D_OX, v['t']]))
    if model.startswith('Cauchy'):
        n = lambda v: v['A'] + v['B'] / (wl / 1000)**2 + 0j
        return (['d', 'A', 'B'], [10.0, 1.7, 0.01], [0.0, 1.3, 0.0], [60.0, 2.4, 0.1],
                lambda v: ([air, n(v)**2, ox, si], [v['d'], C.D_OX]))
    mat = model.split()[0]; ef = eps_mat(mat, wl)
    dox = 0.0 if 'oxide 0' in model else C.D_OX
    if 'oxide free' in model:
        return (['d', 'dox'], [10.0, 3.4], [0.0, 0.0], [60.0, 15.0],
                lambda v: ([air, ef, ox, si], [v['d'], v['dox']]))
    if 'band scale' in model:
        return (['d', 's'], [10.0, 1.0], [0.0, 0.3], [60.0, 2.0],
                lambda v: ([air, eps_mat(mat, wl, v['s']), ox, si], [v['d'], dox]))
    if 'rough' in model:
        return (['d', 'r'], [10.0, 1.0], [0.0, 0.0], [60.0, 20.0],
                lambda v: ([air, rB(ef), ef, ox, si], [v['r'], v['d'], dox]))
    if 'Si damage' in model:
        return (['d', 't', 'fv'], [10.0, 2.0, 0.3], [0.0, 0.0, 0.0], [60.0, 30.0, 0.9],
                lambda v: ([air, ef, ox, B.bruggeman_eps(si, air, 1 - v['fv']), si], [v['d'], C.D_OX, v['t']]))
    return ['d'], [10.0], [0.0], [60.0], lambda v: ([air, ef, ox, si], [v['d'], dox])


def fit(model, sheet=SHEET, n_start=16, lo_wl=245.0, hi_wl=1000.0):
    wl, P, D, ang = L.se(sheet, lo_wl, hi_wl, 2, xl=SRC10 if sheet == SHEET else None); meas = L.ncs(P, D)
    names, x0, lo, hi, st = model_spec(model, wl); lo, hi = np.array(lo, float), np.array(hi, float)
    def res(x):
        v = dict(zip(names, x)); eo, ds = st(v)
        keep = [0] + [i + 1 for i, d in enumerate(ds) if d > 1e-4] + [len(eo) - 1]   # drop vanishing layers
        eo = [eo[i] for i in keep]; ds = [d for d in ds if d > 1e-4]
        N, C_, S_ = B.ncs_regions(wl, ang, [(1.0, eo, eo, ds)], GG.CONV['normalised'])
        return np.concatenate([(N - meas[0]).ravel(), (C_ - meas[1]).ravel(), (S_ - meas[2]).ravel()])
    rng = np.random.default_rng(zlib.crc32((sheet + model).encode())); best = None
    for k in range(n_start):
        s = np.array(x0, float) if k == 0 else lo + rng.random(len(lo)) * (hi - lo)
        r = least_squares(res, np.clip(s, lo + 1e-6, hi - 1e-6), bounds=(lo, hi), max_nfev=600)
        if best is None or r.cost < best.cost: best = r
    rr = best.fun.reshape(3, len(wl), len(ang))
    band = lambda a, b: float(1e3 * np.sqrt(np.mean(rr[:, (wl >= a) & (wl < b), :]**2)))
    J = best.jac; cov = np.linalg.pinv(J.T @ J) * (best.fun @ best.fun) / max(len(best.fun) - len(names), 1)
    return dict(sheet=sheet, model=model, x={k: float(v) for k, v in zip(names, best.x)},
                ci90={k: float(1.645 * np.sqrt(max(cov[i, i], 0))) for i, k in enumerate(names)},
                mse=float(M.mse_ncs(best.fun, len(names))),
                bands={'245-350': band(245, 350), '350-500': band(350, 500), '500-1000': band(500, 1001)})


MODELS = (['bare Si, oxide free', 'Si damage only', 'Cauchy film']
          + ['%s film' % m for m in MATS] + ['%s film, oxide 0' % m for m in MATS]
          + ['%s film + rough' % m for m in MATS] + ['%s film + Si damage' % m for m in MATS])
CHECKS = [(m, sh) for sh in (SHEET, '1_23', '1_21', '1_22')
          for m in ('NPB film', 'NPB film + Si damage', 'NPB film, oxide free', 'NPB film, band scale')]


def main():
    from multiprocessing import Pool
    with Pool(4) as pool:
        res = pool.map(fit, MODELS)
        chk = pool.starmap(fit, CHECKS)
    json.dump(dict(models=res, checks=chk), open(W + 'src10_fits.json', 'w'), indent=1)
    print('%-26s %7s | rms x1e3: %7s %7s %8s | parameters' % ('model', 'MSE', '245-350', '350-500', '500-1000'))
    for r in sorted(res, key=lambda r: r['mse']):
        print('%-26s %7.2f |           %7.1f %7.1f %8.1f | %s' % (r['model'], r['mse'], r['bands']['245-350'], r['bands']['350-500'],
              r['bands']['500-1000'], '  '.join('%s=%.3g' % kv for kv in r['x'].items())))
    print('\nNPB checks: Src10 (%s) next to the NPB tooling pieces' % SHEET)
    for r in chk:
        print('%-7s %-22s MSE %6.2f | UV %5.1f | %s' % ('Src10' if r['sheet'] == SHEET else r['sheet'], r['model'], r['mse'],
              r['bands']['245-350'], '  '.join('%s=%.3g+-%.2g' % (k, v, r['ci90'][k]) for k, v in r['x'].items())))


if __name__ == '__main__':
    main()
