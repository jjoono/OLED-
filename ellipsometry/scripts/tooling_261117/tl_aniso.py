"""Anisotropic films at a fixed 1.5 nm NTVE_JAW: does a biaxial film fit where
the isotropic one does not?

Pieces: Src10 (NPB, SRC10_XLSX), its reference 1_23 (8.3 nm NPB), and the two
nominal 10 nm Bphen Si pieces 2_23 / 2_24 that no film model fitted (section 3).
Oxide NTVE_JAW fixed at 1.5 nm (a bare-wafer native oxide), with the measured
33.68 A beside it.

Film = the material's Gen-Osc (as in its Si .mod files) with the Gaussian band
amplitudes scaled per axis; the UV pole (transitions above the data) stays
isotropic.
  iso        one scale fixed at 1, d free
  uniaxial   a_o (in plane) and a_e (normal) free, d free
  biaxial    a_x, a_y, a_z all free, d free
The data are at one azimuth, so x (in the plane of incidence), y (along s) and
z (normal) are the only axes it can tell apart: s light sees eps_y only, p
light eps_x and eps_z.  Axes are taken along the lab frame (no off-diagonal
terms).  a = 1 on every axis is the isotropic film; (a_x + a_y + a_z) / 3 is
the total band strength relative to the dispersion."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, json, zlib, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_ce_genosc as C, tl_glass_global as GG, tl_bphen10 as B, tl_genosc as GO, tl_ce_thickness as T
import jnk_model as M
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')
SRC10 = _os.environ.get('SRC10_XLSX', _os.path.join(ELLIPS_DATA, 'Src10_10nm.xlsx'))
PIECES = {'Src10': 'NPB', '1_23': 'NPB', '2_23': 'Bphen', '2_24': 'Bphen'}
OXIDES = (1.5, C.D_OX)


def eps_axis(mat, wl, a):
    p = np.array(T.disp(mat)['p'], float); p[2::3] *= a
    return GO.eps_genosc(p, wl)


def ncs_biaxial(wl, ang, ex, ey, ez, ds):
    """ex/ey/ez = per-layer eps lists air..Si (isotropic layers: the same object thrice)."""
    out = []
    for th in ang:
        rp = B.tmm_uni(wl, ex, ez, ds, th)[0]
        rs = B.tmm_uni(wl, ey, ey, ds, th)[1]
        Ip, Is, X = np.abs(rp)**2, np.abs(rs)**2, np.conj(rp) * rs
        I = Ip + Is
        out.append(GG.CONV['normalised']((Is - Ip) / I, 2 * X.real / I, 2 * X.imag / I))
    return tuple(np.stack([o[k] for o in out], 1) for k in range(3))


SPEC = {   # names, x0, lo, hi
    'iso':                 (['d'], [10.0], [0.0], [60.0]),
    'uniaxial':            (['d', 'a_o', 'a_e'], [10.0, 1.0, 1.0], [0.0, 0.0, 0.0], [60.0, 3.0, 3.0]),
    'biaxial':             (['d', 'a_x', 'a_y', 'a_z'], [10.0, 1.0, 1.0, 1.0], [0.0] * 4, [60.0, 3.0, 3.0, 3.0]),
    'no film, Si/void EMA': (['t', 'fv'], [5.0, 0.3], [0.0, 0.0], [30.0, 0.9]),
}


def fit(piece, model, dox, n_start=16):
    mat = PIECES[piece]
    sh, xl = ('Sheet1', SRC10) if piece == 'Src10' else (piece, None)
    wl, P, D, ang = L.se(sh, 245.0, 1000.0, 2, xl=xl); meas = L.ncs(P, D)
    air = np.ones(len(wl), complex); ox = C.N_ce('NTVE_JAW', wl)**2; si = C.N_ce('SI_JAW', wl)**2
    names, x0, lo, hi = SPEC[model]; lo, hi = np.array(lo, float), np.array(hi, float)
    def res(x):
        v = dict(zip(names, x))
        if model == 'no film, Si/void EMA':
            lay = [air, ox, B.bruggeman_eps(si, air, 1 - v['fv']), si]; ds = [dox, v['t']]
            ex = ey = ez = lay
        else:
            a = {'iso': (1, 1, 1), 'uniaxial': (v.get('a_o'), v.get('a_o'), v.get('a_e')),
                 'biaxial': (v.get('a_x'), v.get('a_y'), v.get('a_z'))}[model]
            fx, fy, fz = (eps_axis(mat, wl, ai) for ai in a)
            ex, ey, ez = ([air, f, ox, si] for f in (fx, fy, fz)); ds = [v['d'], dox]
        N, C_, S_ = ncs_biaxial(wl, ang, ex, ey, ez, ds)
        return np.concatenate([(N - meas[0]).ravel(), (C_ - meas[1]).ravel(), (S_ - meas[2]).ravel()])
    rng = np.random.default_rng(zlib.crc32(repr((piece, model, dox)).encode())); best = None
    for k in range(n_start):
        s = np.array(x0, float) if k == 0 else lo + rng.random(len(lo)) * (hi - lo)
        r = least_squares(res, np.clip(s, lo + 1e-6, hi - 1e-6), bounds=(lo, hi), max_nfev=600)
        if best is None or r.cost < best.cost: best = r
    J = best.jac; cov = np.linalg.pinv(J.T @ J) * (best.fun @ best.fun) / max(len(best.fun) - len(names), 1)
    rr = best.fun.reshape(3, len(wl), len(ang))
    band = lambda lo_, hi_: float(1e3 * np.sqrt(np.mean(rr[:, (wl >= lo_) & (wl < hi_), :]**2)))
    return dict(piece=piece, mat=mat, model=model, dox=dox, x={k: float(v) for k, v in zip(names, best.x)},
                ci90={k: float(1.645 * np.sqrt(max(cov[i, i], 0))) for i, k in enumerate(names)},
                mse=float(M.mse_ncs(best.fun, len(names))), uv=band(245, 350), vis=band(350, 1001))


def main():
    from multiprocessing import Pool
    jobs = [(pc, m, dox) for pc in PIECES for dox in OXIDES for m in ('iso', 'uniaxial', 'biaxial')]
    jobs += [(pc, 'no film, Si/void EMA', dox) for pc in ('2_23', '2_24', 'Src10') for dox in OXIDES]
    with Pool(4) as pool:
        res = pool.starmap(fit, jobs)
    json.dump(res, open(W + 'aniso_fits.json', 'w'), indent=1)
    for pc in PIECES:
        print('\n== %s (%s)' % (pc, PIECES[pc]))
        for r in [r for r in res if r['piece'] == pc]:
            print('  NTVE %4.2f nm  %-21s MSE %6.2f | UV %5.1f vis %5.1f | %s'
                  % (r['dox'], r['model'], r['mse'], r['uv'], r['vis'],
                     '  '.join('%s=%.3g+-%.2g' % (k, v, r['ci90'][k]) for k, v in r['x'].items())))


if __name__ == '__main__':
    main()
