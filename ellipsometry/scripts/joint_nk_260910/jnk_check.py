"""Diagnostics before any fitting.

1. How different are the two Ag-free samples (1-5 HATCN, 1-13 MoOx)?
2. Can a bare oxide-on-Si model describe them?  A seed of 4-5 nm should be
   plainly visible above the noise, and the residual tells us where the
   tabulated Si of jnk_ref stops being good enough - that sets the fit window.
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import numpy as np
from scipy.optimize import least_squares
import jnk_data as D, jnk_ref as R, jnk_model as M


def resid_oxide(d_ox, wl, meas, angles, d_film=0.0, n_film=None):
    layers = [R.n_air(wl)]
    ds = []
    if d_film > 0:
        layers.append(n_film); ds.append(d_film)
    layers += [R.n_sio2(wl), R.n_si(wl)]; ds.append(d_ox)
    N, C, S = M.se_ncs(wl, layers, ds, angles)
    return np.concatenate([(N - meas[0]).ravel(), (C - meas[1]).ravel(), (S - meas[2]).ravel()])


def cauchy(wl, A, B):
    return (A + B * 1e4 / wl**2 + 0j) * np.ones(len(wl))


def main():
    wl, P1, D1, _ = D.se('1-5')
    _,  P2, D2, _ = D.se('1-13')
    print('Ag-free samples 1-5 (HATCN 4 nm) vs 1-13 (MoOx 5 nm)')
    print('  max |dPsi| = %.3f deg   max |dDelta| = %.3f deg   (noise level is ~0.02 deg)'
          % (np.abs(P1 - P2).max(), np.abs(D1 - D2).max()))
    for w in (300., 400., 550., 800., 1000.):
        i = np.argmin(abs(wl - w))
        print('    %6.0f nm   dPsi %+7.3f   dDelta %+7.3f' % (w, (P1 - P2)[i, 2], (D1 - D2)[i, 2]))

    print('\nbare oxide-on-Si fit (no seed layer at all), window scan')
    print('%-14s %8s %8s   %s' % ('window (nm)', 'd_ox nm', 'MSE', 'residual by band (x1e3)'))
    for lo, hi in [(260, 1080), (300, 1080), (350, 1080), (400, 1080), (450, 1000)]:
        for sh in ('1-5', '1-13'):
            wl, P, Dl, _ = D.se(sh, lo, hi)
            meas = D.ncs(P, Dl)
            f = lambda p: resid_oxide(p[0], wl, meas, D.ANG)
            r = least_squares(f, [2.0], bounds=([0.0], [50.0]))
            res = r.fun
            n = len(wl)
            bands = []
            for b0, b1 in [(lo, 400), (400, 600), (600, 800), (800, hi)]:
                m = (wl >= b0) & (wl < b1)
                if m.sum() == 0:
                    continue
                sel = np.concatenate([np.tile(m, 5)] * 3)
                bands.append('%d-%d:%5.1f' % (b0, b1, 1e3 * np.sqrt(np.mean(res[sel]**2))))
            print('%-14s %8.2f %8.2f   %s'
                  % ('%s  %d-%d' % (sh, lo, hi), r.x[0], M.mse_ncs(res, 1), '  '.join(bands)))

    print('\nsame, but with a transparent Cauchy seed on top of a 2 nm oxide')
    print('%-8s %8s %8s %8s %8s' % ('sheet', 'd_seed', 'n@633', 'd_ox', 'MSE'))
    for sh in ('1-5', '1-13'):
        wl, P, Dl, _ = D.se(sh, 400, 1080)
        meas = D.ncs(P, Dl)
        def f(p):
            return resid_oxide(p[3], wl, meas, D.ANG, d_film=p[0], n_film=cauchy(wl, p[1], p[2]))
        r = least_squares(f, [4.0, 1.8, 0.02, 2.0],
                          bounds=([0.0, 1.0, -0.5, 0.0], [40.0, 3.5, 0.5, 20.0]))
        n633 = cauchy(np.array([633.0]), r.x[1], r.x[2]).real[0]
        print('%-8s %8.2f %8.3f %8.2f %8.2f' % (sh, r.x[0], n633, r.x[3], M.mse_ncs(r.fun, 4)))


if __name__ == '__main__':
    main()
