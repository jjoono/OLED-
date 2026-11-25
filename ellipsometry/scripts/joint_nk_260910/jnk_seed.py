"""Stage 1 - the two Ag-free samples give the seed dispersion.

Stack: air / seed (Einf + Tauc-Lorentz) / native oxide / Si.

Seed thickness and oxide thickness are optically interchangeable at this
thinness - only their sum is measured - so the oxide is scanned rather than
fitted and every quoted seed thickness carries the assumed oxide with it.
The seed must be allowed to absorb: HATCN and MoOx both have their onset
inside the measured range, and a transparent Cauchy leaves a residual in the
UV that looks exactly like a bad substrate.
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import json, numpy as np
from scipy.optimize import least_squares
import jnk_data as D, jnk_ref as R, jnk_model as M

WL_LO, WL_HI = 260.0, 1080.0
OUT = _os.path.join(ELLIPS_OUT, 'jnk_seed.json')

#        Einf   TL_A   TL_Br  TL_Eo  TL_Eg
P0    = [2.20,  40.0,  1.50,  5.00,  3.00]
LO    = [1.00,   0.0,  0.05,  2.50,  1.50]
HI    = [6.00, 400.0,  8.00,  9.50,  5.50]


def seed_N(p, wl):
    return M.osc_N(p[0], [('tl', p[1], p[2], p[3], p[4])], wl)


def stack(p_seed, d_seed, d_ox, wl):
    return ([R.n_air(wl), seed_N(p_seed, wl), R.n_sio2(wl), R.n_si(wl)],
            [d_seed, d_ox])


def fit(sheet, d_ox, lo=WL_LO, hi=WL_HI, d0=5.0, seed_p0=None):
    wl, P, Dl, _ = D.se(sheet, lo, hi)
    meas = D.ncs(P, Dl)
    def f(p):
        layers, ds = stack(p[1:], p[0], d_ox, wl)
        N, C, S = M.se_ncs(wl, layers, ds, D.ANG)
        return np.concatenate([(N - meas[0]).ravel(), (C - meas[1]).ravel(), (S - meas[2]).ravel()])
    x0 = [d0] + list(seed_p0 if seed_p0 is not None else P0)
    r = least_squares(f, x0, bounds=([0.2] + LO, [40.0] + HI), x_scale='jac')
    return r.x, M.mse_ncs(r.fun, len(x0)), r.fun, wl


def bands(res, wl, edges=(260, 350, 400, 600, 800, 1080)):
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (wl >= a) & (wl < b)
        if m.sum() == 0:
            continue
        sel = np.concatenate([np.tile(m, 5)] * 3)
        out.append('%d-%d:%5.1f' % (a, b, 1e3 * np.sqrt(np.mean(res[sel]**2))))
    return '  '.join(out)


def main():
    res = {}
    print('seed dispersion from the Ag-free samples, %g-%g nm' % (WL_LO, WL_HI))
    print('%-6s %-6s %6s %8s %7s %7s %7s %7s %7s %7s   %s'
          % ('sheet', 'seed', 'd_ox', 'd_seed', 'Einf', 'TL_A', 'TL_Br', 'TL_Eo', 'TL_Eg', 'MSE',
             'residual by band (x1e3)'))
    for seed, sh in D.SEEDS.items():
        best = None
        recs = []
        for d_ox in (0.0, 1.0, 2.0, 3.0, 4.0):
            x, mse, r, wl = fit(sh, d_ox)
            recs.append(dict(d_ox=d_ox, p=list(x), mse=float(mse)))
            print('%-6s %-6s %6.1f %8.2f %7.3f %7.2f %7.3f %7.3f %7.3f %7.2f   %s'
                  % (sh, seed, d_ox, x[0], x[1], x[2], x[3], x[4], x[5], mse, bands(r, wl)))
            if best is None or mse < best['mse']:
                best = recs[-1]
        res[seed] = dict(sheet=sh, scan=recs, best=best)
        print()
    json.dump(res, open(OUT, 'w'), indent=1)
    print('saved ->', OUT)

    wlp = np.array([450., 550., 633., 800.])
    print('\nseed n,k at the 2 nm-oxide solution')
    print('%-6s %s' % ('seed', '   '.join('%6.0f nm' % w for w in wlp)))
    for seed in res:
        rec = [r for r in res[seed]['scan'] if r['d_ox'] == 2.0][0]
        N = seed_N(np.array(rec['p'][1:]), wlp)
        print('%-6s %s   (d_seed = %.2f nm)'
              % (seed, '   '.join('%5.3f+%5.3fj' % (n.real, n.imag) for n in N), rec['p'][0]))


if __name__ == '__main__':
    main()
