"""Si-wafer pieces: air / Cauchy film / native oxide 2 nm / Si over 450-1000 nm.

For an unknown film the danger is the n*d alias, so the thickness comes from an
exhaustive (d, A) grid evaluated directly - no optimiser - and only then is each
of the best distinct minima refined.  The margin to the runner-up branch is
reported so a weak identification is visible."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # tooling.xlsx
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fits, tables, figures
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L
import jnk_ref as R, jnk_model as M
from scipy.optimize import least_squares
from multiprocessing import Pool

D_OX, LO, HI = 2.0, 450.0, 1000.0
OUT = _os.path.join(ELLIPS_OUT, 'cauchy_si.json')

def cauchy(wl, A, B, C=0.0):
    l = wl / 1000.0
    return A + B / l**2 + C / l**4

def resid(p, wl, meas, ang):
    N = cauchy(wl, p[1], p[2], p[3]) + 0j
    n_, c_, s_ = M.se_ncs(wl, [R.n_air(wl), N, R.n_sio2(wl), R.n_si(wl)], [p[0], D_OX], ang)
    return np.concatenate([(n_ - meas[0]).ravel(), (c_ - meas[1]).ravel(), (s_ - meas[2]).ravel()])

def job(sh):
    wl, P, D, ang = L.se(sh, LO, HI, step=4)
    meas = L.ncs(P, D)
    ds = np.arange(1.0, 501.0, 1.0); As = np.arange(1.50, 2.001, 0.025)
    C = np.empty((len(ds), len(As)))
    for i, d in enumerate(ds):
        for j, A in enumerate(As):
            C[i, j] = np.sum(resid([d, A, 0.01, 0.0], wl, meas, ang)**2)
    # local minima in d of the best-over-A profile
    prof = C.min(1)
    mins = [i for i in range(1, len(ds) - 1) if prof[i] <= prof[i - 1] and prof[i] <= prof[i + 1]]
    mins = sorted(mins, key=lambda i: prof[i])[:4]
    wl2, P2, D2, _ = L.se(sh, LO, HI, step=2); meas2 = L.ncs(P2, D2)
    cands = []
    for i in mins:
        A0 = As[np.argmin(C[i])]
        r = least_squares(resid, [ds[i], A0, 0.01, 0.0], args=(wl2, meas2, ang),
                          bounds=([0.5, 1.3, -0.05, -0.005], [800, 2.3, 0.08, 0.005]))
        cands.append((M.mse_ncs(r.fun, 4), [float(v) for v in r.x]))
    cands.sort()
    return sh, cands

if __name__ == '__main__':
    shs = [s for s in L.sheets() if s.split('_')[1].startswith('2')]
    with Pool(4) as pool:
        res = pool.map(job, shs)
    out = {}
    print('%-5s %8s %7s %7s %7s | %-24s' % ('sheet', 'd(nm)', 'n550', 'n633', 'MSE', 'runner-up branch d / MSE'))
    for sh, c in res:
        (mse, x), alt = c[0], (c[1] if len(c) > 1 else (np.nan, [np.nan]))
        n = lambda w: cauchy(np.array([w]), *x[1:])[0]
        out[sh] = dict(d=x[0], cauchy=x[1:], mse=mse, branches=c)
        print('%-5s %8.2f %7.4f %7.4f %7.2f | %8.1f / %6.1f' % (sh, x[0], n(550), n(633), mse, alt[1][0], alt[0]))
    json.dump(out, open(OUT, 'w'), indent=1)
