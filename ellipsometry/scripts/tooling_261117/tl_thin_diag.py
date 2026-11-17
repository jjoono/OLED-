"""Why do 2_23 / 2_24 not fit a 10 nm organic film on 2 nm oxide?"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # tooling.xlsx
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fits, tables, figures
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_osc as O
import jnk_ref as R, jnk_model as M
from scipy.optimize import least_squares
W = _os.path.join(ELLIPS_OUT, '')
OSC = json.load(open(W + 'osc.json'))

def run(sh, mat, mode):
    wl, P, D, ang = L.se(sh, 245, 1000, 2); meas = L.ncs(P, D)
    Nf = O.N_of(np.array(OSC[mat]['q']), wl); air = R.n_air(wl)
    ox, si = R.n_sio2(wl), R.n_si(wl)
    def stack(x):
        if mode == 'oxide free':          # x = d_film, d_ox
            return [air, Nf, ox, si], [x[0], x[1]]
        if mode == 'island (EMA void)':  # x = d_film, film fraction
            return [air, M.bruggeman(Nf, air, x[1]), ox, si], [x[0], 2.0]
        if mode == 'oxide + island':
            return [air, M.bruggeman(Nf, air, x[2]), ox, si], [x[0], x[1]]
        return [air, Nf, ox, si], [x[0], 2.0]
    x0, lo, hi = {'fixed oxide': ([8.0], [0.2], [100]),
                  'oxide free': ([8.0, 2.0], [0.2, 0.0], [100, 30]),
                  'island (EMA void)': ([10.0, 0.7], [0.2, 0.05], [100, 1.0]),
                  'oxide + island': ([10.0, 2.0, 0.7], [0.2, 0.0, 0.05], [100, 30, 1.0])}[mode]
    def res(x):
        ls, ds = stack(x)
        n_, c_, s_ = M.se_ncs(wl, ls, ds, ang)
        return np.concatenate([(n_ - meas[0]).ravel(), (c_ - meas[1]).ravel(), (s_ - meas[2]).ravel()])
    best = None
    for k in range(12):
        s = np.array(x0) * (1 + (0.4 * np.random.default_rng(k).standard_normal(len(x0)) if k else 0))
        r = least_squares(res, np.clip(s, np.array(lo) + 1e-6, np.array(hi) - 1e-6), bounds=(lo, hi))
        if best is None or r.cost < best.cost: best = r
    r = best.fun.reshape(3, len(wl), len(ang))
    band = lambda a, b: 1e3 * np.sqrt(np.mean(r[:, (wl >= a) & (wl < b), :]**2))
    return best.x, M.mse_ncs(best.fun, len(x0)), band(245, 350), band(350, 500), band(500, 1000)

for sh in ('2_23', '2_24', '1_23'):
    for mat in ('Bphen', 'NPB'):
        for mode in ('fixed oxide', 'oxide free', 'island (EMA void)', 'oxide + island'):
            x, mse, b1, b2, b3 = run(sh, mat, mode)
            print('%-5s %-6s %-18s MSE %6.2f  | 245-350 %5.1f 350-500 %5.1f 500-1000 %5.1f | x = %s'
                  % (sh, mat, mode, mse, b1, b2, b3, ', '.join('%.3f' % v for v in x)), flush=True)
    print()
