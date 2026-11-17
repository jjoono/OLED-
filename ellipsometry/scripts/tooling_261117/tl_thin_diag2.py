"""2_23/2_24 round 2: angle offset, a bare oxide, and a free absorbing Cauchy."""
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

def run(sh, mode):
    wl, P, D, ang = L.se(sh, 245, 1000, 2); meas = L.ncs(P, D)
    air, ox, si = R.n_air(wl), R.n_sio2(wl), R.n_si(wl)
    NB = O.N_of(np.array(OSC['Bphen']['q']), wl); NN = O.N_of(np.array(OSC['NPB']['q']), wl)
    l = wl / 1000
    if mode == 'oxide only':                    # d_ox, dth
        f = lambda x: ([air, ox, si], [x[0]], x[1]); x0, lo, hi = [3, 0], [0, -2], [60, 2]
    elif mode == 'Bphen + angle offset':        # d, dth
        f = lambda x: ([air, NB, ox, si], [x[0], 2.0], x[1]); x0, lo, hi = [8, 0], [0.2, -2], [100, 2]
    elif mode == 'NPB + angle offset':
        f = lambda x: ([air, NN, ox, si], [x[0], 2.0], x[1]); x0, lo, hi = [8, 0], [0.2, -2], [100, 2]
    elif mode == 'Cauchy+k (no offset)':        # d, A, B, k
        f = lambda x: ([air, (x[1] + x[2] / l**2) + 1j * x[3] * np.ones(len(wl)), ox, si], [x[0], 2.0], 0.0)
        x0, lo, hi = [8, 1.8, 0.01, 0.01], [0.2, 1.2, -0.05, 0], [100, 4.0, 0.2, 2]
    elif mode == 'Cauchy+k + offset':
        f = lambda x: ([air, (x[1] + x[2] / l**2) + 1j * x[3] * np.ones(len(wl)), ox, si], [x[0], 2.0], x[4])
        x0, lo, hi = [8, 1.8, 0.01, 0.01, 0], [0.2, 1.2, -0.05, 0, -2], [100, 4.0, 0.2, 2, 2]
    def res(x):
        ls, ds, dth = f(x)
        n_, c_, s_ = M.se_ncs(wl, ls, ds, ang + dth)
        return np.concatenate([(n_ - meas[0]).ravel(), (c_ - meas[1]).ravel(), (s_ - meas[2]).ravel()])
    best = None
    for k in range(16):
        s = np.array(x0, float) + (np.random.default_rng(k).standard_normal(len(x0)) * (np.array(hi) - np.array(lo)) * 0.15 if k else 0)
        r = least_squares(res, np.clip(s, np.array(lo) + 1e-6, np.array(hi) - 1e-6), bounds=(lo, hi))
        if best is None or r.cost < best.cost: best = r
    r = best.fun.reshape(3, len(wl), len(ang))
    band = lambda a, b: 1e3 * np.sqrt(np.mean(r[:, (wl >= a) & (wl < b), :]**2))
    return best.x, M.mse_ncs(best.fun, len(x0)), band(245, 350), band(350, 500), band(500, 1000)

for sh in ('2_23', '2_24', '2_21', '1_23'):
    for mode in ('oxide only', 'Bphen + angle offset', 'NPB + angle offset', 'Cauchy+k (no offset)', 'Cauchy+k + offset'):
        x, mse, b1, b2, b3 = run(sh, mode)
        print('%-5s %-22s MSE %6.2f | 245-350 %5.1f 350-500 %5.1f 500-1000 %5.1f | x = %s'
              % (sh, mode, mse, b1, b2, b3, ', '.join('%.3f' % v for v in x)), flush=True)
    print()
