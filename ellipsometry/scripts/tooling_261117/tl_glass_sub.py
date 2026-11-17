"""What does the bare-looking glass piece 2_33 say about the substrate model?
front-only vs incoherent backside (fraction f), 2- vs 3-term Cauchy."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # tooling.xlsx
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fits, tables, figures
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_osc as O, tl_glass as G
import jnk_ref as R, jnk_model as M
from scipy.optimize import least_squares
W = _os.path.join(ELLIPS_OUT, '')
OSC = json.load(open(W + 'osc.json'))

def run(sh, mat, back, three, lo_wl=245.0):
    wl, P, D, ang = L.se(sh, lo_wl, 1000, 2); meas = L.ncs(P, D)
    Nf = O.N_of(np.array(OSC[mat]['q']), wl); l = wl / 1000
    def res(x):
        d, A, B = x[:3]; C = x[3] if three else 0.0; f = x[-1] if back else 0.0
        ng = A + B / l**2 + C / l**4 + 0j
        N, C_, S_ = G.as_reported(*G.ncs_glass(wl, Nf, d, ng, ang, f))
        return np.concatenate([(N - meas[0]).ravel(), (C_ - meas[1]).ravel(), (S_ - meas[2]).ravel()])
    x0 = [1.0, 1.52, 0.004] + ([0.0] if three else []) + ([0.5] if back else [])
    lo = [0.0, 1.40, -0.01] + ([-0.002] if three else []) + ([0.0] if back else [])
    hi = [100., 1.70, 0.03] + ([0.002] if three else []) + ([1.0] if back else [])
    best = None
    for k in range(10):
        s = np.array(x0) + (np.random.default_rng(k).standard_normal(len(x0)) * (np.array(hi) - np.array(lo)) * 0.1 if k else 0)
        r = least_squares(res, np.clip(s, np.array(lo) + 1e-7, np.array(hi) - 1e-7), bounds=(lo, hi))
        if best is None or r.cost < best.cost: best = r
    r = best.fun.reshape(3, len(wl), len(ang))
    band = lambda a, b: 1e3 * np.sqrt(np.mean(r[:, (wl >= a) & (wl < b), :]**2))
    return best.x, M.mse_ncs(best.fun, len(x0)), band(245, 350), band(350, 500), band(500, 1000)

for sh, mat in (('2_33', 'Bphen'), ('2_31', 'Bphen'), ('1_34', 'TAPC')):
    for back in (False, True):
        for three in (False, True):
            x, mse, b1, b2, b3 = run(sh, mat, back, three)
            print('%-5s %-5s backside=%-5s 3-term=%-5s MSE %6.2f | 245-350 %5.1f 350-500 %5.1f 500-1000 %5.1f | x = %s'
                  % (sh, mat, back, three, mse, b1, b2, b3, ', '.join('%.4f' % v for v in x)), flush=True)
    print()
