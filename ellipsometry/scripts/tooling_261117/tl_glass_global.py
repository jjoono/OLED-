"""All 8 glass pieces at once: one glass (2-term Cauchy), one backside collection
fraction f, one thickness per piece, film dispersion fixed to the material each
piece was identified as.  Three readings of how a depolarizing sample's export
maps to N,C,S are tried and the data picks."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # tooling.xlsx
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fits, tables, figures
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_osc as O, tl_glass as G
import jnk_model as M
from scipy.optimize import least_squares
W = _os.path.join(ELLIPS_OUT, '')
OSC = json.load(open(W + 'osc.json'))
ID = {'1_31': 'NPB', '1_32': 'NPB', '1_33': 'NPB', '1_34': 'TAPC',
      '2_31': 'Bphen', '2_32': 'Bphen', '2_33': 'Bphen', '2_34': 'Bphen'}
D0 = {'1_31': 26, '1_32': 46, '1_33': 9, '1_34': 47, '2_31': 44, '2_32': 44, '2_33': 1, '2_34': 3}
SH = list(ID)

def norm(N, C, S):
    P = np.sqrt(N**2 + C**2 + S**2); return N / P, C / P, S / P

CONV = {'psi-from-N': G.as_reported, 'normalised': norm, 'raw NCS': lambda N, C, S: (N, C, S)}

def run(conv, lo_wl=300.0, hi_wl=1000.0):
    data = []
    for sh in SH:
        wl, P, D, ang = L.se(sh, lo_wl, hi_wl, 3)
        data.append((wl, L.ncs(P, D), ang, O.N_of(np.array(OSC[ID[sh]]['q']), wl)))
    cv = CONV[conv]
    def per_piece(x, j):
        A, B, f = x[:3]; wl, meas, ang, Nf = data[j]
        ng = A + B / (wl / 1000)**2 + 0j
        N, C, S = cv(*G.ncs_glass(wl, Nf, x[3 + j], ng, ang, f))
        return np.concatenate([(N - meas[0]).ravel(), (C - meas[1]).ravel(), (S - meas[2]).ravel()])
    res = lambda x: np.concatenate([per_piece(x, j) for j in range(len(SH))])
    x0 = [1.51, 0.005, 0.4] + [D0[s] for s in SH]
    lo = [1.40, -0.01, 0.0] + [0.0] * len(SH); hi = [1.70, 0.03, 1.0] + [150.0] * len(SH)
    best = None
    for f0 in (0.0, 0.3, 0.7):
        x0[2] = f0
        r = least_squares(res, np.clip(x0, np.array(lo) + 1e-7, np.array(hi) - 1e-7), bounds=(lo, hi), x_scale='jac', max_nfev=400)
        if best is None or r.cost < best.cost: best = r
    x = best.x
    per = [M.mse_ncs(per_piece(x, j), 1 + 3 / len(SH)) for j in range(len(SH))]
    return conv, [float(v) for v in x], [float(v) for v in per], float(M.mse_ncs(best.fun, len(x)))

if __name__ == '__main__':
    from multiprocessing import Pool
    with Pool(3) as p:
        out = p.map(run, list(CONV))
    js = {}
    for conv, x, per, tot in out:
        js[conv] = dict(x=x, per=per, total=tot)
        print('%-11s total MSE %6.2f | glass n(550) %.4f  f %.3f' % (conv, tot, x[0] + x[1] / 0.55**2, x[2]))
        print('            ' + '  '.join('%s %s %.2f nm (%.1f)' % (s, ID[s], d, m) for s, d, m in zip(SH, x[3:], per)))
    json.dump(js, open(W + 'glass_global.json', 'w'), indent=1)
