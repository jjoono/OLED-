"""Is the k ~ 0.01 visible tail of TAPC / Bphen required by the data?
Refit with the onset Gaussian (the one below 3.1 eV) forced to zero and compare."""
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
out = {}
for name in ('TAPC', 'Bphen'):
    q0 = np.array(OSC[name]['q']); shs = OSC[name]['sheets']; ds0 = OSC[name]['d']
    ng = (len(q0) - 1) // 3
    kill = [i for i in range(ng) if q0[3 + 3 * i] < 3.1]
    data = [(lambda t: (t[0], L.ncs(t[1], t[2]), t[3]))(L.se(sh, 245, 1000, 2)) for sh in shs]
    def res(x):
        q = x[:len(q0)].copy()
        for i in kill: q[1 + 3 * i] = 0.0
        o = []
        for j, (wl, meas, ang) in enumerate(data):
            n_, c_, s_ = M.se_ncs(wl, [R.n_air(wl), O.N_of(q, wl), R.n_sio2(wl), R.n_si(wl)], [x[len(q0) + j], 2.0], ang)
            o += [(n_ - meas[0]).ravel(), (c_ - meas[1]).ravel(), (s_ - meas[2]).ravel()]
        return np.concatenate(o)
    x0 = np.concatenate([q0, ds0])
    r = least_squares(res, x0, x_scale='jac', max_nfev=3000)
    q = r.x[:len(q0)].copy()
    for i in kill: q[1 + 3 * i] = 0.0
    base = np.concatenate([q0, ds0])
    f0 = lambda x: res(x)
    mse_with = M.mse_ncs(f0(base) if not kill else np.concatenate([
        (lambda wl, meas, ang, d: np.concatenate([(a - b).ravel() for a, b in zip(M.se_ncs(wl, [R.n_air(wl), O.N_of(q0, wl), R.n_sio2(wl), R.n_si(wl)], [d, 2.0], ang), meas)]))(*data[j], ds0[j])
        for j in range(len(shs))]), len(base))
    mse_without = M.mse_ncs(r.fun, len(base) - 3 * len(kill))
    N0 = O.N_of(q0, np.array([450., 550.])); N1 = O.N_of(q, np.array([450., 550.]))
    print('%-6s with tail MSE %.2f  d %s  k450 %.4f n550 %.4f  |  tail removed MSE %.2f  d %s  k450 %.4f n550 %.4f'
          % (name, mse_with, ', '.join('%.2f' % d for d in ds0), N0[0].imag, N0[1].real,
             mse_without, ', '.join('%.2f' % d for d in r.x[len(q0):]), N1[0].imag, N1[1].real))
    out[name] = dict(q=[float(v) for v in q], d=[float(v) for v in r.x[len(q0):]], mse=float(mse_without), sheets=shs)
json.dump(out, open(W + 'osc_notail.json', 'w'), indent=1)
