"""Which structure does layer2 need?  The alternatives that were tried, on one footing.

Transparent range only (570-1000 nm, where the perovskite does not absorb), a
Sellmeier perovskite with its resonance below the window, GraHIL n,k as on
layer1, N,C,S at three angles plus the reflectance from the intensity ratio.
Each model gets random starts over wide bounds; the best is printed.

  one layer on GraHIL 108        the nominal stack
  one layer, GraHIL free         GraHIL thinned / gone?
  one layer + height spread      roughness with a large lateral scale (coherent average)
  one layer + scatter            a homogeneous film that scatters
  split, lower lossless          dense top / low-index lower part, no loss
  split, absorption in GraHIL    the loss put into the GraHIL instead (dedoped PEDOT)
  split, scatter in lower        dense top / lower part = perovskite + void with a lambda^-p extinction
  split, organic matrix          the same with the lower part's host index free
  split, scatter in both         the dense part allowed to scatter too"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import gp_model as P, gp_perov as PV
import ce_osc as osc
from tmm_fix import tmm
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')
N_START = 80
PER = (['d_t', 'einf', 'S', 'lam0', 'd_r'], [40, 1, 0, 380, 0], [460, 4.5, 3.5, 535, 30], [150, 2.6, 1.3, 430, 4])
LOW = (['d_b', 'f_b'], [40, 0.15], [320, 1.0], [160, 0.5])
SCA = (['k_sc', 'p'], [0, 0], [0.5, 8], [0.25, 4])
def join(*parts): return tuple(sum((list(p[i]) for p in parts), []) for i in range(4))
KIND = {
    'one layer on GraHIL 108':     PER,
    'one layer, GraHIL free':      join(PER, (['d_g'], [0], [200], [60])),
    'one layer + height spread':   join(PER, (['sig_s'], [0], [60], [20])),
    'one layer + scatter':         join(PER, SCA),
    'split, lower lossless':       join(PER, LOW),
    'split, absorption in GraHIL': join(PER, LOW, (['A', 'E0', 'Br'], [0, 1.2, 0.15], [6, 2.6, 1.5], [0.5, 2.0, 0.7])),
    'split, scatter in lower':     join(PER, LOW, SCA),
    'split, organic matrix':       join(PER, LOW, SCA, (['n_m'], [1.0], [1.65], [1.4])),
    'split, scatter in both':      join(PER, LOW, SCA, (['k_top'], [0], [0.2], [0.01])),
}
_S = {}
def S():
    if not _S: _S.update(PV.setup(570.0, 1000.0, 3))
    return _S


def model(x, kind):
    s = S(); wl = s['wl']; v = dict(zip(KIND[kind][0], x))
    sc = (633.0 / wl)**v.get('p', 4.0)
    n = np.sqrt(v['einf'] + v['S'] * wl**2 / (wl**2 - v['lam0']**2) + 0j)
    g = s['Ng']
    if 'A' in v:
        g1, g2 = osc.gaussian(s['E'], v['A'], v['Br'], v['E0']); g = np.sqrt(s['Ng']**2 + g1 + 1j * g2); g = np.where(g.imag < 0, -g, g)
    if kind.startswith('one layer'):
        nt = n + 1j * v.get('k_sc', 0.0) * sc
        L = [s['air'], P.bruggeman(nt, s['air'], 0.5), nt, g, s['ox'], s['si']]; D = [v['d_r'], v['d_t'], v.get('d_g', PV.D_G), PV.D_OX]
    else:
        nt = n + 1j * v.get('k_top', 0.0) * sc
        nb = P.bruggeman(n, v.get('n_m', 1.0) * s['air'], v['f_b']) + 1j * v.get('k_sc', 0.0) * sc
        L = [s['air'], P.bruggeman(nt, s['air'], 0.5), nt, nb, g, s['ox'], s['si']]; D = [v['d_r'], v['d_t'], v['d_b'], PV.D_G, PV.D_OX]
    if v.get('sig_s', 0.0) > 0:
        rp, rs = P.rr_coh2(wl, s['ang'], L, D, 1, 2, v['sig_s'], 0.0, nodes=9)
    else:
        rr = [tmm(wl, L, D, th) for th in s['ang']]; rp = np.stack([r[0] for r in rr], 1); rs = np.stack([r[1] for r in rr], 1)
    return P.ncs_from_r(rp, rs), (np.abs(rp)**2 + np.abs(rs)**2) / 2


def job(a):
    kind, seed = a; names, lo, hi, base = KIND[kind]; lo, hi, base = np.array(lo, float), np.array(hi, float), np.array(base, float)
    rng = np.random.default_rng(seed); s = S()
    x0 = base * (1 + 0.6 * (rng.random(len(lo)) - 0.5)) if seed % 2 else lo + rng.random(len(lo)) * (hi - lo)
    def res(x):
        m, R = model(x, kind); return np.r_[P.resid(m, s['meas']), PV.WR * (R / 1.015 - s['Rm']).ravel()]
    try: r = least_squares(res, np.clip(x0, lo + 1e-6, hi - 1e-6), bounds=(lo, hi), x_scale='jac', max_nfev=250)
    except Exception: return kind, 1e9, 1e9, list(x0)
    n3 = 9 * len(s['wl'])
    return kind, float(P.mse_ncs(r.fun[:n3], len(lo))), float(1e3 * np.sqrt(np.mean((r.fun[n3:] / PV.WR)**2))), [float(t) for t in r.x]


def main():
    from multiprocessing import Pool
    jobs = [(k, s) for k in KIND for s in range(N_START if 'spread' not in k else N_START // 3)]
    with Pool(4) as pool: res = pool.map(job, jobs, chunksize=2)
    out = {}
    print('%-30s %7s %7s | best parameters' % ('model (570-1000 nm)', 'MSE', 'rms R'))
    for kind in KIND:
        r = min((r for r in res if r[0] == kind), key=lambda r: r[1] + r[2]); v = dict(zip(KIND[kind][0], r[3]))
        nq = [float(np.sqrt(v['einf'] + v['S'] * w**2 / (w**2 - v['lam0']**2))) for w in (633.0, 1000.0)]
        out[kind] = dict(mse=r[1], rms_R=r[2], x=v, n_633=nq[0], n_1000=nq[1])
        print('%-30s %7.2f %7.1f | n(633) %.2f n(1000) %.2f  ' % (kind, r[1], r[2], *nq) + '  '.join('%s=%.3g' % (k, v[k]) for k in KIND[kind][0] if k not in ('einf', 'S', 'lam0')))
    json.dump(out, open(W + 'explore.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
