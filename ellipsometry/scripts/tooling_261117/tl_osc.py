"""KK-consistent dispersion per material from the thick Si pieces.

film = Einf + 4 Gaussians (CompleteEASE convention, analytic eps1), fitted to the
raw N,C,S of every thick Si piece of that material at once, one shared
dispersion and one thickness per piece.  Started from the model-free k of
tl_pointnk so no band has to be guessed."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # tooling.xlsx
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fits, tables, figures
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L
import jnk_ref as R, jnk_model as M
from scipy.optimize import least_squares
W = _os.path.join(ELLIPS_OUT, '')
HC = 1239.841984
D_OX = 2.0
C = json.load(open(W + 'cauchy_si.json')); PN = json.load(open(W + 'pointnk_si.json'))
MAT = {'NPB': ['1_22', '1_21'], 'TAPC': ['1_24'], 'Bphen': ['2_21', '2_22']}
NG = {'NPB': 7, 'TAPC': 6, 'Bphen': 6}

def N_of(q, wl):
    """q = [Einf, (A, Br, En) x any number of Gaussians]"""
    ng = (len(q) - 1) // 3
    return M.osc_N(q[0], [('gauss', q[1 + 3 * i], q[2 + 3 * i], q[3 + 3 * i]) for i in range(ng)], wl)

def seed_from_pointnk(sh, ng):
    """Fit the Gaussians to the model-free eps2 first - a well-conditioned start.

    ng - 1 absorption Gaussians are spread evenly from the absorption onset to
    the short-wavelength end of the data; the last one sits above the data and
    stands in for everything deeper in the UV (it shapes n, not k)."""
    wl = np.array(PN[sh]['wl']); n = np.array(PN[sh]['n']); k = np.array(PN[sh]['k'])
    m = wl <= 700; wl, n, k = wl[m], n[m], k[m]
    E = HC / wl; e2 = 2 * n * k; e1 = n**2 - k**2
    on = E[e2 > 0.05 * e2.max()].min()
    cen = np.linspace(on + 0.1, E.max() - 0.05, ng - 1)
    q0 = [2.0]
    for c in cen:
        q0 += [max(np.interp(c, E[::-1], e2[::-1]), 0.05), 0.35, c]
    q0 += [2.0, 1.5, 6.5]
    lo = [1.0] + [0, 0.06, 2.6] * (ng - 1) + [0, 0.3, 5.0]
    hi = [4.0] + [20, 1.5, 5.3] * (ng - 1) + [40, 4.0, 9.0]
    def f(q):
        N = N_of(q, wl)
        return np.concatenate([(N**2).imag - e2, 0.5 * ((N**2).real - e1)])
    return least_squares(f, np.clip(q0, np.array(lo) + 1e-6, np.array(hi) - 1e-6), bounds=(lo, hi)).x, lo, hi

def fit_material(name, shs, lo_wl=245.0, hi_wl=1000.0):
    q0, lo, hi = seed_from_pointnk(shs[0], NG[name])
    data = []
    for sh in shs:
        wl, P, D, ang = L.se(sh, lo_wl, hi_wl, step=2)
        data.append((wl, L.ncs(P, D), ang))
    nq = len(q0)
    def res(x):
        q = x[:nq]; out = []
        for j, (wl, meas, ang) in enumerate(data):
            Nf = N_of(q, wl)
            n_, c_, s_ = M.se_ncs(wl, [R.n_air(wl), Nf, R.n_sio2(wl), R.n_si(wl)], [x[nq + j], D_OX], ang)
            out += [(n_ - meas[0]).ravel(), (c_ - meas[1]).ravel(), (s_ - meas[2]).ravel()]
        return np.concatenate(out)
    x0 = np.concatenate([q0, [C[sh]['d'] for sh in shs]])
    LOx = np.concatenate([lo, [1.0] * len(shs)]); HIx = np.concatenate([hi, [200.0] * len(shs)])
    rng = np.random.default_rng(11); best = None
    for k in range(6):
        s = x0 if k == 0 else np.clip(x0 * (1 + 0.15 * rng.standard_normal(len(x0))), LOx + 1e-6, HIx - 1e-6)
        r = least_squares(res, np.clip(s, LOx + 1e-6, HIx - 1e-6), bounds=(LOx, HIx), x_scale='jac', max_nfev=3000)
        if best is None or r.cost < best.cost: best = r
    x = best.x
    per = []
    for j, (wl, meas, ang) in enumerate(data):
        Nf = N_of(x[:nq], wl)
        n_, c_, s_ = M.se_ncs(wl, [R.n_air(wl), Nf, R.n_sio2(wl), R.n_si(wl)], [x[nq + j], D_OX], ang)
        rr = np.concatenate([(n_ - meas[0]).ravel(), (c_ - meas[1]).ravel(), (s_ - meas[2]).ravel()])
        per.append(M.mse_ncs(rr, nq // len(shs) + 1))
    return x[:nq], x[nq:], per

def main():
    out = {}
    wp = np.array([300., 350., 400., 450., 550., 633.])
    from multiprocessing import Pool
    with Pool(3) as pool:
        fits = pool.starmap(fit_material, list(MAT.items()))
    for (name, shs), (q, ds, per) in zip(MAT.items(), fits):
        N = N_of(q, wp)
        out[name] = dict(q=[float(v) for v in q], sheets=shs, d=[float(v) for v in ds], mse=[float(v) for v in per])
        print('%-6s %s' % (name, '  '.join('%s: d=%.2f nm MSE %.2f' % (sh, d, m) for sh, d, m in zip(shs, ds, per))))
        print('       ' + '  '.join('%3.0f: %.3f+%.3fj' % (w, z.real, z.imag) for w, z in zip(wp, N)))
        print('       Gaussians (A, Br eV, En eV -> nm): ' + '; '.join('%.2f %.2f %.2f(%.0f)' % (q[1+3*i], q[2+3*i], q[3+3*i], HC/q[3+3*i]) for i in range((len(q)-1)//3)) + '  Einf %.3f' % q[0])
    json.dump(out, open(W + 'osc.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
