"""Gen-Osc dispersions with parameters an organic film can actually have.

    eps(E) = Einf + A_p / (E_p^2 - E^2) + sum_i Gaussian(A_i, Br_i, En_i)

  * Einf fixed at 1: the transitions above the data are carried by the UV pole,
    so nothing is left for Einf to absorb.  A free Einf trades with the pole
    and ends up meaning nothing.
  * UV pole above the data, 6-12 eV, amplitude >= 0 (CompleteEASE's 'UV Pole').
  * Gaussians only for bands the data shows: centre no lower than the
    absorption onset of the model-free k, amplitude >= 0 (no negative
    oscillators cancelling each other), FWHM 0.10-1.20 eV (vibronic bands of
    small-molecule films).  Nothing below the gap, so k = 0 in the visible by
    construction - these are wide-gap transport materials.

The number of Gaussians is chosen by the data: the smallest count whose MSE is
within 10 % of the best.  Each parameter gets a 90 % confidence interval from
the Jacobian, and anything wider than 50 % of its value is flagged for locking,
the same rule make_mod_v13 applied to the ITO/IZO set."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L
import jnk_ref as R, jnk_model as M, ce_osc as osc
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')
HC = 1239.841984
D_OX = 2.0
EINF = 1.0
MAT = {'NPB': ['1_22', '1_21'], 'TAPC': ['1_24'], 'Bphen': ['2_21', '2_22']}
ONSET_NM = {'NPB': 409.0, 'TAPC': 358.0, 'Bphen': 347.0}     # k > 0.02 in tl_pointnk
NGS = {'NPB': (3, 4, 5, 6), 'TAPC': (2, 3, 4, 5), 'Bphen': (2, 3, 4, 5)}


def eps_genosc(p, wl):
    """p = [A_pole, E_pole, (A, Br, En) x ng]"""
    E = HC / np.asarray(wl, float)
    e1 = EINF + p[0] / (p[1]**2 - E**2); e2 = np.zeros_like(E)
    for i in range((len(p) - 2) // 3):
        a, b = osc.gaussian(E, p[2 + 3 * i], p[3 + 3 * i], p[4 + 3 * i]); e1 += a; e2 += b
    return e1 + 1j * e2


def N_of(p, wl):
    N = np.sqrt(eps_genosc(p, wl).astype(complex))
    return np.where(N.imag < 0, -N, N)


def bounds(name, ng):
    eon = HC / ONSET_NM[name]
    lo = [0.0, 6.0] + [0.0, 0.10, eon - 0.05] * ng
    hi = [500.0, 12.0] + [15.0, 1.20, 5.60] * ng
    return np.array(lo), np.array(hi)


def data_for(shs, step=2):
    out = []
    for sh in shs:
        wl, P, D, ang = L.se(sh, 245.0, 1000.0, step)
        out.append((wl, L.ncs(P, D), ang, R.n_air(wl), R.n_sio2(wl), R.n_si(wl)))
    return out


def residual(x, data, npar):
    p = x[:npar]; out = []
    for j, (wl, meas, ang, air, ox, si) in enumerate(data):
        n_, c_, s_ = M.se_ncs(wl, [air, N_of(p, wl), ox, si], [x[npar + j], D_OX], ang)
        out += [(n_ - meas[0]).ravel(), (c_ - meas[1]).ravel(), (s_ - meas[2]).ravel()]
    return np.concatenate(out)


def starts(name, ng, rng, n=10):
    """Centres spread from the onset into the UV edge, plus random draws."""
    eon = HC / ONSET_NM[name]
    base = np.linspace(eon + 0.15, 5.3, ng)
    s = []
    for k in range(n):
        cen = base if k == 0 else np.sort(rng.uniform(eon, 5.5, ng))
        p = [80.0, 9.0]
        for c in cen:
            p += [1.0, 0.45, c]
        s.append(np.array(p))
    return s


def fit(name, ng, d0):
    shs = MAT[name]; data = data_for(shs)
    lo, hi = bounds(name, ng); npar = len(lo)
    LO = np.concatenate([lo, [1.0] * len(shs)]); HI = np.concatenate([hi, [200.0] * len(shs)])
    rng = np.random.default_rng(zlib_seed(name, ng)); best = None
    for p0 in starts(name, ng, rng):
        x0 = np.clip(np.concatenate([p0, d0]), LO + 1e-6, HI - 1e-6)
        try:
            r = least_squares(residual, x0, args=(data, npar), bounds=(LO, HI), x_scale='jac', max_nfev=2500)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    x = best.x; res = best.fun
    J = best.jac; dof = max(res.size - x.size, 1); s2 = np.sum(res**2) / dof
    cov = np.linalg.pinv(J.T @ J) * s2
    ci = 1.645 * np.sqrt(np.clip(np.diag(cov), 0, None))
    per = []
    for j in range(len(shs)):
        rj = residual(x, [data[j]], npar) if False else None
    # per-piece MSE
    for j in range(len(shs)):
        xx = np.concatenate([x[:npar], [x[npar + j]]])
        per.append(float(M.mse_ncs(residual(xx, [data[j]], npar), npar / len(shs) + 1)))
    at_bound = [i for i in range(npar) if min(abs(x[i] - LO[i]), abs(HI[i] - x[i])) < 1e-4 * max(1, abs(HI[i]))]
    return dict(name=name, ng=ng, p=[float(v) for v in x[:npar]], d=[float(v) for v in x[npar:]],
                ci=[float(v) for v in ci], mse=float(M.mse_ncs(res, x.size)), mse_per=per,
                at_bound=at_bound, sheets=shs)


def zlib_seed(*k):
    import zlib
    return zlib.crc32(repr(k).encode())


def label(i):
    if i == 0: return 'UV Pole Amp'
    if i == 1: return 'UV Pole En'
    g, r = divmod(i - 2, 3)
    return 'G%d %s' % (g + 1, ('Amp', 'Br', 'En')[r])


def main():
    from multiprocessing import Pool
    CS = json.load(open(W + 'cauchy_si.json'))
    jobs = [(name, ng, np.array([CS[s]['d'] for s in MAT[name]])) for name in MAT for ng in NGS[name]]
    with Pool(4) as pool:
        res = pool.starmap(fit, jobs)
    out = {}
    for name in MAT:
        rs = [r for r in res if r['name'] == name]
        best = min(r['mse'] for r in rs)
        pick = min((r for r in rs if r['mse'] <= 1.10 * best), key=lambda r: r['ng'])
        print('%-6s ' % name + '  '.join('ng=%d MSE %.2f%s' % (r['ng'], r['mse'], ' <' if r is pick else '') for r in rs))
        p, ci = np.array(pick['p']), np.array(pick['ci'])
        flags = [label(i) for i in range(len(p)) if p[i] != 0 and ci[i] > 0.5 * abs(p[i])]
        print('       d %s   per-piece MSE %s' % (', '.join('%.2f' % d for d in pick['d']),
                                                   ', '.join('%.2f' % m for m in pick['mse_per'])))
        print('       UV pole %.1f eV^2 at %.2f eV (+-%.1f, +-%.2f)' % (p[0], p[1], ci[0], ci[1]))
        for g in range((len(p) - 2) // 3):
            a, b, e = p[2 + 3 * g: 5 + 3 * g]; ca, cb, ce = ci[2 + 3 * g: 5 + 3 * g]
            print('       G%d  Amp %6.3f +-%5.3f   Br %5.3f +-%5.3f eV   En %5.3f +-%5.3f eV (%3.0f nm)'
                  % (g + 1, a, ca, b, cb, e, ce, HC / e))
        wp = np.array([350., 450., 550., 633.]); N = N_of(p, wp)
        print('       ' + '  '.join('%3.0f: %.3f+%.3fj' % (w, z.real, z.imag) for w, z in zip(wp, N)))
        print('       wide CI (>50 %%): %s   at a bound: %s' % (', '.join(flags) or '-',
              ', '.join(label(i) if i < len(p) else 'd' for i in pick['at_bound']) or '-'))
        out[name] = dict(pick=pick, all=rs, einf=EINF)
    json.dump(out, open(W + 'genosc.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
