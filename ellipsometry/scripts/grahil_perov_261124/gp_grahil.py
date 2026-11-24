"""layer1: GraHIL (PEDOT:PSS : PFI 1:1 w/w, self-organised, PFI-rich at the top)
on Si.  The paper (Nature 611, 688) quotes 75 nm and assumes n = 1.5.

Dispersion: KK-consistent Gaussian basis (gp_model.KKBasis), 200-1000 nm.
Variants, each with the thickness free:
  iso            one homogeneous isotropic layer
  graded         the same material with (eps - 1) running linearly through the
                 depth, top = (1 - g) x, bottom = (1 + g) x mean  (g > 0: the
                 top is optically thinner, as a PFI-rich surface would be)
  uniaxial       eps_e = 1 + s (eps_o - 1): one scale for the out-of-plane axis
  graded+uni     both
  +nonuniform    lateral thickness spread, fitted against N,C,S x sqrt(1 - depol)
Native oxide fixed (OX nm); its effect is scanned at the end.
Beyond 1000 nm the Si wafer is transparent and the detector changes: not fitted."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import gp_load as G, gp_model as P
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')
LO, HI, STEP = 200.0, 1000.0, 2
CEN, WID = P.grid(1.0, 7.4, 0.3, fine=(4.4, 7.4, 0.15))
E_UV = 12.0      # with eps_inf free: a fluorinated ionomer has almost no dispersion, eps_inf = 1 + one pole cannot be that flat
SMOOTH = 0.002
M_GRADE = 8


def setup(lo=LO, hi=HI, step=STEP):
    d = G.load('layer1', lo, hi, step); wl = d['wl']
    return dict(d=d, wl=wl, ang=d['ang'], kb=P.KKBasis(wl, CEN, WID, E_UV), air=P.n_air(wl), ox=P.n_ox(wl), si=P.n_si(wl))


def ncs_uni(wl, ang, eo, ee, ds, sig=0.0):
    """uniaxial stack (axis along z); eo/ee per-layer eps lists air..Si; optional thickness spread."""
    x, w = P._GH[5 if sig > 0 else 1]; out = []
    for th in ang:
        Ip = Is = 0.0; X = 0.0
        for xi, wi in zip(x, w):
            dd = [dj * (1 + sig * xi) for dj in ds[:-1]] + [ds[-1]]          # films vary, the oxide does not
            rp, rs = P.tmm_uni(wl, eo, ee, dd, th)
            Ip = Ip + wi * np.abs(rp)**2; Is = Is + wi * np.abs(rs)**2; X = X + wi * np.conj(rp) * rs
        I = Ip + Is
        out.append(((Is - Ip) / I, 2 * X.real / I, 2 * X.imag / I))
    return tuple(np.stack([o[k] for o in out], 1) for k in range(3))


def unpack(x, variant, nb):
    v = dict(d=x[0], a_uv=x[1], einf=x[2], a=x[3:3 + nb]); k = 3 + nb
    v['g'] = x[k] if 'graded' in variant else 0.0; k += 'graded' in variant
    v['s'] = x[k] if 'uni' in variant else 1.0; k += 'uni' in variant
    v['sig'] = x[k] if 'nonuniform' in variant else 0.0
    return v


def model(x, variant, S, dox):
    v = unpack(x, variant, S['kb'].n); wl = S['wl']
    chi = S['kb'].eps(v['a'], v['a_uv'], v['einf']) - 1.0          # mean susceptibility
    m = M_GRADE if v['g'] != 0 else 1
    eo, ee, ds = [S['air']**2], [S['air']**2], []
    for i in range(m):                                             # top first
        f = 1.0 + v['g'] * (2 * (i + 0.5) / m - 1.0)               # top: 1 - g ... bottom: 1 + g
        eo.append(1 + f * chi); ee.append(1 + v['s'] * f * chi); ds.append(v['d'] / m)
    eo += [S['ox']**2, S['si']**2]; ee += [S['ox']**2, S['si']**2]; ds.append(dox)
    return ncs_uni(wl, S['ang'], eo, ee, ds, v['sig'])


def fit(variant, dox=1.5, x_from=None, lo=LO, hi=HI):
    S = setup(lo, hi); nb = S['kb'].n
    meas = P.data_ncs(S['d'], with_dep='nonuniform' in variant)
    D2 = np.diff(np.eye(nb), 2, axis=0)
    names = ['d', 'a_uv', 'einf'] + ['a%d' % i for i in range(nb)] + (['g'] if 'graded' in variant else []) + (['s'] if 'uni' in variant else []) + (['sig'] if 'nonuniform' in variant else [])
    lo_b = [20.0, 0.0, 1.0] + [0.0] * nb + ([-0.6] if 'graded' in variant else []) + ([0.5] if 'uni' in variant else []) + ([0.0] if 'nonuniform' in variant else [])
    hi_b = [300.0, 300.0, 2.5] + [20.0] * nb + ([0.6] if 'graded' in variant else []) + ([1.5] if 'uni' in variant else []) + ([0.2] if 'nonuniform' in variant else [])
    def res(x):
        return np.r_[P.resid(model(x, variant, S, dox), meas), SMOOTH * (D2 @ x[3:3 + nb])]
    best = None
    starts = [x_from] if x_from is not None else [None]
    for d0 in (107.0, 90.0, 125.0):
        x0 = np.r_[d0, 40.0, 1.6, np.full(nb, 0.01), ([0.0] if 'graded' in variant else []), ([1.0] if 'uni' in variant else []), ([0.02] if 'nonuniform' in variant else [])]
        if x_from is not None:
            x0[:len(x_from[:3 + nb])] = x_from[:3 + nb]
        x0 = np.clip(x0, np.array(lo_b) + 1e-9, np.array(hi_b) - 1e-9)
        r = least_squares(res, x0, bounds=(lo_b, hi_b), x_scale='jac', max_nfev=300)
        if best is None or r.cost < best.cost: best = r
        if x_from is not None: break
    x = best.x; v = unpack(x, variant, nb)
    rd = P.resid(model(x, variant, S, dox), meas)
    J = best.jac; s2 = np.sum(best.fun**2) / max(best.fun.size - x.size, 1)
    ci = 1.645 * np.sqrt(np.clip(np.diag(np.linalg.pinv(J.T @ J) * s2), 0, None))
    rr = rd.reshape(3, len(S['wl']), 3); wl = S['wl']
    band = lambda a, b: float(1e3 * np.sqrt(np.mean(rr[:, (wl >= a) & (wl < b), :]**2)))
    q = np.array([300.0, 450.0, 550.0, 633.0, 800.0]); kq = S['kb'].on(q)
    No = kq.N(v['a'], v['a_uv'], v['einf']); chi = kq.eps(v['a'], v['a_uv'], v['einf']) - 1
    return dict(variant=variant, dox=dox, x=[float(t) for t in x], names=names, d=float(v['d']), d_ci=float(ci[0]),
                g=float(v['g']), s=float(v['s']), sig=float(v['sig']),
                ci={n: float(c) for n, c in zip(names, ci) if n in ('d', 'g', 's', 'sig')},
                mse=float(P.mse_ncs(rd, len(x))), bands={'200-300': band(200, 300), '300-450': band(300, 450), '450-1000': band(450, 1001)},
                n_o={str(int(w)): [float(z.real), float(z.imag)] for w, z in zip(q, No)},
                n_top_bot_550=[float(np.sqrt(1 + (1 - v['g']) * chi[2]).real), float(np.sqrt(1 + (1 + v['g']) * chi[2]).real)],
                n_e_550=float(np.sqrt(1 + v['s'] * chi[2]).real))


def export(r, path):
    """mean ordinary n,k on a 1 nm grid (plus top / bottom of the grade and the extraordinary axis)"""
    wl = np.arange(200.0, 1001.0, 1.0); kb = P.KKBasis(wl, CEN, WID, E_UV)
    v = unpack(np.array(r['x']), r['variant'], kb.n); chi = kb.eps(v['a'], v['a_uv'], v['einf']) - 1
    sq = lambda e: (lambda n: np.where(n.imag < 0, -n, n))(np.sqrt(e))
    cols = [wl, sq(1 + chi).real, sq(1 + chi).imag, sq(1 + (1 - v['g']) * chi).real, sq(1 + (1 + v['g']) * chi).real, sq(1 + v['s'] * chi).real, sq(1 + v['s'] * chi).imag]
    np.savetxt(path, np.c_[tuple(cols)], delimiter=',', header='wavelength_nm,n,k,n_top,n_bottom,n_extraordinary,k_extraordinary', comments='', fmt='%.5f')


def main():
    from multiprocessing import Pool
    variants = ['iso', 'graded', 'uni', 'graded+uni', 'iso+nonuniform', 'graded+uni+nonuniform']
    with Pool(4) as pool:
        res = pool.map(fit, variants)
        best = min(res, key=lambda r: r['mse'] if 'nonuniform' not in r['variant'] else 1e9)
        ox = pool.starmap(fit, [(best['variant'], dox, np.array(best['x'])) for dox in (0.5, 1.0, 2.0, 3.0)])
    json.dump(dict(variants=res, oxide_scan=ox, best=best['variant']), open(W + 'grahil_fits.json', 'w'), indent=1)
    print('%-24s %7s | %7s %7s %8s | %8s %7s %7s %7s | n(550) mean / top / bottom / n_e' % ('variant', 'MSE', '200-300', '300-450', '450-1000', 'd (nm)', 'g', 's', 'sig'))
    for r in res:
        print('%-24s %7.2f | %7.1f %7.1f %8.1f | %6.2f+-%.2f %6.3f %6.3f %6.3f | %.4f / %.4f / %.4f / %.4f   k(300) %.4f k(633) %.4f k(800) %.4f'
              % (r['variant'], r['mse'], r['bands']['200-300'], r['bands']['300-450'], r['bands']['450-1000'], r['d'], r['d_ci'], r['g'], r['s'], r['sig'],
                 r['n_o']['550'][0], *r['n_top_bot_550'], r['n_e_550'], r['n_o']['300'][1], r['n_o']['633'][1], r['n_o']['800'][1]))
    print('oxide scan on %s: ' % best['variant'] + '  '.join('%.1f nm: MSE %.2f d %.2f' % (r['dox'], r['mse'], r['d']) for r in ox))
    export(best, W + 'GraHIL_nk.csv'); print('->', W + 'GraHIL_nk.csv')


if __name__ == '__main__':
    main()
