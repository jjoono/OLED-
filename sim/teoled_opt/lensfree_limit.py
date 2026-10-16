"""Upper limit of a lens-free (planar, emitting straight into air) top-emitting microcavity OLED.
  Ag 100 (McPeak) | ETL d_E (n_t) | EML 2 nm (n_e, dipole sheet at its centre) | HTL d_H (n_t) | top Ag d_Ag | cap d_C (n_c) | air
EML thickness d_eml is a case parameter (default 2 nm sheet; `python3 lensfree_limit.py 25` runs every case with a 25 nm EML,
dipole at its centre, and adds EQE_zone = uniform emission across the EML at the optimum) -> lensfree_limit_eml25.csv.
eta_out = angle-integrated power into air / total dissipated power (CPS), PLQY = 1 so EQE = eta_out (Purcell irrelevant).
Spectrum: delta or Gaussian (FWHM in nm) at lambda0, or the measured Ir(mphmq)2tmd spectrum; organics non-dispersive.
d_E, d_H, d_Ag, d_C optimised globally (differential evolution) per case; 'order1' bounds d_E, d_H <= 150 nm, 'free' <= 400 nm."""
import sys, csv, itertools, numpy as np
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize
sys.path[:0] = ['/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/teoled_opt']
import cps2, tmodel as T

TH = np.arange(0, 90, 0.1) + 0.05
a = np.genfromtxt('/home/user/OLED-/sim/teoled_opt/nk_joint/Ag_HATCN_2-6_8nm_nk.csv', delimiter=',', names=True)
AGJ = np.interp(T.LAMS, a['wavelength_nm'], a['n_joint']) + 1j * np.interp(T.LAMS, a['wavelength_nm'], a['k_joint'])


def eta_air(lam, x, c):
    dE, dH, dAg, dC = x
    i = T.idx(lam); ab = T.AG[i]
    at = AGJ[i] if c['metal'] == 'joint8' else ab
    nt, ne, nc = c['n_t'], c['n_e'], c['n_c']
    above = [(nt, nt, dE), (ab, ab, 100.)]
    below = [(nt, nt, dH), (at, at, dAg)] + ([(nc, nc, dC)] if dC > 0.5 else [])
    d = c.get('d_eml', 2.0)
    St = cps2.Stack(lam, (ne, ne), d, c.get('z0', d / 2), above=above, below=below, n_sub=1.0, n_top=1.0, h=c['theta'])
    P = np.real(cps2.solve(St, npts=3000)['P_tot'])
    return float(np.sum(np.real(cps2.sub_angular(St, TH))) * 0.1 / P)


def spec(c, fine=True):
    """coarse (9 points) for the global search, fine (25 / 41 points) for the polish and the reported EQE"""
    if c['fwhm'] == 'Ir':
        l = np.arange(560, 721, 4.0 if fine else 8.0); w = np.array([T.S[T.idx(x)] for x in l])
    elif c['fwhm'] == 0:
        l = np.array([c['lam0']]); w = np.ones(1)
    else:
        s = c['fwhm'] / 2.3548; h = (2.0 if fine else 1.5) * c['fwhm']
        n = 41 if fine else max(9, int(2 * h / 5.0) + 1)    # search: step <= 5 nm; report: FWHM/10
        l = np.linspace(c['lam0'] - h, c['lam0'] + h, n); w = np.exp(-0.5 * ((l - c['lam0']) / s) ** 2)
    return l, w / w.sum()


def eqe(x, c, l, w):
    return float(sum(wi * eta_air(li, x, c) for li, wi in zip(l, w)))


def run(c):
    l, w = spec(c, fine=False)
    dmax = 150. if c['order'] == 'order1' else 400.
    agb = (7.87, 7.87) if c['metal'] == 'joint8' else (5., 60.)
    bnds = [(15., dmax), (15., dmax), agb, (0., 250.)]
    lam_d = c['lam0'] if c['fwhm'] != 'Ir' else 615
    cd = dict(c, fwhm=0, lam0=lam_d)                     # monochromatic optimum first (cheap, no spectral aliasing)
    rd = differential_evolution(lambda x: -eqe(x, cd, [lam_d], [1.0]), bnds, seed=1, popsize=10, maxiter=40, tol=1e-4, polish=True)
    if c['fwhm'] == 0:
        x, best = rd.x, -rd.fun
    else:                                                # spectral search seeded with it, then multi-start polish on the fine spectrum
        r = differential_evolution(lambda x: -eqe(x, c, l, w), bnds, seed=1, popsize=8, maxiter=30, tol=1e-4, polish=False, x0=rd.x)
        lf, wf = spec(c); x, best = None, -1
        for s0 in (r.x, rd.x):
            rr = minimize(lambda x: -eqe(x, c, lf, wf), s0, method='Nelder-Mead', bounds=bnds, options=dict(maxiter=200, xatol=0.2, fatol=1e-5))
            if -rr.fun > best: x, best = rr.x, -rr.fun
    out = dict(c, EQE=best, dE=x[0], dH=x[1], dAg=x[2], dC=x[3])
    if c['fwhm'] not in (0, 'Ir'):                       # aliasing check on a 0.5 nm grid
        lc = np.arange(c['lam0'] - 2 * c['fwhm'], c['lam0'] + 2 * c['fwhm'] + 0.01, 0.5)
        wc = np.exp(-0.5 * ((lc - c['lam0']) / (c['fwhm'] / 2.3548)) ** 2); out['EQE_check'] = eqe(x, c, lc, wc / wc.sum())
    elif c['fwhm'] == 'Ir':
        lc = np.arange(560, 721, 1.0); wc = np.array([T.S[T.idx(v)] for v in lc]); out['EQE_check'] = eqe(x, c, lc, wc / wc.sum())
    else:
        out['EQE_check'] = best
    d = c.get('d_eml', 2.0)
    if d > 2.0:                                          # uniform emission zone across the EML (PLQY 1: plain average)
        lf, wf = spec(c)
        out['EQE_zone'] = float(np.mean([eqe(x, dict(c, z0=z), lf, wf) for z in np.linspace(d / 14, d - d / 14, 7)]))
    # power budget at lambda0 for the optimum
    lam = c['lam0'] if c['fwhm'] != 'Ir' else 615
    i = T.idx(lam); ab = T.AG[i]; at = AGJ[i] if c['metal'] == 'joint8' else ab
    St = cps2.Stack(lam, (c['n_e'],) * 2, d, d / 2, above=[(c['n_t'],) * 2 + (x[0],), (ab, ab, 100.)],
                    below=[(c['n_t'],) * 2 + (x[1],), (at, at, x[2])] + ([(c['n_c'],) * 2 + (x[3],)] if x[3] > 0.5 else []),
                    n_sub=1.0, n_top=1.0, h=c['theta'])
    r2 = cps2.solve(St, npts=3000)
    out.update(wg=float(np.real(r2['wg'])), spp_evan=float(np.real(r2['spp'])), abs_cone=float(np.real(r2['abs'])))
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in out.items()}, flush=True)
    return out


def cases():
    base = dict(theta=1.0, fwhm=0, lam0=620, metal='bulk', n_t=1.8, n_e=1.8, n_c=2.36, order='free')
    C = []
    for th, fw in itertools.product((0.67, 0.82, 0.9, 1.0), (0, 10, 20, 30, 50, 'Ir')):        # red, main map
        C.append(dict(base, theta=th, fwhm=fw))
    for lam0, fw in itertools.product((460, 530), (0, 20)):                                    # colours
        C.append(dict(base, lam0=lam0, fwhm=fw))
    for fw in (0, 20, 'Ir'):                                                                   # variants at Theta 1 and 0.82
        for th in (1.0, 0.82):
            C.append(dict(base, theta=th, fwhm=fw, metal='joint8'))
            C.append(dict(base, theta=th, fwhm=fw, n_t=1.5))
            C.append(dict(base, theta=th, fwhm=fw, n_c=1.8))
            C.append(dict(base, theta=th, fwhm=fw, order='order1'))
    return C


if __name__ == '__main__':
    D = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0
    C = [dict(c, d_eml=D) for c in cases()] if D != 2.0 else cases()
    with Pool(4) as p:
        res = p.map(run, C, chunksize=1)
    keys = list(res[0].keys())
    with open('lensfree_limit.csv' if D == 2.0 else 'lensfree_limit_eml%g.csv' % D, 'w', newline='') as f:
        wr = csv.DictWriter(f, keys); wr.writeheader(); wr.writerows(res)
