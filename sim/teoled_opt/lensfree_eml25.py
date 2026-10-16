"""EML 25 nm instead of the 2 nm sheet (lensfree_refined.csv cases).
n_t = n_e cases: the EML thickness is optically invisible at a fixed dipole position, so the 2 nm optimum is reused with
d_E, d_H shortened by 11.5 nm (centre-dipole EQE identical, checked) and EQE_zone = uniform emission over the 25 nm EML is added.
n_t = 1.5 cases (EML n 1.8 now a 25 nm high-index slab): re-optimised (DE from the monochromatic optimum + polish from the
shifted 2 nm optimum), then EQE_zone.  -> lensfree_eml25.csv"""
import csv, numpy as np
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize
import lensfree_limit as L

D = 25.0; SH = (D - 2.0) / 2
FIX = lambda v: v if v in ('Ir', 'bulk', 'joint8', 'free', 'order1') else float(v)


def load():
    rows = []
    for r in csv.DictReader(open('lensfree_refined.csv')):
        r = {k: FIX(v) for k, v in r.items()}
        if r['fwhm'] != 'Ir': r['fwhm'] = int(r['fwhm'])
        r['lam0'] = int(r['lam0']); rows.append(r)
    return rows


def zone(x, c):
    l, w = L.spec(c)
    return float(np.mean([L.eqe(x, dict(c, z0=z), l, w) for z in np.linspace(D / 14, D - D / 14, 7)]))


def job(r):
    c = {k: r[k] for k in ('theta', 'fwhm', 'lam0', 'metal', 'n_t', 'n_e', 'n_c', 'order')}
    c['d_eml'] = D
    x2 = np.array([r['dE'] - SH, r['dH'] - SH, r['dAg'], r['dC']])
    l, w = L.spec(c)
    if c['n_t'] == c['n_e']:
        x, e = x2, L.eqe(x2, c, l, w)
    else:
        dmax = 150. if c['order'] == 'order1' else 400.
        bnds = [(15., dmax), (15., dmax), (7.87, 7.87) if c['metal'] == 'joint8' else (5., 60.), (0., 250.)]
        lam_d = c['lam0'] if c['fwhm'] != 'Ir' else 615
        rd = differential_evolution(lambda x: -L.eqe(x, dict(c, fwhm=0, lam0=lam_d), [lam_d], [1.0]), bnds, seed=1,
                                    popsize=10, maxiter=40, tol=1e-4, polish=True)
        x, e = None, -1
        for s0 in (rd.x, np.clip(x2, [b[0] for b in bnds], [b[1] for b in bnds])):
            rr = minimize(lambda v: -L.eqe(v, c, l, w), s0, method='Nelder-Mead', bounds=bnds,
                          options=dict(maxiter=300, xatol=0.2, fatol=1e-5))
            if -rr.fun > e: x, e = rr.x, -rr.fun
    out = dict(c, EQE_2nm=r['EQE'], EQE=e, EQE_zone=zone(x, c), dE=x[0], dH=x[1], dAg=x[2], dC=x[3])
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in out.items()}, flush=True)
    return out


if __name__ == '__main__':
    rows = load()
    rows.sort(key=lambda r: r['n_t'] == r['n_e'])        # slow re-optimisations first
    with Pool(4) as p:
        res = p.map(job, rows, chunksize=1)
    with open('lensfree_eml25.csv', 'w', newline='') as f:
        wr = csv.DictWriter(f, list(res[0].keys())); wr.writeheader(); wr.writerows(res)
