"""Every Si piece against every material's dispersion, thickness the only free
parameter.  The thin pieces cannot carry their own dispersion (n and d trade at
~10 nm), so this is how they are identified; the thick ones are re-scored the
same way so all eight thicknesses come out of one consistent procedure.
The native oxide is the one assumption: thickness is also reported at
d_ox = 1 and 3 nm."""
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
MATS = ('NPB', 'TAPC', 'Bphen')

def fit(sh, mat, d_ox=2.0):
    wl, P, D, ang = L.se(sh, 245, 1000, 2); meas = L.ncs(P, D)
    Nf = O.N_of(np.array(OSC[mat]['q']), wl)
    ox, si, air = R.n_sio2(wl), R.n_si(wl), R.n_air(wl)
    def res(x):
        n_, c_, s_ = M.se_ncs(wl, [air, Nf, ox, si], [x[0], d_ox], ang)
        return np.concatenate([(n_ - meas[0]).ravel(), (c_ - meas[1]).ravel(), (s_ - meas[2]).ravel()])
    best = None
    for d0 in np.arange(1.0, 120.0, 3.0):
        r = least_squares(res, [d0], bounds=([0.2], [300.0]))
        if best is None or r.cost < best.cost: best = r
    return float(best.x[0]), float(M.mse_ncs(best.fun, 1))

if __name__ == '__main__':
    from multiprocessing import Pool
    shs = [s for s in L.sheets() if s.split('_')[1].startswith('2')]
    jobs = [(sh, m) for sh in shs for m in MATS]
    with Pool(4) as p:
        r = p.starmap(fit, jobs)
    tab = {}
    for (sh, m), v in zip(jobs, r):
        tab.setdefault(sh, {})[m] = v
    with Pool(4) as p:
        ox = p.starmap(fit, [(sh, min(tab[sh], key=lambda m: tab[sh][m][1]), dox) for sh in shs for dox in (1.0, 3.0)])
    print('Si pieces, dispersion fixed, thickness only (native oxide 2 nm)')
    print('%-5s | %-17s | %-17s | %-17s | %-6s %s' % ('sheet', 'NPB  d / MSE', 'TAPC  d / MSE', 'Bphen  d / MSE', 'best', 'd at oxide 1 / 3 nm'))
    for i, sh in enumerate(shs):
        t = tab[sh]; b = min(t, key=lambda m: t[m][1])
        tab[sh]['best'] = b; tab[sh]['d_ox1'] = ox[2 * i][0]; tab[sh]['d_ox3'] = ox[2 * i + 1][0]
        print('%-5s | %s | %-6s %.2f / %.2f' % (sh, ' | '.join('%6.2f / %7.2f%s' % (t[m][0], t[m][1], '*' if m == b else ' ') for m in MATS),
              b, ox[2 * i][0], ox[2 * i + 1][0]))
    json.dump(tab, open(W + 'si_id.json', 'w'), indent=1)
