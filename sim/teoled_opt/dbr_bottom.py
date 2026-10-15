"""Bottom mirror: Ag vs IZO-on-chirped-ZnS/LiF-DBR (glass behind, or Ag behind). Top: HATCN10 / Ag 8 nm nominal
(joint-fit nk, 7.87 nm), TFE variants, n1.65 MLA. d_ETL/d_HTL re-optimised (d_ETL up to 300 nm)."""
import csv, numpy as np, tmodel as T
from multiprocessing import Pool
from scipy.optimize import minimize
b = np.load('mla_hemi_1.65.npz'); BT, BR = b['BT'], b['BR']
a = np.genfromtxt('nk_joint/Ag_HATCN_2-6_8nm_nk.csv', delimiter=',', names=True)
AG8 = np.interp(T.LAMS, a['wavelength_nm'], a['n_joint']) + 1j * np.interp(T.LAMS, a['wavelength_nm'], a['k_joint'])
C = {
 'Ag / Al0.2 / Liq (ref)':                dict(al=0.2),
 'Ag / Liq (ref)':                        dict(),
 'DBR(ch) / IZO100 / Liq, glass behind':  dict(mirror='DBR', n_back=1.5, izo=100.),
 'DBR(ch) / IZO50 / Liq, glass behind':   dict(mirror='DBR', n_back=1.5, izo=50.),
 'DBR(ch rev) / IZO50 / Liq, glass':      dict(mirror='DBR', n_back=1.5, izo=50., dbr='chr'),
 'DBR(qw) / IZO50 / Liq, glass':          dict(mirror='DBR', n_back=1.5, izo=50., dbr='qw'),
 'DBR(ch) / IZO50 / Al1.5 / Liq, glass':  dict(mirror='DBR', n_back=1.5, izo=50., al=1.5),
 'Ag / DBR(ch) / IZO50 / Liq':            dict(mirror='DBR+Ag', izo=50.),
}
JOBS = [(n, tfe) for n in C for tfe in (True, 'al2o3')]
def run(j):
    name, tfe = j; kw = dict(C[name], agtop=AG8, tfe=tfe)
    f = lambda x: -T.eqe(x[0], x[1], 10, 7.87, BT, BR, **kw) if (10 <= x[0] <= 300 and 10 <= x[1] <= 200) else 1.0
    starts = ((40, 70), (100, 70), (160, 70), (230, 70))
    r = min((minimize(f, s, method='Nelder-Mead', options=dict(xatol=2, fatol=1e-4, maxiter=50)) for s in starts), key=lambda r: r.fun)
    return (name, 'full TFE' if tfe is True else 'Al2O3->lens', round(r.x[0], 1), round(r.x[1], 1), round(-r.fun, 4))
if __name__ == '__main__':
    with Pool(4) as p: res = p.map(run, JOBS)
    with open('dbr_bottom.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['bottom', 'TFE', 'd_ETL', 'd_HTL', 'EQE']] + res)
    for r in res: print('%-40s %-12s dE %6.1f dH %6.1f EQE %.3f' % r)
