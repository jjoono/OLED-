"""Ag top | ZnS (>=100 nm, ALD H2O barrier) | Al2O3 30 | n1.75 3 um | Al2O3 30 | n1.75 MLA.
(1) EQE map over d_ETL x d_HTL (Ag 8 joint nk, Ag/Al0.2/Liq) for ZnS 0/100/120/150;
(2) Nelder-Mead optimum vs ZnS thickness for Ag 8 and Ag 12."""
import csv, itertools, numpy as np, tmodel as T
import hi175 as H                                 # sets T.NL = 1.75 and builds the effective MLA response
from multiprocessing import Pool
from scipy.optimize import minimize
NK = {8: (H.load(8), 7.87), 12: (H.load(12), 10.71)}
def eqe(dE, dH, ag, zns):
    nk, d = NK[ag]; return H.eqe(dE, dH, d, nk, dict(al=0.2, zns=zns))
DE = range(40, 161, 10); DH = range(40, 241, 20); ZM = (0, 100, 120, 150)
def mapjob(a): z, dE, dH = a; return (z, dE, dH, eqe(dE, dH, 8, z))
def optjob(a):
    ag, z = a
    f = lambda x: -eqe(x[0], x[1], ag, z) if (20 <= x[0] <= 200 and 20 <= x[1] <= 300) else 1.0
    r = min((minimize(f, s, method='Nelder-Mead', options=dict(xatol=1.5, fatol=1e-4, maxiter=60)) for s in ((80, 60), (130, 80), (90, 180))), key=lambda r: r.fun)
    return (ag, z, round(r.x[0], 1), round(r.x[1], 1), round(-r.fun, 4))
if __name__ == '__main__':
    with Pool(4) as p:
        M = p.map(mapjob, list(itertools.product(ZM, DE, DH)), chunksize=8)
        O = p.map(optjob, [(ag, z) for ag in (8, 12) for z in (0, 100, 110, 120, 130, 140, 150, 160, 180)])
    with open('zns_hi_map.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['ZnS', 'd_ETL', 'd_HTL', 'EQE']] + M)
    with open('zns_hi_opt.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['Ag', 'ZnS', 'd_ETL', 'd_HTL', 'EQE']] + O)
    for r in O: print('Ag %2d ZnS %3d  dE %6.1f dH %6.1f EQE %.3f' % r)
