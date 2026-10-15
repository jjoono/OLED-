"""ZnS-capped high-index TFE, re-optimised with dense spectral sampling (27 wavelengths, 560-716 nm; the 6-point
weighting overrates narrow cavity resonances by up to 8 %p). Free: d_ETL (upper bound varied), d_HTL, d_ZnS (>=100).
Ag 8 nm nominal (joint nk), Ag/Al0.2/Liq, HATCN 10."""
import csv, numpy as np, tmodel as T, hi175 as H
from multiprocessing import Pool
from scipy.optimize import minimize
LF = tuple(np.arange(560, 717, 6))
EFF = H.eff(LF)
W = np.array([T.S[T.idx(l)] for l in LF]); W /= W.sum()
NK, DAG = H.load(8), 7.87
def eqe(dE, dH, z):
    eo = [T.eta_out(l, dE, dH, 10, DAG, *EFF[l], agtop=NK, tfe='al2o3', al=0.2, zns=z)[0] for l in LF]
    return T.ETA_RAD * float(W @ np.array(eo))
def run(a):
    emax, zmin = a
    f = lambda x: -eqe(*x) if (20 <= x[0] <= emax and 20 <= x[1] <= 320 and zmin <= x[2] <= 220) else 1.0
    starts = [(min(emax, 90), 80, max(zmin, 110)), (min(emax, 130), 200, max(zmin, 120)), (min(emax, 100), 150, max(zmin, 150))]
    r = min((minimize(f, s, method='Nelder-Mead', options=dict(xatol=2, fatol=2e-4, maxiter=90)) for s in starts), key=lambda r: r.fun)
    return (emax, zmin, *[round(v, 1) for v in r.x], round(-r.fun, 4))
if __name__ == '__main__':
    J = [(e, 100) for e in (80, 100, 130, 160, 200)] + [(130, 0)]
    with Pool(4) as p: res = p.map(run, J)
    with open('zns_hi2.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['dETL_max', 'ZnS_min', 'd_ETL', 'd_HTL', 'd_ZnS', 'EQE']] + res)
    for r in res: print('dETL<=%3d ZnS>=%3d : dE %6.1f dH %6.1f ZnS %6.1f  EQE %.3f' % r)
