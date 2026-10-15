"""Electrically plausible inverted TEOLED candidates: optimal d_ETL / d_HTL (1st and 2nd cavity order)."""
import csv, numpy as np, tmodel as T
from multiprocessing import Pool
from scipy.optimize import minimize
b = np.load('mla_hemi_1.65.npz'); BT, BR = b['BT'], b['BR']
C = {
 'A Al/Liq | TAPC | HATCN50/Ag15 (Batch10 replica)': dict(bottom='Al', dHAT=50, dAg=15),
 'B Ag/Al0.2/Liq | TAPC | HATCN50/Ag15':             dict(al=0.2, dHAT=50, dAg=15),
 'C Ag/Al0.2/Liq | TAPC | HATCN10/Ag12':             dict(al=0.2, dHAT=10, dAg=12),
 'D Ag/Liq | TAPC | HATCN10/Ag12':                   dict(dHAT=10, dAg=12),
 'E Ag/Al0.2/Liq | EML only (RZ at NPB) | NPB | HATCN10/Ag12': dict(al=0.2, dHAT=10, dAg=12, htl='NPB', z0=20.0),
 'F Ag/Liq | TAPC | HATCN5/Ag8 (stretch)':           dict(dHAT=5, dAg=8),
}
def run(a):
    name, order = a; kw = dict(C[name]); dHAT = kw.pop('dHAT'); dAg = kw.pop('dAg')
    lo, hi = (10, 170) if order == 1 else (170, 340)
    def f(x):
        dE, dH = x
        if not (lo <= dE <= hi and 10 <= dH <= 250): return 1.0
        return -T.eqe(dE, dH, dHAT, dAg, BT, BR, **kw)
    starts = ((50, 40), (100, 70), (140, 100)) if order == 1 else ((220, 60), (260, 90))
    best = min((minimize(f, s, method='Nelder-Mead', options=dict(xatol=1.5, fatol=1e-4, maxiter=60)) for s in starts), key=lambda r: r.fun)
    dE, dH = best.x
    # tolerance: EQE drop for +-10 nm in each
    e0 = -best.fun
    tol = [T.eqe(dE + s * 10, dH, dHAT, dAg, BT, BR, **kw) for s in (-1, 1)] + [T.eqe(dE, dH + s * 10, dHAT, dAg, BT, BR, **kw) for s in (-1, 1)]
    return (name, order, round(dE, 1), round(dH, 1), dHAT, dAg, round(e0, 4), round(e0 - min(tol[:2]), 4), round(e0 - min(tol[2:]), 4))
if __name__ == '__main__':
    jobs = [(n, o) for n in C for o in (1, 2)]
    with Pool(4) as p: res = p.map(run, jobs)
    with open('cands.csv', 'w', newline='') as fh:
        csv.writer(fh).writerows([['structure', 'order', 'd_ETL(Bebq2)', 'd_HTL', 'HATCN', 'Ag', 'EQE', 'drop dETL+-10', 'drop dHTL+-10']] + res)
    for r in res: print('%-60s ord%d dE %6.1f dH %6.1f HAT %2d Ag %2d EQE %.3f  (-%.3f / -%.3f)' % r)
