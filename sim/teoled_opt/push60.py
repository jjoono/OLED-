import csv, numpy as np, tmodel as T
from multiprocessing import Pool
from scipy.optimize import minimize
b = np.load('mla_hemi_1.65.npz'); BT, BR = b['BT'], b['BR']
C = {
 'C  Ag/Al0.2 | HATCN10/Ag12 | full TFE (baseline)':        dict(al=0.2, dHAT=10, dAg=12),
 'G  Ag/Al0.2 | HATCN10/Ag8  | full TFE':                   dict(al=0.2, dHAT=10, dAg=8),
 'G1 Ag/Al0.2 | HATCN10/Ag8  | Al2O3 30 -> lens -> ParyC':  dict(al=0.2, dHAT=10, dAg=8, tfe='al2o3'),
 'G2 Ag/Al0.2 | HATCN10/Ag8 film 50% to bulk | full TFE':   dict(al=0.2, dHAT=10, dAg=8, agtop=0.5),
 'G3 Ag/Al0.2 | HATCN10/Ag8 film 50% | Al2O3->lens':        dict(al=0.2, dHAT=10, dAg=8, agtop=0.5, tfe='al2o3'),
 'C1 Ag/Al0.2 | HATCN10/Ag12 | Al2O3->lens':                dict(al=0.2, dHAT=10, dAg=12, tfe='al2o3'),
 'C2 Ag/Al0.2 | HATCN10/Ag12 film 50% | Al2O3->lens':       dict(al=0.2, dHAT=10, dAg=12, agtop=0.5, tfe='al2o3'),
 'H  Ag/Liq   | HATCN10/Ag8  | Al2O3->lens':                dict(dHAT=10, dAg=8, tfe='al2o3'),
 'I  Ag/Al0.2 | HATCN10/Ag10 | Al2O3->lens':                dict(al=0.2, dHAT=10, dAg=10, tfe='al2o3'),
 'J  Ag/Al0.2 | HATCN10/Ag10 film 50% | Al2O3->lens':       dict(al=0.2, dHAT=10, dAg=10, agtop=0.5, tfe='al2o3'),
}
def run(name):
    kw = dict(C[name]); dHAT = kw.pop('dHAT'); dAg = kw.pop('dAg'); kw.setdefault('agtop', 'meas')
    def f(x):
        if not (10 <= x[0] <= 170 and 10 <= x[1] <= 200): return 1.0
        return -T.eqe(x[0], x[1], dHAT, dAg, BT, BR, **kw)
    r = min((minimize(f, s, method='Nelder-Mead', options=dict(xatol=1.5, fatol=1e-4, maxiter=60)) for s in ((95, 70), (110, 50))), key=lambda r: r.fun)
    return (name, round(r.x[0], 1), round(r.x[1], 1), round(-r.fun, 4))
if __name__ == '__main__':
    with Pool(4) as p: res = p.map(run, list(C))
    with open('push60.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['structure', 'd_ETL', 'd_HTL', 'EQE']] + res)
    for r in res: print('%-58s dE %6.1f dH %6.1f EQE %.3f' % r)
