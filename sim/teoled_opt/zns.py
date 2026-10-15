import csv, numpy as np, tmodel as T
from multiprocessing import Pool
from scipy.optimize import minimize
b = np.load('mla_hemi_1.65.npz'); BT, BR = b['BT'], b['BR']
def opt(a):
    agt, dAg, zns, tfe = a
    f = lambda x: -T.eqe(max(x[0], 5), max(x[1], 10), 10, dAg, BT, BR, agtop=agt, zns=zns, tfe=tfe)
    best = None
    for x0 in ((100, 60), (110, 80)):
        r = minimize(f, x0, method='Nelder-Mead', options=dict(xatol=2, fatol=1e-4, maxiter=50))
        if best is None or r.fun < best.fun: best = r
    return (agt, dAg, zns, tfe, round(best.x[0], 1), round(best.x[1], 1), -best.fun)
jobs = [(agt, dAg, z, tfe) for agt, dAg in (('meas', 6), ('meas', 8), ('bulk', 8)) for z in (0, 20, 40, 60, 80) for tfe in (True,)]
jobs += [('meas', 8, z, False) for z in (0, 40, 60)]
if __name__ == '__main__':
    with Pool(4) as p: res = p.map(opt, jobs)
    with open('zns.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['agtop', 'dAg', 'ZnS', 'TFE', 'dE', 'dH', 'EQE']] + res)
    for r in res: print('%-5s Ag %d ZnS %3d TFE %-5s dE %6.1f dH %6.1f EQE %.3f' % r)
