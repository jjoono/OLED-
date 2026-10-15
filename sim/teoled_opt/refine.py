import itertools, csv, numpy as np, tmodel as T
from multiprocessing import Pool
from scipy.optimize import minimize
b = np.load('mla_hemi_1.65.npz'); BT, BR = b['BT'], b['BR']
def opt(a):
    tag, dAg, x0, kw = a
    f = lambda x: -T.eqe(max(x[0], 5), max(x[1], 10), kw.pop('dHAT', 10) if False else kw.get('dHAT', 10), dAg, BT, BR,
                          **{k: v for k, v in kw.items() if k != 'dHAT'})
    r = minimize(f, x0, method='Nelder-Mead', options=dict(xatol=2, fatol=1e-4, maxiter=60))
    return (tag, dAg, round(r.x[0], 1), round(r.x[1], 1), -r.fun)
jobs = []
for agt in ('meas', 'bulk'):
    for dAg in (6, 7, 8, 10):
        for x0 in ((100, 70), (220, 80)):
            jobs.append((f'Ag top {agt}', dAg, x0, dict(agtop=agt)))
# variants at Ag 8 meas
V = {'Ag/Al0.2/Liq': dict(al=0.2), 'Al/Liq (Batch)': dict(bottom='Al'), 'no TFE (lens on Ag)': dict(tfe=False),
     'z0 = 2.5 (ETL side)': dict(z0=2.5), 'z0 = 22.5 (HTL side)': dict(z0=22.5)}
for name, kw in V.items():
    jobs.append((name, 8, (100, 70), dict(agtop='meas', **kw)))
# index study (Ag 8 meas): (n_EML, n_ETL)
for ne, nt in ((1.81, 1.81), (1.65, 1.81), (1.81, 1.65), (1.65, 1.65), (1.55, 1.55)):
    for x0 in ((100, 70), (220, 80)):
        jobs.append((f'n_EML {ne} / n_ETL {nt}', 8, x0, dict(agtop='meas', n_eml=ne, n_etl=nt)))
if __name__ == '__main__':
    with Pool(4) as p: res = p.map(opt, jobs)
    with open('refine.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['case', 'dAg', 'dE', 'dH', 'EQE']] + res)
    for r in res: print('%-28s Ag %2d  dE %6.1f  dH %6.1f  EQE %.3f' % r)
