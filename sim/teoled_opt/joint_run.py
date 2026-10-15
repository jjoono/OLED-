"""Re-run the best candidates with the 2026-08 thin-Ag-on-HATCN optical constants.
TRinv  = T/R point inversion (dft_seedlayer_screen/data/nk, branch ellipsometry-tmm-fix-ag-seed-260820)
jointA = TRinv spectrum rescaled so that n, k at 550 nm equal the SE+T/R joint fit (README_260910_joint.md)
         -- a stand-in until the joint CSVs (jnk_nk/Ag_HATCN_*_nk.csv) are available."""
import csv, numpy as np, tmodel as T
from multiprocessing import Pool
from scipy.optimize import minimize
b = np.load('mla_hemi_1.65.npz'); BT, BR = b['BT'], b['BR']
JOINT550 = {7: (0.184, 3.360), 8: (0.153, 3.450)}
def load(d, kind):
    a = np.genfromtxt(f'nk_joint/Ag{d}nm_on_HATCN5_TRinv.csv', delimiter=',', names=True)
    n = np.interp(T.LAMS, a['wavelength_nm'], a['n']); k = np.interp(T.LAMS, a['wavelength_nm'], a['k'])
    if kind == 'jointA':
        n *= JOINT550[d][0] / np.interp(550, a['wavelength_nm'], a['n']); k *= JOINT550[d][1] / np.interp(550, a['wavelength_nm'], a['k'])
    return n + 1j * k
CASES = [(d, kind, bot, tfe) for d in (7, 8) for kind in ('TRinv', 'jointA') for bot in ('Al0.2', 'Liq') for tfe in (True, 'al2o3')]
def run(c):
    d, kind, bot, tfe = c; nkAg = load(d, kind); kw = dict(agtop=nkAg, tfe=tfe, al=0.2 if bot == 'Al0.2' else 0.0)
    f = lambda x: -T.eqe(x[0], x[1], 5, d, BT, BR, **kw) if (10 <= x[0] <= 170 and 10 <= x[1] <= 200) else 1.0
    r = min((minimize(f, s, method='Nelder-Mead', options=dict(xatol=1.5, fatol=1e-4, maxiter=60)) for s in ((100, 65), (110, 50))), key=lambda r: r.fun)
    return (d, kind, bot, 'full TFE' if tfe is True else 'Al2O3->lens', round(r.x[0], 1), round(r.x[1], 1), round(-r.fun, 4))
if __name__ == '__main__':
    with Pool(4) as p: res = p.map(run, CASES)
    with open('joint_run.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['Ag_nm', 'nk', 'bottom', 'TFE', 'd_ETL', 'd_HTL', 'EQE']] + res)
    for r in res: print('Ag %d %-6s Ag/%-5s %-12s dE %6.1f dH %6.1f EQE %.3f' % r)
