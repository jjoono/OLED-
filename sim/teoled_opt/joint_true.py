"""Candidates with the SE+T/R joint-fit n,k of Ag on HATCN (thin_Ag_nk_260910, 'joint' columns).
Ag thickness = SE-fitted thickness of each piece (README_260910_joint.md sec.5), nominal in the label."""
import csv, numpy as np, tmodel as T
from multiprocessing import Pool
from scipy.optimize import minimize
b = np.load('mla_hemi_1.65.npz'); BT, BR = b['BT'], b['BR']
PIECES = {6: ('1-8', 6.18), 7: ('2-5', 6.84), 8: ('2-6', 7.87), 10: ('2-7', 8.97), 12: ('2-8', 10.71)}
def load(nom):
    sh, _ = PIECES[nom]
    a = np.genfromtxt(f'nk_joint/Ag_HATCN_{sh}_{nom}nm_nk.csv', delimiter=',', names=True)
    return np.interp(T.LAMS, a['wavelength_nm'], a['n_joint']) + 1j * np.interp(T.LAMS, a['wavelength_nm'], a['k_joint'])
CASES = [(nom, bot, tfe) for nom in PIECES for bot in ('Al0.2', 'Liq') for tfe in (True, 'al2o3')]
def run(c):
    nom, bot, tfe = c; d = PIECES[nom][1]
    kw = dict(agtop=load(nom), tfe=tfe, al=0.2 if bot == 'Al0.2' else 0.0)
    f = lambda x: -T.eqe(x[0], x[1], 10, d, BT, BR, **kw) if (10 <= x[0] <= 170 and 10 <= x[1] <= 200) else 1.0
    r = min((minimize(f, s, method='Nelder-Mead', options=dict(xatol=1.5, fatol=1e-4, maxiter=60)) for s in ((100, 70), (110, 50))), key=lambda r: r.fun)
    return (nom, d, bot, 'full TFE' if tfe is True else 'Al2O3->lens', round(r.x[0], 1), round(r.x[1], 1), round(-r.fun, 4))
if __name__ == '__main__':
    with Pool(4) as p: res = p.map(run, CASES)
    with open('joint_true.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['Ag_nominal', 'Ag_fitted', 'bottom', 'TFE', 'd_ETL', 'd_HTL', 'EQE']] + res)
    for r in res: print('Ag %2d (%.2f) Ag/%-5s %-12s dE %6.1f dH %6.1f EQE %.3f' % r)
