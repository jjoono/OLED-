"""Ag top | Al2O3 30 | n=1.75 layer 3 um | Al2O3 30 | n=1.75 hemisphere MLA (1.5-dyad TFE with a high-index organic).
The 3 um layer is treated as incoherent (round-trip path >> coherence length) and index-matched to the lens body;
the upper Al2O3 30 nm is an angle-dependent lossless interface between them, folded into the MLA response:
  BT_eff = BT (I - R_I BR)^-1 T_I,   BR_eff = R_I + T_I BR (I - R_I BR)^-1 T_I."""
import csv, numpy as np, tmodel as T
from multiprocessing import Pool
from scipy.optimize import minimize
T.NL = 1.75
b = np.load('mla_hemi_1.75.npz'); BT0, BR0 = b['BT'], b['BR']
I90 = np.eye(90)
def eff(lam_list=(580, 600, 615, 630, 650, 680)):
    out = {}
    for lam in lam_list:
        R = T.R_from_lens([(T.AL2O3[T.idx(lam)], 30.0)], lam, 1.75)    # Al2O3 30 nm between two n=1.75 media
        Tt = 1 - R
        X = np.linalg.solve(I90 - R[:, None] * BR0, np.diag(Tt))
        out[lam] = (BT0 @ X, np.diag(R) + Tt[:, None] * (BR0 @ X))
    return out
EFF = eff()
def eqe(dE, dH, dAg, nkAg, kw, eff_on=True, lams=(580, 600, 615, 630, 650, 680)):
    w = np.array([T.S[T.idx(l)] for l in lams]); w = w / w.sum(); eo = []
    for l in lams:
        BT, BR = EFF[l] if eff_on else (BT0, BR0)
        eo.append(T.eta_out(l, dE, dH, 10, dAg, BT, BR, agtop=nkAg, tfe='al2o3', **kw)[0])
    return T.ETA_RAD * float(w @ np.array(eo))
PIECES = {7: ('2-5', 6.84), 8: ('2-6', 7.87), 10: ('2-7', 8.97), 12: ('2-8', 10.71)}
def load(nom):
    sh, _ = PIECES[nom]; a = np.genfromtxt(f'nk_joint/Ag_HATCN_{sh}_{nom}nm_nk.csv', delimiter=',', names=True)
    return np.interp(T.LAMS, a['wavelength_nm'], a['n_joint']) + 1j * np.interp(T.LAMS, a['wavelength_nm'], a['k_joint'])
JOBS = [(nom, bot, e) for nom in PIECES for bot in ('Al0.2', 'Liq') for e in (True,)] + [(8, 'Al0.2', False), (8, 'Liq', False)]
def run(j):
    nom, bot, e = j; kw = dict(al=0.2 if bot == 'Al0.2' else 0.0); nk = load(nom); d = PIECES[nom][1]
    f = lambda x: -eqe(x[0], x[1], d, nk, kw, e) if (10 <= x[0] <= 170 and 10 <= x[1] <= 200) else 1.0
    r = min((minimize(f, s, method='Nelder-Mead', options=dict(xatol=1.5, fatol=1e-4, maxiter=60)) for s in ((100, 70), (110, 50))), key=lambda r: r.fun)
    return (nom, d, bot, 'with upper Al2O3' if e else 'no upper Al2O3', round(r.x[0], 1), round(r.x[1], 1), round(-r.fun, 4))
if __name__ == '__main__':
    with Pool(4) as p: res = p.map(run, JOBS)
    with open('hi175.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['Ag_nominal', 'Ag_fit', 'bottom', 'interface', 'd_ETL', 'd_HTL', 'EQE']] + res)
    for r in res: print('Ag %2d (%.2f) Ag/%-5s %-17s dE %6.1f dH %6.1f EQE %.3f' % r)
