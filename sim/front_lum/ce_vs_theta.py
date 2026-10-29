"""On-axis CE and EQE versus the horizontal dipole ratio Theta, structures fixed at their Theta = 0.76 optimum
(results_ce/*.json, back diffuser mu_a d = 0.01, best S/g); the glass reference is also re-optimised at every Theta.
-> ce_vs_theta.csv"""
import json, csv, numpy as np
from scipy.optimize import minimize
import ce_opt as C
def struct(f):
    s = json.load(open(f))['structure']; return [s[k] for k in ('d_ETL', 'd_HTL', 'lam_c', 'd_IZO_front', 'd_IZO_back', 'n_s')]
CASES = [('B N0', 0, struct('results_ce/B_N0_mua0.01.json')), ('B N4', 4, struct('results_ce/B_N4_mua0.01.json')),
         ('B N6', 6, struct('results_ce/B_N6_mua0.01.json'))]
x_ref = [json.load(open('results_ce/REF.json'))['structure'][k] for k in ('d_ETL', 'd_HTL')]
rows = []
for th in (2 / 3, 0.76, 0.85, 0.95, 1.0):
    C.THETA = th
    r_fix = C.evaluate('REF', x_ref, step=10)
    m = minimize(lambda x: -C.evaluate('REF', list(x), step=10)['CE'], x_ref, method='Nelder-Mead', options=dict(maxiter=80, xatol=1, fatol=1e-3))
    r_opt = C.evaluate('REF', list(m.x), step=10)
    row = dict(Theta=round(th, 3), REF_fixed_CE=round(r_fix['CE'], 1), REF_fixed_EQE=round(100 * r_fix['EQE'], 1),
               REF_reopt_CE=round(r_opt['CE'], 1), REF_reopt_EQE=round(100 * r_opt['EQE'], 1))
    for name, N, x in CASES:
        r = C.evaluate('B', x, N, 0.01, step=10); row[name + '_CE'] = round(r['CE'], 1); row[name + '_EQE'] = round(100 * r['EQE'], 1)
    rows.append(row); print(row, flush=True)
with open('ce_vs_theta.csv', 'w', newline='') as fo:
    w = csv.DictWriter(fo, list(rows[0])); w.writeheader(); w.writerows(rows)
