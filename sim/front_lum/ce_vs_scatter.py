"""On-axis CE of the optimised DBR-4 structure (results_ce/B_N4_mua0.01.json) for every back-diffuser setting
(S, g, mu_a d) -> how weak (thin) can the scattering layer be?   -> ce_vs_scatter.csv"""
import json, csv, numpy as np
import ce_opt as C
d = json.load(open('results_ce/B_N4_mua0.01.json')); s = d['structure']
x = [s['d_ETL'], s['d_HTL'], s['lam_c'], s['d_IZO_front'], s['d_IZO_back'], s['n_s']]
js = list(range(0, 301, 10)); res = [C.point(('B', j, tuple(x), 4, 1500)) for j in js]
ref = json.load(open('results_ce/REF.json'))['results']['CE']
rows = []
for mua in (0.0, 0.01, 0.03):
    for S, g in C.SG:
        B = C.BZ['S%g_g%g_a%g' % (S, g, mua)]
        lam = np.array(js) + 400.0; p = C.I[js] * lam; p /= p.sum(); F = np.array([r[0] for r in res]); k = C.Q * p * F / (1 - C.Q + C.Q * (p @ F))
        conv = 683 * C.V(lam) * 1239.84 / lam; Ld = np.array([r[1] for r in res]); Lr = []
        for r in res:
            Ps, Tf, Rf, ns = r[3]; U = np.linalg.solve(np.eye(C.NB) - B * Rf[None, :], B @ Ps); Lr.append((U[:3] * Tf[:3]).sum() / C.OM / ns ** 2)
        ce = float(conv * k @ (Ld + np.array(Lr))); rows.append(dict(mua=mua, S=S, g=g, CE=round(ce, 1), ratio=round(ce / ref, 2)))
        print(rows[-1], flush=True)
with open('ce_vs_scatter.csv', 'w', newline='') as fo:
    w = csv.DictWriter(fo, list(rows[0])); w.writeheader(); w.writerows(rows)
