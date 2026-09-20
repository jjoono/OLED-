#!/usr/bin/env python3
"""Is n_SPP the governing variable, or is it n_e?

Same stack as run_fig3c.py, at two fixed ETL thicknesses, for three families of
ETL ordinary index.  If n_SPP is what matters, the three families must collapse
onto one curve when plotted against n_SPP and must not when plotted against n_e.
"""
import numpy as np, csv
import cps2, nspp, materials as M

LAM, NPTS = 550.0, 12000
D_LIST = [60.0, 100.0, 150.0]
NO_LIST = [1.70, 1.80, 1.90]


def stack(no, ne, d):
    return cps2.Stack(LAM, (1.8, 1.8), 20.0, 10.0,
                      above=[(no, ne, d), (M.AG, M.AG, 100.0)],
                      below=[(1.8, 1.8, 50.0), (M.ITO, M.ITO, 50.0)], n_sub=1.8)


rows = []
for no in NO_LIST:
    for ne in np.round(np.arange(1.36, no + 1e-9, 0.02), 3):
        for d in D_LIST:
            r = cps2.solve(stack(no, float(ne), d), npts=NPTS)
            rows.append(['n_o=%.2f' % no, no, float(ne),
                         nspp.n_spp(M.AG, no, float(ne)), d,
                         r['air'] + r['sub'], r['spp'], r['abs'], r['wg']])
for name, (no, ne) in M.ETL.items():
    for d in D_LIST:
        r = cps2.solve(stack(no, ne, d), npts=NPTS)
        rows.append([name, no, ne, nspp.n_spp(M.AG, no, ne), d,
                     r['air'] + r['sub'], r['spp'], r['abs'], r['wg']])

with open('fig3c_collapse.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['series', 'n_o', 'n_e', 'n_SPP_Ag', 'd_ETL_nm', 'eta_sub',
                'spp', 'abs', 'wg'])
    for r in rows:
        w.writerow([r[0]] + ['%.4f' % x for x in r[1:4]] + ['%.1f' % r[4]]
                   + ['%.6f' % x for x in r[5:]])
print('fig3c_collapse.csv : %d rows' % len(rows))

a = np.array([[r[2], r[3], r[4], r[5]] for r in rows])
lab = np.array([r[0] for r in rows])
print('\nspread of eta_sub between the three n_o families at d = 150 nm')
for x, name in ((0, 'same n_e  '), (1, 'same n_SPP')):
    g = np.linspace(1.45, 1.68, 12) if x == 0 else np.linspace(1.62, 1.95, 12)
    dev = []
    for v in g:
        vals = []
        for no in NO_LIST:
            m = (lab == 'n_o=%.2f' % no) & (a[:, 2] == 150.0)
            if m.sum() < 2:
                continue
            xs, ys = a[m][:, x], a[m][:, 3]
            o = np.argsort(xs)
            if xs[o][0] <= v <= xs[o][-1]:
                vals.append(np.interp(v, xs[o], ys[o]))
        if len(vals) == 3:
            dev.append(max(vals) - min(vals))
    print('   binned by %s : mean spread %.4f, max spread %.4f  (n=%d bins)'
          % (name, np.mean(dev), np.max(dev), len(dev)))
