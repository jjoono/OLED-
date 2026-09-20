#!/usr/bin/env python3
"""Fig 3(c): what a low out-of-plane index in the ETL buys, and when.

    air | Ag 100 nm | ETL (n_o, n_e, d_ETL) | EML 1.8, 20 nm, dipole at centre
        | HTL 1.8, 50 nm | ITO 50 nm (Koenig) | substrate 1.8

550 nm, PLQY = 1, isotropic dipole orientation.  n_sub = n_EML = 1.8, so the
substrate light line and the organic light line coincide at u = 1: everything
below it reaches the substrate, everything above it is SPP / near field.
"""
import numpy as np, csv, sys
import cps2, nspp, materials as M

LAM, NPTS = 550.0, 12000
D = np.concatenate([np.arange(10., 200., 2.), np.arange(200., 401., 5.)])
NE = np.round(np.arange(1.40, 1.8001, 0.02), 3)


def stack(no, ne, d):
    return cps2.Stack(LAM, (1.8, 1.8), 20.0, 10.0,
                      above=[(no, ne, d), (M.AG, M.AG, 100.0)],
                      below=[(1.8, 1.8, 50.0), (M.ITO, M.ITO, 50.0)], n_sub=1.8)


def sweep(no, ne, label):
    rows = []
    for d in D:
        r = cps2.solve(stack(no, ne, float(d)), npts=NPTS)
        rows.append([label, no, ne, d, r['air'], r['sub'], r['wg'], r['spp'],
                     r['abs'], r['air'] + r['sub']])
    return rows


def threshold(rows, col, target, above):
    """Smallest d at which the quantity first reaches the target and stays there."""
    a = np.array([[r[3], r[col]] for r in rows], dtype=float)
    ok = a[:, 1] >= target if above else a[:, 1] <= target
    if not ok.any():
        return float('nan')
    i = len(ok) - 1
    while i > 0 and ok[i - 1]:
        i -= 1
    return a[i, 0]


rows = []
for ne in NE:
    rows += sweep(1.8, float(ne), 'no1.80')
    print('  n_e = %.2f done' % ne, file=sys.stderr)
for name, (no, ne) in M.ETL.items():
    rows += sweep(no, ne, name)
    print('  %s done' % name, file=sys.stderr)

with open('fig3c_sweep.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['series', 'n_o', 'n_e', 'd_ETL_nm', 'air', 'sub_confined',
                'wg', 'spp', 'abs', 'eta_sub'])
    for r in rows:
        w.writerow([r[0], '%.4f' % r[1], '%.4f' % r[2], '%.1f' % r[3]]
                   + ['%.6f' % x for x in r[4:]])

out = []
for lab in ['no1.80'] + list(M.ETL):
    for ne in (NE if lab == 'no1.80' else [M.ETL[lab][1]]):
        sub = [r for r in rows if r[0] == lab and abs(r[2] - ne) < 1e-9]
        no = sub[0][1]
        out.append([lab, no, ne, nspp.n_spp(M.AG, no, ne),
                    threshold(sub, 7, 0.10, False),     # SPP < 10 %
                    threshold(sub, 9, 0.90, True),      # eta_sub > 90 %
                    max(r[9] for r in sub)])
with open('fig3c_threshold.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['series', 'n_o', 'n_e', 'n_SPP_Ag', 'd_SPP_below_10pct_nm',
                'd_eta_sub_above_90pct_nm', 'max_eta_sub'])
    for r in out:
        w.writerow([r[0], '%.4f' % r[1], '%.4f' % r[2], '%.4f' % r[3],
                    '%.1f' % r[4], '%.1f' % r[5], '%.4f' % r[6]])

print('\n  n_e   n_SPP    d(SPP<10%)  d(eta_sub>90%)  max eta_sub')
for r in out:
    tag = '' if r[0] == 'no1.80' else '   <- ' + r[0]
    print('  %.2f  %.4f   %6.0f nm    %6.0f nm       %.3f%s'
          % (r[2], r[3], r[4], r[5], r[6], tag))
