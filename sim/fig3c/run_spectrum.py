#!/usr/bin/env python3
"""Power-dissipation spectrum resolved in the in-plane wavevector
n_eff = k_x/k0, for the Fig 3(c) stack.  TM and TE are separated because the
surface plasmon lives entirely in TM."""
import numpy as np, csv, sys
import cps2, nspp, materials as M

LAM = 550.0
GRID = np.linspace(1.00, 2.60, 3201)
GRID = GRID[np.abs(GRID - 1.8) > 1.5e-3]        # skip the u = 1 branch point


def stack(no, ne, d):
    return cps2.Stack(LAM, (1.8, 1.8), 20.0, 10.0,
                      above=[(no, ne, d), (M.AG, M.AG, 100.0)],
                      below=[(1.8, 1.8, 50.0), (M.ITO, M.ITO, 50.0)], n_sub=1.8)


def run(d_list, ne_list):
    rows = []
    for d in d_list:
        for no, ne, lab in ne_list:
            with np.errstate(all='ignore'):
                tm, te = cps2.spectrum(stack(no, ne, d), GRID)
            for g, a, b in zip(GRID, tm, te):
                rows.append([lab, no, ne, d, g, a, b])
    return rows


FAM = [(1.8, x, 'n_e=%.2f' % x) for x in (1.80, 1.70, 1.60, 1.50)]
REAL = [(M.ETL[k][0], M.ETL[k][1], k) for k in ('B3PyMPM', 'B4PyMPM')]
rows = run([60.0, 100.0, 150.0], FAM + REAL)

with open('fig3c_spectrum.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['series', 'n_o', 'n_e', 'd_ETL_nm', 'n_eff', 'TM', 'TE'])
    for r in rows:
        w.writerow([r[0], '%.4f' % r[1], '%.4f' % r[2], '%.0f' % r[3],
                    '%.5f' % r[4], '%.6e' % r[5], '%.6e' % r[6]])

a = np.array([[r[3], r[4], r[5]] for r in rows], dtype=float)
lab = np.array([r[0] for r in rows])
print('TM peak position (k_x/k0)')
print('  series      n_SPP(inf)   d=60     d=100    d=150')
for no, ne, L in FAM + REAL:
    line = []
    for d in (60., 100., 150.):
        m = (lab == L) & (a[:, 0] == d) & (a[:, 1] > 1.55)
        line.append(a[m][np.argmax(a[m][:, 2]), 1])
    print('  %-10s   %.4f     ' % (L, nspp.n_spp(M.AG, no, ne))
          + '   '.join('%.4f' % x for x in line))
