#!/usr/bin/env python3
"""EQE against the horizontal dipole ratio, on a generic stack.

    reflector 100 nm | ETL n=1.8 (d_ETL) | EML n=1.8, 20 nm, dipole at the centre
                     | HTL n=1.8 (d_HTL) | TCO 1.864 + 0.0032i, 50 nm | substrate
                     + hemispherical microlens array of the same index

550 nm, PLQY = 1.  No material is named: all three transport/emissive layers
share one organic index, and the TCO is the Koenig ITO.  The two devices differ
only in the reflector and the substrate/MLA index.

  reference  Al, n_sub = 1.50, thicknesses at their own cavity optimum
  proposed   Ag, n_sub = 1.80, deliberately un-tuned at 200/200 nm

EQE = eta_sub * eta_ext with eta_ext = p/[p + (1-p)A'] (Eq. 2), A' = 1 - <R_LED>
flux-averaged over the substrate angles, p the single-pass escape probability of
the MLA taken from the Fig. 2 ray trace.  Only eta_sub depends on the dipole
orientation; A' and p do not.
"""
import numpy as np, csv
import cps2, materials as M

AL = 0.958336812040867 + 6.68678039868499j
P_MLA = {1.50: 0.380, 1.80: 0.305}
DEV = {'reference': dict(refl=AL, n_sub=1.50, d_ETL=80.0, d_HTL=230.0),
       'proposed':  dict(refl=M.AG, n_sub=1.80, d_ETL=200.0, d_HTL=200.0)}
THETA = np.round(np.arange(0.50, 1.0001, 0.01), 3)
TH = np.radians(np.arange(0.25, 90.0, 0.5))
W = np.cos(TH) * np.sin(TH)


def stack(d, h):
    return cps2.Stack(550.0, (1.8, 1.8), 20.0, 10.0,
                      above=[(1.8, 1.8, d['d_ETL']), (d['refl'], d['refl'], 100.0)],
                      below=[(1.8, 1.8, d['d_HTL']), (M.ITO, M.ITO, 50.0)],
                      n_sub=d['n_sub'], h=h)


rows = []
summary = {}
for name, d in DEV.items():
    A = 1.0 - np.sum(cps2.stack_reflectance(stack(d, 2 / 3), np.degrees(TH)) * W) / np.sum(W)
    p = P_MLA[d['n_sub']]
    eta_ext = p / (p + (1 - p) * A)
    summary[name] = (A, p, eta_ext)
    for t in THETA:
        r = cps2.solve_pol(stack(d, float(t)), npts=16000)
        es = r['air'] + r['sub']
        rows.append([name, d['n_sub'], float(t), es, r['spp'], r['wg'], r['abs'],
                     A, p, eta_ext, es * eta_ext])

with open('fig3_theta.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['device', 'n_substrate', 'theta_horizontal', 'eta_sub', 'spp',
                'wg', 'abs', 'Aprime', 'p_escape', 'eta_ext', 'EQE'])
    for r in rows:
        w.writerow([r[0], '%.2f' % r[1], '%.2f' % r[2]] + ['%.6f' % x for x in r[3:]])

a = np.array([[r[2], r[10]] for r in rows])
lab = np.array([r[0] for r in rows])
print('device      A_prime   p      eta_ext | EQE at theta = 0.50 / 0.67 / 0.83 / 1.00')
for name in DEV:
    A, p, e = summary[name]
    m = a[lab == name]
    v = [np.interp(t, m[:, 0], m[:, 1]) for t in (0.50, 2 / 3, 0.83, 1.00)]
    print('%-10s  %.4f  %.3f  %.4f  | ' % (name, A, p, e)
          + '  '.join('%.4f' % x for x in v)
          + '   (swing %.1f %%p)' % (100 * (v[3] - v[0])))
