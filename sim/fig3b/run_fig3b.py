#!/usr/bin/env python3
"""Fig 3(b): mode-resolved power budget vs organic thickness and substrate index.

        air | Ag 100 nm | organic d_org (isotropic n=1.8, emitter at centre)
            | ITO 50 nm (1.86 + 0.003i) | glass n_sub

lambda = 550 nm, PLQY = 1, isotropic dipole orientation (horizontal ratio 2/3).
Same stack and same five channels as Planar_sweep22_preprint.m; only the
u-quadrature is different (see cps.py).
"""
import numpy as np, csv, sys
import cps

LAM = 550.0
N_SUBS = [1.50, 1.65, 1.80]
D = np.arange(10.0, 500.0 + 1e-9, 2.0)        # 2 nm steps -> smooth curves
NPTS = 12000
P_FREE = 0.5                                   # analytic, in this normalisation

rows = []
for ns in N_SUBS:
    for d in D:
        r = cps.solve(float(d), ns, lam=LAM, npts=NPTS)
        rows.append([d, ns, r['air'], r['sub'], r['wg'], r['spp'], r['abs'],
                     r['air'] + r['sub'], r['P_tot'] / P_FREE])
    print('  n_sub = %.2f done' % ns, file=sys.stderr)

hdr = ['d_organic_nm', 'n_substrate', 'EQE_air', 'EQE_sub_confined',
       'EQE_wg', 'EQE_spp', 'EQE_abs', 'eta_sub_total', 'Purcell_factor']
with open('fig3b_modes.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(hdr)
    for r in rows:
        w.writerow(['%.1f' % r[0], '%.2f' % r[1]] + ['%.6f' % x for x in r[2:]])
print('fig3b_modes.csv : %d rows' % len(rows))

a = np.array(rows)
print('\n closure  max|sum-1| = %.2e' % np.max(np.abs(a[:, 2:7].sum(1) - 1)))
COL = dict(d=0, ns=1, air=2, sub=3, wg=4, spp=5, ab=6, eta=7, F=8)
print('\n   n_sub | d where SPP<10%  SPP@500nm  wg@500nm  eta_sub@500nm   max eta_sub')
for ns in N_SUBS:
    m = a[a[:, COL['ns']] == ns]
    below = m[m[:, COL['spp']] < 0.10]
    dthr = below[0, 0] if len(below) else float('nan')
    i = int(np.argmax(m[:, COL['eta']]))
    print('   %.2f  |     %6.0f nm       %6.3f    %6.3f       %6.3f      %6.3f (%3.0f nm)'
          % (ns, dthr, m[-1, COL['spp']], m[-1, COL['wg']], m[-1, COL['eta']],
             m[i, COL['eta']], m[i, COL['d']]))
