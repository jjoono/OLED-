#!/usr/bin/env python3
"""EQE against microlens aspect ratio, generic stack, full recycling series.

    reflector 100 nm | ETL n=1.8 200 nm | EML n=1.8 20 nm, dipole at the centre
                     | HTL n=1.8 200 nm | TCO 1.8636+0.0032i 50 nm | substrate 1.80
                     + hexagonally close-packed spherical caps, index-matched,
                       infinite planar substrate and infinite array

550 nm, PLQY = 1, isotropic emitter.

The closed form eta_ext = p/[p+(1-p)A'] is NOT used here: it assumes the surface
randomises the angles at every bounce, and a shallow lens does not, which is
exactly what the low-aspect-ratio end of this sweep is about.  The angle-resolved
B_T and B_R come from the ray trace and the series is summed term by term.
"""
import numpy as np, csv, sys
sys.path.insert(0, '/home/user/OLED-/sim/fig3c')
import cps2, materials as M
import raytrace as rt, series

AL = 0.958336812040867 + 6.68678039868499j
N_SUB = 1.80
AR = np.round(np.concatenate([np.arange(0.10, 0.51, 0.05),
                              np.arange(0.60, 1.51, 0.10)]), 2)
TH = np.arange(90) + 0.5
W = np.cos(np.radians(TH)) * np.sin(np.radians(TH))
NRAY = 30000


def dev(refl):
    return cps2.Stack(550.0, (1.8, 1.8), 20.0, 10.0,
                      above=[(1.8, 1.8, 200.0), (refl, refl, 100.0)],
                      below=[(1.8, 1.8, 200.0), (M.ITO, M.ITO, 50.0)], n_sub=N_SUB)


base = {}
for name, refl in (('Al', AL), ('Ag', M.AG)):
    S = dev(refl)
    r = cps2.solve_pol(S, npts=16000)
    base[name] = (r['air'] + r['sub'], cps2.stack_reflectance(S, TH),
                  cps2.sub_angular(S, TH))
    A = 1 - np.sum(base[name][1] * W) / np.sum(W)
    print('  %s : eta_sub = %.4f   A_prime = %.4f' % (name, base[name][0], A),
          file=sys.stderr, flush=True)

rows = []
BT_store, BR_store = {}, {}
for x in np.concatenate([[0.02], AR]):
    BT, BR = rt.bsdf(float(x), N_SUB, N=NRAY, max_events=80, seed=17)
    BT_store[float(x)] = BT
    BR_store[float(x)] = BR
    p_rand = float(np.sum(BT * W) / np.sum(W))
    row = [x, p_rand]
    for name in ('Al', 'Ag'):
        es, R, P = base[name]
        e, _ = series.eta_ext(BT, BR, R, P, n_term=80)
        row += [e, es * e]
    rows.append(row)
    print('  AR = %.2f done  (p_rand %.4f, EQE_Ag %.4f)' % (x, p_rand, row[5]),
          file=sys.stderr, flush=True)

np.savez_compressed('mla_bsdf.npz', aspect=np.array(list(BT_store)),
                    BT=np.array([BT_store[k] for k in BT_store]),
                    BR=np.array([BR_store[k] for k in BR_store]),
                    theta_deg=TH, n_mla=N_SUB)

with open('mla_aspect.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['aspect_ratio', 'p_single_pass', 'eta_ext_Al', 'EQE_Al',
                'eta_ext_Ag', 'EQE_Ag'])
    for r in rows:
        w.writerow(['%.2f' % r[0]] + ['%.6f' % x for x in r[1:]])

a = np.array(rows)
print('\n  AR     p      eta_ext(Al)  EQE(Al)   eta_ext(Ag)  EQE(Ag)')
for r in rows:
    print('  %.2f  %.4f    %.4f     %.4f     %.4f     %.4f' % tuple(r))
m = a[a[:, 0] >= 0.1]
for k, name in ((5, 'Ag'), (3, 'Al')):
    i = int(np.argmax(m[:, k]))
    thr = m[m[:, k] >= 0.95 * m[i, k], 0][0]
    print('  %s : best EQE %.4f at AR = %.2f, within 5 %% of it from AR = %.2f'
          % (name, m[i, k], m[i, 0], thr))
