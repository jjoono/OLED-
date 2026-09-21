#!/usr/bin/env python3
"""The lens-free reference for the aspect-ratio panel.

A flat substrate/air interface does not randomise anything: light inside the
escape cone leaves with the Fresnel transmittance and everything else comes
straight back at the same polar angle.  So B_T is the unpolarised Fresnel
transmittance and B_R is diagonal, and the same recycling series gives the
reference EQE.  It cannot be read off the AR -> 0 end of the sweep, because even
a 2 %-aspect-ratio lens still redirects light.
"""
import numpy as np, csv, sys, os
sys.path.insert(0, '/home/user/OLED-/sim/fig3c')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cps2, materials as M, series

AL = 0.958336812040867 + 6.68678039868499j
N_SUB = 1.80
TH = np.arange(90) + 0.5
W = np.cos(np.radians(TH)) * np.sin(np.radians(TH))
HERE = os.path.dirname(os.path.abspath(__file__))


def dev(refl):
    return cps2.Stack(550.0, (1.8, 1.8), 20.0, 10.0,
                      above=[(1.8, 1.8, 200.0), (refl, refl, 100.0)],
                      below=[(1.8, 1.8, 200.0), (M.ITO, M.ITO, 50.0)], n_sub=N_SUB)


def fresnel_bsdf(n):
    ci = np.cos(np.radians(TH))
    s2 = n ** 2 * (1.0 - ci ** 2)
    ct = np.sqrt(np.clip(1.0 - s2, 0.0, None))
    rs = (n * ci - ct) / (n * ci + ct)
    rp = (n * ct - ci) / (n * ct + ci)
    T = np.where(s2 >= 1.0, 0.0, np.clip(1.0 - 0.5 * (rs ** 2 + rp ** 2), 0.0, 1.0))
    return T, np.diag(1.0 - T)


BT, BR = fresnel_bsdf(N_SUB)
out = {}
for name, refl in (('Al', AL), ('Ag', M.AG)):
    S = dev(refl)
    r = cps2.solve_pol(S, npts=16000)
    es = r['air'] + r['sub']
    R, P = cps2.stack_reflectance(S, TH), cps2.sub_angular(S, TH)
    e, _ = series.eta_ext(BT, BR, R, P, n_term=400)
    out[name] = (es, 1 - np.sum(R * W) / np.sum(W), float(np.sum(BT * W) / np.sum(W)),
                 e, es * e)
    print('  %s : eta_sub %.4f  A_prime %.4f  p_flat %.4f  eta_ext %.4f  EQE %.4f'
          % ((name,) + out[name]))

with open(os.path.join(HERE, 'mla_flat_reference.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['reflector', 'eta_sub', 'A_prime', 'p_flat', 'eta_ext_flat', 'EQE_flat'])
    for k, v in out.items():
        w.writerow([k] + ['%.6f' % x for x in v])
print('mla_flat_reference.csv written')
