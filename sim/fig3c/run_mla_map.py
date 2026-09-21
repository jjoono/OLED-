#!/usr/bin/env python3
"""Converged planar map for the MLA device, on the d_ETL x d_HTL grid.

    air | Ag 100 | B3PyMPM (d_ETL) | TCTA:B3PyMPM 25 (dipole at the centre)
        | TAPC (d_HTL) | ITO 50 | glass 1.5        @ 550 nm, PLQY = 1

Same stack as Planar_sweep22_preprint.m with EML_position = 4, but with the
substitution quadrature instead of the uniform u grid, and with the five
channels binned per polarisation.  Exports everything the MLA recycling step
needs, so that EQE_MLA is a small matrix product once the BSDF is in hand.
"""
import numpy as np, csv, sys, scipy.io
import cps2, materials as M

LAM, NPTS = 550.0, 16000
B3o, B3e = 1.820533808087541, 1.608639620557389      # material.l_B3_o_JO / l_B3_e_JO
Eo, Ee   = 1.832730000000000, 1.671220000000000      # material.l_TCTA_B3_o_JO / _e_
To, Te   = 1.690750000000000, 1.663750000000000      # material.l_TAPC_o_JO / _e_
D_ETL = np.arange(10.0, 500.1, 10.0)
D_HTL = np.arange(10.0, 500.1, 10.0)
TH = np.arange(90) + 0.5                              # substrate angle, degrees


def stack(dE, dH):
    return cps2.Stack(LAM, (Eo, Ee), 25.0, 12.5,
                      above=[(B3o, B3e, dE), (M.AG, M.AG, 100.0)],
                      below=[(To, Te, dH), (M.ITO, M.ITO, 50.0)], n_sub=1.5)


nE, nH = len(D_ETL), len(D_HTL)
ch = np.zeros((nE, nH, 6))
P0 = np.zeros((nE, nH))
Psub = np.zeros((nE, nH, 90))
RLED = np.zeros((nE, nH, 90))

for a, dE in enumerate(D_ETL):
    for b, dH in enumerate(D_HTL):
        S = stack(float(dE), float(dH))
        r = cps2.solve_pol(S, npts=NPTS)
        ch[a, b] = [r['air'], r['sub'], r['wg'], r['spp'], r['abs'],
                    r['P_tot'] / 0.5]
        P0[a, b] = r['air'] + r['sub']
        Psub[a, b] = cps2.sub_angular(S, TH) / r['P_tot']
        RLED[a, b] = cps2.stack_reflectance(S, TH)
    print('  d_ETL = %3.0f nm done' % dE, file=sys.stderr, flush=True)

with open('mla_planar_map.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['d_ETL_nm', 'd_HTL_nm', 'air', 'sub_confined', 'wg', 'spp',
                'abs', 'eta_sub', 'Purcell', 'sum_of_channels'])
    for a, dE in enumerate(D_ETL):
        for b, dH in enumerate(D_HTL):
            c = ch[a, b]
            w.writerow(['%.0f' % dE, '%.0f' % dH]
                       + ['%.6f' % x for x in c[:5]]
                       + ['%.6f' % P0[a, b], '%.6f' % c[5],
                          '%.6f' % c[:5].sum()])

scipy.io.savemat('mla_inputs.mat', {
    'd_ETL': D_ETL, 'd_HTL': D_HTL, 'theta_deg': TH,
    'P0': P0, 'P_sub': Psub, 'R_LED': RLED,
    'channels': ch[:, :, :5], 'Purcell': ch[:, :, 5],
    'readme': ('P0(i,j) is the power delivered into the substrate mode; '
               'P_sub(i,j,:) is its angular distribution per degree, already '
               'normalised by the total dissipated power; R_LED(i,j,:) is the '
               'unpolarised stack reflectance seen from the substrate. '
               'channels = [air sub_confined wg spp abs].')}, do_compression=True)

print('mla_planar_map.csv and mla_inputs.mat written: %d x %d grid' % (nE, nH))
print('closure: max |sum-1| = %.2e' % np.max(np.abs(ch[:, :, :5].sum(2) - 1)))
print('eta_sub: min %.4f  max %.4f  (at d_ETL=%.0f, d_HTL=%.0f)'
      % (P0.min(), P0.max(), D_ETL[np.argmax(P0) // nH], D_HTL[np.argmax(P0) % nH]))


def rough(y):
    return float(np.sqrt(np.mean(np.diff(y, 2) ** 2)))


print('roughness of eta_sub along d_ETL (mean over d_HTL) = %.2e'
      % np.mean([rough(P0[:, b]) for b in range(nH)]))
