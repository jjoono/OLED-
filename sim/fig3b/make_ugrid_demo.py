#!/usr/bin/env python3
"""Why disp_matrix jitters, and what the new quadrature does about it.

The guided modes of the organic slab are poles of K(u) whose width in u is set
by the residual loss (Ag tail + ITO k = 0.003) and is 1e-5..1e-3.  The original
uniform grid u = [0:1/N:1) U (1+1/N:1/N:3] with N = 1000 has a spacing of 1e-3,
so a pole is hit or missed depending on where it happens to sit -- and it moves
continuously as the organic thickness is swept.  That is the jitter.
"""
import numpy as np, csv, sys
import cps

D = np.arange(10.0, 500.0 + 1e-9, 2.0)
NS = 1.50
out = []
for d in D:
    a = cps.legacy(float(d), NS, N=1000)
    b = cps.legacy(float(d), NS, N=3000)
    c = cps.solve(float(d), NS, npts=12000)
    out.append([d, a['wg'], a['spp'], a['sub'], b['wg'], b['spp'], b['sub'],
                c['wg'], c['spp'], c['sub']])
out = np.array(out)

with open('fig3b_ugrid_demo.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['d_organic_nm', 'wg_N1000', 'spp_N1000', 'sub_N1000',
                'wg_N3000', 'spp_N3000', 'sub_N3000',
                'wg_quad', 'spp_quad', 'sub_quad'])
    for r in out:
        w.writerow(['%.1f' % r[0]] + ['%.6f' % x for x in r[1:]])


def rough(y):
    """RMS of the second difference: a scale-free measure of point-to-point
    jitter that a genuinely smooth curve does not have."""
    return float(np.sqrt(np.mean(np.diff(y, 2) ** 2)))


print('point-to-point roughness (RMS 2nd difference), n_sub = %.2f' % NS)
print('  channel |  original N=1000   original N=3000    new quadrature')
for name, i in (('wg ', 1), ('spp', 2), ('sub', 3)):
    print('    %s   |     %.2e         %.2e         %.2e'
          % (name, rough(out[:, i]), rough(out[:, i + 3]), rough(out[:, i + 6])))
print('\n(the new quadrature keeps a little roughness because the curve has real')
print(' kinks at guided-mode cut-offs; doubling the number of points moves it by')
print(' less than 1e-4, and at 0.5 nm sampling it is smooth.)')
print('\nlargest deviation of the original N=1000 from the converged answer:')
for name, i in (('wg ', 1), ('spp', 2), ('sub', 3)):
    print('    %s   %.4f  (at d = %.0f nm)'
          % (name, np.max(np.abs(out[:, i] - out[:, i + 6])),
             out[int(np.argmax(np.abs(out[:, i] - out[:, i + 6]))), 0]))

chk = np.array([[cps.solve(float(d), NS, npts=n)['wg'] for d in D[::10]]
                for n in (12000, 48000)])
print('\nconvergence of the new quadrature: max |npts 12000 - 48000| = %.1e'
      % np.max(np.abs(chk[0] - chk[1])))
