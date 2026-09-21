#!/usr/bin/env python3
"""EQE against microlens aspect ratio, generic stack, full recycling series.

    reflector 100 nm | ETL n=1.8 200 nm | EML n=1.8 20 nm, dipole at the centre
                     | HTL n=1.8 200 nm | TCO 1.8636+0.0032i 50 nm | substrate 1.80
                     + hexagonally close-packed lenses, index-matched to the
                       substrate, infinite planar substrate and infinite array

550 nm, PLQY = 1, isotropic emitter.

Lens shape.  Up to AR = 1 the lens is a spherical cap of base radius r and
height h -- the shape a reflowed lens takes, and the one that is close-packed
without overlap because its sphere centre, z = h - R with R = (r^2+h^2)/2h,
sits below the base plane.  At AR = 1 the cap is exactly a hemisphere.  Beyond
AR = 1 the centre rises above the base plane, the widest circle of the sphere
becomes R > r, and the caps would cut into their neighbours, so a taller lens
is modelled as a half-ellipsoid (x^2+y^2)/r^2 + z^2/h^2 = 1, whose widest
circle is the base circle for every height.  The two families meet exactly at
the hemisphere, so the curve is continuous there by construction.

mla_aspect_alt.csv keeps both families outside their own range, for the record:
the overlapping cap above AR = 1 (it puts a spurious 3 %p dip at AR ~ 1.1) and
the half-ellipsoid below it (its rim meets the substrate vertically, unlike a
reflowed lens).

The closed form eta_ext = p/[p+(1-p)A'] is NOT used here: it assumes the surface
randomises the angles at every bounce, and a shallow lens does not, which is
exactly what the low-aspect-ratio end of this sweep is about.  The angle-resolved
B_T and B_R come from the ray trace and the series is summed term by term.
"""
import numpy as np, csv, sys, os
sys.path.insert(0, '/home/user/OLED-/sim/fig3c')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concurrent.futures import ProcessPoolExecutor
import cps2, materials as M
import raytrace as rt, series

AL = 0.958336812040867 + 6.68678039868499j
N_SUB = 1.80
TH = np.arange(90) + 0.5
W = np.cos(np.radians(TH)) * np.sin(np.radians(TH))
NRAY = 120000          # rays per 1-degree incidence bin
NPROC = 4
HERE = os.path.dirname(os.path.abspath(__file__))

AR_CAP = np.round(np.arange(0.10, 1.0001, 0.05), 2)      # spherical cap branch
AR_ELL = np.round(np.arange(1.05, 1.5001, 0.05), 2)      # half-ellipsoid branch
ALT_CAP = (1.05, 1.10, 1.15, 1.20, 1.30, 1.40, 1.50)     # overlapping cap, for the record
ALT_ELL = (0.20, 0.40, 0.60, 0.80, 0.90)                 # ellipsoid below AR = 1


def shape_of(ar):
    return 'cap' if ar <= 1.0 else 'ellipsoid'


def dev(refl):
    return cps2.Stack(550.0, (1.8, 1.8), 20.0, 10.0,
                      above=[(1.8, 1.8, 200.0), (refl, refl, 100.0)],
                      below=[(1.8, 1.8, 200.0), (M.ITO, M.ITO, 50.0)], n_sub=N_SUB)


def _one_bin(a):
    """One incidence bin of one aspect ratio, so the 90 bins can go out to the pool."""
    ar, i, th, shape = a
    # deterministic, and the same stream at a given bin for every aspect ratio:
    # common random numbers, so neighbouring aspect ratios move together
    rng = np.random.default_rng([17, i, 0 if shape == 'ellipsoid' else 1])
    esc, rc = rt._trace_chunk(ar, N_SUB, NRAY, 80, rng, float(th), shape)
    good = np.isfinite(rc)
    ang = np.degrees(np.arccos(np.clip(rc[good], 0.0, 1.0)))
    h, _ = np.histogram(ang, bins=np.linspace(0.0, 90.0, 91))
    return i, esc.mean(), h / NRAY


def bsdf_par(ar, ex, shape):
    BT = np.zeros(90)
    BR = np.zeros((90, 90))
    for i, t, col in ex.map(_one_bin, [(ar, i, TH[i], shape) for i in range(90)]):
        BT[i] = t
        BR[:, i] = col
    return BT, BR


def main():
    base = {}
    for name, refl in (('Al', AL), ('Ag', M.AG)):
        S = dev(refl)
        r = cps2.solve_pol(S, npts=16000)
        base[name] = (r['air'] + r['sub'], cps2.stack_reflectance(S, TH),
                      cps2.sub_angular(S, TH))
        A = 1 - np.sum(base[name][1] * W) / np.sum(W)
        print('  %s : eta_sub = %.4f   A_prime = %.4f' % (name, base[name][0], A),
              file=sys.stderr, flush=True)

    def evaluate(ar, shape, ex):
        BT, BR = bsdf_par(float(ar), ex, shape)
        row = [ar, float(np.sum(BT * W) / np.sum(W))]
        for name in ('Al', 'Ag'):
            es, R, P = base[name]
            e, _ = series.eta_ext(BT, BR, R, P, n_term=80)
            row += [e, es * e]
        print('  AR %.2f %-9s p %.4f  EQE_Al %.4f  EQE_Ag %.4f  closure %.6f'
              % (ar, shape, row[1], row[3], row[5], float(np.mean(BT + BR.sum(0)))),
              file=sys.stderr, flush=True)
        return row, BT, BR

    rows, alt_rows = [], []
    BT_store, BR_store = {}, {}
    with ProcessPoolExecutor(NPROC) as ex:
        for x in np.concatenate([[0.02], AR_CAP, AR_ELL]):
            row, BT, BR = evaluate(float(x), shape_of(float(x)), ex)
            rows.append(row + [shape_of(float(x))])
            BT_store[float(x)], BR_store[float(x)] = BT, BR
        for x in ALT_CAP:
            alt_rows.append(evaluate(float(x), 'cap', ex)[0] + ['cap (overlapping)'])
        for x in ALT_ELL:
            alt_rows.append(evaluate(float(x), 'ellipsoid', ex)[0] + ['ellipsoid'])

    np.savez_compressed(os.path.join(HERE, 'mla_bsdf.npz'),
                        aspect=np.array(list(BT_store)),
                        BT=np.array([BT_store[k] for k in BT_store]),
                        BR=np.array([BR_store[k] for k in BR_store]),
                        theta_deg=TH, n_mla=N_SUB, n_ray=NRAY)

    head = ['aspect_ratio', 'p_single_pass', 'eta_ext_Al', 'EQE_Al',
            'eta_ext_Ag', 'EQE_Ag', 'lens_shape']
    for fn, data in (('mla_aspect.csv', rows), ('mla_aspect_alt.csv', alt_rows)):
        with open(os.path.join(HERE, fn), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(head)
            for r in data:
                w.writerow(['%.2f' % r[0]] + ['%.6f' % x for x in r[1:6]] + [r[6]])

    a = np.array([r[:6] for r in rows])
    m = a[a[:, 0] >= 0.1]
    print('\n  AR     p      eta_ext(Al)  EQE(Al)   eta_ext(Ag)  EQE(Ag)')
    for r in rows:
        print('  %.2f  %.4f    %.4f     %.4f     %.4f     %.4f   %s'
              % tuple(r[:6] + [r[6]]))
    for k, name in ((5, 'Ag'), (3, 'Al')):
        i = int(np.argmax(m[:, k]))
        thr = m[m[:, k] >= 0.95 * m[i, k], 0][0]
        d2 = np.diff(m[:, k], 2)
        print('  %s : best EQE %.4f at AR = %.2f, within 5 %% of it from AR = %.2f,'
              '  RMS 2nd difference %.2e'
              % (name, m[i, k], m[i, 0], thr, np.sqrt(np.mean(d2 ** 2))))
    print('\n  the other shape family outside its own range, for the record')
    for r in alt_rows:
        print('   AR %.2f %-18s p %.4f  EQE_Al %.4f  EQE_Ag %.4f'
              % (r[0], r[6], r[1], r[3], r[5]))


if __name__ == '__main__':
    main()
