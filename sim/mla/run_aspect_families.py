#!/usr/bin/env python3
"""How the EQE-vs-aspect-ratio curve depends on how the lens shape is defined.

Three constructions, all hexagonal with pitch 2r and index-matched to the substrate:

  A  cap/ellipsoid   spherical cap of FIXED base radius r (h = AR r), close-packed at
                     every AR; half-ellipsoid above AR = 1.   (mla_aspect.csv, high
                     statistics -- the panel as drawn)
  B  trunc/bullet    sphere of FIXED radius r = half pitch, cut at height h, so the base
                     circle shrinks below the hemisphere (r_base = 2AR r/(1+AR^2)) and the
                     flat gap grows; above AR = 1 a cylinder of radius r topped by a
                     hemisphere.  This is the usual way of building an MLA of variable
                     height in a ray-tracing package.
  C  ellipsoid       a hemisphere scaled in z by AR at every AR.

Same generic stack as run_aspect.py; 30k rays per incidence bin.
Usage: run_aspect_families.py [n_sub] [families, e.g. A,B] [output csv]
"""
import numpy as np, csv, sys, os
sys.path.insert(0, '/home/user/OLED-/sim/fig3c')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concurrent.futures import ProcessPoolExecutor
import cps2, materials as M
import raytrace as rt, series

AL = 0.958336812040867 + 6.68678039868499j
N_SUB = float(sys.argv[1]) if len(sys.argv) > 1 else 1.80
FAMS = sys.argv[2].split(',') if len(sys.argv) > 2 else ['B', 'C']
OUT = sys.argv[3] if len(sys.argv) > 3 else 'mla_aspect_families.csv'
TH = np.arange(90) + 0.5
W = np.cos(np.radians(TH)) * np.sin(np.radians(TH))
NRAY = 30000
HERE = os.path.dirname(os.path.abspath(__file__))
AR = np.round(np.arange(0.2, 1.5001, 0.1), 2)


def dev(refl):
    return cps2.Stack(550.0, (1.8, 1.8), 20.0, 10.0,
                      above=[(1.8, 1.8, 200.0), (refl, refl, 100.0)],
                      below=[(1.8, 1.8, 200.0), (M.ITO, M.ITO, 50.0)], n_sub=N_SUB)


def _one_bin(a):
    ar, i, th, shape = a
    rng = np.random.default_rng([17, i, 7])
    esc, rc = rt._trace_chunk(ar, N_SUB, NRAY, 80, rng, float(th), shape)
    good = np.isfinite(rc)
    ang = np.degrees(np.arccos(np.clip(rc[good], 0.0, 1.0)))
    h, _ = np.histogram(ang, bins=np.linspace(0.0, 90.0, 91))
    return i, esc.mean(), h / NRAY


def bsdf_par(ar, ex, shape):
    BT = np.zeros(90); BR = np.zeros((90, 90))
    for i, t, col in ex.map(_one_bin, [(ar, i, TH[i], shape) for i in range(90)]):
        BT[i] = t; BR[:, i] = col
    return BT, BR


def main():
    base = {}
    for name, refl in (('Al', AL), ('Ag', M.AG)):
        S = dev(refl); r = cps2.solve_pol(S, npts=16000)
        base[name] = (r['air'] + r['sub'], cps2.stack_reflectance(S, TH), cps2.sub_angular(S, TH))
    allf = {'A': lambda ar: 'cap' if ar <= 1.0 else 'ellipsoid',
            'B': lambda ar: 'trunc' if ar <= 1.0 else 'bullet',
            'C': lambda ar: 'ellipsoid'}
    fams = {k: allf[k] for k in FAMS}
    rows = []
    with ProcessPoolExecutor(4) as ex:
        for fam, pick in fams.items():
            for x in AR:
                shape = pick(float(x))
                BT, BR = bsdf_par(float(x), ex, shape)
                row = [fam, shape, x, float(np.sum(BT * W) / np.sum(W))]
                for name in ('Al', 'Ag'):
                    es, R, P = base[name]
                    e, _ = series.eta_ext(BT, BR, R, P, n_term=80)
                    row += [e, es * e]
                rows.append(row)
                print('  %s %-9s AR %.2f  p %.4f  EQE_Al %.4f  EQE_Ag %.4f' % (fam, shape, x, row[3], row[5], row[7]),
                      file=sys.stderr, flush=True)
    with open(os.path.join(HERE, OUT), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['family', 'lens_shape', 'aspect_ratio', 'p_single_pass', 'eta_ext_Al', 'EQE_Al', 'eta_ext_Ag', 'EQE_Ag'])
        for r in rows:
            w.writerow(r[:2] + ['%.2f' % r[2]] + ['%.6f' % v for v in r[3:]])
    print('done')


if __name__ == '__main__':
    main()
