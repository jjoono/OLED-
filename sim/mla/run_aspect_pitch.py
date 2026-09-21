#!/usr/bin/env python3
"""Family A (fixed-base cap / ellipsoid) with the lenses pulled apart: pitch = 2.4 r and
2.8 r (fill factors 0.63 and 0.46) instead of touching (0.907).  Usage: run_aspect_pitch.py <pitch> <out.csv>"""
import numpy as np, csv, sys, os
sys.path.insert(0, '/home/user/OLED-/sim/fig3c'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concurrent.futures import ProcessPoolExecutor
import cps2, materials as M, raytrace as rt, series
PITCH = float(sys.argv[1]); OUT = sys.argv[2]
AL = 0.958336812040867 + 6.68678039868499j
N_SUB = 1.80; TH = np.arange(90) + 0.5; W = np.cos(np.radians(TH)) * np.sin(np.radians(TH)); NRAY = 30000
AR = np.round(np.arange(0.2, 1.5001, 0.1), 2)
def dev(refl):
    return cps2.Stack(550.0, (1.8, 1.8), 20.0, 10.0, above=[(1.8, 1.8, 200.0), (refl, refl, 100.0)],
                      below=[(1.8, 1.8, 200.0), (M.ITO, M.ITO, 50.0)], n_sub=N_SUB)
def _one_bin(a):
    ar, i, th, shape = a
    rt.PITCH = PITCH
    rng = np.random.default_rng([17, i, 9])
    esc, rc = rt._trace_chunk(ar, N_SUB, NRAY, 80, rng, float(th), shape)
    good = np.isfinite(rc); ang = np.degrees(np.arccos(np.clip(rc[good], 0.0, 1.0)))
    h, _ = np.histogram(ang, bins=np.linspace(0.0, 90.0, 91)); return i, esc.mean(), h / NRAY
def main():
    base = {}
    for name, refl in (('Al', AL), ('Ag', M.AG)):
        S = dev(refl); r = cps2.solve_pol(S, npts=16000)
        base[name] = (r['air'] + r['sub'], cps2.stack_reflectance(S, TH), cps2.sub_angular(S, TH))
    rows = []
    with ProcessPoolExecutor(4) as ex:
        for x in AR:
            shape = 'cap' if x <= 1.0 else 'ellipsoid'
            BT = np.zeros(90); BR = np.zeros((90, 90))
            for i, t, col in ex.map(_one_bin, [(float(x), i, TH[i], shape) for i in range(90)]):
                BT[i] = t; BR[:, i] = col
            row = [PITCH, shape, x, float(np.sum(BT * W) / np.sum(W))]
            for name in ('Al', 'Ag'):
                es, R, P = base[name]; e, _ = series.eta_ext(BT, BR, R, P, n_term=80); row += [e, es * e]
            rows.append(row)
            print('  pitch %.1f %-9s AR %.2f  p %.4f  EQE_Al %.4f  EQE_Ag %.4f  closure %.5f' % (PITCH, shape, x, row[3], row[5], row[7], float(np.mean(BT + BR.sum(0)))), file=sys.stderr, flush=True)
    with open(OUT, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['pitch_over_r', 'lens_shape', 'aspect_ratio', 'p_single_pass', 'eta_ext_Al', 'EQE_Al', 'eta_ext_Ag', 'EQE_Ag'])
        for r in rows: w.writerow([r[0], r[1], '%.2f' % r[2]] + ['%.6f' % v for v in r[3:]])
    print('done')
if __name__ == '__main__':
    main()
