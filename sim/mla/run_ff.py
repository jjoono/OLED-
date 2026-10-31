#!/usr/bin/env python3
"""Hexagonal MLA fill factor sweep at aspect ratio 0.5 and 1.0 (Family A: spherical cap of fixed base radius r).
Fill factor FF = (pi / (2 sqrt 3)) (2 r / pitch)^2, i.e. pitch = 2 r sqrt(0.9069 / FF); FF 0.907 = touching.
Same device as run_aspect_pitch.py (550 nm, n_sub 1.8, Al or Ag reflector, laterally infinite, 80 round trips).
Usage: run_ff.py <out.csv>"""
import numpy as np, csv, sys, os
sys.path.insert(0, '/home/user/OLED-/sim/fig3c'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concurrent.futures import ProcessPoolExecutor
import cps2, materials as M, raytrace as rt, series
OUT = sys.argv[1]
AL = 0.958336812040867 + 6.68678039868499j
N_SUB = 1.80; TH = np.arange(90) + 0.5; W = np.cos(np.radians(TH)) * np.sin(np.radians(TH)); NRAY = 30000
FFS = (0.9069, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
def dev(refl):
    return cps2.Stack(550.0, (1.8, 1.8), 20.0, 10.0, above=[(1.8, 1.8, 200.0), (refl, refl, 100.0)],
                      below=[(1.8, 1.8, 200.0), (M.ITO, M.ITO, 50.0)], n_sub=N_SUB)
def _one_bin(a):
    ar, i, th, pitch = a
    rt.PITCH = pitch
    rng = np.random.default_rng([17, i, 9])
    esc, rc = rt._trace_chunk(ar, N_SUB, NRAY, 80, rng, float(th), 'cap')
    good = np.isfinite(rc); ang = np.degrees(np.arccos(np.clip(rc[good], 0.0, 1.0)))
    h, _ = np.histogram(ang, bins=np.linspace(0.0, 90.0, 91)); return i, esc.mean(), h / NRAY
def main():
    base = {}
    for name, refl in (('Al', AL), ('Ag', M.AG)):
        S = dev(refl); r = cps2.solve_pol(S, npts=16000)
        base[name] = (r['air'] + r['sub'], cps2.stack_reflectance(S, TH), cps2.sub_angular(S, TH))
    rows = []
    with ProcessPoolExecutor(4) as ex:
        for ff in FFS:
            pitch = 2.0 * np.sqrt(0.9069 / ff)
            for ar in (0.5, 1.0):
                BT = np.zeros(90); BR = np.zeros((90, 90))
                for i, t, col in ex.map(_one_bin, [(ar, i, TH[i], pitch) for i in range(90)]):
                    BT[i] = t; BR[:, i] = col
                row = [ff, pitch, ar, float(np.sum(BT * W) / np.sum(W))]
                for name in ('Al', 'Ag'):
                    es, R, P = base[name]; e, _ = series.eta_ext(BT, BR, R, P, n_term=80); row += [e, es * e]
                rows.append(row)
                print('FF %.3f pitch %.3f AR %.1f  p %.4f  EQE_Al %.4f  EQE_Ag %.4f' % (ff, pitch, ar, row[3], row[5], row[7]), flush=True)
    with open(OUT, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['fill_factor', 'pitch_over_r', 'aspect_ratio', 'p_single_pass', 'eta_ext_Al', 'EQE_Al', 'eta_ext_Ag', 'EQE_Ag'])
        w.writerows(rows)
if __name__ == '__main__':
    main()
