#!/usr/bin/env python3
"""Does a small gap between neighbouring lenses change the hemisphere?  Full BSDF -> EQE at
AR 0.55 and 1.0 for pitch 2.0 (touching), 2.05, 2.1, 2.2 r.  30k rays per bin."""
import numpy as np, csv, sys, os
sys.path.insert(0, '/home/user/OLED-/sim/fig3c'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concurrent.futures import ProcessPoolExecutor
import cps2, materials as M, raytrace as rt, series
AL = 0.958336812040867 + 6.68678039868499j
N_SUB = 1.80; TH = np.arange(90) + 0.5; W = np.cos(np.radians(TH)) * np.sin(np.radians(TH)); NRAY = 30000
def dev(refl):
    return cps2.Stack(550.0, (1.8, 1.8), 20.0, 10.0, above=[(1.8, 1.8, 200.0), (refl, refl, 100.0)],
                      below=[(1.8, 1.8, 200.0), (M.ITO, M.ITO, 50.0)], n_sub=N_SUB)
def _one_bin(a):
    ar, i, th, pitch = a
    rt.PITCH = pitch
    rng = np.random.default_rng([17, i, 13])
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
        for pitch in (2.0, 2.05, 2.1, 2.2):
            for ar in (0.55, 1.0):
                BT = np.zeros(90); BR = np.zeros((90, 90))
                for i, t, col in ex.map(_one_bin, [(ar, i, TH[i], pitch) for i in range(90)]):
                    BT[i] = t; BR[:, i] = col
                row = [pitch, ar, float(np.sum(BT * W) / np.sum(W))]
                for name in ('Al', 'Ag'):
                    es, R, P = base[name]; e, _ = series.eta_ext(BT, BR, R, P, n_term=80); row += [e, es * e]
                rows.append(row)
                print('  pitch %.2f AR %.2f  p %.4f  EQE_Al %.4f  EQE_Ag %.4f' % (pitch, ar, row[2], row[4], row[6]), file=sys.stderr, flush=True)
    with open('mla_gap_test.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['pitch_over_r', 'aspect_ratio', 'p_single_pass', 'eta_ext_Al', 'EQE_Al', 'eta_ext_Ag', 'EQE_Ag'])
        for r in rows: w.writerow(['%.2f' % r[0], '%.2f' % r[1]] + ['%.6f' % v for v in r[2:]])
    print('done')
if __name__ == '__main__':
    main()
