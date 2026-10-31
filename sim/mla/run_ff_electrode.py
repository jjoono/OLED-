#!/usr/bin/env python3
"""Fill-factor sweep (run_ff.py geometry, same seeds) with the transparent electrode changed:
ITO 50 nm (reference, = run_ff.py) | Ag 8 nm with bulk Ag n,k (continuous-film limit) | Ag 8 nm with the measured
HATCN-seeded 8-nm Ag n,k (teoled_opt/nk_joint).  Ag reflector 100 nm; 550 nm, n_sub 1.8, laterally infinite, 80 round trips.
The BSDFs are stored in mla_fill_factor_bsdf.npz.  Usage: run_ff_electrode.py <out.csv>"""
import numpy as np, csv, sys, os
sys.path.insert(0, '/home/user/OLED-/sim/fig3c'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concurrent.futures import ProcessPoolExecutor
import cps2, materials as M, series
from run_ff import _one_bin, FFS, TH, W, N_SUB
OUT = sys.argv[1]
a = np.genfromtxt('/home/user/OLED-/sim/teoled_opt/nk_joint/Ag_HATCN_2-6_8nm_nk.csv', delimiter=',', names=True)
AG8 = complex(np.interp(550, a['wavelength_nm'], a['n_joint']), np.interp(550, a['wavelength_nm'], a['k_joint']))
ELEC = {'ITO50': (M.ITO, 50.0), 'Ag8_bulk': (M.AG, 8.0), 'Ag8_meas': (AG8, 8.0)}
def dev(tc, d):
    return cps2.Stack(550.0, (1.8, 1.8), 20.0, 10.0, above=[(1.8, 1.8, 200.0), (M.AG, M.AG, 100.0)],
                      below=[(1.8, 1.8, 200.0), (tc, tc, d)], n_sub=N_SUB)
def main():
    base = {}
    for name, (tc, d) in ELEC.items():
        S = dev(tc, d); r = cps2.solve_pol(S, npts=16000)
        base[name] = (float(np.real(r['air'] + r['sub'])), np.real(cps2.stack_reflectance(S, TH)), np.real(cps2.sub_angular(S, TH)))
        Pn = base[name][2] / base[name][2].sum()
        print('%-9s n,k %s  eta_sub %.3f  A\' %.3f' % (name, np.round(tc, 3), base[name][0], 1 - float((Pn * base[name][1]).sum())), flush=True)
    rows, store = [], {}
    with ProcessPoolExecutor(4) as ex:
        for ff in FFS:
            pitch = 2.0 * np.sqrt(0.9069 / ff)
            for ar in (0.5, 1.0):
                BT = np.zeros(90); BR = np.zeros((90, 90))
                for i, t, col in ex.map(_one_bin, [(ar, i, TH[i], pitch) for i in range(90)]):
                    BT[i] = t; BR[:, i] = col
                store['BT_ff%.3f_ar%.1f' % (ff, ar)] = BT; store['BR_ff%.3f_ar%.1f' % (ff, ar)] = BR
                row = [ff, pitch, ar, float(np.sum(BT * W) / np.sum(W))]
                for name in ELEC:
                    es, R, P = base[name]; e, _ = series.eta_ext(BT, BR, R, P, n_term=80); row += [es, e, es * e]
                rows.append(row)
                print('FF %.3f AR %.1f  ' % (ff, ar) + '  '.join('%s %.4f' % (n, row[6 + 3 * k]) for k, n in enumerate(ELEC)), flush=True)
    np.savez_compressed(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mla_fill_factor_bsdf.npz'), **store)
    with open(OUT, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['fill_factor', 'pitch_over_r', 'aspect_ratio', 'p_single_pass'] +
                                      [c % n for n in ELEC for c in ('eta_sub_%s', 'eta_ext_%s', 'EQE_%s')])
        w.writerows(rows)
if __name__ == '__main__':
    main()
