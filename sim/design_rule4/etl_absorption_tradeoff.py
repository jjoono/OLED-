"""Rule 4 with an absorbing ETL: EQE vs ETL thickness for k_ETL = 0 ... 3e-3 (k added to n_o and n_e of B3PyMPM).
Stack of Fig. 4c type: Ag 100 | B3PyMPM d_ETL | TCTA:B3PyMPM 25 (emitter at centre) | TAPC 200 | IZO 50 | substrate n 1.8
with a close-packed hemispherical MLA (n 1.8), laterally infinite device (matrix series, 100 round trips).
Theta 2/3, PLQY 1, emission spectrum of Ir(dmppy-ph)2tmd, 400-700 nm in 10-nm steps."""
import numpy as np, sys, csv
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/orange']
import orange_reflectors as O

NS = 1.8
C = np.load('/home/user/OLED-/sim/mla/cone_bsdf.npz'); B = C['hemi_cp_n1.8']; BT, BR = B[:90].sum(0), B[90:180][::-1]
DS = (30, 50, 70, 100, 150, 200, 300, 400)
KS = (0.0, 3e-4, 1e-3, 3e-3)


def point(args):
    k, d, j = args
    S = O.cps2.Stack(float(O.LAM[j]), (O.EMo[j], O.EMe[j]), 25.0, 12.5,
                     above=[(O.B3o[j] + 1j * k, O.B3e[j] + 1j * k, float(d)), (O.AG[j], O.AG[j], 100.0)],
                     below=[(O.TPo[j], O.TPe[j], 200.0), (O.IZO[j], O.IZO[j], 50.0)], n_sub=NS, h=2 / 3)
    r = O.cps2.solve_pol(S, npts=3000)
    return float(np.real(r['air'] + r['sub'])), np.real(O.cps2.sub_angular(S, O.TH)), np.real(O.cps2.stack_reflectance(S, O.TH))


if __name__ == '__main__':
    wp = (O.LAM * O.SPEC)[O.LAMJ]; wp = wp / wp.sum()
    rows = []
    for k in KS:
        for d in DS:
            with Pool(4) as p: res = p.map(point, [(k, d, j) for j in O.LAMJ])
            ES = np.array([r[0] for r in res]); PS = np.array([r[1] for r in res]); RS = np.array([r[2] for r in res])
            es = float(wp @ ES); WL = wp * ES / (wp @ ES); Pn = PS / PS.sum(1, keepdims=True)
            Ap = float(WL @ (1 - (Pn * RS).sum(1)))
            ext = float(sum(WL[i] * O.series.eta_ext(BT, BR, RS[i], PS[i], 100)[0] for i in range(len(O.LAMJ))))
            rows.append(dict(k_ETL=k, d_ETL=d, eta_sub=es, A_prime=Ap, eta_sta=ext, EQE=es * ext))
            print('k %.0e  d %3d | eta_sub %.3f  A\' %.3f  eta_sta %.3f  EQE %.3f' % (k, d, es, Ap, ext, es * ext), flush=True)
    with open('/home/user/OLED-/sim/design_rule4/etl_absorption_tradeoff.csv', 'w', newline='') as fh:
        wr = csv.DictWriter(fh, list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
