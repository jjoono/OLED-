"""Sensitivity of the green BE device (device_eqe.py stack, cone array, Theta 0.76, PLQY 1) to a 5 nm Cs-rich layer
(80 wt% Cs2CO3 in B3PyMPM) at the Ag interface: B3PyMPM 165 + Cs-rich 5 instead of B3PyMPM 170 (ETL total unchanged).
Cs-rich layer: n 1.7 (assumed), k swept (unknown; Cs suboxide / free-carrier absorption)."""
import numpy as np, sys
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/be_cone', '/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/mla']
import device_eqe as D, cps2, series
BT, BR = D.EXT['cone']


def one(a):
    k, j = a
    cs = 1.7 + 1j * k
    S = cps2.Stack(float(D.LAM[j]), (D.EMo[j], D.EMe[j]), 25.0, 12.5,
                   above=[(D.B3o[j], D.B3e[j], 10.0), (D.B3o[j], D.B3e[j], 165.0), (cs, cs, 5.0), (D.AG[j], D.AG[j], 100.0)],
                   below=[(D.TCo[j], D.TCe[j], 10.0), (D.TPo[j], D.TPe[j], 200.0), (D.HAT[j], D.HAT[j], 5.0), (D.ITO[j], D.ITO[j], 50.0)],
                   n_sub=D.NS, h=0.76)
    r = cps2.solve_pol(S, npts=3000); R = np.real(cps2.stack_reflectance(S, D.TH)); P = np.real(cps2.sub_angular(S, D.TH))
    es = float(np.real(r['air'] + r['sub'])); Pn = P / P.sum()
    return es, float(1 - Pn @ R), series.eta_ext(BT, BR, R, P, 300)[0]


if __name__ == '__main__':
    js = list(range(0, 301, 4)); I = np.clip(np.asarray(D.spec.l_I_Irppy2acac, float)[:301], 0, None)
    w = (D.LAM * I)[js]; w /= w.sum(); base = None
    for k in (0.0, 0.02, 0.05, 0.1, 0.2, 0.5):
        with Pool(4) as p: r = np.array(p.map(one, [(k, j) for j in js]))
        es, Ap = w @ r[:, 0], w @ r[:, 1]; eqe = w @ (r[:, 0] * r[:, 2])
        base = base or eqe
        print('k %.2f: eta_sub %.4f  A\' %.4f  eta_ext %.4f  EQE %.4f  (change %+.2f %%p)' % (k, es, Ap, eqe / es, eqe, 100 * (eqe - base)), flush=True)
