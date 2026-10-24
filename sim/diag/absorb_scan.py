"""Can extra parasitic absorption explain ref 20 % / MLA 51 % / MLA+reflector 58 % ?
For each absorption case the reference angular EL is re-fitted (fit_ref_angles: ETL scale, TAPC scale, Theta, IZO free,
dipole 1 nm from the ETL/EML interface), then the three EQEs are computed at that stack (dev3.three, PLQY 1).
  kd: k of the 60 nm B3PyMPM:Cs2CO3 10 wt% layer;  kh: k of the three 5 nm HATCN layers (HATCN/TAPC interface CT absorption)
Measured (device_eval, -90..90 deg): ref 0.200, MLA+reflector 0.579; MLA without reflector 0.51 (GUI value, not re-evaluated)."""
import sys, csv, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
sys.path[:0] = ['/home/user/OLED-/sim/diag']
import fit_ref_angles as F, dev3
Z = 1.0
CASES = [(0, 0)] + [(kd, 0) for kd in (0.005, 0.01, 0.02, 0.04)] + [(0, kh) for kh in (0.1, 0.2, 0.4, 0.8)]
MEAS = np.array([0.200, 0.51, 0.579])


def fit(c):
    kd, kh = c; S = F.measured()[0]; Rm = S / S[0]; w = (S[0] / S[0].max() > 0.15).astype(float)
    def cost(x):
        f, g, t, di = x
        M = F.model((f, g, Z, di, 1.0, kd, kh))[1]; I = (1 - t) * M[:, 0] + t * M[:, 1]; R = (I / I[:, :1]).T
        e = np.log(R[1:8] / Rm[1:8]) * w
        return float(np.sqrt((e ** 2).sum() / (w.sum() * 7)))
    best = None
    for x0 in ([0.94, 0.78, 0.71, 10.0], [0.92, 0.75, 0.72, 30.0]):
        r = minimize(cost, x0, method='Nelder-Mead', bounds=[(0.7, 1.2), (0.4, 1.5), (0.5, 1.0), (0.0, 90.0)],
                     options=dict(maxiter=300, xatol=1e-3, fatol=1e-6))
        if best is None or r.fun < best.fun: best = r
    return c, best.fun, best.x


if __name__ == '__main__':
    with Pool(4) as p:
        fits = p.map(fit, CASES)
    rows = []
    for (kd, kh), rms, (f, g, t, di) in fits:
        with Pool(4) as p:
            o = dev3.three(p, f=f, kd=kd, theta=t, geo=(g, Z, di, 1.0, kh), N=150000, nrun=8)
        m = np.array([o['ref'], o['cathode'], o['full']])
        rows.append(dict(kd=kd, kh=kh, rms=round(rms, 4), ETL=round(185 * f, 1), TAPC=round(210 * g, 1), Theta=round(t, 3), IZO=round(di, 1),
                         eta_sub177=round(o['es'], 3), ref=round(m[0], 3), MLA=round(m[1], 3), MLA_refl=round(m[2], 3),
                         meas_over_model=' / '.join('%.2f' % v for v in MEAS / m), MLA_ref=round(m[1] / m[0], 2), refl_ref=round(m[2] / m[0], 2)))
        print(rows[-1], flush=True)
    with open('absorb_scan.csv', 'w', newline='') as fo:
        wr = csv.DictWriter(fo, list(rows[0])); wr.writeheader(); wr.writerows(rows)
    print('measured MLA/ref %.2f, refl/ref %.2f' % (MEAS[1] / MEAS[0], MEAS[2] / MEAS[0]))
