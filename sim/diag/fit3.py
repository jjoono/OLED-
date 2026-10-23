"""Which single change (with a common internal-efficiency factor s) reproduces the measured 16 / 51 / 59 %?
Scans: ETL tooling f, absorption of the 50 wt% layer (k_cs), absorption of the 10 wt% layer (k_d), Theta."""
import numpy as np, sys, csv
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/diag']
import dev3

MEAS = np.array([0.16, 0.51, 0.59])
def fit(m):
    s = float(MEAS @ m / (m @ m)); res = MEAS - s * m
    return s, float(np.sqrt((res ** 2).mean())), s * m

rows = []
with Pool(4) as p:
    scans = [('ETL tooling f', 'f', np.arange(0.6, 1.41, 0.05)),
             ('k of 50 wt% layer (5 nm)', 'kcs', np.array([0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8])),
             ('k of 10 wt% layer (60 nm)', 'kd', np.array([0, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05])),
             ('Theta', 'theta', np.array([0.6, 0.67, 0.76, 0.85, 0.95]))]
    for lab, key, vals in scans:
        print(lab)
        for v in vals:
            o = dev3.three(p, **{key: float(v)}, N=40000, nrun=4)
            m = np.array([o['ref'], o['cathode'], o['full']]); s, rms, pred = fit(m)
            rows.append(dict(scan=lab, value=float(v), ref=m[0], mla=m[1], mla_refl=m[2], s=s, rms=rms,
                             pred_ref=pred[0], pred_mla=pred[1], pred_refl=pred[2]))
            print('  %6.3f: model %.3f %.3f %.3f | ratios %.2f %.2f | best s %.3f -> %.3f %.3f %.3f (rms %.4f)'
                  % (v, *m, m[1] / m[0], m[2] / m[1], s, *pred, rms), flush=True)
with open('fit3.csv', 'w', newline='') as fh:
    wr = csv.DictWriter(fh, list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
print('measured ratios: MLA/ref %.2f, refl/MLA %.2f' % (MEAS[1] / MEAS[0], MEAS[2] / MEAS[1]))
