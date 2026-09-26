"""For each n_sub (= n_MLA), find the ETL/HTL thicknesses that maximise the EQE, separately for
eq. (3) and for the matrix series, and compare the optima (Fig. 2(d)-(f) stack otherwise:
ITO 50 nm, 20 nm EML with centred isotropic dipole, all organics n = 1.8, 100 nm reflector, 550 nm).
Exhaustive grid (global on the grid) followed by a 1 nm local refinement around the best point."""
import numpy as np, scipy.io as sio, sys, csv, os
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/mla']
import cps2, series, materials as M
TH = np.arange(90) + 0.5; W = np.cos(np.radians(TH)) * np.sin(np.radians(TH))
B = sio.loadmat('/home/user/OLED-/sim/mla/lt_hemisphere_bsdf.mat')['BSDF_MLA']; NS = np.round(np.arange(1.30, 2.001, 0.05), 2)
def bsdf(n):
    b = B[:, :, int(np.argmin(abs(NS - n)))]; return b[:90].sum(0), b[90:180][::-1]
def ev(refl, n, de, dh):
    S = cps2.Stack(550.0, (1.8, 1.8), 20.0, 10.0, above=[(1.8, 1.8, de), (refl, refl, 100.0)],
                   below=[(1.8, 1.8, dh), (M.ITO, M.ITO, 50.0)], n_sub=n)
    BT, BR = bsdf(n); r = cps2.solve_pol(S, npts=4000); es = float(np.real(r['air'] + r['sub']))
    R = np.real(cps2.stack_reflectance(S, TH)); P = np.real(cps2.sub_angular(S, TH))
    A = 1 - (R @ W) / W.sum(); p = (BT @ W) / W.sum(); e3 = p / (p + (1 - p) * A)
    es_, _ = series.eta_ext(BT, BR, R, P, 300)
    return es, e3, es_, float(np.sum(P[TH > 70]) / P.sum())
def job(a):
    name, n, de, dh = a; return (name, n, de, dh) + ev(M.AG if name == 'Ag' else M.AL, n, de, dh)
if __name__ == '__main__':
    ETL = np.arange(20, 401, 10.0); HTL = np.arange(20, 401, 20.0)
    jobs = [(m, float(n), de, dh) for m in ('Ag', 'Al') for n in NS for de in ETL for dh in HTL]
    with Pool(4) as p: res = p.map(job, jobs, chunksize=50)
    out = [['reflector', 'n_sub', 'method', 'd_ETL', 'd_HTL', 'eta_sub', 'eta_ext', 'EQE', 'eta_ext_other_method_at_same_d', 'Psub_beyond70', 'max_eta_ext_any_d']]
    for m in ('Ag', 'Al'):
        for n in NS:
            rs = [r for r in res if r[0] == m and r[1] == n]
            for meth, i in (('eq3', 5), ('series', 6)):
                b = max(rs, key=lambda r: r[4] * r[i])
                # local refinement +-10 nm in 2 nm steps
                loc = [(b[2] + x, b[3] + y) for x in range(-10, 11, 2) for y in range(-20, 21, 4) if b[2] + x > 5 and b[3] + y > 5]
                with Pool(4) as p: rl = p.map(job, [(m, n, x, y) for x, y in loc])
                b = max(rl + [b], key=lambda r: r[4] * r[i])
                mx = max(r[i] for r in rs)
                out.append([m, n, meth, b[2], b[3], b[4], b[i], b[4] * b[i], b[11 - i], b[7], mx])
                print(out[-1], flush=True)
    with open('/home/user/OLED-/sim/opt_nsub/opt_nsub.csv', 'w', newline='') as f: csv.writer(f).writerows(out)
