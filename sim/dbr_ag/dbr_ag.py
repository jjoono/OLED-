"""Fig. 2(j)-style comparison: ITO 50 nm + chirped ZnS/LiF DBR with air behind it, the same DBR
backed by 100 nm Ag, and a plain 100 nm Ag reflector.  Common part as in Fig. 2(g)-(j):
glass n = 1.5 / ITO 150 nm / 420 nm organics (n = 1.8) / reflector; green spectrum 430-700 nm."""
import numpy as np, sys, os, csv
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'design_rule4'))
import fig2d as F
LAM = np.arange(430.0, 701.0); sel = np.isin(F.LAM, LAM)
common = [(F.ITO, 150.0), (F.ORG, F.ORG_D)]
def S(extra, exit_=F.AIR): return dict(layers=common + extra, exit=exit_)
chirp = F.dbr_chirp(10, 55.0, 91.0, 1.40)
T = {'ITO50 + chirped DBR (air)': S([(F.ITO, 50.0)] + chirp),
     'ITO50 + chirped DBR + Ag': S([(F.ITO, 50.0)] + chirp + [(F.AG, 100.0)]),
     'Ag': S([(F.AG, 100.0)])}
# DBR+Ag re-optimised: fewer pairs and other thicknesses, since Ag now closes the escape cone
LC = np.arange(430.0, 701.0, 10.0); selc = np.isin(F.LAM, LC)
thc = np.linspace(0, np.pi/2, 91); thc = 0.5*(thc[1:]+thc[:-1]); wc = np.cos(thc)*np.sin(thc); wc /= wc.sum()
spc = np.clip(np.interp(LC, F.LAM, F.GREEN), 0, None); spc /= spc.sum()
def loss(args):
    n, dz, dl, fac = args
    lay = [(F.ITO, 50.0)] + [l for l in F.dbr_chirp(n, dz, dl, fac) if l[1] > 0] + [(F.AG, 100.0)]
    R, _ = F.RT([(m[selc], d) for m, d in common + lay], np.full(LC.shape, 1.5+0j), F.AIR[selc], thc, LC)
    return (float(spc @ ((1-R) @ wc)), n, dz, dl, fac)
from multiprocessing import Pool
if True:
    grid = [(n, dz, dl, fac) for n in (1, 2, 3, 4, 6) for dz in range(0, 121, 10) for dl in range(0, 201, 20) for fac in (1.0, 1.4)]
    with Pool(4) as P: best = min(P.map(loss, grid, chunksize=20))
_, n, dz, dl, fac = best
T['ITO50 + DBR + Ag (re-opt)'] = S([(F.ITO, 50.0)] + [l for l in F.dbr_chirp(n, dz, dl, fac) if l[1] > 0] + [(F.AG, 100.0)])
T['LiF only + Ag (no ITO50)'] = S([(F.LIF, 100.0), (F.AG, 100.0)])
p = 0.43
th = np.linspace(0, np.pi/2, 181); th = 0.5*(th[1:]+th[:-1])
spec = np.clip(np.interp(LAM, F.LAM, F.GREEN), 0, None); spec /= spec.sum()
rows = []
for k, s in T.items():
    W = F.weighted(k, 1.5, F.GREEN, LAM, T)
    R, Tr = F.RT([(m[sel], d) for m, d in s['layers']], np.full(LAM.shape, 1.5+0j), s['exit'][sel], th, LAM)
    rows.append((k, W, spec @ (1-R), spec @ Tr))
    A = W['loss']; print(f"{k:32s} absorbed {100*W['absorbed']:.2f}  transmitted {100*W['transmitted']:.2f}  A' {100*A:.2f}  eta_ext(eq3,p=0.43) {p/(p+(1-p)*A):.3f}")
print('re-opt DBR+Ag:', best)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dbr_ag_curves.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['theta_deg'] + [f'{k}: 1-R (%)' for k, *_ in rows] + [f'{k}: T (%)' for k, *_ in rows] + ['cos.sin (norm)'])
    cs = np.cos(th)*np.sin(th); cs /= cs.max()
    for i in range(len(th)): w.writerow([np.degrees(th[i])] + [100*r[2][i] for r in rows] + [100*r[3][i] for r in rows] + [cs[i]])
