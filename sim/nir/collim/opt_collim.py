"""DE optimisation of the top DBR (free TiO2/SiO2 thicknesses, TiO2 first next to the slab) and slab choice for a
forward-directed NIR OLED.  Usage: python3 opt_collim.py <dev.npz> <lam0> <fwhm> <q> <npairs> <objective> <tag>
objective: Leq10 | C20 | E_z3_a1 | EQE ; result -> results/<tag>.json"""
import sys, json, os, time, numpy as np
from scipy.optimize import differential_evolution
import collim as C
dev = C.load_dev(sys.argv[1]); lam0, fwhm, q, npr, obj, tag = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), sys.argv[6], sys.argv[7]
slabs = [C.load_slab(f) for f in C.SLABS]; names = [os.path.basename(f)[:-4] for f in C.SLABS]
bd = [(0, 250) if k % 2 == 0 else (0, 350) for k in range(2 * npr)] + [(0, len(slabs) - 0.001)]
def f(y, step=10):
    D = C.evaluate(dev, 'slabdbr', lam0, fwhm, q, step=step, x=y[:-1], slab=slabs[int(y[-1])]); return C.metrics(D)
t0 = time.time()
r = differential_evolution(lambda y: -f(y)[obj], bd, seed=2, popsize=10, maxiter=120, tol=1e-5, polish=False, init='latinhypercube')
# local polish (continuous thicknesses, slab fixed)
from scipy.optimize import minimize
s = int(r.x[-1]); lo = np.array([b[0] for b in bd[:-1]]); hi = np.array([b[1] for b in bd[:-1]])
m = minimize(lambda x: -f(np.r_[np.clip(x, lo, hi), s + 0.5])[obj], r.x[:-1], method='Nelder-Mead', options=dict(maxiter=1500, xatol=0.3, fatol=1e-6))
xb = np.clip(m.x, lo, hi); fin = f(np.r_[xb, s + 0.5], step=2)
D = C.evaluate(dev, 'slabdbr', lam0, fwhm, q, step=2, x=xb, slab=slabs[s])
out = dict(tag=tag, dev=sys.argv[1], lam0=lam0, fwhm=fwhm, q=q, npairs=npr, objective=obj, slab=names[s],
           dbr=[round(float(v), 1) for v in xb], metrics=fin, D=[float(v) for v in D], minutes=round((time.time() - t0) / 60, 1))
os.makedirs('results', exist_ok=True); json.dump(out, open('results/%s.json' % tag, 'w'), indent=1)
print(tag, obj, names[s], {k: round(v, 4) for k, v in fin.items()}, flush=True)
