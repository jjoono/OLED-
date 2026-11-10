"""Lens-only optimization at a FIXED OLED stack (Al, ETL 61.3 / HTL 97.5 nm, 550 nm) in the multi-pass model:
hemisphere (1 variable: height) versus axisymmetric freeform (11 variables: 5 radial + 5 height control points, height
scale), same objective, same tracer.  Objectives: substrate-to-air extraction (total) and the power reaching air in
40-60 and 60-80 deg.  Search: (1+8)-ES with step-size adaptation (10 generations), started at the hemispherical optimum (restart-control
protocol), 4,000 rays per bin; winners re-evaluated at 20,000 rays.  Output: fixedcav/results.json"""
import os, json, time, numpy as np
from scipy.interpolate import PchipInterpolator
import hftrace as T

OUT = 'fixedcav'; os.makedirs(OUT, exist_ok=True)
DV = np.load('../selectivity/data/devices.npz')
P0 = DV['Al_P'][75].reshape(90, 2).sum(1); R = DV['Al_R'][75].reshape(90, 2).mean(1)
r = np.linspace(0, 1, 20001)
HX = np.array([0.25881905, 0.5, 0.70710678, 0.8660254, 0.96592583]); HY = HX[::-1].copy()


def profile(x):
    px = np.concatenate(([0.0], np.sort(x[0:5]), [1.0])); py = np.concatenate(([1.0], x[5:10], [0.0]))
    return PchipInterpolator(px, py)(r) * x[10]


def feasible(x):
    xs = np.sort(x[0:5])
    if np.any(x[0:5] < 0) or np.any(x[0:5] > 1) or np.any(x[5:10] < 0) or np.any(x[5:10] > 1.5) or not (0.1 <= x[10] <= 3):
        return False
    q = np.stack([np.concatenate(([0.0], xs, [1.0])), np.concatenate(([1.0], x[5:10], [0.0]))], 1)
    gaps = np.hypot(*np.diff(q, axis=0).T)                               # control-point spacing (lens radii)
    if gaps.min() < 0.05 or gaps.max() > 1.0:
        return False
    h = profile(x)
    if h.min() < -1e-9 or np.abs(np.gradient(h, r)).max() > 60:          # height field: no (near-)vertical walls
        return False
    pts = np.stack([np.concatenate(([0.0], xs, [1.0])), np.concatenate(([1.0], x[5:10], [0.0])) * x[10]], 1)
    v = np.diff(pts, axis=0); ang = np.degrees(np.abs(np.arctan2(v[1:, 0] * v[:-1, 1] - v[1:, 1] * v[:-1, 0], (v[1:] * v[:-1]).sum(1))))
    return bool(ang.max() <= 120)


def evaluate(x, N):
    BT, BR, lost = T.bsdf(r, profile(x), n=1.51, N=N)
    sol = np.linalg.solve(np.eye(90) - R[:, None] * BR.T, P0); air = sol @ BT
    return dict(total=float(air.sum()), b4060=float(air[40:60].sum()), b6080=float(air[60:80].sum()), BT=BT, BR=BR)


CACHE = {}


def f(x, obj, N=4000):
    k = (tuple(np.round(x, 6)), N)
    if k not in CACHE:
        CACHE[k] = evaluate(x, N)
    return CACHE[k][obj]


def hemi_opt(obj):
    best = None
    for s in np.linspace(0.4, 2.6, 12):
        x = np.concatenate((HX, HY, [s])); v = f(x, obj)
        if best is None or v > best[0]:
            best = (v, s)
    lo, hi = best[1] - 0.2, best[1] + 0.2
    for s in np.linspace(lo, hi, 7):
        x = np.concatenate((HX, HY, [s])); v = f(x, obj)
        if v > best[0]:
            best = (v, s)
    return np.concatenate((HX, HY, [best[1]])), best[0]


def es(obj, x0, gens=10, lam=8, sig=0.08, seed=0):
    rng = np.random.default_rng(seed); x = x0.copy(); fx = f(x, obj); hist = [fx]; ne = 0
    scale = np.r_[np.ones(10), 2.0]
    for g in range(gens):
        succ = 0; bestc = None
        for _ in range(lam):
            for _try in range(50):
                c = x + sig * scale * rng.normal(size=11)
                if feasible(c):
                    break
            else:
                continue
            fc = f(c, obj); ne += 1
            if bestc is None or fc > bestc[1]:
                bestc = (c, fc)
        if bestc is not None and bestc[1] > fx:
            x, fx = bestc; succ = 1
        sig *= 1.5 if succ else 0.75
        hist.append(fx)
        print(f'  {obj} gen {g} best {fx:.5f} sig {sig:.3f}', flush=True)
    return x, fx, ne, hist


if __name__ == '__main__':
    res = json.load(open(f'{OUT}/results.json')) if os.path.exists(f'{OUT}/results.json') else {}
    from scipy.io import loadmat
    MB = np.load('bsdfs.npz')
    lt = {'total': MB['free_total_x'], 'b4060': MB['free_40-60_x'], 'b6080': MB['free_60-80_x']}
    for obj in ('total', 'b4060'):
        if obj in res:
            continue
        t0 = time.time()
        xh, fh = hemi_opt(obj)
        xl = np.concatenate((lt[obj][0:10], [lt[obj][12]]))
        xf, ff, ne, hist = es(obj, xh, seed=1)
        fin = {k: evaluate(x, 20000) for k, x in (('hemi', xh), ('free_es', xf), ('free_lighttools', xl))}
        res[obj] = dict(hemi_x=xh.tolist(), free_x=xf.tolist(), lt_x=xl.tolist(), evals=ne, hist=hist,
                        **{f'{k}_{o}': v[o] for k, v in fin.items() for o in ('total', 'b4060', 'b6080')}, sec=round(time.time() - t0))
        np.savez(f'{OUT}/{obj}_bsdfs.npz', **{f'{k}_{m}': v[m] for k, v in fin.items() for m in ('BT', 'BR')})
        json.dump(res, open(f'{OUT}/results.json', 'w'), indent=1)
        print(obj, 'hemi %.5f  freeform(ES, fixed cavity) %.5f  gain %.4f  | LightTools freeform at this cavity %.5f  [%ds]' %
              (fin['hemi'][obj], fin['free_es'][obj], fin['free_es'][obj] / fin['hemi'][obj], fin['free_lighttools'][obj], time.time() - t0), flush=True)
