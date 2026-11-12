"""
refine.py -- local re-optimisation (Powell, bounded) of a design variant, seeded from a
previously optimised design (opt_<base>.json).  Used for material / A' variants, other
emitters (thicknesses scaled by the wavelength ratio), DBR pair count etc.

usage: python3 refine.py '<json spec>'
spec keys: tag, base (tag of seed result), job (ref|refdbr|prop|prophyb), emitter, preset, N, Nb,
           scale (thickness scale for seeding), fixed (dict), maxfev, x0 (optional explicit seed)
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import optimize as opt  # noqa: E402
from devices import run_proposed, run_reference  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def bounds_for(job, ag_b, d_ag_min=10.0, n_m_bounds=(1.6, 2.0)):
    if job == "ref":
        return [(10, 500), (20, 500), (d_ag_min, 60), (0, 160)]
    if job == "refdbr":
        return [(10, 500), (20, 500), (d_ag_min, 60), (0, 220), (20, 220), (20, 260)]
    if job == "prop":
        return [(10, 500), (20, 600), ag_b or (d_ag_min, 25), (0, 220), (20, 220), (20, 300), n_m_bounds]
    if job == "prophyb":
        return [(10, 500), (20, 600), ag_b or (d_ag_min, 25), (0, 220), (20, 220), (20, 300), n_m_bounds,
                (20, 150), (20, 160), (40, 260)]
    raise ValueError(job)


def params_for(job, x, N, Nb):
    if job == "ref":
        return opt._ref_params(x)
    if job == "refdbr":
        return opt._refdbr_params(x, N)
    if job == "prop":
        return opt._prop_params(x, N)
    return opt._prophyb_params(x, N, Nb)


def main(spec):
    t0 = time.time()
    job = spec["job"]
    N = spec.get("N", 3)
    Nb = spec.get("Nb", 2)
    pm, pf, ag_b = opt.preset_mats(spec.get("preset"))
    fixed = {**pf, **spec.get("fixed", {})}
    ctx, s0 = opt.make_ctx(spec["emitter"], mats=pm, step=spec.get("step", 4.0))
    a_sl = spec.get("a_sl", 0.0)
    if "x0" in spec:
        x0 = np.array(spec["x0"], float)
    else:
        with open(os.path.join(HERE, "results", f"opt_{spec['base']}.json")) as fh:
            base = json.load(fh)
        x0 = np.array(base["x"], float)
        sc = spec.get("scale", 1.0)
        if job in ("prop", "prophyb"):
            idx_len = [0, 1, 3, 4, 5] + ([7, 8, 9] if job == "prophyb" else [])
        elif job == "refdbr":
            idx_len = [0, 1, 3, 4, 5]
        else:
            idx_len = [0, 1, 3]
        for i in idx_len:
            if i < len(x0):
                x0[i] *= sc
        if spec.get("x0_override"):
            for k, v in spec["x0_override"].items():
                x0[int(k)] = v
    b = bounds_for(job, ag_b, d_ag_min=spec.get("d_ag_min", 10.0))
    x0 = np.clip(x0, [lo for lo, hi in b], [hi for lo, hi in b])

    def f(x):
        P = {**fixed, **params_for(job, x, N, Nb)}
        du = spec.get("du", 1e-3)
        if job.startswith("ref"):
            return -run_reference(ctx, P, s0, du=du)["summary"]["CE"]
        return -run_proposed(ctx, P, s0, n_th=91, du=du, a_sl_pass=a_sl)["summary"]["CE"]

    f0 = f(x0)
    # Nelder-Mead (bounded) from an initial simplex of +/-8 % steps; Powell with bounds returned
    # points worse than the seed for these aliasing-prone objectives, so it is not used.
    lo = np.array([l for l, h in b]); hi = np.array([h for l, h in b])
    sim = [x0]
    for i in range(len(x0)):
        xi = x0.copy()
        step = 0.08 * max(abs(x0[i]), 0.05 * (hi[i] - lo[i]))
        xi[i] = x0[i] + step if x0[i] + step <= hi[i] else x0[i] - step
        sim.append(xi)
    res = minimize(f, x0, method="Nelder-Mead", bounds=b,
                   options=dict(maxfev=spec.get("maxfev", 400), xatol=0.3, fatol=0.05, initial_simplex=np.array(sim), adaptive=True))
    if res.fun > f0:  # never return worse than the seed
        res.x, res.fun = x0, f0
    ctxf, s0f = opt.make_ctx(spec["emitter"], mats=pm, full=True)
    P = {**fixed, **params_for(job, res.x, N, Nb)}
    if job.startswith("ref"):
        fin = run_reference(ctxf, P, s0f, du=2.5e-4)["summary"]
    else:
        fin = run_proposed(ctxf, P, s0f, du=2.5e-4, a_sl_pass=a_sl)["summary"]
    out = dict(job=job, emitter=spec["emitter"], N=N, Nb=Nb, x=list(map(float, res.x)), fun=float(res.fun),
               f_seed=float(f0), params=opt._jsonable(P), final=fin, fixed=opt._jsonable(fixed),
               preset=spec.get("preset"), a_sl=a_sl, nfev=int(res.nfev), seconds=time.time() - t0,
               method="powell-refine", base=spec.get("base"))
    with open(os.path.join(HERE, "results", f"opt_{spec['tag']}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(spec["tag"], "seed", -f0, "->", -res.fun, "final CE", fin["CE"], "EQE", fin["eqe"], f"{time.time()-t0:.0f}s")


if __name__ == "__main__":
    main(json.loads(sys.argv[1]))
