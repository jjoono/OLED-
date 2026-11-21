"""
freeform.py -- aperiodic (free-form) top filter for the best architecture found by arch.py.

The periodic DBR of arch.py (first layer + N identical H/L pairs) is relaxed: every filter layer thickness is
free, a thick low-index top spacer is added (its far interface is a weak, phase-tunable reflector), and the
organic thicknesses, IZO thickness and SL index are re-optimised with it.  Objective: on-axis cd/A on the
arch.py objective grid; L-BFGS-B with a parallel central-difference gradient, multistart over the spacer
fringe; final evaluation on the 0.5-nm grid.  Robustness: Monte Carlo over random layer-thickness errors
(sigma = 1 % and 2 % of each thickness, independent layers).

usage: python3 freeform.py <base case> [spacer_lo spacer_hi]      (writes results/freeform_<base>.json)
"""
from __future__ import annotations

import copy
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import arch as A  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
G = {}


def unpack(z):
    """z = [d_etl, d_htl, d_ag, 100*n_m, filter layers (H, L, H, ..., H, spacer L)]"""
    P = copy.deepcopy(G["P0"])
    P.update(d_etl=z[0], d_htl=z[1], d_ag=z[2], n_m=z[3] / 100.0)
    P["dbr"] = [float(t) for t in z[4:]]
    return P


def ce(z, ctx=None, s0=None, **kw):
    c = G["case"]
    return A.run(ctx or G["ctx"], s0 if s0 is not None else G["s0"], unpack(z), c["arch"], c["exit"], c["mode"],
                 n_th=kw.get("n_th", 91), dbr=c.get("dbr"), du=kw.get("du", 1e-3), chunk=kw.get("chunk"))["CE"]


def _f(z):
    return -ce(z)


def grad(z, pool, h=0.5):
    pts = []
    for i in range(len(z)):
        e = np.zeros(len(z)); e[i] = h
        pts += [z + e, z - e]
    v = np.array(pool.map(_f, pts))
    return (v[0::2] - v[1::2]) / (2 * h)


def main(base, lo=300.0, hi=2400.0):
    d = json.load(open(os.path.join(A.OUT, f"{base}.json")))
    c = A.CASES[base]
    ctx, s0, ag_b = A.make(c["emitter"], c["anode"])
    G.update(ctx=ctx, s0=s0, case=c, P0=d["params"])
    D = d["params"]["dbr"]
    filt = [D["dH1"]] + [D["dL"], D["dH"]] * int(D["N"])
    z0 = np.array([d["params"]["d_etl"], d["params"]["d_htl"], d["params"]["d_ag"], 100 * d["params"]["n_m"]] + filt + [0.0])
    bnd = [(10, 500), (20, 600), ag_b, (160, 200)] + [(5, 400)] * len(filt) + [(0, 2600)]
    t0 = time.time()
    with Pool(int(os.environ.get("FF_WORKERS", "4"))) as pool:
        base_ce = -_f(z0)
        # spacer scan (one fringe per ~ lambda/2n; fine enough to land on every maximum)
        ts = np.arange(lo, hi, 12.0)
        zs = [np.r_[z0[:-1], t] for t in ts]
        v = -np.array(pool.map(_f, zs))
        order = np.argsort(v)[::-1]
        starts = []
        for k in order:                                    # two best, separated maxima
            if all(abs(ts[k] - ts[j]) > 120 for j in starts):
                starts.append(k)
            if len(starts) == 2:
                break
        runs = []
        for k in starts:
            z = np.r_[z0[:-1], ts[k]]
            r = minimize(_f, z, jac=lambda zz: grad(zz, pool), method="L-BFGS-B", bounds=bnd,
                         options=dict(maxiter=40, ftol=1e-7, gtol=1e-3))
            runs.append(dict(spacer_start=float(ts[k]), CE_start=float(v[k]), CE=float(-r.fun), nit=int(r.nit), z=r.x.tolist()))
            print(f"{base}: start spacer {ts[k]:.0f} nm ({v[k]:.0f}) -> {-r.fun:.0f} cd/A after {r.nit} iterations", flush=True)
        best = max(runs, key=lambda q: q["CE"])
        zb = np.array(best["z"])
        # robustness: random thickness errors on every layer (organics, IZO, filter, spacer)
        rng = np.random.default_rng(7)
        mc = {}
        for sig in (0.01, 0.02):
            zz = []
            for _ in range(24):
                e = np.ones(len(zb)); e[[0, 1, 2]] += sig * rng.standard_normal(3); e[4:] += sig * rng.standard_normal(len(zb) - 4)
                zz.append(zb * e)
            vals = -np.array(pool.map(_f, zz))
            mc[str(sig)] = dict(mean=float(vals.mean()), p10=float(np.percentile(vals, 10)), min=float(vals.min()))
            zz = []
            for _ in range(24):                            # same for the periodic design (reference)
                e = np.ones(len(z0)); e[[0, 1, 2]] += sig * rng.standard_normal(3); e[4:] += sig * rng.standard_normal(len(z0) - 4)
                zz.append(z0 * e)
            vals = -np.array(pool.map(_f, zz))
            mc[str(sig) + "_periodic"] = dict(mean=float(vals.mean()), p10=float(np.percentile(vals, 10)), min=float(vals.min()))
    ctxf, s0f, _ = A.make(c["emitter"], c["anode"], final=True)
    fin = A.run(ctxf, s0f, unpack(zb), c["arch"], c["exit"], c["mode"], n_th=181, du=2.5e-4, dbr=c.get("dbr"), chunk=40)
    fin0 = d["final"]
    out = dict(base=base, case=c, periodic_objective_CE=base_ce, periodic_final_CE=fin0["CE"], runs=runs,
               best_objective_CE=best["CE"], final=fin, params=unpack(zb), robustness=mc, seconds=time.time() - t0)
    with open(os.path.join(HERE, "results", f"freeform_{base}.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"{base}: periodic {fin0['CE']:.0f} -> free-form {fin['CE']:.0f} cd/A (0.5-nm grid); EQE {fin['eqe']:.3f}; "
          f"A' {100 * fin['Aprime']:.2f}%; robustness {json.dumps({k: round(v['mean']) for k, v in mc.items()})}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], *map(float, sys.argv[2:4]))
