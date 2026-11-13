"""Ag_n (n = 1-4) on HATCN and on MoOx (Mo3O9): charge growth and 2D-vs-3D.

    python scripts/170_agn_wetting.py --gfn2            # build + relax Ag (substrate frozen)
    QC_THREADS=4 python scripts/170_agn_wetting.py --pbe0  # PBE0-D3BJ single points (resumable)
    -> runs/agn_wetting.json

Why: the literature survey (notes/lit_survey/REVIEW.md) refuted the sign-only
form of the charge veto (Ag is delta+ on HATCN too) and found no data on how
the Ag charge evolves with cluster size on an organic acceptor vs an oxide
(gap 6). The question the veto really asks is dynamic: does the substrate keep
pulling charge as the cluster grows, and does that favour flat (wetting) or
3D (balling-up) clusters? So for each substrate and n:

  q_tot(n)      total Mulliken charge on Ag_n
  E_adh(n)      E(substrate) + E(Ag_n, same geometry, free) - E(complex)
  dE_2D3D(n)    E(complex, flat cluster) - E(complex, upright/3D cluster)
                negative = substrate prefers the flat (wetting) cluster

Geometries: idealised Ag_n (Ag-Ag 2.80 A) anchored on the PBE0 deepest Ag1
site (structures/*_Ag1_deep.xyz, package F), "flat" = in the local surface
plane, "up" = pointing away from it (n=2: upright dimer, n=3: upright
triangle, n=4: tetrahedron). Ag relaxed at GFN2 with the substrate frozen
(GFN2 overbinds Ag-Ag, so the relaxation is short and only removes clashes);
energies and charges at PBE0-D3BJ/def2-SVP.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "runs", "agn_wetting.json")
CK = os.path.join(ROOT, "runs", "agn_wetting_ckpt")
SUBS = {"HATCN": "HATCN_Ag1_deep.xyz", "Mo3O9": "Mo3O9_Ag1_deep.xyz"}
D = 2.80
EV = 27.211386


def load():
    return json.load(open(OUT)) if os.path.exists(OUT) else {}


def put(path, val):
    d = load()
    cur = d
    for k in path[:-1]:
        cur = cur.setdefault(k, {})
    cur[path[-1]] = val
    json.dump(d, open(OUT, "w"), indent=1)


def frame(X, ag1):
    """Local surface normal n (away from the substrate) and two tangents."""
    near = X[np.argsort(np.linalg.norm(X - ag1, axis=1))[:8]]
    c = near.mean(0)
    u, s, vt = np.linalg.svd(near - c)
    n = vt[-1]
    if np.dot(ag1 - c, n) < 0:
        n = -n
    t1 = np.cross(n, [0, 0, 1.0])
    if np.linalg.norm(t1) < 1e-3:
        t1 = np.cross(n, [1.0, 0, 0])
    t1 /= np.linalg.norm(t1)
    t2 = np.cross(n, t1)
    return n, t1, t2


def clusters(X, ag1):
    n, t1, t2 = frame(X, ag1)
    h = D * np.sqrt(3) / 2
    out = {1: {"flat": [ag1]}}
    out[2] = {"flat": [ag1, ag1 + D * t1], "up": [ag1, ag1 + D * n]}
    out[3] = {"flat": [ag1, ag1 + D * t1, ag1 + D / 2 * t1 + h * t2],
              "up": [ag1, ag1 + D * t1, ag1 + D / 2 * t1 + h * n]}
    # rhombus flat; tetrahedron standing on an edge-sharing base
    out[4] = {"flat": [ag1, ag1 + D * t1, ag1 + D / 2 * t1 + h * t2, ag1 + 1.5 * D * t1 + h * t2],
              "up": [ag1, ag1 + D * t1, ag1 + D / 2 * t1 + h * t2,
                     ag1 + D / 2 * t1 + h / 3 * t2 + D * np.sqrt(2 / 3) * n]}
    return {k: {s: np.array(v) for s, v in d.items()} for k, d in out.items()}


def gfn2_relax(S, X, A):
    from ase import Atoms
    from ase.constraints import FixAtoms
    from ase.optimize import BFGS
    from tblite.ase import TBLite
    nag = len(A)
    a = Atoms(S + ["Ag"] * nag, positions=np.vstack([X, A]))
    a.set_constraint(FixAtoms(indices=list(range(len(S)))))
    # short: remove clashes, keep the shape class (GFN2 overbinds Ag-Ag).
    # Oxide clusters with several Ag sometimes do not converge at 300 K
    # electronic temperature; retry hotter, then keep the start geometry.
    for et in (300.0, 1500.0, 3000.0):
        a.positions = np.vstack([X, A])
        a.calc = TBLite(method="GFN2-xTB", verbosity=0, multiplicity=1 + (nag % 2),
                        max_iterations=400, mixer_damping=0.1, electronic_temperature=et)
        try:
            BFGS(a, logfile=None, maxstep=0.1).run(fmax=0.15, steps=60)
            return a.positions[len(S):].copy(), float(a.get_potential_energy()), a.get_charges()[len(S):]
        except Exception:
            continue
    return np.array(A), float("nan"), np.full(nag, np.nan)


def stage_gfn2():
    from pathgeom import read_xyz
    for name, fn in SUBS.items():
        S, X = read_xyz(os.path.join(ROOT, "structures", fn))
        X = np.asarray(X, float)
        k = S.index("Ag")
        ag1 = X[k]
        S = [s for i, s in enumerate(S) if i != k]
        X = np.delete(X, k, axis=0)
        put([name, "substrate"], {"S": S, "X": X.round(5).tolist()})
        for n, shapes in clusters(X, ag1).items():
            for shape, A in shapes.items():
                if name in load() and f"{n}_{shape}" in load()[name].get("gfn2", {}):
                    continue
                t0 = time.time()
                P, e, q = gfn2_relax(S, X, A)
                put([name, "gfn2", f"{n}_{shape}"], {"Ag": P.round(5).tolist(), "E": e,
                                                     "q_tot": round(float(q.sum()), 4),
                                                     "sec": round(time.time() - t0)})
                print(name, n, shape, "q_tot", round(float(q.sum()), 3), flush=True)


def scf_ck(m, tag):
    from pyscf import dft
    os.makedirs(CK, exist_ok=True)
    ck = os.path.join(CK, tag + ".chk")
    mf = (dft.UKS(m) if m.spin else dft.RKS(m)).density_fit()
    mf.xc = "PBE0"
    mf.grids.level = 4
    mf.conv_tol = 1e-9
    mf.level_shift = 0.2
    if os.path.exists(ck):
        dm0 = mf.from_chk(ck)
    else:
        g = (dft.UKS(m) if m.spin else dft.RKS(m)).density_fit()
        g.xc, g.grids.level, g.conv_tol, g.max_cycle, g.level_shift, g.damp = "PBE", 2, 1e-6, 150, 0.3, 0.3
        g.kernel()
        dm0 = g.make_rdm1()
    mf.chkfile = ck
    mf.max_cycle = 40
    e = mf.kernel(dm0)
    if not mf.converged:
        mf = mf.newton()
        mf.chkfile = ck
        mf.max_cycle = 60
        e = mf.kernel(mf.make_rdm1())
    return mf, e


def stage_pbe0():
    import qclocal as q
    from pyscf import lib
    lib.num_threads(int(os.environ.get("QC_THREADS", "4")))
    jobs = [(n, s) for n in (1, 2, 4, 3) for s in (("flat",) if n == 1 else ("flat", "up"))]
    for name in SUBS:
        sub = load()[name]["substrate"]
        S, X = sub["S"], np.array(sub["X"])
        if "E_sub" not in load()[name]:
            mf, e = scf_ck(q.mol(S, X, 1), f"{name}_sub")
            put([name, "E_sub"], float(e))
    for n, shape in jobs:
        for name in SUBS:
            key = f"{n}_{shape}"
            if key in load()[name].get("pbe0", {}):
                continue
            sub = load()[name]["substrate"]
            S, X = sub["S"], np.array(sub["X"])
            A = np.array(load()[name]["gfn2"][key]["Ag"])
            mult = 1 + (n % 2)
            t0 = time.time()
            mfc, ec = scf_ck(q.mol(S + ["Ag"] * n, np.vstack([X, A]), mult), f"{name}_{key}")
            mff, ef = scf_ck(q.mol(["Ag"] * n, A, mult), f"{name}_{key}_free")
            qq = q.mulliken(mfc)[len(S):]
            d3 = (q.d3bj(S + ["Ag"] * n, np.vstack([X, A])) - q.d3bj(S, X)
                  - (q.d3bj(["Ag"] * n, A) if n > 1 else 0.0)) * EV
            e_adh = (load()[name]["E_sub"] + ef - ec) * EV - d3
            put([name, "pbe0", key], {"E_cx": float(ec), "E_free": float(ef),
                                      "conv": bool(mfc.converged and mff.converged),
                                      "q_Ag": [round(float(v), 4) for v in qq],
                                      "q_tot": round(float(qq.sum()), 4),
                                      "E_adh": round(float(e_adh), 4),
                                      "E_adh_per_Ag": round(float(e_adh) / n, 4),
                                      "sec": round(time.time() - t0)})
            print(name, key, "q_tot", round(float(qq.sum()), 3), "E_adh", round(float(e_adh), 3),
                  f"{time.time()-t0:.0f}s", flush=True)
    # summary
    d = load()
    summ = {}
    for name in SUBS:
        p = d[name].get("pbe0", {})
        row = {}
        for n in (1, 2, 3, 4):
            f_, u_ = p.get(f"{n}_flat"), p.get(f"{n}_up")
            row[n] = {"q_tot_flat": f_ and f_["q_tot"], "q_tot_up": u_ and u_["q_tot"],
                      "E_adh_flat": f_ and f_["E_adh"], "E_adh_up": u_ and u_["E_adh"],
                      "dE_2D3D": round((f_["E_cx"] - u_["E_cx"]) * EV, 4) if f_ and u_ else None}
        summ[name] = row
    put(["summary"], summ)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    if "--gfn2" in sys.argv:
        stage_gfn2()
    if "--pbe0" in sys.argv:
        stage_pbe0()
