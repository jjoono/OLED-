"""Ag on Bphen(-)Cs(+) with the whole Bphen, at the cost of a fragment.

    QC_THREADS=4 python scripts/155_bphencs_fragment.py   (resumable)
    -> runs/bphencs_fragment.json, checkpoints in runs/bphencs_frag_ckpt/

Replaces the Gaussian package (bphencs_jobs*), which did not run on the
workstation, and the whole-molecule PySCF run, which did not finish one SCF
in 2 h on 4 cores. Same scheme as 149/151, adapted to a radical:

  locate   Ag position relaxed at PBE-D3 on a 4 A fragment (molecule frozen),
           checkpointed per step
  score    E_b = [E(frag) + E(Ag) - min(E_cx singlet, E_cx triplet)]_PBE0
                 - [D3(BphenCs+Ag) - D3(BphenCs)]  (whole molecule)
  Ag2      PBE0 single point of Ag2 at the Cs site (GFN2 geometry): charge
           on two Ag against one, the veto descriptor at PBE0

Fragments always keep the Cs: the extra electron sits on the phenanthroline
and its donor must stay with it (the fragment is then a doublet). The 4,7-
phenyls enter the electronic part only if they lie within the cut radius;
they always enter D3. Sites: the deepest of each kind from the GFN2 survey
(runs/bphen_cs_seed.json). Relaxation spin: singlet at the Cs site (0.93 eV
below the triplet in the phen model), triplet elsewhere (lower at the N side).
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("f149", os.path.join(HERE, "149_fragment_pbe0.py"))
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)
q = f.q
from pyscf import dft  # noqa: E402

ROOT = os.path.join(HERE, "..")
XYZ = os.path.join(ROOT, "structures", "common", "BphenCs_opt.xyz")
SURVEY = os.path.join(ROOT, "runs", "bphen_cs_seed.json")
OUT = os.path.join(ROOT, "runs", "bphencs_fragment.json")
CKPT = os.path.join(ROOT, "runs", "bphencs_frag_ckpt")
EV = f.EV


def load():
    return json.load(open(OUT)) if os.path.exists(OUT) else {}


def save(d):
    json.dump(d, open(OUT, "w"), indent=1)


def frag_with_cs(S, X, centre, rcut):
    old = f.R_CUT
    f.R_CUT = rcut
    FS, FX, ncap = f.fragment(S, X, centre)
    f.R_CUT = old
    if "Cs" not in FS:
        ics = S.index("Cs")
        FS = FS + ["Cs"]
        FX = np.vstack([FX, X[ics]])
    return FS, np.asarray(FX), ncap


def relax(S, X, ag0, lab, mult):
    ck_path = os.path.join(CKPT, f"{lab}.json")
    ck = json.load(open(ck_path)) if os.path.exists(ck_path) else None
    if ck and ck.get("relaxed"):
        return np.array(ck["best"]["pos"]), ck["steps"]
    os.makedirs(CKPT, exist_ok=True)
    FS, FX, _ = frag_with_cs(S, X, ag0, 4.0)
    Zs = np.array([q.Z[s] for s in S + ["Ag"]])
    from dftd3.interface import RationalDampingParam, DispersionModel
    st = {"dm": None, "n": ck["steps"] if ck else 0, "best": ck["best"] if ck else None}

    def fg(p):
        m = q.mol(FS + ["Ag"], np.vstack([FX, p]), mult)
        mf = dft.UKS(m).density_fit()
        mf.xc = "PBE"
        mf.grids.level = 2
        mf.conv_tol = 1e-7
        mf.max_cycle = 200
        mf.level_shift = 0.2
        e = mf.kernel(st["dm"]) if st["dm"] is not None else mf.kernel()
        if not mf.converged:
            mf = mf.newton()
            e = mf.kernel(mf.make_rdm1())
        st["dm"] = mf.make_rdm1()
        g = mf.nuc_grad_method().kernel()[-1] / q.BOHR
        r = DispersionModel(Zs, np.vstack([X, p]) / q.BOHR).get_dispersion(
            RationalDampingParam(method="pbe"), grad=True)
        e = (e + r["energy"]) * EV
        g = (g + r["gradient"][-1] / q.BOHR) * EV
        st["n"] += 1
        print(f"    {lab} step {st['n']:2d}  E {e:.4f} eV  |g| {np.linalg.norm(g):.3f} eV/A",
              flush=True)
        if st["best"] is None or e < st["best"]["E"]:
            st["best"] = {"E": e, "pos": [float(v) for v in p], "g": float(np.linalg.norm(g))}
        json.dump({"start": [float(v) for v in ag0], "steps": st["n"], "best": st["best"],
                   "relaxed": False}, open(ck_path, "w"), indent=1)
        return e, g

    if not (st["best"] and st["best"]["g"] < 0.05):
        x0 = np.array(st["best"]["pos"]) if st["best"] else ag0
        minimize(fg, x0, jac=True, method="L-BFGS-B",
                 bounds=[(v - 0.8, v + 0.8) for v in ag0],
                 options={"gtol": 0.05, "maxiter": max(3, 12 - st["n"])})
    json.dump({"start": [float(v) for v in ag0], "steps": st["n"], "best": st["best"],
               "relaxed": True}, open(ck_path, "w"), indent=1)
    return np.array(st["best"]["pos"]), st["n"]


def score(d, S, X, lab, ag):
    row = d.setdefault("sites", {}).setdefault(lab, {})
    FS, FX, ncap = frag_with_cs(S, X, ag, 5.0)
    if "E_ag" not in d:
        d["E_ag"] = float(q.scf(q.mol(["Ag"], [[0, 0, 0]], 2))[1])
        save(d)
    if "E_frag" not in row:
        mf, e = q.scf(q.mol(FS, FX, 2))
        row.update(E_frag=float(e), frag_atoms=len(FS), caps=ncap, frag_conv=bool(mf.converged))
        save(d)
    for mult in (1, 3):
        if f"E_m{mult}" in row:
            continue
        t0 = time.time()
        mf, e = q.scf(q.mol(FS + ["Ag"], np.vstack([FX, ag]), mult))
        mq = q.mulliken(mf)
        row.update({f"E_m{mult}": float(e), f"conv_m{mult}": bool(mf.converged),
                    f"q_Ag_m{mult}": round(float(mq[-1]), 4),
                    f"q_Cs_m{mult}": round(float(mq[FS.index("Cs")]), 4),
                    f"sec_m{mult}": round(time.time() - t0)})
        save(d)
        print(f"  {lab} mult {mult}: {row[f'sec_m{mult}']}s", flush=True)
    e_cx = min(row["E_m1"], row["E_m3"])
    elec = (row["E_frag"] + d["E_ag"] - e_cx) * EV
    d3 = (q.d3bj(S + ["Ag"], np.vstack([X, ag])) - q.d3bj(S, X)) * EV
    dd = np.linalg.norm(X - ag, axis=1)
    k = int(dd.argmin())
    row.update(E_b=round(elec - d3, 4), elec=round(elec, 4), disp=round(-d3, 4),
               state="singlet" if row["E_m1"] <= row["E_m3"] else "triplet",
               singlet_minus_triplet_eV=round((row["E_m1"] - row["E_m3"]) * EV, 4),
               pos=[round(float(v), 4) for v in ag], near=f"{S[k]}{k}",
               d_near=round(float(dd.min()), 3),
               d_Cs=round(float(np.linalg.norm(X[S.index("Cs")] - ag)), 3),
               d_H_min=round(float(min(dd[i] for i, s in enumerate(S) if s == "H")), 3))
    save(d)
    print(f"E_b {lab:<10} {row['E_b']:.3f} eV ({row['state']}; elec {row['elec']:.3f}, "
          f"disp {row['disp']:.3f}; q_Ag {row['q_Ag_m1' if row['state'] == 'singlet' else 'q_Ag_m3']:+.3f})"
          f"  near {row['near']} d {row['d_near']:.2f}, Ag-Cs {row['d_Cs']:.2f}", flush=True)


def ag2(d, S, X):
    if "ag2" in d and "q_Ag2_total" in d["ag2"]:
        return
    s154 = importlib.import_module("154_bphencs_gaussian")
    sv = json.load(open(SURVEY))
    deep = sorted(sv["sites"], key=lambda k: -sv["sites"][k]["E_b"])[0]
    P = s154.ag2_start(S, X, sv, deep)
    FS, FX, _ = frag_with_cs(S, X, P.mean(axis=0), 5.0)
    t0 = time.time()
    mf, e = q.scf(q.mol(FS + ["Ag", "Ag"], np.vstack([FX, P]), 2))
    mq = q.mulliken(mf)
    d["ag2"] = {"site": deep, "pos": [[round(float(v), 4) for v in p] for p in P],
                "q_Ag": [round(float(v), 4) for v in mq[-2:]],
                "q_Ag2_total": round(float(mq[-2:].sum()), 4),
                "q_Cs": round(float(mq[FS.index("Cs")]), 4), "conv": bool(mf.converged),
                "Ag_Ag_A": round(float(np.linalg.norm(P[0] - P[1])), 3),
                "sec": round(time.time() - t0)}
    save(d)
    print("Ag2", d["ag2"], flush=True)


def main():
    from pyscf import lib
    lib.num_threads(int(os.environ.get("QC_THREADS", "4")))
    S, X = f.read_xyz(XYZ)
    X = np.asarray(X, float)
    sv = json.load(open(SURVEY))
    by_el = {}
    for lab, v in sorted(sv["sites"].items(), key=lambda kv: -kv[1]["E_b"]):
        el = "".join(c for c in v["near"] if c.isalpha())
        by_el.setdefault(el, lab)
    d = load()
    order = [("Cs", 1), ("ag2", None), ("N", 3), ("C", 3)]
    for el, mult in order:
        if el == "ag2":
            ag2(d, S, X)
            continue
        if el not in by_el:
            continue
        lab = by_el[el]
        if "E_b" in d.get("sites", {}).get(lab, {}):
            continue
        ag, n = relax(S, X, np.array(sv["sites"][lab]["pos"]), lab, mult)
        d.setdefault("sites", {}).setdefault(lab, {})["relax_steps"] = n
        save(d)
        score(d, S, X, lab, ag)
    print("done", flush=True)


if __name__ == "__main__":
    main()
