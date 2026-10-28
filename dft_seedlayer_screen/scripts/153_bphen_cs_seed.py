"""Ag on n-doped Bphen: the Bphen(-)-Cs(+) complex as the doped-film model.

    python scripts/153_bphen_cs_seed.py            # all stages, resumable
    -> runs/bphen_cs_seed.json, structures/common/BphenCs_opt.xyz

Bphen:Cs2CO3 films are n-doped by Cs: Cs gives its 6s electron to Bphen and
sits as Cs+ in the phenanthroline pocket (the [Cs(Bphen)] complex). At 20 wt%
Cs2CO3 that is about one Cs per two Bphen, so the doped surface mixes this
complex with neat Bphen (chelate 0.62 eV, already known).

Stages
  1. GFN2: optimise BphenCs (doublet); Cs charge and Cs-N distances.
  2. GFN2: Ag sites, molecule frozen -- the script-148 site set (lone pairs,
     hollows, heteroatom atop/bridge, pockets) plus starts around Cs. The
     complex is a doublet, so Ag1 is an even-electron system: singlet.
  3. GFN2: Ag2 at the deepest site (three starts) -> dq(Ag1->Ag2), the
     same descriptor as for every other material.
  4. PBE0-D3BJ/def2-SVP, whole molecule (the fragment scheme of 149 would cut
     the radical's pi system away from its Cs donor): E_b at the deepest
     distinct sites. Ag* + Bphen(-)Cs(+) are two radicals, so the complex is
     run as a closed-shell singlet and as a triplet and the lower is used --
     a weakly coupled pair is not a closed-shell singlet.
"""
import importlib
import json
import os
import sys
import time

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from tblite.ase import TBLite

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pathgeom import read_xyz, _bonded  # noqa: E402

s148 = importlib.import_module("148_common_materials_survey")
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "runs", "bphen_cs_seed.json")
XYZ = os.path.join(ROOT, "structures", "common", "BphenCs_opt.xyz")
EV = 27.211386


def load():
    return json.load(open(OUT)) if os.path.exists(OUT) else {}


def save(d):
    json.dump(d, open(OUT, "w"), indent=1)


def calc(mult, et=300.0):
    return TBLite(method="GFN2-xTB", verbosity=0, multiplicity=mult, max_iterations=400,
                  mixer_damping=0.2, electronic_temperature=et)


def single(S, X, mult):
    a = Atoms(S, positions=X)
    a.calc = calc(mult)
    return a.get_potential_energy(), a.get_charges()


def relax(S, X, n_fixed, mult, fmax=0.05, steps=200):
    for et in (300.0, 1500.0):
        a = Atoms(S, positions=X)
        if n_fixed:
            a.set_constraint(FixAtoms(indices=list(range(n_fixed))))
        a.calc = calc(mult, et)
        try:
            BFGS(a, logfile=None, maxstep=0.2).run(fmax=fmax, steps=steps)
            return a.get_potential_energy(), a.get_charges(), a.positions.copy()
        except Exception:
            continue
    return None


def write_xyz(p, S, X, comment=""):
    with open(p, "w") as f:
        f.write(f"{len(S)}\n{comment}\n")
        for s, x in zip(S, X):
            f.write(f"{s} {x[0]:.6f} {x[1]:.6f} {x[2]:.6f}\n")


def stage1(d):
    if "complex" in d:
        return read_xyz(XYZ)
    S, X = read_xyz(os.path.join(ROOT, "structures", "BphenCs.xyz"))
    r = relax(S, np.asarray(X, float), 0, 2, fmax=0.03, steps=600)
    e, q, P = r
    write_xyz(XYZ, S, P, "BphenCs GFN2 doublet")
    ics = S.index("Cs")
    nN = [i for i, s in enumerate(S) if s == "N"]
    d["complex"] = {"E_gfn2_eV": round(float(e), 4), "q_Cs": round(float(q[ics]), 3),
                    "Cs_N_A": [round(float(np.linalg.norm(P[ics] - P[i])), 3) for i in nN]}
    save(d)
    print("complex", d["complex"], flush=True)
    return S, P


def cs_starts(S, X, n=60, r=3.2):
    """Directions around Cs with room for Ag (H >= 2.3 A, heavy >= 2.8 A)."""
    ics = S.index("Cs")
    k = np.arange(n) + 0.5
    th = np.arccos(1 - 2 * k / n)
    ph = np.pi * (1 + 5 ** 0.5) * k
    dirs = np.stack([np.cos(ph) * np.sin(th), np.sin(ph) * np.sin(th), np.cos(th)], 1)
    out = []
    for j, u in enumerate(dirs):
        p = X[ics] + r * u
        dd = np.linalg.norm(X - p, axis=1)
        okH = all(dd[i] >= 2.3 for i, s in enumerate(S) if s == "H")
        okX = all(dd[i] >= 2.8 for i, s in enumerate(S) if s not in ("H", "Cs"))
        if okH and okX and all(np.linalg.norm(p - q) > 2.0 for _, q in out):
            out.append((f"Cs{j}", p))
    return out


def stage2(d, S, X):
    if "sites" in d and d.get("sites_done"):
        return
    n = len(S)
    e_mol, _ = single(S, X, 2)
    e_ag, _ = single(["Ag"], [[0, 0, 0]], 2)
    res = d.get("sites", {})
    todo = s148.sites(S, X) + cs_starts(S, X)
    for lab, pos in todo:
        if lab in res:
            continue
        r = relax(S + ["Ag"], np.vstack([X, pos]), n, 1)
        if r is None:
            continue
        e, q, P = r
        dd = np.linalg.norm(X - P[n], axis=1)
        k = int(dd.argmin())
        hmin = min(dd[i] for i, s in enumerate(S) if s == "H")
        res[lab] = {"E_b": round(float(e_mol + e_ag - e), 4), "q": round(float(q[n]), 4),
                    "q_Cs": round(float(q[S.index("Cs")]), 4), "near": f"{S[k]}{k}",
                    "d": round(float(dd.min()), 3), "H_min": round(float(hmin), 3),
                    "pos": [round(float(v), 4) for v in P[n]]}
        d["sites"] = res
        save(d)
    d["sites_done"] = True
    save(d)
    top = sorted(res.items(), key=lambda kv: -kv[1]["E_b"])[:8]
    for lab, v in top:
        print(f"site {lab:<16} E_b(GFN2) {v['E_b']:.3f}  q {v['q']:+.3f}  near {v['near']}"
              f"  H_min {v['H_min']:.2f}", flush=True)


def picks(d, nmax=3):
    ranked = sorted(d["sites"].items(), key=lambda kv: -kv[1]["E_b"])
    out = []
    for lab, v in ranked:
        if any(np.linalg.norm(np.array(v["pos"]) - np.array(d["sites"][p]["pos"])) < 2.0
               for p in out):
            continue
        out.append(lab)
        if len(out) == nmax:
            break
    return out


def stage3(d, S, X):
    if "ag2" in d:
        return
    n = len(S)
    deep = picks(d, 1)[0]
    ag1 = np.array(d["sites"][deep]["pos"])
    q1 = d["sites"][deep]["q"]
    dd = np.linalg.norm(X - ag1, axis=1)
    nv = ag1 - X[int(dd.argmin())]
    nv /= np.linalg.norm(nv)
    t1 = np.cross(nv, [0, 0, 1.0])
    if np.linalg.norm(t1) < 1e-6:
        t1 = np.cross(nv, [1.0, 0, 0])
    t1 /= np.linalg.norm(t1)
    t2 = np.cross(nv, t1)
    best = None
    for th in (0, 2 * np.pi / 3, 4 * np.pi / 3):
        p2 = ag1 + 2.6 * (np.cos(th) * t1 + np.sin(th) * t2)
        r = relax(S + ["Ag", "Ag"], np.vstack([X, ag1, p2]), n, 2)
        if r and (best is None or r[0] < best[0]):
            best = r
    q2 = best[1][n:]
    d["ag2"] = {"site": deep, "q_Ag1": q1, "q_Ag2_total": round(float(q2.sum()), 4),
                "dq_Ag1_Ag2": round(float(q2.sum() - q1), 4),
                "q_Cs_Ag2": round(float(best[1][S.index("Cs")]), 4),
                "Ag_Ag_A": round(float(np.linalg.norm(best[2][n] - best[2][n + 1])), 3)}
    save(d)
    print("ag2", d["ag2"], flush=True)


def stage4(d, S, X):
    import qclocal as q
    from pyscf import lib
    lib.num_threads(int(os.environ.get("QC_THREADS", "4")))
    pb = d.setdefault("pbe0", {})
    if "E_mol" not in pb:
        t0 = time.time()
        mf, e = q.scf(q.mol(S, X, 2))
        pb["E_mol"] = float(e)
        pb["mol_converged"] = bool(mf.converged)
        mfa, ea = q.scf(q.mol(["Ag"], [[0, 0, 0]], 2))
        pb["E_ag"] = float(ea)
        pb["D3_mol"] = float(q.d3bj(S, X))
        save(d)
        print(f"pbe0 reference energies done {time.time() - t0:.0f}s", flush=True)
    # the deepest site of each kind (nearest element): GFN2 ranks within an
    # element, not across (it overbinds N 4-5x, carbon 1.2x)
    by_el = {}
    for lab, v in sorted(d["sites"].items(), key=lambda kv: -kv[1]["E_b"]):
        el = "".join(c for c in v["near"] if c.isalpha())
        by_el.setdefault(el, lab)
    for lab in by_el.values():
        if lab in pb and "E_b" in pb[lab]:
            continue
        ag = np.array(d["sites"][lab]["pos"])
        S2, X2 = S + ["Ag"], np.vstack([X, ag])
        row = pb.get(lab, {})
        for mult in (1, 3):
            key = f"E_m{mult}"
            if key in row:
                continue
            t0 = time.time()
            mf, e = q.scf(q.mol(S2, X2, mult))
            row[key] = float(e)
            row[f"conv_m{mult}"] = bool(mf.converged)
            row[f"q_Ag_m{mult}"] = round(float(q.mulliken(mf)[-1]), 4)
            row[f"q_Cs_m{mult}"] = round(float(q.mulliken(mf)[S.index("Cs")]), 4)
            row[f"sec_m{mult}"] = round(time.time() - t0)
            pb[lab] = row
            save(d)
            print(f"pbe0 {lab} mult {mult}: {row[f'sec_m{mult}']}s", flush=True)
        d3 = (q.d3bj(S2, X2) - pb["D3_mol"]) * EV
        e_cx = min(row["E_m1"], row["E_m3"])
        elec = (pb["E_mol"] + pb["E_ag"] - e_cx) * EV
        row.update(E_b=round(elec - d3, 4), elec=round(elec, 4), D3=round(d3, 4),
                   state="singlet" if row["E_m1"] <= row["E_m3"] else "triplet",
                   singlet_minus_triplet_eV=round((row["E_m1"] - row["E_m3"]) * EV, 4),
                   near=d["sites"][lab]["near"], H_min=d["sites"][lab]["H_min"])
        pb[lab] = row
        save(d)
        print(f"E_b {lab:<16} {row['E_b']:.3f} eV ({row['state']}, elec {row['elec']:.3f}, "
              f"disp {-row['D3']:.3f})  near {row['near']}", flush=True)


if __name__ == "__main__":
    d = load()
    S, X = stage1(d)
    X = np.asarray(X, float)
    stage2(d, S, X)
    stage3(d, S, X)
    if "--no-pbe0" not in sys.argv:
        stage4(d, S, X)
