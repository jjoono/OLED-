"""Re-place Ag at PBE0 where the GFN2 site is crowded, then score it as in 149.

    python scripts/151_crowded_sites_pbe0.py T2T DPEPO PO-T2T NBphen

Script 149 takes the Ag position from the GFN2 survey (molecule frozen, Ag
relaxed). On open sites that position is within 0.03 eV of PBE0. On sites
squeezed between ortho hydrogens it is not: GFN2 pulls Ag to 1.7-2.2 A from H
(NBphen pocket 1.70, DPEPO P=O 1.69, PO-T2T 1.97, T2T 2.20), a contact PBE0
finds repulsive -- T2T scored 0.07 eV there, below benzene.

Here, for each such site:
  1. push Ag out along (site - nearest heteroatom) until no H is closer than
     2.4 A;
  2. relax the Ag position only, molecule frozen, at PBE0/def2-SVP on a small
     fragment (R_RELAX, default 4 A) plus D3BJ of the whole molecule;
  3. score the relaxed position exactly as script 149 does (FRAG_RCUT fragment,
     whole-molecule D3).
Results go into runs/common_pbe0_fragment.json under "<site>_pbe0relax",
alongside the GFN2-position value, which is kept for the record.
"""
import os, sys, json, time
import importlib.util

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("f149", os.path.join(HERE, "149_fragment_pbe0.py"))
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)
q = f.q

ROOT = os.path.join(HERE, "..")
R_RELAX = float(os.environ.get("R_RELAX", "4.0"))
H_MIN = 2.4
HET = {"N", "O", "P", "Al", "B", "Si"}


def push_out(S, X, ag):
    het = [i for i, s in enumerate(S) if s in HET]
    k = min(het, key=lambda i: np.linalg.norm(X[i] - ag))
    u = ag - X[k]
    u /= np.linalg.norm(u)
    H = np.array([X[i] for i, s in enumerate(S) if s == "H"])
    p = ag.copy()
    for _ in range(40):
        if np.linalg.norm(H - p, axis=1).min() >= H_MIN:
            break
        p = p + 0.05 * u
    return p, S[k] + str(k)


def relax_ag(S, X, ag0):
    """Ag-only minimisation of E_PBE0(fragment+Ag) + D3(molecule+Ag).

    The fragment is cut once, around the start; Ag moves < 1 A, so it stays
    inside it."""
    old = f.R_CUT
    f.R_CUT = R_RELAX
    FS, FX, _ = f.fragment(S, X, ag0)
    f.R_CUT = old
    Zs = np.array([q.Z[s] for s in S + ["Ag"]])
    from dftd3.interface import RationalDampingParam, DispersionModel
    state = {"dm": None, "n": 0}

    def fg(p):
        m = q.mol(FS + ["Ag"], np.vstack([FX, p]), 2)
        mf, e = q.scf(m, conv=1e-8, dm0=state["dm"])
        state["dm"] = mf.make_rdm1()
        g = mf.nuc_grad_method().kernel()[-1] / q.BOHR          # Ha/A
        model = DispersionModel(Zs, np.vstack([X, p]) / q.BOHR)
        r = model.get_dispersion(RationalDampingParam(method="pbe0"), grad=True)
        e += r["energy"]
        g = g + r["gradient"][-1] / q.BOHR
        state["n"] += 1
        print(f"    step {state['n']:2d}  E {e:.6f}  |g| {np.linalg.norm(g) * f.EV:.3f} eV/A",
              flush=True)
        return e, g

    res = minimize(fg, ag0, jac=True, method="BFGS",
                   options={"gtol": 0.03 / f.EV, "maxiter": 25})
    return res.x, state["n"]


def main(names):
    sv = json.load(open(os.path.join(ROOT, "runs", "common_survey_xtb.json")))
    for name in names:
        S, X = read(name)
        site = sv[name]["pbe0_sites"][0]
        key = f"{site}_pbe0relax"
        done = json.load(open(f.OUT))
        if key in done.get(name, {}):
            continue
        ag_g = np.array(sv[name]["sites"][site]["pos"])
        t0 = time.time()
        ag0, het = push_out(S, X, ag_g)
        print(f"{name}: {site} pushed {np.linalg.norm(ag0 - ag_g):.2f} A out from {het}", flush=True)
        ag, nstep = relax_ag(S, X, ag0)
        d = np.linalg.norm(X - ag, axis=1)
        r = f.eb(S, X, ag)
        r.update(seconds=round(time.time() - t0), relax_steps=nstep, R_cut=f.R_CUT,
                 R_relax=R_RELAX, pos=[round(float(v), 4) for v in ag],
                 near=f"{S[int(d.argmin())]}{int(d.argmin())}",
                 d_near=round(float(d.min()), 3),
                 d_H_min=round(float(min(d[i] for i, s in enumerate(S) if s == "H")), 3))
        done = json.load(open(f.OUT))
        row = done.setdefault(name, {})
        row[key] = r
        good = {k: v["E_b"] for k, v in row.items() if isinstance(v, dict) and "E_b" in v}
        row["_deepest"] = max(good, key=good.get)
        row["_E_b_deepest"] = good[row["_deepest"]]
        json.dump(done, open(f.OUT, "w"), indent=1)
        print(f"{name:<9} {key:<26} near {r['near']:<5} E_b {r['E_b']:.3f}  "
              f"d {r['d_near']:.2f}  H_min {r['d_H_min']:.2f}  {r['seconds']}s", flush=True)


def read(name):
    S, X = f.read_xyz(os.path.join(ROOT, "structures", "common", f"{name}.xyz"))
    return S, np.asarray(X, float)


if __name__ == "__main__":
    main(sys.argv[1:])
