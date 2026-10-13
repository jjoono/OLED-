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
ANCHOR = {"N", "O"}           # atoms Ag can bind; P, B, Si, Al here are saturated


def push_out(S, X, ag):
    """The open position closest to the binding heteroatom(s).

    Two cheaper rules failed: pushing along (GFN2 site - atom) sent the DPEPO
    start 2 A into the molecule (|g| 39 eV/A), and pushing along the lone
    pair runs into neighbouring rings on twisted molecules. So: the
    N and O atoms within 3 A of the GFN2 site are the anchors; candidate points
    lie on shells of 2.3-3.1 A around each anchor (600 directions); a point
    is open if no H is within 2.3 A, no other heavy atom within 2.8 A and no
    anchor within 2.2 A; the
    open point with the smallest summed distance to the anchors is the start
    (for a chelate pocket this is between both N). None open -> the
    heteroatom is buried and not a site."""
    het = [i for i, s in enumerate(S) if s in ANCHOR and np.linalg.norm(X[i] - ag) < 3.0]
    if not het:
        het = [min((i for i, s in enumerate(S) if s in ANCHOR),
                   key=lambda i: np.linalg.norm(X[i] - ag))]
    n = 600
    k = np.arange(n) + 0.5
    th = np.arccos(1 - 2 * k / n)
    ph = np.pi * (1 + 5 ** 0.5) * k
    dirs = np.stack([np.cos(ph) * np.sin(th), np.sin(ph) * np.sin(th), np.cos(th)], 1)
    Hs = np.array([X[i] for i, s in enumerate(S) if s == "H"])
    heavy = np.array([X[i] for i, s in enumerate(S) if s != "H" and i not in het])
    best = None
    for a in het:
        for r in (2.3, 2.5, 2.7, 2.9, 3.1):
            P = X[a] + r * dirs
            dH = np.linalg.norm(P[:, None] - Hs[None], axis=2).min(1)
            dX = np.linalg.norm(P[:, None] - heavy[None], axis=2).min(1)
            dA = np.linalg.norm(P[:, None] - X[het][None], axis=2).min(1)
            ok = (dH >= 2.3) & (dX >= 2.8) & (dA >= 2.2)
            for p in P[ok]:
                s = sum(np.linalg.norm(p - X[b]) for b in het)
                if best is None or s < best[0]:
                    best = (s, p)
    label = ",".join(S[i] + str(i) for i in het)
    return (None if best is None else best[1]), label


# one file per molecule: two runners write at the same time
CKPT = os.path.join(ROOT, "runs", "crowded_relax_ckpt")


def ckpt_load(name):
    p = os.path.join(CKPT, f"{name}.json")
    return json.load(open(p)) if os.path.exists(p) else None


def ckpt_save(name, rec):
    os.makedirs(CKPT, exist_ok=True)
    json.dump(rec, open(os.path.join(CKPT, f"{name}.json"), "w"), indent=1)


def relax_ag(S, X, ag0, name, x0=None, steps0=0, best0=None):
    """Ag-only minimisation of E_PBE0(fragment+Ag) + D3(molecule+Ag).

    The fragment is cut once, around the start ag0; Ag moves < 1 A, so it
    stays inside it. Every step writes the lowest point so far to CKPT: the
    container restarts every several hours and a 5-step relaxation of a
    40-atom fragment takes 6-10 h. A resumed run starts from that point with a
    fresh Hessian, keeping the same fragment and box (both tied to ag0)."""
    old = f.R_CUT
    f.R_CUT = R_RELAX
    FS, FX, _ = f.fragment(S, X, ag0)
    f.R_CUT = old
    Zs = np.array([q.Z[s] for s in S + ["Ag"]])
    from dftd3.interface import RationalDampingParam, DispersionModel
    state = {"dm": None, "n": steps0, "best": best0}

    def fg(p):
        m = q.mol(FS + ["Ag"], np.vstack([FX, p]), 2)
        # grid 2 and 1e-7 are enough to locate the minimum; the score at the
        # end uses the standard settings of 149
        mf, e = q.scf(m, conv=1e-7, grid=2, dm0=state["dm"])
        state["dm"] = mf.make_rdm1()
        g = mf.nuc_grad_method().kernel()[-1] / q.BOHR          # Ha/A
        model = DispersionModel(Zs, np.vstack([X, p]) / q.BOHR)
        r = model.get_dispersion(RationalDampingParam(method="pbe0"), grad=True)
        e = (e + r["energy"]) * f.EV                             # eV
        g = (g + r["gradient"][-1] / q.BOHR) * f.EV              # eV/A
        state["n"] += 1
        print(f"    step {state['n']:2d}  E {e:.4f} eV  |g| {np.linalg.norm(g):.3f} eV/A",
              flush=True)
        if state["best"] is None or e < state["best"]["E"]:
            state["best"] = {"E": e, "pos": [float(v) for v in p],
                             "g": float(np.linalg.norm(g))}
        ckpt_save(name, {"start": [float(v) for v in ag0], "steps": state["n"],
                         "best": state["best"], "relaxed": False})
        return e, g

    # eV and A: the identity starting Hessian then means a first step of
    # g x 1 A^2/eV, about right for an Ag-N stretch (k ~ 1-3 eV/A^2); in
    # hartree units the first step was 0.01 A
    # a box of +-0.8 A around the start caps any wild first step
    x_start = ag0 if x0 is None else np.asarray(x0, float)
    if state["best"] is not None and state["best"]["g"] < 0.05:
        return np.array(state["best"]["pos"]), state["n"]
    minimize(fg, x_start, jac=True, method="L-BFGS-B",
             bounds=[(v - 0.8, v + 0.8) for v in ag0],
             options={"gtol": 0.05, "maxiter": max(3, 12 - steps0)})
    ag = np.array(state["best"]["pos"])
    ckpt_save(name, {"start": [float(v) for v in ag0], "steps": state["n"],
                     "best": state["best"], "relaxed": True})
    return ag, state["n"]


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
        ck = ckpt_load(name)
        if ck:
            ag0 = np.array(ck["start"])
            if ck["relaxed"]:
                ag, nstep = np.array(ck["best"]["pos"]), ck["steps"]
            else:
                print(f"{name}: resuming after {ck['steps']} steps, "
                      f"|g| {ck['best']['g']:.3f} eV/A", flush=True)
                ag, nstep = relax_ag(S, X, ag0, name, x0=ck["best"]["pos"],
                                     steps0=ck["steps"], best0=ck["best"])
            score(name, key, S, X, ag, nstep, t0)
            continue
        ag0, het = push_out(S, X, ag_g)
        if ag0 is None:
            print(f"{name}: {site} anchors {het} buried, no open start", flush=True)
            done = json.load(open(f.OUT))
            done.setdefault(name, {})[key] = {"buried": het}
            json.dump(done, open(f.OUT, "w"), indent=1)
            continue
        dh = {S[i] + str(i): round(float(np.linalg.norm(X[i] - ag0)), 2)
              for i in range(len(S)) if S[i] in HET and np.linalg.norm(X[i] - ag0) < 3.5}
        print(f"{name}: {site} start {np.linalg.norm(ag0 - ag_g):.2f} A from the GFN2 site, "
              f"anchors {dh}", flush=True)
        ag, nstep = relax_ag(S, X, ag0, name)
        score(name, key, S, X, ag, nstep, t0)


def score(name, key, S, X, ag, nstep, t0):
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
