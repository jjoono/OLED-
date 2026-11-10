"""Does HATCN's charge hand-back survive outside a single molecule?

    python scripts/145_hatcn_raft_xtb.py

The periodic slab (script 137) answered the oxide's version of the cluster
objection: on extended MoO3(010) silver stays ionised at every size, more so
than on Mo3O9. The symmetric objection to the organic is that one molecule in
vacuum is not a film -- neighbouring molecules add intermolecular sites, and a
neighbour's acceptor orbital could take the charge a lone molecule gives back.

Three HATCN discs, flat and coplanar, at the centre-to-centre distance where the
closest intermolecular heavy-atom contact is 3.2 A (van der Waals contact for
N...N), in a triangle -- the smallest raft with an interior intermolecular
pocket. Ag at the aza pocket of one molecule facing the interior, at the aza
pocket facing outward, and in the three-molecule hollow; then Ag2 at the
deepest. GFN2-xTB, molecules frozen, silver relaxed, the level validated on
the package-G clusters (it reproduces the PBE0 ordering and HATCN's hand-back,
+0.31 -> +0.06 per atom).
"""
import os, sys, json

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from tblite.ase import TBLite

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathgeom import STRUCT, read_xyz

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "runs", "hatcn_raft_xtb.json")
GAP = 3.2


def raft():
    S, X = read_xyz(os.path.join(STRUCT, "HATCN.xyz"))
    X = np.asarray(X, float)
    X -= X.mean(axis=0)
    n = np.linalg.svd(X)[2][-1]
    # rotate the molecular plane onto xy
    z = np.array([0, 0, 1.0])
    v = np.cross(n, z)
    if np.linalg.norm(v) > 1e-8:
        ang = np.arccos(np.clip(n @ z, -1, 1))
        v /= np.linalg.norm(v)
        K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K
        X = X @ R.T
    X[:, 2] = 0.0

    def place(D):
        mols = []
        for k in range(3):
            phi = 2 * np.pi * k / 3 + np.pi / 2
            c = D / np.sqrt(3) * np.array([np.cos(phi), np.sin(phi), 0])
            mols.append(X + c)
        return mols

    lo, hi = 8.0, 20.0
    for _ in range(50):
        D = (lo + hi) / 2
        m = place(D)
        dmin = min(np.min(np.linalg.norm(m[a][:, None] - m[b][None], axis=2))
                   for a in range(3) for b in range(a + 1, 3))
        lo, hi = (D, hi) if dmin < GAP else (lo, D)
    m = place(hi)
    return S * 3, np.vstack(m), len(S), hi


def relax(S, X, n_sub, nag):
    a = Atoms(S, positions=X)
    a.set_constraint(FixAtoms(indices=list(range(n_sub))))
    for et in (300.0, 1500.0):
        a.calc = TBLite(method="GFN2-xTB", verbosity=0, multiplicity=1 + nag % 2,
                        max_iterations=400, mixer_damping=0.2, electronic_temperature=et)
        try:
            BFGS(a, logfile=None, maxstep=0.2).run(fmax=0.05, steps=250)
            return a.get_potential_energy(), a.get_charges(), a.positions.copy()
        except Exception:
            a.positions = np.asarray(X, float)
    return None


def aza_pockets(S, X, nper):
    """Midpoints of ring-N pairs 2.6-2.9 A apart, per molecule."""
    out = []
    for m in range(3):
        idx = [m * nper + i for i in range(nper) if S[i] == "N"]
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                i, j = idx[a], idx[b]
                if 2.6 < np.linalg.norm(X[i] - X[j]) < 2.9:
                    out.append((m, (X[i] + X[j]) / 2))
    return out


def main():
    S, X, nper, D = raft()
    n_sub = len(S)
    cen = X.mean(axis=0)
    pk = aza_pockets(S, X, nper)
    # molecule 0's pockets: the one nearest the raft centre faces inward, the
    # furthest faces out
    p0 = sorted([p for m, p in pk if m == 0], key=lambda p: np.linalg.norm(p - cen))
    starts = {"aza_inward": p0[0] + [0, 0, 2.2], "aza_outward": p0[-1] + [0, 0, 2.2],
              "raft_hollow": cen + [0, 0, 2.6]}
    res = {"center_distance_A": round(D, 3), "sites": {}}
    best = None
    for lab, p in starts.items():
        r = relax(S + ["Ag"], np.vstack([X, p]), n_sub, 1)
        if r is None:
            continue
        e, q, P = r
        res["sites"][lab] = {"E": round(float(e), 4), "q_Ag": round(float(q[-1]), 4)}
        print(f"  Ag1 {lab:<12} q {q[-1]:+.3f}  E {e:.3f}", flush=True)
        if best is None or e < best[0]:
            best = (e, lab, P[-1].copy(), float(q[-1]))
    e1, lab, ag1, q1 = best
    best2 = None
    for th in (0, 2 * np.pi / 3, 4 * np.pi / 3):
        p2 = ag1 + 2.6 * np.array([np.cos(th), np.sin(th), 0.0])
        r = relax(S + ["Ag", "Ag"], np.vstack([X, ag1, p2]), n_sub, 2)
        if r and (best2 is None or r[0] < best2[0]):
            best2 = r
    res["deepest"] = lab
    res["q_Ag1"] = round(q1, 4)
    if best2:
        res["q_Ag2_per_atom"] = round(float(best2[1][n_sub:].mean()), 4)
        res["Ag2_Ag_Ag_A"] = round(float(np.linalg.norm(best2[2][-1] - best2[2][-2])), 3)
    res["single_molecule_reference"] = {"q_Ag1": 0.308, "q_Ag2_per_atom": 0.055,
                                        "source": "GFN2 on structures/clusters/HATCN_Ag1,2"}
    json.dump(res, open(OUT, "w"), indent=1)
    print(f"raft: deepest {lab}, q(Ag1) {q1:+.3f}, q(Ag2)/atom "
          f"{res.get('q_Ag2_per_atom', float('nan')):+.3f}   "
          f"(single molecule: +0.308, +0.055)")


if __name__ == "__main__":
    main()
