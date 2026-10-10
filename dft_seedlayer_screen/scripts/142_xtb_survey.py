"""Close the screening table's biggest gap: a multi-site survey of every
candidate, at GFN2-xTB, validated against the PBE0 surveys that exist.

    python scripts/142_xtb_survey.py              # all candidates
    python scripts/142_xtb_survey.py HATCN Bphen  # a subset

Nineteen of the 25 rows in the screening table are single-site lower bounds, and
HATCN showed how wrong one site can be (0.604 -> 1.631 eV). Re-surveying them at
PBE0 is weeks of workstation time; at GFN2-xTB it is minutes per candidate. The
price is the level, so the survey is only allowed to do two things:

  1. FIND the deepest site -- which a cheap method can do if it ranks sites the
     way PBE0 does, checked here on the five candidates PBE0 surveyed
  2. MEASURE the charge on Ag1 and the per-atom charge on Ag2 at that site --
     the observable package G kept, whose ordering GFN2 reproduces on all
     sixteen package-G geometries (script 142's sibling check in the notes)

GFN2's absolute E_b is not used for ranking against the PBE0 rows.

Sites: atop / bridge / hollow from script 131, no cap on how many, PLUS
heteroatom pockets -- pairs of N, O, S or F 2.4-3.6 A apart and not bonded,
which is the chelate script 131 missed on Bphen. Ag relaxed with the molecule
frozen, as everywhere else in the project.
"""
import os, sys, json, importlib

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from tblite.ase import TBLite

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathgeom import CANDIDATES, STRUCT, read_xyz, relaxed_file, contact, _bonded
s131 = importlib.import_module("131_site_spread")
s131.MAX_SITES = 1000

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "runs", "xtb_survey.json")
SKIP = {"TPBi"}
HETERO = {"N", "O", "S", "F"}
D_AG = 2.60


def calc(nag):
    return TBLite(method="GFN2-xTB", verbosity=0, multiplicity=1 + (nag % 2),
                  max_iterations=400, mixer_damping=0.2)


def energy(S, X, nag):
    a = Atoms(S, positions=X)
    a.calc = calc(nag)
    return a.get_potential_energy(), a


def relax(S, X, n_sub, nag):
    a = Atoms(S, positions=X)
    a.set_constraint(FixAtoms(indices=list(range(n_sub))))
    for et in (300.0, 1500.0):
        a.calc = TBLite(method="GFN2-xTB", verbosity=0, multiplicity=1 + (nag % 2),
                        max_iterations=400, mixer_damping=0.2,
                        electronic_temperature=et)
        try:
            BFGS(a, logfile=None, maxstep=0.2).run(fmax=0.05, steps=200)
            return a.get_potential_energy(), a.get_charges(), a.positions.copy()
        except Exception:
            a.positions = np.asarray(X, float)
    return None


def pockets(S, X):
    B = _bonded(S, X)
    cen = X.mean(axis=0)
    nrm = np.linalg.svd(X - cen)[2][-1]
    out = []
    het = [i for i, s in enumerate(S) if s in HETERO]
    for a in range(len(het)):
        for b in range(a + 1, len(het)):
            i, j = het[a], het[b]
            d = float(np.linalg.norm(X[i] - X[j]))
            if B[i][j] or not 2.4 < d < 3.6:
                continue
            mid = (X[i] + X[j]) / 2
            v = mid - cen
            # both the in-plane outward direction and the face normal
            for u in (v / (np.linalg.norm(v) + 1e-9), nrm, -nrm):
                if np.linalg.norm(u) < 0.5:
                    continue
                p = mid + 2.0 * u
                lim = np.array([contact(x) for x in S])
                for _ in range(40):
                    if float((np.linalg.norm(X - p, axis=1) - lim).min()) >= 0:
                        break
                    p = p + 0.1 * u
                out.append((f"pocket{S[i]}{i}{S[j]}{j}", p))
    return out


def survey(tag, fn):
    S, X = read_xyz(os.path.join(STRUCT, relaxed_file(fn)))
    k = S.index("Ag")
    sub_s = [s for i, s in enumerate(S) if i != k]
    sub_x = np.delete(X, k, axis=0)
    n = len(sub_s)
    e_mol, _ = energy(sub_s, sub_x, 0)
    e_ag, _ = energy(["Ag"], [[0, 0, 0]], 1)
    sites = [(f"{kind}{i}", pos) for kind, i, pos, _ in s131.unique_sites(sub_s, sub_x)]
    sites += pockets(sub_s, sub_x)
    sites.append(("stage1", X[k]))
    res = {}
    best = None
    for lab, pos in sites:
        r = relax(sub_s + ["Ag"], np.vstack([sub_x, pos]), n, 1)
        if r is None:
            continue
        e, q, P = r
        eb = e_mol + e_ag - e
        d = np.linalg.norm(sub_x - P[n], axis=1)
        res[lab] = {"E_b": round(float(eb), 4), "q": round(float(q[n]), 4),
                    "near": f"{sub_s[int(d.argmin())]}{int(d.argmin())}",
                    "d": round(float(d.min()), 3)}
        if best is None or eb > best[0]:
            best = (eb, lab, P[n].copy(), float(q[n]))
    # Ag2 at the deepest site: three orientations in the local tangent plane
    eb, lab, ag1, q1 = best
    d = np.linalg.norm(sub_x - ag1, axis=1)
    nv = ag1 - sub_x[int(d.argmin())]
    nv /= np.linalg.norm(nv)
    t1 = np.cross(nv, [0, 0, 1.0])
    if np.linalg.norm(t1) < 1e-6:
        t1 = np.cross(nv, [1.0, 0, 0])
    t1 /= np.linalg.norm(t1)
    t2 = np.cross(nv, t1)
    best2 = None
    for th in (0, 2 * np.pi / 3, 4 * np.pi / 3):
        p2 = ag1 + D_AG * (np.cos(th) * t1 + np.sin(th) * t2)
        r = relax(sub_s + ["Ag", "Ag"], np.vstack([sub_x, ag1, p2]), n, 2)
        if r and (best2 is None or r[0] < best2[0]):
            best2 = r
    out = {"sites": res, "deepest_site": lab, "E_b_deepest": round(float(eb), 4),
           "q_Ag1": round(q1, 4), "n_sites": len(res)}
    if best2:
        out["q_Ag2_per_atom"] = round(float(best2[1][n:].mean()), 4)
        out["hand_back"] = round(q1 - out["q_Ag2_per_atom"], 4)
        out["E_bind_Ag2"] = round(float(e_mol + 2 * e_ag - best2[0]), 4)
    return out


def main():
    only = set(sys.argv[1:])
    done = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for tag, fn, rule, mult in CANDIDATES:
        if tag in SKIP or (only and tag not in only) or tag in done:
            continue
        try:
            done[tag] = survey(tag, fn)
        except Exception as ex:
            print(f"{tag}: FAILED {ex!r}", flush=True)
            continue
        json.dump(done, open(OUT, "w"), indent=1)
        r = done[tag]
        print(f"{tag:<11} {r['n_sites']:>3} sites  deepest {r['deepest_site']:<20} "
              f"E_b {r['E_b_deepest']:.3f}  q1 {r['q_Ag1']:+.3f}  "
              f"q2/atom {r.get('q_Ag2_per_atom', float('nan')):+.3f}", flush=True)


if __name__ == "__main__":
    main()
