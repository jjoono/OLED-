"""GFN2-xTB site survey and Ag2 charge for the common OLED materials.

    python scripts/148_common_materials_survey.py [names...]

Same survey as script 142 -- atop/bridge/hollow sites by Morgan class plus
heteroatom pockets, Ag relaxed on the frozen molecule -- restricted to the sites
that can plausibly be deepest on a 50-105 atom molecule: anything whose anchor
is a heteroatom, every ring hollow, and every heteroatom pocket. Pure C-C atop
and bridge sites on hydrocarbon rims never won on any of the 25 candidates
already surveyed, and on molecules this size there are dozens of them.

GFN2 is used only to FIND sites (validated 4/4 against the PBE0 surveys on
organics) and for the Ag2 charge (validated ordering on all sixteen package-G
geometries). Its E_b is not used for ranking -- its rank correlation with PBE0
on organics is only 0.66, because it overbinds nitrile nitrogens 4-5x. The
three deepest chemically distinct sites per molecule go on to PBE0 in script
149.
"""
import os, sys, json, importlib

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathgeom import read_xyz
s142 = importlib.import_module("142_xtb_survey")
s131 = s142.s131

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
COMMON = os.path.join(ROOT, "structures", "common")
OUT = os.path.join(ROOT, "runs", "common_survey_xtb.json")
HET = {"N", "O", "P", "S", "F", "Si", "Al", "B"}
D_AG = 2.60


def lone_pair_sites(S, X):
    """A start along the lone pair of every N, O and F, 2.3 A out.

    The atop starts of script 131 aim along the radial direction from the
    molecular centroid, which on a twisted molecule need not be where a
    pyridyl nitrogen's lone pair points: on 3TPYMB every atop start relaxed
    into the methyl hydrogens and the survey missed the N site entirely (GFN2
    0.414 against 0.467 at the nitrogen). The lone pair is the opposite of the
    sum of the bond directions.
    """
    from pathgeom import _bonded
    B = _bonded(S, X)
    out = []
    for i, s in enumerate(S):
        if s not in ("N", "O", "F"):
            continue
        nb = [j for j in range(len(S)) if B[i][j]]
        if not nb or (s == "N" and len(nb) >= 3):
            continue
        u = -sum((X[j] - X[i]) / np.linalg.norm(X[j] - X[i]) for j in nb)
        if np.linalg.norm(u) < 1e-6:
            continue
        u /= np.linalg.norm(u)
        out.append((f"lp{s}{i}", X[i] + 2.3 * u))
    return out


def sites(S, X):
    out = lone_pair_sites(S, X)
    for kind, i, pos, clear in s131.unique_sites(S, X):
        if kind == "hollow" or S[i] in HET:
            out.append((f"{kind}{i}", pos))
    out += s142.pockets(S, X)
    return out


def survey(name, previous=None):
    S, X = read_xyz(os.path.join(COMMON, f"{name}.xyz"))
    X = np.asarray(X, float)
    n = len(S)
    e_mol, _ = s142.energy(S, X, 0)
    e_ag, _ = s142.energy(["Ag"], [[0, 0, 0]], 1)
    res = dict(previous["sites"]) if previous else {}
    todo = lone_pair_sites(S, X) if previous else sites(S, X)
    for lab, pos in todo:
        if lab in res:
            continue
        r = s142.relax(S + ["Ag"], np.vstack([X, pos]), n, 1)
        if r is None:
            continue
        e, q, P = r
        d = np.linalg.norm(X - P[n], axis=1)
        k = int(d.argmin())
        res[lab] = {"E_b": round(float(e_mol + e_ag - e), 4), "q": round(float(q[n]), 4),
                    "near": f"{S[k]}{k}", "near_el": S[k], "d": round(float(d.min()), 3),
                    "pos": [round(float(v), 4) for v in P[n]]}
    ranked = sorted(res.items(), key=lambda kv: -kv[1]["E_b"])
    # three deepest chemically distinct sites: distinct nearest atom, and at
    # least 1.5 A apart, so that equivalent relaxed endpoints are not repeated
    pick = []
    for lab, v in ranked:
        if any(np.linalg.norm(np.array(v["pos"]) - np.array(res[p]["pos"])) < 1.5 for p in pick):
            continue
        pick.append(lab)
        if len(pick) == 3:
            break
    deep = ranked[0][1]
    ag1 = np.array(deep["pos"])
    d = np.linalg.norm(X - ag1, axis=1)
    nv = ag1 - X[int(d.argmin())]
    nv /= np.linalg.norm(nv)
    t1 = np.cross(nv, [0, 0, 1.0])
    if np.linalg.norm(t1) < 1e-6:
        t1 = np.cross(nv, [1.0, 0, 0])
    t1 /= np.linalg.norm(t1)
    t2 = np.cross(nv, t1)
    best2 = None
    for th in (0, 2 * np.pi / 3, 4 * np.pi / 3):
        p2 = ag1 + D_AG * (np.cos(th) * t1 + np.sin(th) * t2)
        r = s142.relax(S + ["Ag", "Ag"], np.vstack([X, ag1, p2]), n, 2)
        if r and (best2 is None or r[0] < best2[0]):
            best2 = r
    out = {"n_atoms": n, "n_sites": len(res), "sites": res, "pbe0_sites": pick,
           "deepest_site": ranked[0][0], "q_Ag1": deep["q"]}
    if best2:
        q2 = best2[1][n:]
        out["q_Ag2_total"] = round(float(q2.sum()), 4)
        out["delta_q_Ag1_Ag2"] = round(float(q2.sum() - deep["q"]), 4)
    return out


def main():
    args = sys.argv[1:]
    supplement = "--lone-pairs" in args
    args = [a for a in args if not a.startswith("--")]
    names = args or sorted(f[:-4] for f in os.listdir(COMMON) if f.endswith(".xyz"))
    for name in names:
        done = json.load(open(OUT)) if os.path.exists(OUT) else {}
        prev = done.get(name)
        if supplement:
            if prev is None or prev.get("lone_pairs_added"):
                continue
        elif prev is not None:
            continue
        try:
            r = survey(name, prev if supplement else None)
            r["lone_pairs_added"] = True
        except Exception as ex:
            print(f"{name}: FAILED {ex!r}", flush=True)
            continue
        done = json.load(open(OUT)) if os.path.exists(OUT) else {}
        done[name] = r
        json.dump(done, open(OUT, "w"), indent=1)
        print(f"{name:<9} {r['n_atoms']:>4} atoms {r['n_sites']:>3} sites  deepest "
              f"{r['deepest_site']} ({r['sites'][r['deepest_site']]['near']})  "
              f"q1 {r['q_Ag1']:+.3f}  dq(Ag1->Ag2) {r.get('delta_q_Ag1_Ag2', float('nan')):+.3f}",
              flush=True)


if __name__ == "__main__":
    main()
