"""GFN2-xTB Ag site survey + Ag2 charge for the design library.

    python scripts/162_design_sites.py [--workers 4] [--max N]  ->  runs/design/sites.jsonl

Same survey as script 148 (lone-pair starts, heteroatom pockets, ring hollows
and heteroatom atop sites; Ag relaxed on the frozen GFN2 molecule), with two
changes for a library of symmetric molecules:
  * starts are de-duplicated by RDKit canonical atom class (the xyz written by
    script 161 keeps RDKit atom order), so HATCN gets 1 lone-pair start per
    N class rather than 12
  * one face of a pocket only (the molecules are planar or C2/C3 symmetric)

Kept per molecule: GFN2 E_b and q(Ag1) at the deepest site, the nearest atom,
the shortest Ag-H distance (a crowded start, Ag-H < 2.4 A, is the case where
GFN2 and PBE0 disagree -- script 151), the top-3 distinct sites, and
dq(Ag1->Ag2) at the deepest site. GFN2's absolute E_b overbinds N-rich
acceptors (HATCN 3.29 vs PBE0 1.63); script 163 maps it to PBE0.
Order: highest EA first (from props.jsonl), so the interesting part of the
library finishes first.
"""
import importlib
import json
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
LIB = os.path.join(ROOT, "runs", "design", "library.json")
PROPS = os.path.join(ROOT, "runs", "design", "props.jsonl")
OUT = os.path.join(ROOT, "runs", "design", "sites.jsonl")
XYZ = os.path.join(ROOT, "structures", "design")
D_AG = 2.60
HET = {"N", "O", "F", "Cl"}


def ranks(smiles):
    from rdkit import Chem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return list(Chem.CanonicalRankAtoms(m, breakTies=False))


def starts(S, X, rk):
    from pathgeom import _bonded, contact
    s131 = importlib.import_module("131_site_spread")
    s131.MAX_SITES = 1000
    B = _bonded(S, X)
    cen = X.mean(0)
    nrm = np.linalg.svd(X - cen)[2][-1]
    out, seen = [], set()
    lim = np.array([contact(x) for x in S])

    def push(p, u):
        for _ in range(40):
            if float((np.linalg.norm(X - p, axis=1) - lim).min()) >= 0:
                break
            p = p + 0.1 * u
        return p

    for i, s in enumerate(S):  # lone pairs
        if s not in ("N", "O", "F"):
            continue
        nb = [j for j in range(len(S)) if B[i][j]]
        if not nb or (s == "N" and len(nb) >= 3) or ("lp", rk[i]) in seen:
            continue
        u = -sum((X[j] - X[i]) / np.linalg.norm(X[j] - X[i]) for j in nb)
        if np.linalg.norm(u) < 1e-6:
            continue
        seen.add(("lp", rk[i]))
        u /= np.linalg.norm(u)
        out.append((f"lp{s}{i}", X[i] + 2.3 * u))
    het = [i for i, s in enumerate(S) if s in ("N", "O", "F")]
    for a in range(len(het)):  # pockets
        for b in range(a + 1, len(het)):
            i, j = het[a], het[b]
            d = float(np.linalg.norm(X[i] - X[j]))
            key = ("pk", tuple(sorted((rk[i], rk[j]))))
            if B[i][j] or not 2.4 < d < 3.6 or key in seen:
                continue
            seen.add(key)
            mid = (X[i] + X[j]) / 2
            v = mid - cen
            for tag, u in (("in", v / (np.linalg.norm(v) + 1e-9)), ("face", nrm)):
                if np.linalg.norm(u) < 0.5:
                    continue
                out.append((f"pocket{S[i]}{i}{S[j]}{j}{tag}", push(mid + 2.0 * u, u)))
    for kind, i, pos, clear in s131.unique_sites(S, X):  # hollows, heteroatom atop
        if not (kind == "hollow" or S[i] in HET):
            continue
        key = (kind, rk[i])
        if key in seen:
            continue
        seen.add(key)
        out.append((f"{kind}{i}", pos))
    return out


def work(item):
    name, smiles = item
    t0 = time.time()
    s142 = importlib.import_module("142_xtb_survey")
    from pathgeom import read_xyz
    S, X = read_xyz(os.path.join(XYZ, f"{name}.xyz"))
    X = np.asarray(X, float)
    n = len(S)
    rk = ranks(smiles)
    if len(rk) != n:
        return {"name": name, "error": "atom order mismatch"}
    e_mol, _ = s142.energy(S, X, 0)
    e_ag, _ = s142.energy(["Ag"], [[0, 0, 0]], 1)
    H = np.array([s == "H" for s in S])
    res = {}
    for lab, pos in starts(S, X, rk):
        r = s142.relax(S + ["Ag"], np.vstack([X, pos]), n, 1)
        if r is None:
            continue
        e, q, P = r
        d = np.linalg.norm(X - P[n], axis=1)
        k = int(d.argmin())
        res[lab] = {"E_b": round(float(e_mol + e_ag - e), 4), "q": round(float(q[n]), 4),
                    "near": f"{S[k]}{k}", "d": round(float(d.min()), 3),
                    "dH": round(float(d[H].min()), 3) if H.any() else None,
                    "pos": [round(float(v), 4) for v in P[n]]}
    if not res:
        return {"name": name, "error": "no site converged"}
    ranked = sorted(res.items(), key=lambda kv: -kv[1]["E_b"])
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
    out = {"name": name, "n_sites": len(res), "deepest": ranked[0][0],
           "E_b_gfn2": deep["E_b"], "q_Ag1": deep["q"], "near": deep["near"],
           "near_el": "".join(c for c in deep["near"] if c.isalpha()),
           "dH": deep["dH"], "crowded": bool(deep["dH"] is not None and deep["dH"] < 2.4),
           "top3": {p: res[p] for p in pick}, "sec": round(time.time() - t0)}
    if best2:
        q2 = best2[1][n:]
        out["q_Ag2_total"] = round(float(q2.sum()), 4)
        out["dq"] = round(float(q2.sum() - deep["q"]), 4)
    return out


def main():
    from multiprocessing import Pool
    a = sys.argv
    workers = int(a[a.index("--workers") + 1]) if "--workers" in a else 4
    nmax = int(a[a.index("--max") + 1]) if "--max" in a else None
    lib = json.load(open(LIB))
    props = {json.loads(l)["name"]: json.loads(l) for l in open(PROPS) if l.strip()}
    done = set()
    if os.path.exists(OUT):
        done = {json.loads(l)["name"] for l in open(OUT) if l.strip()}
    ea_min = float(a[a.index("--ea-min") + 1]) if "--ea-min" in a else -99
    refs = {k for k, v in lib.items() if v["family"] == "reference"}
    todo = [k for k, p in props.items() if "error" not in p and k not in done
            and (p["EA"] >= ea_min or k in refs)]
    todo.sort(key=lambda k: -props[k]["EA"])
    if "--names" in a:
        want = a[a.index("--names") + 1].split(",")
        todo = [k for k in want if k not in done]
    if nmax:
        todo = todo[:nmax]
    print(f"{len(todo)} to do, {len(done)} done", flush=True)
    with Pool(workers, maxtasksperchild=5) as p, open(OUT, "a") as f:
        for r in p.imap_unordered(work, [(k, lib[k]["smiles"]) for k in todo]):
            f.write(json.dumps(r) + "\n")
            f.flush()
            print(r["name"], r.get("E_b_gfn2"), r.get("q_Ag1"), r.get("dq"), r.get("near"),
                  r.get("n_sites"), r.get("sec"), r.get("error", ""), flush=True)


if __name__ == "__main__":
    main()
