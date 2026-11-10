"""Is the 0.53 eV gap between Mo3O9 and HATCN a basis-set artefact?

    python scripts/143_eb_bsse_basis.py

E_b is the first axis of the two-dimensional criterion (bind hard enough to
nucleate densely, leave silver metallic), and every E_b in the project is
PBE0-D3/def2-SVP without counterpoise. def2-SVP has no diffuse functions and a
small basis overbinds through basis-set superposition, by an amount that depends
on how many basis functions crowd the adatom -- so it need not cancel between an
oxide and an organic.

For the deepest site of each package-G substrate (geometries in
structures/clusters/<tag>_Ag1.xyz, the same ones whose E_b is quoted):

    E_b      = E(mol) + E(Ag) - E(complex)                       as published
    E_b^CP   = E(mol in full basis) + E(Ag in full basis) - E(complex)
    def2-TZVP, with and without counterpoise, for HATCN and Mo3O9

The D3 energy depends only on geometry and is identical in every variant, so it
is added once. The def2-SVP uncorrected value is the check that PySCF reproduces
the published E_b.

Restartable; results in runs/eb_bsse_basis.json.
"""
import os, sys, json, time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qclocal as q
from pathgeom import STRUCT, read_xyz
from pyscf import lib

lib.num_threads(int(os.environ.get("QC_THREADS", "12")))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "runs",
                   "eb_bsse_basis.json")
EV = 27.211386
PUBLISHED = {"HATCN": 1.631, "Mo3O9": 2.164, "F4TCNQ": 1.208, "benzene": 0.205}
PLAN = [("benzene", "def2-svp"), ("Mo3O9", "def2-svp"), ("F4TCNQ", "def2-svp"),
        ("HATCN", "def2-svp"), ("Mo3O9", "def2-tzvp"), ("HATCN", "def2-tzvp")]


def e(S, X, mult, basis, ghost=()):
    mf, en = q.scf(q.mol(S, X, mult, basis, ghost=ghost), conv=1e-9, grid=3)
    if not mf.converged:
        raise RuntimeError("SCF not converged")
    return float(en)


def main():
    done = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for tag, basis in PLAN:
        key = f"{tag}_{basis}"
        if key in done:
            continue
        t0 = time.time()
        S, X = read_xyz(os.path.join(STRUCT, "clusters", f"{tag}_Ag1.xyz"))
        X = np.asarray(X, float)
        k = S.index("Ag")
        mol_idx = [i for i in range(len(S)) if i != k]
        ms, mx = [S[i] for i in mol_idx], X[mol_idx]
        e_cx = e(S, X, 2, basis)
        e_mol = e(ms, mx, 1, basis)
        e_ag = e(["Ag"], X[[k]], 2, basis)
        e_mol_g = e(S, X, 1, basis, ghost={k})
        e_ag_g = e(S, X, 2, basis, ghost=set(mol_idx))
        d3 = q.d3bj(S, X) - q.d3bj(ms, mx)
        eb = (e_mol + e_ag - e_cx) * EV - d3 * EV
        eb_cp = (e_mol_g + e_ag_g - e_cx) * EV - d3 * EV
        done[key] = {"E_b": round(eb, 4), "E_b_CP": round(eb_cp, 4),
                     "BSSE": round(eb - eb_cp, 4), "D3_part": round(-d3 * EV, 4),
                     "published_svp": PUBLISHED[tag],
                     "seconds": round(time.time() - t0)}
        json.dump(done, open(OUT, "w"), indent=1)
        print(f"{key:<18} E_b {eb:.3f}  CP {eb_cp:.3f}  BSSE {eb - eb_cp:.3f}  "
              f"(published SVP {PUBLISHED[tag]:.3f})  {done[key]['seconds']}s", flush=True)
    report(done)


def report(done):
    print("\nGap Mo3O9 - HATCN in each variant:")
    for b in ("def2-svp", "def2-tzvp"):
        h, m = done.get(f"HATCN_{b}"), done.get(f"Mo3O9_{b}")
        if h and m:
            print(f"  {b:<10} raw {m['E_b'] - h['E_b']:+.3f} eV   "
                  f"CP {m['E_b_CP'] - h['E_b_CP']:+.3f} eV")


if __name__ == "__main__":
    main()
