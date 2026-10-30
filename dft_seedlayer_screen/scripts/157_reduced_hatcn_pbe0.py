"""PBE0 check of script 156: Ag on HATCN that its underlayer has reduced.

    QC_THREADS=4 python scripts/157_reduced_hatcn_pbe0.py  ->  runs/reduced_hatcn_pbe0.json

Models and Ag positions as in 156 (GFN2, molecule frozen):
  neutral   HATCN + Ag at the aza pocket             (glass)
  ion_pair  HATCN(-) Cs(+), Cs 3.5 A below the ring centre on the far face
Whole molecule, PBE0-D3BJ/def2-SVP (31-32 atoms; no fragment needed). The
ion-pair complex is run as singlet and triplet, the lower is used. The
neutral case is recomputed at its own GFN2 Ag position so the two values are
a matched pair (the published 1.631 eV used the PBE0-optimised position).
"""
import importlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qclocal as q  # noqa: E402
from pathgeom import read_xyz  # noqa: E402

s156 = importlib.import_module("156_reduced_hatcn_gfn2")
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "runs", "reduced_hatcn_pbe0.json")
EV = 27.211386


def geometries():
    S0, X0 = read_xyz(os.path.join(ROOT, "structures", "clusters", "HATCN_Ag1.xyz"))
    X0 = np.asarray(X0, float)
    k = S0.index("Ag")
    ag0 = X0[k]
    S = [s for i, s in enumerate(S0) if i != k]
    X = np.delete(X0, k, axis=0)
    cen = X.mean(axis=0)
    nrm = np.linalg.svd(X - cen)[2][-1]
    if np.dot(ag0 - cen, nrm) < 0:
        nrm = -nrm
    _, _, Xc = s156.relax(S + ["Cs"], np.vstack([X, cen - 3.5 * nrm]), len(S), 0, 2)
    out = {}
    for name, Sm, Xm, mult in (("neutral", S, X, 1), ("ion_pair", S + ["Cs"], Xc, 2)):
        n = len(Sm)
        _, _, P = s156.relax(Sm + ["Ag"], np.vstack([Xm, ag0]), n, 0, 1 if mult == 2 else 2)
        out[name] = (Sm, Xm, mult, P[n])
    return out


def main():
    from pyscf import lib
    lib.num_threads(int(os.environ.get("QC_THREADS", "4")))
    d = json.load(open(OUT)) if os.path.exists(OUT) else {}
    geo = geometries()
    if "E_ag" not in d:
        d["E_ag"] = float(q.scf(q.mol(["Ag"], [[0, 0, 0]], 2))[1])
    for name, (Sm, Xm, mult, ag) in geo.items():
        row = d.setdefault(name, {})
        if "E_b" in row:
            continue
        t0 = time.time()
        if "E_mol" not in row:
            row["E_mol"] = float(q.scf(q.mol(Sm, Xm, mult))[1])
            json.dump(d, open(OUT, "w"), indent=1)
        S2, X2 = Sm + ["Ag"], np.vstack([Xm, ag])
        states = (2,) if mult == 1 else (1, 3)
        for m in states:
            if f"E_m{m}" in row:
                continue
            mf, e = q.scf(q.mol(S2, X2, m))
            mq = q.mulliken(mf)
            row[f"E_m{m}"] = float(e)
            row[f"q_Ag_m{m}"] = round(float(mq[-1]), 4)
            row[f"q_HATCN_m{m}"] = round(float(mq[:30].sum()), 4)
            if "Cs" in Sm:
                row[f"q_Cs_m{m}"] = round(float(mq[Sm.index("Cs")]), 4)
            json.dump(d, open(OUT, "w"), indent=1)
            print(f"{name} mult {m} done", flush=True)
        e_cx = min(row[f"E_m{m}"] for m in states)
        best = min(states, key=lambda m: row[f"E_m{m}"])
        d3 = (q.d3bj(S2, X2) - q.d3bj(Sm, Xm)) * EV
        elec = (row["E_mol"] + d["E_ag"] - e_cx) * EV
        row.update(E_b=round(elec - d3, 4), elec=round(elec, 4), disp=round(-d3, 4),
                   state={1: "singlet", 2: "doublet", 3: "triplet"}[best],
                   q_Ag=row[f"q_Ag_m{best}"], sec=round(time.time() - t0),
                   Ag_N_min_A=round(float(min(np.linalg.norm(ag - Xm[i])
                                              for i, s in enumerate(Sm) if s == "N")), 3))
        json.dump(d, open(OUT, "w"), indent=1)
        print(f"E_b {name:<9} {row['E_b']:.3f} eV ({row['state']}, elec {row['elec']:.3f}, "
              f"q_Ag {row['q_Ag']:+.3f})", flush=True)


if __name__ == "__main__":
    main()
