"""Scouting (GFN2): does HATCN stay a good Ag seed once its underlayer has given it an electron?

    python scripts/156_reduced_hatcn_gfn2.py  ->  runs/reduced_hatcn_gfn2.json

On glass HATCN is neutral and takes charge from Ag. On an n-doped ETL (the
charge-generation junction) or on an HTL, electrons flow from the layer
below into HATCN. If the top molecules are reduced, the one-electron
capacity that made HATCN a good seed is already used. Three models, Ag in
the same aza pocket, molecule frozen:

  neutral   HATCN                                  (glass)
  ion_pair  HATCN(-) with Cs(+) 3.5 A below the ring centre, on the face away
            from Ag -- an electron supplied from the layer underneath
  anion     HATCN(-) alone, charge -1 (upper bound: no counter-charge)

For each: Ag1 E_b and charge, and Ag2 at the same pocket (dq Ag1->Ag2).
GFN2 overbinds N sites 2-5x; only the direction between models is read.
"""
import json
import os
import sys

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from tblite.ase import TBLite

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pathgeom import read_xyz  # noqa: E402

ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "runs", "reduced_hatcn_gfn2.json")


def calc(charge, mult, et=300.0):
    return TBLite(method="GFN2-xTB", verbosity=0, charge=charge, multiplicity=mult,
                  max_iterations=400, mixer_damping=0.2, electronic_temperature=et)


def energy(S, X, charge, mult):
    a = Atoms(S, positions=X)
    a.calc = calc(charge, mult)
    return a.get_potential_energy(), a.get_charges()


def relax(S, X, nfix, charge, mult):
    for et in (300.0, 1500.0):
        a = Atoms(S, positions=X)
        a.set_constraint(FixAtoms(indices=list(range(nfix))))
        a.calc = calc(charge, mult, et)
        try:
            BFGS(a, logfile=None, maxstep=0.2).run(fmax=0.05, steps=300)
            return a.get_potential_energy(), a.get_charges(), a.positions.copy()
        except Exception:
            continue
    return None


def main():
    S0, X0 = read_xyz(os.path.join(ROOT, "structures", "clusters", "HATCN_Ag1.xyz"))
    X0 = np.asarray(X0, float)
    k = S0.index("Ag")
    ag_start = X0[k]
    S = [s for i, s in enumerate(S0) if i != k]
    X = np.delete(X0, k, axis=0)
    cen = X.mean(axis=0)
    nrm = np.linalg.svd(X - cen)[2][-1]
    if np.dot(ag_start - cen, nrm) < 0:
        nrm = -nrm
    cs = cen - 3.5 * nrm          # under the ring centre, away from the Ag side

    models = {
        "neutral":  (S, X, 0, 1),
        "ion_pair": (S + ["Cs"], np.vstack([X, cs]), 0, 2),
        "anion":    (S, X, -1, 2),
    }
    _, _, Xc = relax(S + ["Cs"], np.vstack([X, cs]), len(S), 0, 2)  # Cs settles, HATCN frozen
    models["ion_pair"] = (S + ["Cs"], Xc, 0, 2)

    e_ag, _ = energy(["Ag"], [[0, 0, 0]], 0, 2)
    res = {}
    for name, (Sm, Xm, ch, mult) in models.items():
        n = len(Sm)
        e_m, q_m = energy(Sm, Xm, ch, mult)
        mult1 = 1 if mult == 2 else 2
        r1 = relax(Sm + ["Ag"], np.vstack([Xm, ag_start]), n, ch, mult1)
        e1, q1, P1 = r1
        # Ag2: second atom 2.6 A from the first, along the pocket edge
        u = np.cross(nrm, P1[n] - cen)
        u /= np.linalg.norm(u)
        r2 = relax(Sm + ["Ag", "Ag"], np.vstack([Xm, P1[n], P1[n] + 2.6 * u]), n, ch, mult)
        e2, q2, P2 = r2
        row = {"E_b": round(float(e_m + e_ag - e1), 3), "q_Ag1": round(float(q1[n]), 3),
               "q_HATCN_before": round(float(q_m[:30].sum()), 3),
               "q_HATCN_with_Ag": round(float(q1[:30].sum()), 3),
               "q_Ag2_total": round(float(q2[n:].sum()), 3),
               "dq_Ag1_Ag2": round(float(q2[n:].sum() - q1[n]), 3),
               "Ag_N_min_A": round(float(min(np.linalg.norm(P1[n] - P1[i])
                                             for i in range(30) if Sm[i] == "N")), 3)}
        if "Cs" in Sm:
            row["q_Cs"] = round(float(q_m[Sm.index("Cs")]), 3)
            row["Ag_Cs_A"] = round(float(np.linalg.norm(P1[n] - P1[Sm.index("Cs")])), 3)
        res[name] = row
        print(name, row, flush=True)
    json.dump(res, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
