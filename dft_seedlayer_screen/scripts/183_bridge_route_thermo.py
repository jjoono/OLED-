"""Thermodynamic checks for the bridge-forming route to Hz-C(CN)2-Hz (and Hz-NCN-Hz).

    python scripts/183_bridge_route_thermo.py     -> runs/design/bridge_route_thermo.json

GFN2-xTB + ALPB(acetonitrile), RDKit/MMFF start, one conformer. Reaction
energies only (no barriers, no entropy); screening level, sign and size.

Route assumed (from cyameluric chloride, see notes/SYNTHESIS_HZ_CCN_HZ.md):
  1  Cl3Hz + CH(CN)2(-)          -> Cl2Hz-C(CN)2(-) + HCl-equivalent
  2  Cl2Hz-C(CN)2(-) + Cl3Hz     -> Cl2Hz-C(CN)2-HzCl2 + Cl(-)
     competing: attack through a nitrile N -> ketenimine Cl2Hz-C(CN)=C=N-HzCl2
  3  Cl -> CN on the dimer with CN(-)
     competing: CN(-) attacks the bridgehead C of a heptazine and expels the
     stabilised Cl2Hz-C(CN)2(-) anion (bridge cleavage)
Same questions for the NCN bridge (N-cyano amine vs carbodiimide isomer).
"""
import json
import os

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "runs", "design", "bridge_route_thermo.json")
HA = 27.211386
GEO = os.path.join(ROOT, "structures", "bridge_route")


def core(s1, s2, s3):
    return f"{s1}c1nc2nc({s2})nc3nc({s3})nc(n1)n23"


def unit(s2, s3):   # second heptazine as a substituent
    return f"c4nc5nc({s2})nc6nc({s3})nc(n4)n56"


SP = {
    # cyano end products
    "HzCN3": core("N#C", "C#N", "C#N"),
    "CCN": core("N#C", "C(C#N)(C#N)" + unit("C#N", "C#N"), "C#N"),
    "CCN_ketenimine": core("N#C", "C(C#N)=C=N" + unit("C#N", "C#N"), "C#N"),
    "CCN_anion": core("N#C", "[C-](C#N)C#N", "C#N"),
    "NCN": core("N#C", "N(C#N)" + unit("C#N", "C#N"), "C#N"),
    "NCN_carbodiimide": core("N#C", "N=C=N" + unit("C#N", "C#N"), "C#N"),
    "NCN_anion": core("N#C", "[N-]C#N", "C#N"),
    # chloro intermediates
    "HzCl3": core("Cl", "Cl", "Cl"),
    "Cl2Hz_CN": core("Cl", "C#N", "Cl"),
    "Cl_CCN_anion": core("Cl", "[C-](C#N)C#N", "Cl"),
    "Cl_CCN": core("Cl", "C(C#N)(C#N)" + unit("Cl", "Cl"), "Cl"),
    "Cl_CCN_ketenimine": core("Cl", "C(C#N)=C=N" + unit("Cl", "Cl"), "Cl"),
    "Cl_CCN_1CN": core("N#C", "C(C#N)(C#N)" + unit("Cl", "Cl"), "Cl"),
    "Cl_NCN_anion": core("Cl", "[N-]C#N", "Cl"),
    "Cl_NCN": core("Cl", "N(C#N)" + unit("Cl", "Cl"), "Cl"),
    "Cl_NCN_carbodiimide": core("Cl", "N=C=N" + unit("Cl", "Cl"), "Cl"),
    "Cl_NCN_1CN": core("N#C", "N(C#N)" + unit("Cl", "Cl"), "Cl"),
    "CN-": "[C-]#N",
    "Cl-": "[Cl-]",
    "malononitrile_anion": "[CH-](C#N)C#N",
    "cyanamide_anion": "[NH-]C#N",
    "HCN": "C#N",
    "HCl": "Cl",
}

RX = {
    "step2_C_attack (target)": ({"Cl_CCN_anion": 1, "HzCl3": 1}, {"Cl_CCN": 1, "Cl-": 1}),
    "step2_N_attack (ketenimine)": ({"Cl_CCN_anion": 1, "HzCl3": 1}, {"Cl_CCN_ketenimine": 1, "Cl-": 1}),
    "step2_NCN_N_attack (target)": ({"Cl_NCN_anion": 1, "HzCl3": 1}, {"Cl_NCN": 1, "Cl-": 1}),
    "step2_NCN_carbodiimide": ({"Cl_NCN_anion": 1, "HzCl3": 1}, {"Cl_NCN_carbodiimide": 1, "Cl-": 1}),
    "step3_CCN_Cl_to_CN (wanted)": ({"Cl_CCN": 1, "CN-": 1}, {"Cl_CCN_1CN": 1, "Cl-": 1}),
    "step3_CCN_bridge_cleavage": ({"Cl_CCN": 1, "CN-": 1}, {"Cl2Hz_CN": 1, "Cl_CCN_anion": 1}),
    "step3_NCN_Cl_to_CN (wanted)": ({"Cl_NCN": 1, "CN-": 1}, {"Cl_NCN_1CN": 1, "Cl-": 1}),
    "step3_NCN_bridge_cleavage": ({"Cl_NCN": 1, "CN-": 1}, {"Cl2Hz_CN": 1, "Cl_NCN_anion": 1}),
    "final_CCN_isomer (ketenimine - target)": ({"CCN": 1}, {"CCN_ketenimine": 1}),
    "final_NCN_isomer (carbodiimide - target)": ({"NCN": 1}, {"NCN_carbodiimide": 1}),
    "final_CCN_cleavage_by_CN-": ({"CCN": 1, "CN-": 1}, {"HzCN3": 1, "CCN_anion": 1}),
    "final_NCN_cleavage_by_CN-": ({"NCN": 1, "CN-": 1}, {"HzCN3": 1, "NCN_anion": 1}),
}


def embed(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    ch = Chem.GetFormalCharge(m)
    if m.GetNumAtoms() == 1:
        return [a.GetSymbol() for a in m.GetAtoms()], np.zeros((1, 3)), ch
    ids = list(AllChem.EmbedMultipleConfs(m, numConfs=8, randomSeed=7))
    if not ids:
        p = AllChem.ETKDGv3(); p.useRandomCoords = True; p.randomSeed = 7
        ids = list(AllChem.EmbedMultipleConfs(m, numConfs=8, params=p))
    try:
        res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=2000)
    except Exception:
        res = AllChem.UFFOptimizeMoleculeConfs(m, maxIters=2000)
    best = min(zip(ids, res), key=lambda r: r[1][1])[0]
    return [a.GetSymbol() for a in m.GetAtoms()], m.GetConformer(best).GetPositions(), ch


def gfn2(S, X, ch):
    from ase import Atoms
    from ase.optimize import BFGS
    from tblite.ase import TBLite
    a = Atoms(S, positions=X)
    a.calc = TBLite(method="GFN2-xTB", charge=ch, multiplicity=1, verbosity=0,
                    solvation=("alpb", "acetonitrile"))
    if len(S) > 1:
        BFGS(a, logfile=None, maxstep=0.1).run(fmax=0.03, steps=2000)
    return float(a.get_potential_energy()), a.positions.copy()


def main():
    out = json.load(open(OUT)) if os.path.exists(OUT) else {"E_eV": {}}
    E = out["E_eV"]
    for name, smi in SP.items():
        if name in E:
            continue
        S, X, ch = embed(smi)
        E[name], XO = gfn2(S, X, ch)
        os.makedirs(GEO, exist_ok=True)
        with open(os.path.join(GEO, f"{name}.xyz"), "w") as f:
            f.write(f"{len(S)}\n{name} GFN2/ALPB(MeCN) charge {ch}\n")
            f.writelines(f"{a} {x:.6f} {y:.6f} {z:.6f}\n" for a, (x, y, z) in zip(S, XO))
        json.dump(out, open(OUT, "w"), indent=1)
        print(f"{name:22s} q={ch:+d} E={E[name]:.3f} eV", flush=True)
    out["reactions_eV"] = {k: round(sum(n * E[s] for s, n in p.items()) - sum(n * E[s] for s, n in r.items()), 2)
                           for k, (r, p) in RX.items()}
    json.dump(out, open(OUT, "w"), indent=1)
    for k, v in out["reactions_eV"].items():
        print(f"{k:45s} {v:+.2f} eV")


if __name__ == "__main__":
    main()
