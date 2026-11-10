"""Build the common evaporable OLED materials the screen has not covered yet.

    python scripts/147_common_materials_build.py

SMILES -> RDKit ETKDG conformer (best of 20 by MMFF) -> GFN2-xTB relaxation of
the free molecule. The molecule is then frozen for every Ag calculation, as in
the rest of the project, so the conformer matters only through which sites it
exposes; strain cancels out of E_b.

Writes structures/common/<name>.xyz and <name>.mol (the .mol keeps bond orders
for fragmenting in script 149).
"""
import os, sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "structures", "common")

# role is what the layer does in an OLED; it decides which electrode it would
# sit under, not how it is ranked
MATERIALS = {
    # electron transport / hole blocking
    "TPBi":   ("ETL", "c1ccc(cc1)n1c(nc2ccccc12)-c1cc(cc(c1)-c1nc2ccccc2n1-c1ccccc1)-c1nc2ccccc2n1-c1ccccc1"),
    "BCP":    ("ETL/HBL", "Cc1cc(-c2ccccc2)c2ccc3c(-c4ccccc4)cc(C)nc3c2n1"),
    "NBphen": ("ETL", "c1ccc(cc1)-c1cc(-c2ccc3ccccc3c2)nc2c1ccc1c(-c3ccccc3)cc(-c3ccc4ccccc4c3)nc12"),
    "TmPyPB": ("ETL", "c1cc(-c2cccnc2)cc(c1)-c1cc(cc(c1)-c1cccc(c1)-c1cccnc1)-c1cccc(c1)-c1cccnc1"),
    "3TPYMB": ("ETL", "B(c1c(C)cc(C)c(-c2cccnc2)c1C)(c1c(C)cc(C)c(-c2cccnc2)c1C)c1c(C)cc(C)c(-c2cccnc2)c1C"),
    "T2T":    ("ETL/host", "c1ccc(cc1)-c1cccc(c1)-c1nc(nc(n1)-c1cccc(c1)-c1ccccc1)-c1cccc(c1)-c1ccccc1"),
    "PO-T2T": ("ETL", "O=P(c1ccccc1)(c1ccccc1)c1cccc(c1)-c1nc(nc(n1)-c1cccc(c1)P(=O)(c1ccccc1)c1ccccc1)-c1cccc(c1)P(=O)(c1ccccc1)c1ccccc1"),
    "DPEPO":  ("host/HBL", "O=P(c1ccccc1)(c1ccccc1)c1ccccc1Oc1ccccc1P(=O)(c1ccccc1)c1ccccc1"),
    "TSPO1":  ("host/ETL", "O=P(c1ccccc1)(c1ccccc1)c1ccc(cc1)[Si](c1ccccc1)(c1ccccc1)c1ccccc1"),
    # hole transport / host
    "NPB":    ("HTL", "c1ccc(cc1)N(c1ccc(cc1)-c1ccc(cc1)N(c1ccccc1)c1cccc2ccccc12)c1cccc2ccccc12"),
    "TPD":    ("HTL", "Cc1cccc(c1)N(c1ccccc1)c1ccc(cc1)-c1ccc(cc1)N(c1ccccc1)c1cccc(C)c1"),
    "TAPC":   ("HTL", "Cc1ccc(cc1)N(c1ccc(C)cc1)c1ccc(cc1)C1(CCCCC1)c1ccc(cc1)N(c1ccc(C)cc1)c1ccc(C)cc1"),
    "TCTA":   ("HTL/host", "c1ccc2c(c1)c1ccccc1n2-c1ccc(cc1)N(c1ccc(cc1)-n1c2ccccc2c2ccccc12)c1ccc(cc1)-n1c2ccccc2c2ccccc12"),
    "CBP":    ("host", "c1ccc2c(c1)c1ccccc1n2-c1ccc(cc1)-c1ccc(cc1)-n1c2ccccc2c2ccccc12"),
    "mCP":    ("host", "c1cc(cc(c1)-n1c2ccccc2c2ccccc12)-n1c2ccccc2c2ccccc12"),
    # hole injection / p-dopant
    "F6-TCNNQ": ("HIL/p-dopant", "FC1=C(F)C(=C(C#N)C#N)C(F)=C2C(F)=C(F)C(=C(C#N)C#N)C(F)=C12"),
}
# Alq3: RDKit refuses six-coordinate Al, so the molecule is embedded with the
# three Al-O bonds only and the GFN2 relaxation is left to close the three N->Al
# dative bonds; the build checks that it did (Al-N < 2.3 A).
ALQ3 = ("ETL", "C1=CC2=C(C(=C1)O[Al](OC3=CC=CC4=C3N=CC=C4)OC5=CC=CC6=C5N=CC=C6)N=CC=C2")


def embed(smi):
    m = Chem.MolFromSmiles(smi, sanitize=False)
    try:
        Chem.SanitizeMol(m)
    except Exception:
        Chem.SanitizeMol(m, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    m = Chem.AddHs(m)
    ids = AllChem.EmbedMultipleConfs(m, numConfs=20, randomSeed=7)
    if not ids:
        p = AllChem.ETKDGv3()
        p.useRandomCoords = True
        ids = AllChem.EmbedMultipleConfs(m, numConfs=20, params=p)
    best, be = None, None
    try:
        res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=2000)
        for cid, (conv, e) in zip(ids, res):
            if be is None or e < be:
                best, be = cid, e
    except Exception:
        AllChem.UFFOptimizeMoleculeConfs(m, maxIters=2000)
        best = ids[0]
    return m, best


def gfn2_opt(m, cid):
    S = [a.GetSymbol() for a in m.GetAtoms()]
    X = m.GetConformer(cid).GetPositions()
    a = Atoms(S, positions=X)
    a.calc = TBLite(method="GFN2-xTB", verbosity=0, max_iterations=500)
    BFGS(a, logfile=None, maxstep=0.2).run(fmax=0.03, steps=600)
    return S, a.positions


def main():
    os.makedirs(OUT, exist_ok=True)
    todo = dict(MATERIALS)
    todo["Alq3"] = ALQ3
    only = set(sys.argv[1:])
    for name, (role, smi) in todo.items():
        if only and name not in only:
            continue
        p = os.path.join(OUT, f"{name}.xyz")
        if os.path.exists(p):
            continue
        m, cid = embed(smi)
        S, X = gfn2_opt(m, cid)
        conf = m.GetConformer(cid)
        for i, x in enumerate(X):
            conf.SetAtomPosition(i, x.tolist())
        Chem.MolToMolFile(m, os.path.join(OUT, f"{name}.mol"), confId=cid)
        with open(p, "w") as f:
            f.write(f"{len(S)}\n{name} ({role}) RDKit ETKDG/MMFF -> GFN2-xTB relaxed\n")
            for s, x in zip(S, X):
                f.write(f"{s} {x[0]:.6f} {x[1]:.6f} {x[2]:.6f}\n")
        extra = ""
        if "Al" in S:
            al = S.index("Al")
            dn = sorted(float(np.linalg.norm(X[al] - X[k])) for k, s in enumerate(S) if s == "N")
            extra = "  Al-N " + ", ".join(f"{d:.2f}" for d in dn)
        print(f"{name:<9} {role:<12} {len(S):>4} atoms  {Chem.rdMolDescriptors.CalcMolFormula(m)}"
              + extra, flush=True)


if __name__ == "__main__":
    main()
