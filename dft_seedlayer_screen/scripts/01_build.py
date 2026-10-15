"""Build seed-layer models + Ag-adsorbed starting structures.

Descriptor: Ag single-atom binding energy E_b on each seed surface model,
compared against Ag-Ag cohesion (Ag2 dimer at same level; bulk E_coh = 2.95 eV).
Strong E_b (>~ Ag-Ag) -> nucleation-dense, 2D-like growth; weak -> Volmer-Weber islands.
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import os

OUT = os.path.join(os.path.dirname(__file__), "..", "structures")
os.makedirs(OUT, exist_ok=True)

SMILES = {
    "HATCN":     "N#Cc1nc2c3nc(C#N)c(C#N)nc3c3nc(C#N)c(C#N)nc3c2nc1C#N",
    "TPBi":      "c1ccc(-n2c(-c3cc(-c4nc5ccccc5n4-c4ccccc4)cc(-c4nc5ccccc5n4-c4ccccc4)c3)nc3ccccc32)cc1",
    "pbPPhenB":  "c1ccc(-c2cc(-c3ccc(-c4cc(-c5ccccc5)nc5c4ccc4cccnc45)cc3)c3ccc4cccnc4c3n2)cc1",
}
EXPECTED = {"HATCN": "C18N12", "TPBi": "C45H30N6", "pbPPhenB": "C42H26N4"}

def write_xyz(path, symbols, coords, comment=""):
    with open(path, "w") as f:
        f.write(f"{len(symbols)}\n{comment}\n")
        for s, c in zip(symbols, coords):
            f.write(f"{s} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n")

def mol_to_arrays(mol):
    conf = mol.GetConformer()
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                     conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
    return syms, xyz

# ---------- organic molecules ----------
for name, smi in SMILES.items():
    mol = Chem.MolFromSmiles(smi)
    assert mol is not None, f"bad SMILES {name}"
    formula = CalcMolFormula(mol)
    print(name, "formula:", formula, "(expected", EXPECTED[name] + ")")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    syms, xyz = mol_to_arrays(mol)
    xyz -= xyz.mean(axis=0)
    write_xyz(os.path.join(OUT, f"{name}.xyz"), syms, xyz, name)

    # Ag placements
    if name == "HATCN":
        # (a) Ag above core center; core = ring N atoms (aromatic N)
        core_idx = [a.GetIdx() for a in mol.GetAtoms()
                    if a.GetSymbol() == "N" and a.GetIsAromatic()]
        core = xyz[core_idx]
        center = core.mean(axis=0)
        # normal of best-fit plane
        u, s, vt = np.linalg.svd(core - center)
        n = vt[2] / np.linalg.norm(vt[2])
        ag = center + 2.6 * n
        write_xyz(os.path.join(OUT, f"{name}_Ag_face.xyz"),
                  syms + ["Ag"], np.vstack([xyz, ag]), f"{name}+Ag face")
        # (b) Ag near a nitrile N, in plane, along C#N axis
        for a in mol.GetAtoms():
            if a.GetSymbol() == "N" and not a.GetIsAromatic():
                nidx = a.GetIdx()
                cidx = a.GetNeighbors()[0].GetIdx()
                break
        d = xyz[nidx] - xyz[cidx]; d /= np.linalg.norm(d)
        ag = xyz[nidx] + 2.3 * d
        write_xyz(os.path.join(OUT, f"{name}_Ag_CN.xyz"),
                  syms + ["Ag"], np.vstack([xyz, ag]), f"{name}+Ag nitrile")
    elif name == "pbPPhenB":
        # Ag chelated in one phenanthroline N,N pocket
        ri = mol.GetRingInfo()
        Ns = [a.GetIdx() for a in mol.GetAtoms()
              if a.GetSymbol() == "N" and a.GetIsAromatic()]
        # find N pair < 3.2 A apart (chelate pocket)
        best = None
        for i in range(len(Ns)):
            for j in range(i + 1, len(Ns)):
                r = np.linalg.norm(xyz[Ns[i]] - xyz[Ns[j]])
                if r < 3.4 and (best is None or r < best[0]):
                    best = (r, Ns[i], Ns[j])
        assert best, "no chelate pocket found"
        _, i1, i2 = best
        mid = 0.5 * (xyz[i1] + xyz[i2])
        away = mid - xyz.mean(axis=0); away /= np.linalg.norm(away)
        ag = mid + 1.9 * away
        write_xyz(os.path.join(OUT, f"{name}_Ag_chelate.xyz"),
                  syms + ["Ag"], np.vstack([xyz, ag]), f"{name}+Ag chelate")
    elif name == "TPBi":
        # Ag at a benzimidazole pyridine-type N (sp2 N with 2 neighbors)
        for a in mol.GetAtoms():
            if a.GetSymbol() == "N" and a.GetIsAromatic() and a.GetDegree() == 2:
                nidx = a.GetIdx(); break
        nb = [n.GetIdx() for n in mol.GetAtomWithIdx(nidx).GetNeighbors()]
        bis = xyz[nidx] - 0.5 * (xyz[nb[0]] + xyz[nb[1]])
        bis /= np.linalg.norm(bis)
        ag = xyz[nidx] + 2.3 * bis
        write_xyz(os.path.join(OUT, f"{name}_Ag_N.xyz"),
                  syms + ["Ag"], np.vstack([xyz, ag]), f"{name}+Ag N")

# ---------- LiF (001) cluster: 4x4x2 rock salt fragment, a = 4.026 A ----------
a = 4.026; d = a / 2
syms, xyz = [], []
for ix in range(4):
    for iy in range(4):
        for iz in range(2):
            for (ox, oy, oz, s) in [(0, 0, 0, "Li"), (1, 0, 0, "F")]:
                pass
# build proper rocksalt: site parity decides element
syms, xyz = [], []
for ix in range(4):
    for iy in range(4):
        for iz in range(2):
            s = "Li" if (ix + iy + iz) % 2 == 0 else "F"
            syms.append(s); xyz.append([ix * d, iy * d, iz * d])
xyz = np.array(xyz, float)
write_xyz(os.path.join(OUT, "LiF32.xyz"), syms, xyz, "LiF 4x4x2 cluster (bulk positions)")
# Ag on top-F site nearest the surface center (z max layer)
top = xyz[:, 2].max()
cands = [i for i in range(len(syms)) if syms[i] == "F" and abs(xyz[i, 2] - top) < 1e-3]
cxy = xyz[:, :2].mean(axis=0)
fi = min(cands, key=lambda i: np.linalg.norm(xyz[i, :2] - cxy))
ag = xyz[fi] + np.array([0, 0, 2.5])
write_xyz(os.path.join(OUT, "LiF32_Ag.xyz"), syms + ["Ag"], np.vstack([xyz, ag]), "LiF+Ag on F-top")
print("LiF cluster:", len(syms), "atoms; Ag over F index", fi)

# ---------- Mo3O9 cluster (standard MoO3 gas-phase model) ----------
# ring of 3 MoO2 units bridged by O; start from idealized geometry, will relax
mo_r = 1.9
sy, xy = [], []
for k in range(3):
    th = 2 * np.pi * k / 3
    mo = np.array([2.0 * np.cos(th), 2.0 * np.sin(th), 0.0])
    sy.append("Mo"); xy.append(mo)
    # two terminal O per Mo (up/down, outward)
    outw = mo / np.linalg.norm(mo)
    for zs in (+1, -1):
        sy.append("O"); xy.append(mo + 1.7 * (0.6 * outw + zs * 0.8 * np.array([0, 0, 1])))
    # bridging O between this Mo and next
    th2 = 2 * np.pi * (k + 0.5) / 3
    sy.append("O"); xy.append(np.array([2.3 * np.cos(th2), 2.3 * np.sin(th2), 0.0]))
xy = np.array(xy)
write_xyz(os.path.join(OUT, "Mo3O9.xyz"), sy, xy, "Mo3O9 idealized")
ag = np.array([0.0, 0.0, 2.3])  # above ring center
write_xyz(os.path.join(OUT, "Mo3O9_Ag.xyz"), sy + ["Ag"], np.vstack([xy, ag]), "Mo3O9+Ag top")
print("Mo3O9 built")

# ---------- Ag references ----------
write_xyz(os.path.join(OUT, "Ag1.xyz"), ["Ag"], np.zeros((1, 3)), "Ag atom")
write_xyz(os.path.join(OUT, "Ag2.xyz"), ["Ag", "Ag"],
          np.array([[0, 0, 0], [0, 0, 2.53]]), "Ag2 dimer")
print("done")
