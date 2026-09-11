"""Quick screen of new candidates: F4TCNQ (nitrile/CT anchor) and Liq (O,N chelate).
RDKit -> GFN2-xTB opt -> PBE-D3BJ/def2-SVP CP single points (same protocol).
"""
import numpy as np, os, subprocess, shutil, re, json
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STR, RUNS = os.path.join(BASE, "structures"), os.path.join(BASE, "runs")

def write_xyz(path, syms, xyz, cm=""):
    with open(path, "w") as f:
        f.write(f"{len(syms)}\n{cm}\n")
        for s, c in zip(syms, xyz):
            f.write(f"{s} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n")

def run_xtb(xyz_path, tag, uhf=0):
    wd = os.path.join(RUNS, tag); os.makedirs(wd, exist_ok=True)
    shutil.copy(xyz_path, os.path.join(wd, "in.xyz"))
    with open(os.path.join(wd, "xtb.log"), "w") as log:
        subprocess.run(["xtb", "in.xyz", "--gfn", "2", "--chrg", "0",
                        "--uhf", str(uhf), "--opt", "tight"],
                       cwd=wd, stdout=log, stderr=subprocess.STDOUT)
    txt = open(os.path.join(wd, "xtb.log"), errors="ignore").read()
    m = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", txt)
    print(tag, "xtb:", m[-1] if m else "FAILED", flush=True)

# ---- build ----
mols = {
    "F4TCNQ": "N#CC(C#N)=C1C(F)=C(F)C(=C(C#N)C#N)C(F)=C1F",
    "Liq": "[Li+].[O-]c1cccc2cccnc12",
}
for name, smi in mols.items():
    mol = Chem.MolFromSmiles(smi)
    print(name, CalcMolFormula(mol), flush=True)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000) if name != "Liq" else None
    conf = mol.GetConformer()
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                     conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
    write_xyz(os.path.join(STR, f"{name}.xyz"), syms, xyz, name)
    run_xtb(os.path.join(STR, f"{name}.xyz"), name)
    # place Ag
    s2, x2 = [], []
    lines = open(os.path.join(RUNS, name, "xtbopt.xyz")).read().strip().splitlines()
    for l in lines[2:2 + int(lines[0])]:
        p = l.split(); s2.append(p[0]); x2.append(list(map(float, p[1:4])))
    x2 = np.array(x2)
    if name == "F4TCNQ":
        # Ag at a nitrile N along C#N axis (find N with 1 neighbor C at ~1.16 A)
        ni = next(i for i, s in enumerate(s2) if s == "N")
        d = np.linalg.norm(x2 - x2[ni], axis=1); d[ni] = 9e9
        ci = int(np.argmin(d))
        ax = x2[ni] - x2[ci]; ax /= np.linalg.norm(ax)
        ag = x2[ni] + 2.3 * ax
    else:
        # Liq: Ag near O and N side, in-plane outward from ring centroid
        oi = next(i for i, s in enumerate(s2) if s == "O")
        ni = next(i for i, s in enumerate(s2) if s == "N")
        mid = 0.5 * (x2[oi] + x2[ni])
        away = mid - x2.mean(axis=0); away /= np.linalg.norm(away)
        ag = mid + 2.0 * away
    write_xyz(os.path.join(STR, f"{name}_Ag.xyz"), s2 + ["Ag"], np.vstack([x2, ag]))
    run_xtb(os.path.join(STR, f"{name}_Ag.xyz"), f"{name}_Ag", uhf=1)

# ---- DFT CP single points ----
import psi4
psi4.set_memory("5 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_newcand.out"), False)
psi4.set_options({"basis": "def2-svp", "scf_type": "df", "reference": "uks",
                  "maxiter": 200, "guess": "sad"})
H2EV = 27.211386

def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0]); syms, xyz = [], []
    for l in lines[2:2+n]:
        p = l.split(); syms.append(p[0]); xyz.append(p[1:4])
    return syms, xyz

def gstr(syms, xyz, ghost=None, mult=1):
    s = f"0 {mult}\n"
    for i, (sym, c) in enumerate(zip(syms, xyz)):
        t = f"Gh({sym})" if ghost and i in ghost else sym
        s += f"{t} {c[0]} {c[1]} {c[2]}\n"
    return s + "symmetry c1\nno_reorient\nno_com\n"

res = {}
for tag in ["F4TCNQ_Ag", "Liq_Ag"]:
    syms, xyz = read_xyz(os.path.join(RUNS, tag, "xtbopt.xyz"))
    agi = len(syms) - 1
    e_cx = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, None, 2)))
    e_sub = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, {agi}, 1)))
    e_ag = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, set(range(agi)), 2)))
    eb = (e_sub + e_ag - e_cx) * H2EV
    res[tag] = eb
    print(f"{tag}: E_b(CP) = {eb:.3f} eV", flush=True)
    psi4.core.clean()

json.dump(res, open(os.path.join(RUNS, "newcand_binding_eV.json"), "w"), indent=2)
