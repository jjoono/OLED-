"""B3PyMPM + Ag (pyridyl N site) and pyridine + Ag reference. Same protocol."""
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

def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0]); syms, xyz = [], []
    for l in lines[2:2+n]:
        p = l.split(); syms.append(p[0]); xyz.append(list(map(float, p[1:4])))
    return syms, np.array(xyz)

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

def build(smi, name, seed=5):
    mol0 = Chem.MolFromSmiles(smi)
    print(name, "formula:", CalcMolFormula(mol0), flush=True)
    mol = Chem.AddHs(mol0)
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=3000)
    conf = mol.GetConformer()
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                     conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
    write_xyz(os.path.join(STR, f"{name}.xyz"), syms, xyz, name)
    run_xtb(os.path.join(STR, f"{name}.xyz"), name)
    return mol, read_xyz(os.path.join(RUNS, name, "xtbopt.xyz"))

def place_ag_at_pyN(mol, syms, xyz, name, uhf=1):
    # pyridine-type N: aromatic, degree 2
    for a in mol.GetAtoms():
        if a.GetSymbol() == "N" and a.GetIsAromatic() and a.GetDegree() == 2:
            ni = a.GetIdx(); nb = [n.GetIdx() for n in a.GetNeighbors()]; break
    bis = xyz[ni] - 0.5 * (xyz[nb[0]] + xyz[nb[1]])
    bis /= np.linalg.norm(bis)
    ag = xyz[ni] + 2.3 * bis
    write_xyz(os.path.join(STR, f"{name}_Ag.xyz"), syms + ["Ag"], np.vstack([xyz, ag]))
    run_xtb(os.path.join(STR, f"{name}_Ag.xyz"), f"{name}_Ag", uhf=uhf)

# B3PyMPM
smi = "Cc1nc(-c2cc(-c3cccnc3)cc(-c3cccnc3)c2)cc(-c2cc(-c3cccnc3)cc(-c3cccnc3)c2)n1"
mol, (s0, x0) = build(smi, "B3PyMPM")
place_ag_at_pyN(mol, s0, x0, "B3PyMPM")

# pyridine reference
molp, (s1, x1) = build("c1ccncc1", "pyridine")
place_ag_at_pyN(molp, s1, x1, "pyridine")

# DFT CP
import psi4
psi4.set_memory("6 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_b3pympm.out"), False)
psi4.set_options({"basis": "def2-svp", "scf_type": "df", "reference": "uks",
                  "maxiter": 300, "guess": "sad"})
H2EV = 27.211386

def gstr(syms, xyz, ghost=None, mult=1):
    st = f"0 {mult}\n"
    for i, (sym, c) in enumerate(zip(syms, xyz)):
        t = f"Gh({sym})" if ghost and i in ghost else sym
        st += f"{t} {c[0]} {c[1]} {c[2]}\n"
    return st + "symmetry c1\nno_reorient\nno_com\n"

res = {}
for tag in ["pyridine_Ag", "B3PyMPM_Ag"]:
    syms, xyz = read_xyz(os.path.join(RUNS, tag, "xtbopt.xyz"))
    n = len(syms); agi = n - 1
    e_cx = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, None, 2))); psi4.core.clean()
    e_sub = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, {agi}, 1))); psi4.core.clean()
    e_ag = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, set(range(agi)), 2))); psi4.core.clean()
    res[tag] = (e_sub + e_ag - e_cx) * H2EV
    print(f"{tag}: E_b(CP) = {res[tag]:.3f} eV", flush=True)
    json.dump(res, open(os.path.join(RUNS, "b3pympm_binding_eV.json"), "w"), indent=2)
