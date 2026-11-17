"""Proxy for HATCN nitrile site: benzonitrile (PhCN) + ZnS molecule / Zn atom.
GFN2 opt -> PBE-D3BJ/def2-SVP CP binding (RKS singlets).
"""
import numpy as np, os, subprocess, shutil, re, json
from rdkit import Chem
from rdkit.Chem import AllChem

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

def run_xtb(xyz_path, tag):
    wd = os.path.join(RUNS, tag); os.makedirs(wd, exist_ok=True)
    shutil.copy(xyz_path, os.path.join(wd, "in.xyz"))
    with open(os.path.join(wd, "xtb.log"), "w") as log:
        subprocess.run(["xtb", "in.xyz", "--gfn", "2", "--chrg", "0",
                        "--uhf", "0", "--opt", "tight"],
                       cwd=wd, stdout=log, stderr=subprocess.STDOUT)

mol = Chem.AddHs(Chem.MolFromSmiles("N#Cc1ccccc1"))
AllChem.EmbedMolecule(mol, randomSeed=1); AllChem.MMFFOptimizeMolecule(mol)
conf = mol.GetConformer()
syms = [a.GetSymbol() for a in mol.GetAtoms()]
xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                 conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
write_xyz(os.path.join(STR, "PhCN.xyz"), syms, xyz)
run_xtb(os.path.join(STR, "PhCN.xyz"), "PhCN")

s0, x0 = read_xyz(os.path.join(RUNS, "PhCN", "xtbopt.xyz"))
ni = next(i for i, s in enumerate(s0) if s == "N")
d = np.linalg.norm(x0 - x0[ni], axis=1); d[ni] = 9e9
ci = int(np.argmin(d))
ax = x0[ni] - x0[ci]; ax /= np.linalg.norm(ax)

# PhCN + ZnS (N...Zn-S linear)
zn = x0[ni] + 2.2 * ax; s_at = x0[ni] + 4.3 * ax
write_xyz(os.path.join(STR, "PhCN_ZnS.xyz"), s0 + ["Zn", "S"], np.vstack([x0, zn, s_at]))
run_xtb(os.path.join(STR, "PhCN_ZnS.xyz"), "PhCN_ZnS")
# PhCN + Zn
write_xyz(os.path.join(STR, "PhCN_Zn.xyz"), s0 + ["Zn"], np.vstack([x0, zn]))
run_xtb(os.path.join(STR, "PhCN_Zn.xyz"), "PhCN_Zn")

import psi4
psi4.set_memory("5 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_phcn.out"), False)
psi4.set_options({"basis": "def2-svp", "scf_type": "df", "reference": "rks",
                  "maxiter": 300, "guess": "sad"})
H2EV = 27.211386

def gstr(syms, xyz, ghost=None):
    s = "0 1\n"
    for i, (sym, c) in enumerate(zip(syms, xyz)):
        t = f"Gh({sym})" if ghost and i in ghost else sym
        s += f"{t} {c[0]} {c[1]} {c[2]}\n"
    return s + "symmetry c1\nno_reorient\nno_com\n"

jp = os.path.join(RUNS, "zns_binding_eV.json")
res = json.load(open(jp)) if os.path.exists(jp) else {"Zn_on_HATCN": 0.226}
for tag, nad in [("PhCN_ZnS", 2), ("PhCN_Zn", 1)]:
    syms, xyz = read_xyz(os.path.join(RUNS, tag, "xtbopt.xyz"))
    n = len(syms)
    ad = set(range(n - nad, n))
    e_cx = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz))); psi4.core.clean()
    e_sub = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, ad))); psi4.core.clean()
    e_ad = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, set(range(n - nad))))); psi4.core.clean()
    eb = (e_sub + e_ad - e_cx) * H2EV
    res[tag] = eb
    print(f"{tag}: E_b(CP) = {eb:.3f} eV", flush=True)

json.dump(res, open(os.path.join(RUNS, "zns_binding_eV.json"), "w"), indent=2)
