"""HTL site chemistry: Ag on triphenylamine (TAPC proxy), N-phenylcarbazole (TCTA proxy),
and benzene pi-face (generic hydrocarbon reference). Same protocol.
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

def build(smi, name, seed=3):
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    conf = mol.GetConformer()
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                     conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
    write_xyz(os.path.join(STR, f"{name}.xyz"), syms, xyz, name)
    run_xtb(os.path.join(STR, f"{name}.xyz"), name)
    return read_xyz(os.path.join(RUNS, name, "xtbopt.xyz"))

# --- triphenylamine: Ag above N along local C3 axis ---
s, x = build("c1ccc(N(c2ccccc2)c3ccccc3)cc1", "TPA")
ni = next(i for i, a in enumerate(s) if a == "N")
# neighbors: 3 C
d = np.linalg.norm(x - x[ni], axis=1); d[ni] = 9e9
nbrs = np.argsort(d)[:3]
normal = np.cross(x[nbrs[1]] - x[nbrs[0]], x[nbrs[2]] - x[nbrs[0]])
normal /= np.linalg.norm(normal)
ag = x[ni] + 2.5 * normal
write_xyz(os.path.join(STR, "TPA_Ag.xyz"), s + ["Ag"], np.vstack([x, ag]))
run_xtb(os.path.join(STR, "TPA_Ag.xyz"), "TPA_Ag", uhf=1)

# --- N-phenylcarbazole: Ag above carbazole ring face ---
s2, x2 = build("c1ccc(-n2c3ccccc3c3ccccc32)cc1", "PhCz")
# carbazole plane: use the 8 aromatic C of the two fused rings + N
ni2 = next(i for i, a in enumerate(s2) if a == "N")
heavy = [i for i, a in enumerate(s2) if a != "H"]
c0 = x2[heavy].mean(axis=0)
u, sv, vt = np.linalg.svd(x2[heavy] - c0)
nrm = vt[2] / np.linalg.norm(vt[2])
ag2 = x2[ni2] + 3.0 * nrm
write_xyz(os.path.join(STR, "PhCz_Ag.xyz"), s2 + ["Ag"], np.vstack([x2, ag2]))
run_xtb(os.path.join(STR, "PhCz_Ag.xyz"), "PhCz_Ag", uhf=1)

# --- benzene pi reference ---
s3, x3 = build("c1ccccc1", "benzene")
c3 = x3.mean(axis=0)
heavy = [i for i, a in enumerate(s3) if a != "H"]
u, sv, vt = np.linalg.svd(x3[heavy] - x3[heavy].mean(axis=0))
nrm3 = vt[2] / np.linalg.norm(vt[2])
ag3 = x3[heavy].mean(axis=0) + 2.6 * nrm3
write_xyz(os.path.join(STR, "benzene_Ag.xyz"), s3 + ["Ag"], np.vstack([x3, ag3]))
run_xtb(os.path.join(STR, "benzene_Ag.xyz"), "benzene_Ag", uhf=1)

# --- DFT CP binding ---
import psi4
psi4.set_memory("5 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_htl.out"), False)
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
for tag in ["TPA_Ag", "PhCz_Ag", "benzene_Ag"]:
    syms, xyz = read_xyz(os.path.join(RUNS, tag, "xtbopt.xyz"))
    n = len(syms); agi = n - 1
    e_cx = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, None, 2))); psi4.core.clean()
    e_sub = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, {agi}, 1))); psi4.core.clean()
    e_ag = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, set(range(agi)), 2))); psi4.core.clean()
    eb = (e_sub + e_ag - e_cx) * H2EV
    res[tag] = eb
    print(f"{tag}: E_b(CP) = {eb:.3f} eV", flush=True)

json.dump(res, open(os.path.join(RUNS, "htl_binding_eV.json"), "w"), indent=2)
