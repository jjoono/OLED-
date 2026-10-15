"""Bphen, Cs2CO3, Bphen:Cs (n-doped model) as Ag seeds. Same protocol.
- Bphen_Ag: Ag in phen N,N pocket (doublet)
- Cs2CO3_Ag: Ag on carbonate O (doublet)
- BphenCs_Ag: Cs in pocket (charge-transfer complex), Ag on phen pi face (singlet)
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

# --- Bphen: 4,7-diphenyl-1,10-phenanthroline ---
smi = "c1ccc(-c2ccnc3c2ccc2c(-c4ccccc4)ccnc23)cc1"
mol = Chem.AddHs(Chem.MolFromSmiles(smi))
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
print("Bphen formula:", CalcMolFormula(Chem.MolFromSmiles(smi)), "(expect C24H16N2)", flush=True)
AllChem.EmbedMolecule(mol, randomSeed=11)
AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
conf = mol.GetConformer()
syms = [a.GetSymbol() for a in mol.GetAtoms()]
xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                 conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
write_xyz(os.path.join(STR, "Bphen.xyz"), syms, xyz)
run_xtb(os.path.join(STR, "Bphen.xyz"), "Bphen")
s0, x0 = read_xyz(os.path.join(RUNS, "Bphen", "xtbopt.xyz"))

Ns = [i for i, a in enumerate(s0) if a == "N"]
assert len(Ns) == 2
mid = 0.5 * (x0[Ns[0]] + x0[Ns[1]])
away = mid - x0.mean(axis=0); away /= np.linalg.norm(away)

# Bphen + Ag (pocket)
ag = mid + 1.9 * away
write_xyz(os.path.join(STR, "Bphen_Ag.xyz"), s0 + ["Ag"], np.vstack([x0, ag]))
run_xtb(os.path.join(STR, "Bphen_Ag.xyz"), "Bphen_Ag", uhf=1)

# Bphen + Cs (pocket, CT complex) -> then Ag on phen pi face
cs = mid + 2.9 * away
write_xyz(os.path.join(STR, "BphenCs.xyz"), s0 + ["Cs"], np.vstack([x0, cs]))
run_xtb(os.path.join(STR, "BphenCs.xyz"), "BphenCs", uhf=1)
s1, x1 = read_xyz(os.path.join(RUNS, "BphenCs", "xtbopt.xyz"))
# phen core plane normal (use ring N + nearby aromatic C within 3 A of N midpoint)
midc = 0.5 * (x1[Ns[0]] + x1[Ns[1]])
core = [i for i, a in enumerate(s1) if a in ("C", "N")
        and np.linalg.norm(x1[i] - midc) < 4.0]
c0 = x1[core].mean(axis=0)
u, sv, vt = np.linalg.svd(x1[core] - c0)
nrm = vt[2] / np.linalg.norm(vt[2])
ag1 = c0 + 2.8 * nrm
write_xyz(os.path.join(STR, "BphenCs_Ag.xyz"), s1 + ["Ag"], np.vstack([x1, ag1]))
run_xtb(os.path.join(STR, "BphenCs_Ag.xyz"), "BphenCs_Ag", uhf=0)  # radical pair -> singlet

# --- Cs2CO3 monomer + Ag ---
# rough start: planar CO3 with two Cs above O-O bridges
co3 = np.array([[0,0,0],[1.29,0,0],[-0.645,1.117,0],[-0.645,-1.117,0]])  # C,O,O,O
cs1 = np.array([2.6, 1.6, 0.8]); cs2 = np.array([-2.6, 0.0, 0.8])
sy = ["C","O","O","O","Cs","Cs"]
xy = np.vstack([co3, cs1, cs2])
write_xyz(os.path.join(STR, "Cs2CO3.xyz"), sy, xy)
run_xtb(os.path.join(STR, "Cs2CO3.xyz"), "Cs2CO3")
s2, x2 = read_xyz(os.path.join(RUNS, "Cs2CO3", "xtbopt.xyz"))
# Ag near an O away from both Cs
oi = [i for i, a in enumerate(s2) if a == "O"]
csi = [i for i, a in enumerate(s2) if a == "Cs"]
best_o = max(oi, key=lambda i: min(np.linalg.norm(x2[i]-x2[j]) for j in csi))
away2 = x2[best_o] - x2[[i for i,a in enumerate(s2) if a=="C"][0]]
away2 /= np.linalg.norm(away2)
ag2 = x2[best_o] + 2.4 * away2
write_xyz(os.path.join(STR, "Cs2CO3_Ag.xyz"), s2 + ["Ag"], np.vstack([x2, ag2]))
run_xtb(os.path.join(STR, "Cs2CO3_Ag.xyz"), "Cs2CO3_Ag", uhf=1)

# --- DFT CP binding ---
import psi4
psi4.set_memory("5 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_bphen.out"), False)
H2EV = 27.211386

def gstr(syms, xyz, ghost=None, mult=1):
    st = f"0 {mult}\n"
    for i, (sym, c) in enumerate(zip(syms, xyz)):
        t = f"Gh({sym})" if ghost and i in ghost else sym
        st += f"{t} {c[0]} {c[1]} {c[2]}\n"
    return st + "symmetry c1\nno_reorient\nno_com\n"

res = {}
def cp_binding(tag, cx_mult, sub_mult, ad_mult, extra=None):
    syms, xyz = read_xyz(os.path.join(RUNS, tag, "xtbopt.xyz"))
    n = len(syms); agi = n - 1
    opts = {"basis": "def2-svp", "scf_type": "df", "reference": "uks",
            "maxiter": 300, "guess": "sad"}
    if extra: opts.update(extra)
    psi4.set_options(opts)
    e_cx = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, None, cx_mult))); psi4.core.clean()
    e_sub = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, {agi}, sub_mult))); psi4.core.clean()
    e_ag = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, set(range(agi)), ad_mult))); psi4.core.clean()
    eb = (e_sub + e_ag - e_cx) * H2EV
    res[tag] = eb
    print(f"{tag}: E_b(CP) = {eb:.3f} eV", flush=True)
    json.dump(res, open(os.path.join(RUNS, "bphen_cs_binding_eV.json"), "w"), indent=2)

cp_binding("Bphen_Ag", 2, 1, 2)
cp_binding("Cs2CO3_Ag", 2, 1, 2)
cp_binding("BphenCs_Ag", 1, 2, 2)  # complex singlet; sub = BphenCs doublet; Ag doublet
