"""n-doping effect isolated: Ag binding to Bphen radical anion (chrg -1).
Complex [Bphen-Ag]^- singlet; fragments Bphen^- (doublet) + Ag (doublet).
Geometry: xtb opt at chrg -1 starting from neutral Bphen_Ag pocket geometry.
"""
import os, subprocess, shutil, re, json
import numpy as np

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")

wd = os.path.join(RUNS, "BphenAg_anion"); os.makedirs(wd, exist_ok=True)
shutil.copy(os.path.join(RUNS, "Bphen_Ag", "xtbopt.xyz"), os.path.join(wd, "in.xyz"))
with open(os.path.join(wd, "xtb.log"), "w") as log:
    subprocess.run(["xtb", "in.xyz", "--gfn", "2", "--chrg", "-1", "--uhf", "0",
                    "--opt", "tight"], cwd=wd, stdout=log, stderr=subprocess.STDOUT)
txt = open(os.path.join(wd, "xtb.log"), errors="ignore").read()
m = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", txt)
print("BphenAg_anion xtb:", m[-1] if m else "FAILED", flush=True)

import psi4
psi4.set_memory("5 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_bphen_anion.out"), False)
psi4.set_options({"basis": "def2-svp", "scf_type": "df", "reference": "uks",
                  "maxiter": 300, "guess": "sad"})
H2EV = 27.211386

def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0]); syms, xyz = [], []
    for l in lines[2:2+n]:
        p = l.split(); syms.append(p[0]); xyz.append(p[1:4])
    return syms, xyz

def gstr(syms, xyz, ghost=None, charge=0, mult=1):
    st = f"{charge} {mult}\n"
    for i, (sym, c) in enumerate(zip(syms, xyz)):
        t = f"Gh({sym})" if ghost and i in ghost else sym
        st += f"{t} {c[0]} {c[1]} {c[2]}\n"
    return st + "symmetry c1\nno_reorient\nno_com\n"

syms, xyz = read_xyz(os.path.join(wd, "xtbopt.xyz"))
n = len(syms); agi = n - 1
e_cx = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, None, -1, 1))); psi4.core.clean()
e_sub = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, {agi}, -1, 2))); psi4.core.clean()
e_ag = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, set(range(agi)), 0, 2))); psi4.core.clean()
eb = (e_sub + e_ag - e_cx) * H2EV
print(f"BphenAg_anion: E_b(CP) = {eb:.3f} eV", flush=True)

j = os.path.join(RUNS, "bphen_cs_binding_eV.json")
d = json.load(open(j)) if os.path.exists(j) else {}
d["Bphen_anion_Ag"] = eb
json.dump(d, open(j, "w"), indent=2)
