"""Zn atom and ZnS molecule binding on HATCN (nitrile N site) — same protocol.
Zn: closed-shell singlet atom. ZnS molecule: singlet, Zn-end to N.
"""
import numpy as np, os, subprocess, shutil, re, json

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STR, RUNS = os.path.join(BASE, "structures"), os.path.join(BASE, "runs")

def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0]); syms, xyz = [], []
    for l in lines[2:2+n]:
        p = l.split(); syms.append(p[0]); xyz.append(list(map(float, p[1:4])))
    return syms, np.array(xyz)

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

# HATCN optimized geometry; find a nitrile N (1 bonded C) like before
syms, xyz = read_xyz(os.path.join(RUNS, "HATCN", "xtbopt.xyz"))
# nitrile N: N whose nearest neighbor C is at ~1.16 A
ni, ci = None, None
for i, s in enumerate(syms):
    if s != "N": continue
    d = np.linalg.norm(xyz - xyz[i], axis=1); d[i] = 9e9
    j = int(np.argmin(d))
    if d[j] < 1.25:  # triple bond C#N
        ni, ci = i, j; break
ax = xyz[ni] - xyz[ci]; ax /= np.linalg.norm(ax)

# Zn on nitrile N
zn = xyz[ni] + 2.2 * ax
write_xyz(os.path.join(STR, "HATCN_Zn.xyz"), syms + ["Zn"], np.vstack([xyz, zn]))
run_xtb(os.path.join(STR, "HATCN_Zn.xyz"), "HATCN_Zn", uhf=0)

# ZnS molecule, Zn-end toward N (linear N...Zn-S)
zn = xyz[ni] + 2.2 * ax
s_at = xyz[ni] + (2.2 + 2.1) * ax
write_xyz(os.path.join(STR, "HATCN_ZnS.xyz"), syms + ["Zn", "S"],
          np.vstack([xyz, zn, s_at]))
run_xtb(os.path.join(STR, "HATCN_ZnS.xyz"), "HATCN_ZnS", uhf=0)

# free ZnS molecule
write_xyz(os.path.join(STR, "ZnS.xyz"), ["Zn", "S"],
          np.array([[0,0,0],[0,0,2.1]]))
run_xtb(os.path.join(STR, "ZnS.xyz"), "ZnS", uhf=0)

# ---- DFT CP binding ----
import psi4
psi4.set_memory("5 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_zns.out"), False)
psi4.set_options({"basis": "def2-svp", "scf_type": "df", "reference": "uks",
                  "maxiter": 200, "guess": "sad"})
H2EV = 27.211386

def gstr(syms, xyz, ghost=None, mult=1):
    s = f"0 {mult}\n"
    for i, (sym, c) in enumerate(zip(syms, xyz)):
        t = f"Gh({sym})" if ghost and i in ghost else sym
        s += f"{t} {c[0]} {c[1]} {c[2]}\n"
    return s + "symmetry c1\nno_reorient\nno_com\n"

def read_run(tag):
    return read_xyz(os.path.join(RUNS, tag, "xtbopt.xyz"))

res = {}
# Zn on HATCN: fragments = HATCN (singlet) + Zn (singlet); complex singlet
s1, x1 = read_run("HATCN_Zn")
n = len(s1)
e_cx = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(s1, x1, None, 1)))
e_sub = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(s1, x1, {n-1}, 1)))
e_ad = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(s1, x1, set(range(n-1)), 1)))
res["Zn_on_HATCN"] = (e_sub + e_ad - e_cx) * H2EV
print(f"Zn on HATCN: {res['Zn_on_HATCN']:.3f} eV", flush=True)
psi4.core.clean()

# ZnS molecule on HATCN: fragments = HATCN + ZnS(molecule, singlet)
s2, x2 = read_run("HATCN_ZnS")
n = len(s2)
e_cx = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(s2, x2, None, 1)))
e_sub = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(s2, x2, {n-2, n-1}, 1)))
e_ad = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(s2, x2, set(range(n-2)), 1)))
res["ZnS_on_HATCN"] = (e_sub + e_ad - e_cx) * H2EV
print(f"ZnS(molecule) on HATCN: {res['ZnS_on_HATCN']:.3f} eV", flush=True)
psi4.core.clean()

json.dump(res, open(os.path.join(RUNS, "zns_binding_eV.json"), "w"), indent=2)
