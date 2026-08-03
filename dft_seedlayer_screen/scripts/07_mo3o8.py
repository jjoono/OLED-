"""Reduced MoOx model: remove one terminal O from Mo3O9 -> Mo3O8 (O vacancy site),
GFN2 reopt of bare cluster, place Ag at vacancy, GFN2 opt, write structures.
"""
import numpy as np, os, subprocess, shutil, re, sys

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")
STR = os.path.join(BASE, "structures")

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

def run_xtb(xyz_path, tag, uhf=0, opt=True):
    wd = os.path.join(RUNS, tag); os.makedirs(wd, exist_ok=True)
    shutil.copy(xyz_path, os.path.join(wd, "in.xyz"))
    args = ["xtb", "in.xyz", "--gfn", "2", "--chrg", "0", "--uhf", str(uhf)]
    if opt: args += ["--opt", "tight"]
    with open(os.path.join(wd, "xtb.log"), "w") as log:
        subprocess.run(args, cwd=wd, stdout=log, stderr=subprocess.STDOUT)
    txt = open(os.path.join(wd, "xtb.log"), errors="ignore").read()
    m = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", txt)
    print(tag, m[-1] if m else "FAILED", flush=True)

syms, xyz = read_xyz(os.path.join(RUNS, "Mo3O9", "xtbopt.xyz"))
# remove first terminal O (index 1: attached to Mo0, z>0)
del_i = 1
s2 = [s for i, s in enumerate(syms) if i != del_i]
x2 = np.delete(xyz, del_i, axis=0)
write_xyz(os.path.join(STR, "Mo3O8.xyz"), s2, x2, "Mo3O8 (O vacancy)")
run_xtb(os.path.join(STR, "Mo3O8.xyz"), "Mo3O8", uhf=0, opt=True)  # closed-shell Mo(IV) site

# Ag near the vacancy site (where O was)
s3, x3 = read_xyz(os.path.join(RUNS, "Mo3O8", "xtbopt.xyz"))
mo0 = x3[0]
# direction of removed O: up/outward from Mo0
d = np.array([xyz[del_i][0]-xyz[0][0], xyz[del_i][1]-xyz[0][1], xyz[del_i][2]-xyz[0][2]])
d /= np.linalg.norm(d)
ag = mo0 + 2.6 * d
write_xyz(os.path.join(STR, "Mo3O8_Ag.xyz"), s3 + ["Ag"], np.vstack([x3, ag]), "Mo3O8+Ag at vacancy")
run_xtb(os.path.join(STR, "Mo3O8_Ag.xyz"), "Mo3O8_Ag", uhf=1, opt=True)
print("done")
