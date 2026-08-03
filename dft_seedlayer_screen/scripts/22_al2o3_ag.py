"""Ag adatom binding on Al2O3 surface model + clean Al(111) cluster reference.
Al2O3: (Al2O3)4 cluster (relaxed) as amorphous native-oxide proxy, Ag on top-O.
Al(111): small Al13 cuboctahedral-ish cluster, Ag on surface hollow.
Same protocol: GFN2 geom (Al2O3, Al cluster fixed for Ag) -> PBE-D3BJ/def2-SVP CP.
"""
import numpy as np, os, subprocess, shutil, re, json

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

def run_xtb(xyz_path, tag, uhf=0, fix=None):
    wd = os.path.join(RUNS, tag); os.makedirs(wd, exist_ok=True)
    shutil.copy(xyz_path, os.path.join(wd, "in.xyz"))
    args = ["xtb", "in.xyz", "--gfn", "2", "--chrg", "0", "--uhf", str(uhf), "--opt", "tight"]
    if fix:
        open(os.path.join(wd, "x.inp"), "w").write("$fix\n  atoms: " + ",".join(map(str, fix)) + "\n$end\n")
        args += ["--input", "x.inp"]
    with open(os.path.join(wd, "xtb.log"), "w") as log:
        subprocess.run(args, cwd=wd, stdout=log, stderr=subprocess.STDOUT)
    txt = open(os.path.join(wd, "xtb.log"), errors="ignore").read()
    m = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", txt)
    print(tag, "xtb:", m[-1] if m else "FAILED", flush=True)
    return m[-1] if m else None

# --- (Al2O3)4 cluster: start from stoichiometric arrangement, relax with GFN2 ---
# build a rough Al8O12 cluster: alternate Al and O on a distorted grid
rng = np.random.default_rng(3)
al = np.array([[0,0,0],[2.8,0,0],[1.4,2.4,0],[1.4,0.8,2.3],
               [0,2.4,2.3],[2.8,2.4,2.3],[1.4,4.0,1.1],[4.2,1.2,1.1]], float)
o = np.array([[1.4,0,0],[0.7,1.2,0],[2.1,1.2,0],[0.7,1.2,2.3],[2.1,1.2,2.3],
              [1.4,2.4,1.1],[0,1.2,1.1],[2.8,1.2,1.1],[1.4,0.8,0.9],
              [1.4,3.2,2.0],[3.5,0.6,0.6],[3.5,1.8,1.6]], float)
syms = ["Al"]*8 + ["O"]*12
xyz = np.vstack([al, o])
write_xyz(os.path.join(STR, "Al2O3cl.xyz"), syms, xyz, "Al8O12 start")
run_xtb(os.path.join(STR, "Al2O3cl.xyz"), "Al2O3cl", uhf=0)

s0, x0 = read_xyz(os.path.join(RUNS, "Al2O3cl", "xtbopt.xyz"))
# Ag on a surface O with lowest coordination (fewest Al within 2.2 A), placed outward
oi = [i for i, s in enumerate(s0) if s == "O"]
def coord(i):
    return sum(1 for j, s in enumerate(s0) if s == "Al" and np.linalg.norm(x0[i]-x0[j]) < 2.2)
best_o = min(oi, key=coord)
com = x0.mean(axis=0)
outw = x0[best_o] - com; outw /= np.linalg.norm(outw)
ag = x0[best_o] + 2.2 * outw
write_xyz(os.path.join(STR, "Al2O3cl_Ag.xyz"), s0 + ["Ag"], np.vstack([x0, ag]), "Al2O3+Ag")
run_xtb(os.path.join(STR, "Al2O3cl_Ag.xyz"), "Al2O3cl_Ag", uhf=1, fix=list(range(1, 21)))

# --- clean Al cluster (Al13) reference ---
# icosahedral-ish Al13: center + 12 around at ~2.86 A
c = np.array([0,0,0.])
shell = []
phi = (1+5**0.5)/2
verts = []
for a in [-1,1]:
    for b in [-1,1]:
        verts += [(0,a,b*phi),(a,b*phi,0),(b*phi,0,a)]
verts = np.array(verts); verts = verts/np.linalg.norm(verts[0])*2.86
al13 = np.vstack([c, verts])
write_xyz(os.path.join(STR, "Al13.xyz"), ["Al"]*13, al13, "Al13")
run_xtb(os.path.join(STR, "Al13.xyz"), "Al13", uhf=1)
s1, x1 = read_xyz(os.path.join(RUNS, "Al13", "xtbopt.xyz"))
top = x1[np.argmax(x1[:,2])]
ag1 = top + np.array([0,0,2.6])
write_xyz(os.path.join(STR, "Al13_Ag.xyz"), ["Al"]*13+["Ag"], np.vstack([x1, ag1]), "Al13+Ag")
run_xtb(os.path.join(STR, "Al13_Ag.xyz"), "Al13_Ag", uhf=0, fix=list(range(1,14)))

# --- DFT CP ---
import psi4
psi4.set_memory("5 GB"); psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_al2o3.out"), False)
H2EV = 27.211386
ROB = {"basis":"def2-svp","scf_type":"df","reference":"uks","maxiter":400,"guess":"sad",
       "level_shift":1.0,"level_shift_cutoff":1e-3,"damping_percentage":15.0}

def gstr(syms, xyz, ghost=None, mult=1):
    st=f"0 {mult}\n"
    for i,(s,c) in enumerate(zip(syms,xyz)):
        t=f"Gh({s})" if ghost and i in ghost else s
        st+=f"{t} {c[0]} {c[1]} {c[2]}\n"
    return st+"symmetry c1\nno_reorient\nno_com\n"

def en(syms,xyz,ghost,mult):
    psi4.set_options(ROB)
    e=psi4.energy("pbe-d3bj",molecule=psi4.geometry(gstr(syms,xyz,ghost,mult))); psi4.core.clean(); return e

res={}
for tag,sub_mult in [("Al2O3cl_Ag",1),("Al13_Ag",2)]:
    syms,xyz=read_xyz(os.path.join(RUNS,tag,"xtbopt.xyz"))
    agi=len(syms)-1
    try:
        e_cx=en(syms,xyz,None,2)
        e_sub=en(syms,xyz,{agi},sub_mult)
        e_ag=en(syms,xyz,set(range(agi)),2)
        res[tag]=(e_sub+e_ag-e_cx)*H2EV
        print(f"{tag}: E_b(CP) = {res[tag]:.3f} eV",flush=True)
    except Exception as ex:
        print(f"{tag} FAILED: {ex}",flush=True)
    json.dump(res,open(os.path.join(RUNS,"al2o3_binding_eV.json"),"w"),indent=2)
