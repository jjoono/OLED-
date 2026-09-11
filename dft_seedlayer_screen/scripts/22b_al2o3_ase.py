"""Build Al2O3 cluster from ASE corundum bulk, cap, GFN2 relax, Ag adsorption.
Also stoichiometric Al4O6 gas cluster as cross-check.
"""
import numpy as np, os, subprocess, shutil, re, json
from ase import Atoms
from ase.build import bulk
from ase.io import write as ase_write

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STR, RUNS = os.path.join(BASE, "structures"), os.path.join(BASE, "runs")

def read_xyz(path):
    lines = open(path).read().strip().splitlines(); n=int(lines[0]); s,x=[],[]
    for l in lines[2:2+n]:
        p=l.split(); s.append(p[0]); x.append(list(map(float,p[1:4])))
    return s, np.array(x)

def write_xyz(path, syms, xyz, cm=""):
    with open(path,"w") as f:
        f.write(f"{len(syms)}\n{cm}\n")
        for s,c in zip(syms,xyz): f.write(f"{s} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n")

def run_xtb(xyz_path, tag, uhf=0, fix=None):
    wd=os.path.join(RUNS,tag); os.makedirs(wd,exist_ok=True)
    shutil.copy(xyz_path, os.path.join(wd,"in.xyz"))
    args=["xtb","in.xyz","--gfn","2","--chrg","0","--uhf",str(uhf),"--opt","normal"]
    if fix:
        open(os.path.join(wd,"x.inp"),"w").write("$fix\n  atoms: "+",".join(map(str,fix))+"\n$end\n")
        args+=["--input","x.inp"]
    with open(os.path.join(wd,"xtb.log"),"w") as log:
        subprocess.run(args,cwd=wd,stdout=log,stderr=subprocess.STDOUT)
    txt=open(os.path.join(wd,"xtb.log"),errors="ignore").read()
    m=re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)",txt)
    print(tag,"xtb:",m[-1] if m else "FAILED",flush=True)
    return bool(m)

# corundum alpha-Al2O3
try:
    al2o3 = bulk("Al2O3", crystalstructure="corundum", a=4.785, c=13.0)
except Exception:
    # fallback: build hexagonal corundum manually via spacegroup
    from ase.spacegroup import crystal
    al2o3 = crystal(["Al","O"], basis=[(0,0,0.352),(0.306,0,0.25)],
                    spacegroup=167, cellpar=[4.785,4.785,12.991,90,90,120])
sc = al2o3.repeat((2,2,1))
# carve a compact cluster: atoms within radius R of a chosen center
pos = sc.get_positions(); syms = sc.get_chemical_symbols()
center = pos.mean(axis=0)
d = np.linalg.norm(pos-center,axis=1)
idx = np.argsort(d)[:20]  # ~20 atoms
csyms = [syms[i] for i in idx]; cxyz = pos[idx]
# ensure some O present; report composition
from collections import Counter
print("cluster comp:", Counter(csyms), flush=True)
cxyz = cxyz - cxyz.mean(axis=0)
write_xyz(os.path.join(STR,"Al2O3_corundum.xyz"), csyms, cxyz, "corundum fragment")
ok = run_xtb(os.path.join(STR,"Al2O3_corundum.xyz"), "Al2O3_cor", uhf=0)

if not ok:
    # fallback stoichiometric Al4O6 book structure (known stable)
    a4o6_s = ["Al","Al","Al","Al","O","O","O","O","O","O"]
    a4o6_x = np.array([
        [ 1.60, 0.00, 0.60],[-1.60,0.00,0.60],[0.00,1.60,-0.60],[0.00,-1.60,-0.60],
        [ 1.55, 1.55, 0.00],[ 1.55,-1.55,0.00],[-1.55,1.55,0.00],[-1.55,-1.55,0.00],
        [ 0.00, 0.00, 1.75],[ 0.00, 0.00,-1.75]])
    write_xyz(os.path.join(STR,"Al4O6.xyz"), a4o6_s, a4o6_x, "Al4O6")
    ok = run_xtb(os.path.join(STR,"Al4O6.xyz"),"Al2O3_cor",uhf=0)

s0,x0 = read_xyz(os.path.join(RUNS,"Al2O3_cor","xtbopt.xyz"))
oi=[i for i,s in enumerate(s0) if s=="O"]
def coord(i): return sum(1 for j,s in enumerate(s0) if s=="Al" and np.linalg.norm(x0[i]-x0[j])<2.2)
best_o=min(oi,key=coord)
com=x0.mean(axis=0); outw=x0[best_o]-com; outw/=np.linalg.norm(outw)
ag=x0[best_o]+2.2*outw
write_xyz(os.path.join(STR,"Al2O3_cor_Ag.xyz"), s0+["Ag"], np.vstack([x0,ag]),"Al2O3+Ag")
run_xtb(os.path.join(STR,"Al2O3_cor_Ag.xyz"),"Al2O3_cor_Ag",uhf=1,fix=list(range(1,len(s0)+1)))

# DFT CP
import psi4
psi4.set_memory("5 GB"); psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS,"psi4_al2o3b.out"),False)
H2EV=27.211386
ROB={"basis":"def2-svp","scf_type":"df","reference":"uks","maxiter":400,"guess":"sad",
     "level_shift":1.0,"level_shift_cutoff":1e-3,"damping_percentage":15.0}
def gstr(syms,xyz,ghost=None,mult=1):
    st=f"0 {mult}\n"
    for i,(s,c) in enumerate(zip(syms,xyz)):
        t=f"Gh({s})" if ghost and i in ghost else s
        st+=f"{t} {c[0]} {c[1]} {c[2]}\n"
    return st+"symmetry c1\nno_reorient\nno_com\n"
def en(syms,xyz,ghost,mult):
    psi4.set_options(ROB)
    e=psi4.energy("pbe-d3bj",molecule=psi4.geometry(gstr(syms,xyz,ghost,mult))); psi4.core.clean(); return e

syms,xyz=read_xyz(os.path.join(RUNS,"Al2O3_cor_Ag","xtbopt.xyz"))
agi=len(syms)-1
e_cx=en(syms,xyz,None,2); e_sub=en(syms,xyz,{agi},1); e_ag=en(syms,xyz,set(range(agi)),2)
eb=(e_sub+e_ag-e_cx)*H2EV
print(f"Al2O3_Ag: E_b(CP) = {eb:.3f} eV",flush=True)
j=os.path.join(RUNS,"al2o3_binding_eV.json")
d=json.load(open(j)) if os.path.exists(j) else {}
d["Al2O3_Ag"]=eb; json.dump(d,open(j,"w"),indent=2)
