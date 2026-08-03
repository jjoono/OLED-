"""Stoichiometric Al2O3: Al4O6 gas cluster (Al:O = 2:3) -> proper AlOx proxy.
Cross-check against Al-rich Al10O10 (0.91 eV, overbinding).
"""
import numpy as np, os, subprocess, shutil, re, json

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STR, RUNS = os.path.join(BASE, "structures"), os.path.join(BASE, "runs")

def read_xyz(p):
    L=open(p).read().strip().splitlines(); n=int(L[0]); s,x=[],[]
    for l in L[2:2+n]:
        q=l.split(); s.append(q[0]); x.append(list(map(float,q[1:4])))
    return s,np.array(x)
def write_xyz(p,s,x,c=""):
    with open(p,"w") as f:
        f.write(f"{len(s)}\n{c}\n")
        for a,b in zip(s,x): f.write(f"{a} {b[0]:.6f} {b[1]:.6f} {b[2]:.6f}\n")
def run_xtb(p,tag,uhf=0,fix=None):
    wd=os.path.join(RUNS,tag); os.makedirs(wd,exist_ok=True)
    shutil.copy(p,os.path.join(wd,"in.xyz"))
    args=["xtb","in.xyz","--gfn","2","--chrg","0","--uhf",str(uhf),"--opt","normal"]
    if fix:
        open(os.path.join(wd,"x.inp"),"w").write("$fix\n  atoms: "+",".join(map(str,fix))+"\n$end\n")
        args+=["--input","x.inp"]
    with open(os.path.join(wd,"xtb.log"),"w") as log:
        subprocess.run(args,cwd=wd,stdout=log,stderr=subprocess.STDOUT)
    m=re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)",open(os.path.join(wd,"xtb.log"),errors="ignore").read())
    print(tag,"xtb:",m[-1] if m else "FAILED",flush=True); return bool(m)

# Al4O6: two Al2O3 units, cage. Build with Al tetrahedral-ish, O bridging.
s=["Al","Al","Al","Al","O","O","O","O","O","O"]
x=np.array([
 [ 1.70, 1.00, 0.00],[-1.70, 1.00, 0.00],[ 0.00,-1.60, 1.10],[ 0.00,-1.60,-1.10],
 [ 0.00, 1.75, 0.00],   # bridging O between Al1,Al2
 [ 1.10,-0.40, 1.05],   # O between Al1,Al3
 [ 1.10,-0.40,-1.05],   # O between Al1,Al4
 [-1.10,-0.40, 1.05],   # O between Al2,Al3
 [-1.10,-0.40,-1.05],   # O between Al2,Al4
 [ 0.00,-2.90, 0.00]])  # terminal-ish O on Al3/Al4
write_xyz(os.path.join(STR,"Al4O6.xyz"),s,x,"Al4O6")
ok=run_xtb(os.path.join(STR,"Al4O6.xyz"),"Al4O6",uhf=0)
if not ok:
    print("Al4O6 relax failed"); raise SystemExit

s0,x0=read_xyz(os.path.join(RUNS,"Al4O6","xtbopt.xyz"))
# Ag on the most exposed O (max distance from centroid among O)
oi=[i for i,a in enumerate(s0) if a=="O"]
com=x0.mean(axis=0)
best_o=max(oi,key=lambda i: np.linalg.norm(x0[i]-com))
outw=x0[best_o]-com; outw/=np.linalg.norm(outw)
ag=x0[best_o]+2.2*outw
write_xyz(os.path.join(STR,"Al4O6_Ag.xyz"),s0+["Ag"],np.vstack([x0,ag]),"Al4O6+Ag")
run_xtb(os.path.join(STR,"Al4O6_Ag.xyz"),"Al4O6_Ag",uhf=1,fix=list(range(1,len(s0)+1)))

import psi4
psi4.set_memory("5 GB"); psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS,"psi4_al4o6.out"),False)
H2EV=27.211386
ROB={"basis":"def2-svp","scf_type":"df","reference":"uks","maxiter":400,"guess":"sad",
     "level_shift":1.0,"level_shift_cutoff":1e-3,"damping_percentage":15.0}
def gstr(s,x,gh=None,m=1):
    st=f"0 {m}\n"
    for i,(a,c) in enumerate(zip(s,x)):
        t=f"Gh({a})" if gh and i in gh else a
        st+=f"{t} {c[0]} {c[1]} {c[2]}\n"
    return st+"symmetry c1\nno_reorient\nno_com\n"
def en(s,x,gh,m):
    psi4.set_options(ROB); e=psi4.energy("pbe-d3bj",molecule=psi4.geometry(gstr(s,x,gh,m))); psi4.core.clean(); return e
s1,x1=read_xyz(os.path.join(RUNS,"Al4O6_Ag","xtbopt.xyz"))
agi=len(s1)-1
e_cx=en(s1,x1,None,2); e_sub=en(s1,x1,{agi},1); e_ag=en(s1,x1,set(range(agi)),2)
eb=(e_sub+e_ag-e_cx)*H2EV
print(f"Al4O6_Ag (stoichiometric AlOx): E_b(CP) = {eb:.3f} eV",flush=True)
j=os.path.join(RUNS,"al2o3_binding_eV.json"); d=json.load(open(j)) if os.path.exists(j) else {}
d["Al4O6_Ag_stoich"]=eb; json.dump(d,open(j,"w"),indent=2)
