"""Clean DFT diffusion barrier for TPBi and B3PyMPM.
Reference: optimized on-site complex (robust SCF).
Bridge: Ag lifted above the midpoint of the two nearest aromatic N atoms, at the
same height Ag sits above its bound N. z-scan a few heights; reject SCF states that
deviate >5 eV from on-site (failed excited states). Barrier = min-valid-bridge - onsite.
"""
import numpy as np, os, json, psi4

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")
psi4.set_memory("5 GB"); psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS,"psi4_tpbi_b3.out"), False)
H2EV=27.211386
PLAIN={"basis":"def2-svp","scf_type":"df","reference":"uks","maxiter":300,"guess":"sad"}
ROB={"basis":"def2-svp","scf_type":"df","reference":"uks","maxiter":400,"guess":"sad",
     "level_shift":1.0,"level_shift_cutoff":1e-3,"damping_percentage":15.0}

def rx(p):
    L=open(p).read().strip().splitlines();n=int(L[0]);s=[];x=[]
    for l in L[2:2+n]:
        q=l.split();s.append(q[0]);x.append(list(map(float,q[1:4])))
    return s,np.array(x)
def gstr(s,x,m):
    return f"0 {m}\n"+"".join(f"{a} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n" for a,c in zip(s,x))+"symmetry c1\nno_reorient\nno_com\n"
def SP(s,x,m):
    for opts in (PLAIN, ROB):
        psi4.set_options(opts)
        try:
            e=psi4.energy("pbe-d3bj",molecule=psi4.geometry(gstr(s,x,m)))
            psi4.core.clean(); return e
        except Exception:
            psi4.core.clean(); continue
    return None

def run(tag, agfile):
    s,x=rx(agfile); sub_s,sub_x,ag=s[:-1],x[:-1],x[-1]
    e_site=SP(s,x,2)
    if e_site is None: print(f"{tag}: on-site FAIL"); return None
    # bound N = aromatic N nearest Ag; neighbor N = next nearest N
    Ns=[i for i,a in enumerate(sub_s) if a=="N"]
    o=sorted(Ns,key=lambda i:np.linalg.norm(sub_x[i]-ag))
    Nb, Nn = o[0], o[1]
    hAg=np.linalg.norm(ag-sub_x[Nb])   # Ag height above bound N
    mid=0.5*(sub_x[Nb]+sub_x[Nn])
    outw=ag-sub_x[Nb]; outw/=np.linalg.norm(outw)  # approx surface-outward dir
    valids=[]
    for h in (hAg-0.3, hAg, hAg+0.3, hAg+0.6):
        e=SP(sub_s+["Ag"], np.vstack([sub_x, mid+h*outw]), 2)
        if e is None: continue
        if (e - e_site)*H2EV > 5.0: continue   # reject failed excited state
        valids.append(e)
        print(f"  {tag} bridge h={h:.2f}: {e:.6f} (rel {(e-e_site)*H2EV:.3f} eV)",flush=True)
    if not valids: print(f"{tag}: bridge FAIL"); return None
    Eb=(min(valids)-e_site)*H2EV
    print(f"{tag}: on-site={e_site:.6f}, diffusion barrier ~ {Eb:.3f} eV",flush=True)
    return Eb

jp=os.path.join(RUNS,"diffusion_barrier_eV.json")
res=json.load(open(jp)) if os.path.exists(jp) else {}
res.pop("B3PyMPM",None)  # remove garbage
for tag,run_dir in [("TPBi","TPBi_Ag"),("B3PyMPM","B3PyMPM_Ag")]:
    eb=run(tag, os.path.join(RUNS,run_dir,"xtbopt.xyz"))
    if eb is not None: res[tag+"_diff"]=eb
    json.dump(res,open(jp,"w"),indent=2)
print(json.dumps(res,indent=2))
