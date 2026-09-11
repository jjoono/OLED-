"""Ag diffusion barrier on TAPC (TPA proxy), TPBi, B3PyMPM.
Same frozen-substrate lateral-drag method. Reference on-site energy computed with
robust SCF at the optimized geometry; bridge points scanned.
"""
import numpy as np, os, json, sys, psi4

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")
psi4.set_memory("5 GB"); psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS,"psi4_diff_org.out"), False)
H2EV=27.211386
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
    psi4.set_options(ROB)
    try:
        e=psi4.energy("pbe-d3bj",molecule=psi4.geometry(gstr(s,x,m)))
    except Exception:
        psi4.core.clean(); return None
    psi4.core.clean(); return e
def relaxed(sub_s,sub_x,ag_xy,normal,m,zs):
    best=None
    for dz in zs:
        e=SP(sub_s+["Ag"], np.vstack([sub_x, ag_xy+dz*normal]), m)
        if e is not None and (best is None or e<best): best=e
    return best

def onsite_energy(agfile, mult):
    s,x=rx(agfile)
    return SP(s,x,mult), s, x

# each: (tag, ag_geom_file, find-two-sites function)
def sites_TPBi(s,x,ag):
    # pyridine-type aromatic N deg2; pick nearest two to Ag
    Ns=[i for i,a in enumerate(s) if a=="N"]
    o=sorted(Ns,key=lambda i:np.linalg.norm(x[i]-ag))
    def site(i):
        nb=[j for j in range(len(s)) if np.linalg.norm(x[i]-x[j])<1.5 and j!=i]
        bis=x[i]-0.5*(x[nb[0]]+x[nb[1]]) if len(nb)>=2 else x[i]-x.mean(0)
        bis/=np.linalg.norm(bis); return x[i]+2.3*bis, bis
    return site(o[0]), site(o[1])
def sites_TPA(s,x,ag):
    # amine N + a ring; Ag sits over N/pi. two nearest heavy sites
    com=x.mean(0)
    N=[i for i,a in enumerate(s) if a=="N"][0]
    # site A: over N; site B: over an adjacent ring centroid
    outw=ag-x[N]; outw/=np.linalg.norm(outw) if np.linalg.norm(outw)>0 else 1
    A=(x[N]+2.5*outw, outw)
    # ring centroid ~ nearest 6 C to Ag
    C=[i for i,a in enumerate(s) if a=="C"]
    near=sorted(C,key=lambda i:np.linalg.norm(x[i]-ag))[:6]
    cen=x[near].mean(0); nB=cen-com; nB/=np.linalg.norm(nB)
    B=(cen+2.6*nB, nB)
    return A,B
def sites_B3(s,x,ag):
    Ns=[i for i,a in enumerate(s) if a=="N"]
    o=sorted(Ns,key=lambda i:np.linalg.norm(x[i]-ag))
    def site(i):
        nb=[j for j in range(len(s)) if np.linalg.norm(x[i]-x[j])<1.5 and j!=i]
        bis=x[i]-0.5*(x[nb[0]]+x[nb[1]]) if len(nb)>=2 else x[i]-x.mean(0)
        bis/=np.linalg.norm(bis); return x[i]+2.3*bis, bis
    return site(o[0]), site(o[1])

CASES={
 "TAPC":("TPA_Ag", sites_TPA),
 "TPBi":("TPBi_Ag", sites_TPBi),
 "B3PyMPM":("B3PyMPM_Ag", sites_B3),
}
jp=os.path.join(RUNS,"diffusion_barrier_eV.json")
res=json.load(open(jp)) if os.path.exists(jp) else {}
only=sys.argv[1:] if len(sys.argv)>1 else list(CASES)

for tag in only:
    if tag in res: print(f"{tag} cached {res[tag]:.3f}"); continue
    run,sfun=CASES[tag]
    agfile=os.path.join(RUNS,run,"xtbopt.xyz")
    s,x=rx(agfile); sub_s=s[:-1]; sub_x=x[:-1]; ag=x[-1]
    e_site=SP(s,x,2)  # on-site complex, doublet
    if e_site is None:
        print(f"{tag}: on-site FAILED"); continue
    (A,nA),(B,nB)=sfun(sub_s,sub_x,ag)
    nrm=(nA+nB); nrm/=np.linalg.norm(nrm)
    Emax=None
    for t in (0.4,0.5,0.6):
        lat=(1-t)*A+t*B
        e=relaxed(sub_s,sub_x,lat,nrm,2,[-0.2,0.2,0.5])
        print(f"  {tag} t={t:.1f}: {'FAIL' if e is None else f'{e:.6f}'}",flush=True)
        if e is not None and (Emax is None or e>Emax): Emax=e
    if Emax is None:
        print(f"{tag}: bridge FAILED"); continue
    Eb=(Emax-e_site)*H2EV
    res[tag]=Eb
    print(f"{tag}: on-site={e_site:.6f}, barrier ~ {Eb:.3f} eV",flush=True)
    json.dump(res,open(jp,"w"),indent=2)

print(json.dumps(res,indent=2))
