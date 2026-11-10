"""Faster diffusion-barrier estimate. Order: clean Al, AlOx (small, fast),
then HATCN with reduced sampling (npath=4, 3 z-points). Frozen substrate, DFT.
Barrier = max(E_path) - E(start site).
"""
import numpy as np, os, json, sys, psi4

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")
psi4.set_memory("5 GB"); psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS,"psi4_diff2.out"), False)
H2EV=27.211386
PLAIN={"basis":"def2-svp","scf_type":"df","reference":"uks","maxiter":250,"guess":"sad"}
ROB=dict(PLAIN, level_shift=1.0, level_shift_cutoff=1e-3, damping_percentage=15.0, maxiter=350)

def read_xyz(p):
    L=open(p).read().strip().splitlines();n=int(L[0]);s=[];x=[]
    for l in L[2:2+n]:
        q=l.split();s.append(q[0]);x.append(list(map(float,q[1:4])))
    return s,np.array(x)
def gstr(s,x,m):
    st=f"0 {m}\n"+ "".join(f"{a} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n" for a,c in zip(s,x))
    return st+"symmetry c1\nno_reorient\nno_com\n"
def SP(s,x,m,opts):
    psi4.set_options(opts)
    try:
        e=psi4.energy("pbe-d3bj",molecule=psi4.geometry(gstr(s,x,m)))
    except Exception:
        psi4.core.clean(); return None
    psi4.core.clean(); return e
def relaxed(sub_s,sub_x,ag_xy,normal,m,opts,zs):
    best=None
    for dz in zs:
        e=SP(sub_s+["Ag"], np.vstack([sub_x, ag_xy+dz*normal]), m, opts)
        if e is not None and (best is None or e<best): best=e
    return best
def barrier(tag, agfile, siteA, siteB, normal, mult, opts, npath, zs):
    s,x=read_xyz(agfile)
    if s[-1]=="Ag": s,x=s[:-1],x[:-1]
    E=[]
    for t in np.linspace(0,1,npath):
        lat=(1-t)*siteA+t*siteB
        e=relaxed(s,x,lat,normal,mult,opts,zs)
        E.append(e); print(f"  {tag} t={t:.2f}: {'FAIL' if e is None else f'{e:.6f}'}",flush=True)
    v=[e for e in E if e is not None]
    if len(v)<2 or E[0] is None:
        print(f"{tag}: barrier FAILED",flush=True); return None
    Eb=(max(v)-E[0])*H2EV
    print(f"{tag}: diffusion barrier ~ {Eb:.3f} eV",flush=True); return Eb

jp=os.path.join(RUNS,"diffusion_barrier_eV.json")
res=json.load(open(jp)) if os.path.exists(jp) else {}
only=sys.argv[1:] if len(sys.argv)>1 else ["cleanAl","AlOx","HATCN"]

# clean Al (Al13)
if "cleanAl" in only and "cleanAl" not in res:
    s,x=read_xyz(os.path.join(RUNS,"Al13_Ag","xtbopt.xyz")); sub_x=x[:-1]; ag=x[-1]; com=sub_x.mean(0)
    o=sorted(range(len(sub_x)),key=lambda i:np.linalg.norm(sub_x[i]-ag))
    def site(i,d=2.6):
        w=sub_x[i]-com;w/=np.linalg.norm(w);return sub_x[i]+d*w,w
    A,nA=site(o[0]); B,nB=site(o[1]); nrm=(nA+nB);nrm/=np.linalg.norm(nrm)
    res["cleanAl"]=barrier("cleanAl",os.path.join(RUNS,"Al13_Ag","xtbopt.xyz"),A,B,nrm,1,ROB,4,[-0.3,0,0.4])
    json.dump(res,open(jp,"w"),indent=2)

# AlOx (Al4O6)
if "AlOx" in only and "AlOx" not in res:
    s,x=read_xyz(os.path.join(RUNS,"Al4O6_Ag","xtbopt.xyz")); sub_s=s[:-1];sub_x=x[:-1];ag=x[-1];com=sub_x.mean(0)
    O=[i for i,a in enumerate(sub_s) if a=="O"]; o=sorted(O,key=lambda i:np.linalg.norm(sub_x[i]-ag))
    def osite(i,d=2.2):
        w=sub_x[i]-com;w/=np.linalg.norm(w);return sub_x[i]+d*w,w
    A,nA=osite(o[0]);B,nB=osite(o[1]);nrm=(nA+nB);nrm/=np.linalg.norm(nrm)
    res["AlOx"]=barrier("AlOx",os.path.join(RUNS,"Al4O6_Ag","xtbopt.xyz"),A,B,nrm,2,ROB,4,[-0.3,0,0.4])
    json.dump(res,open(jp,"w"),indent=2)

# HATCN reduced
if "HATCN" in only and "HATCN" not in res:
    s,x=read_xyz(os.path.join(RUNS,"HATCN_Ag_CN","xtbopt.xyz")); sub_s=s[:-1];sub_x=x[:-1];ag=x[-1]
    N=[i for i,a in enumerate(sub_s) if a=="N"]
    nit=[i for i in N if (lambda d:(d.__setitem__(i,9),d.min())[1])(np.linalg.norm(sub_x-sub_x[i],axis=1))<1.25]
    o=sorted(nit,key=lambda i:np.linalg.norm(sub_x[i]-ag))
    def cn(i):
        d=np.linalg.norm(sub_x-sub_x[i],axis=1);d[i]=9;ci=int(np.argmin(d));ax=sub_x[i]-sub_x[ci];ax/=np.linalg.norm(ax);return sub_x[i]+2.3*ax,ax
    A,axA=cn(o[0]);B,axB=cn(o[1]);nrm=(axA+axB);nrm/=np.linalg.norm(nrm)
    res["HATCN"]=barrier("HATCN",os.path.join(RUNS,"HATCN_Ag_CN","xtbopt.xyz"),A,B,nrm,2,PLAIN,4,[-0.2,0.2,0.5])
    json.dump(res,open(jp,"w"),indent=2)

print(json.dumps(res,indent=2))
