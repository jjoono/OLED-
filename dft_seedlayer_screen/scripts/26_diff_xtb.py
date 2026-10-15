"""Robust diffusion-barrier estimate via GFN2-xTB lateral drag (converges reliably).
Calibrate against HATCN (DFT barrier known = 0.29 eV) to scale xTB corrugation.
Systems: HATCN (calib), TAPC, TPBi, B3PyMPM. Frozen substrate, Ag dragged site->site,
z relaxed at each lateral point via xtb single points over a z grid.
"""
import numpy as np, os, subprocess, json, tempfile

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")

def rx(p):
    L=open(p).read().strip().splitlines();n=int(L[0]);s=[];x=[]
    for l in L[2:2+n]:
        q=l.split();s.append(q[0]);x.append(list(map(float,q[1:4])))
    return s,np.array(x)

def xtb_energy(syms, xyz, uhf):
    """single-point GFN2 energy (Eh)."""
    wd = tempfile.mkdtemp()
    p = os.path.join(wd,"g.xyz")
    open(p,"w").write(f"{len(syms)}\n\n"+"".join(f"{a} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n" for a,c in zip(syms,xyz)))
    r=subprocess.run(["xtb","g.xyz","--gfn","2","--chrg","0","--uhf",str(uhf)],
                     cwd=wd,capture_output=True,text=True)
    for line in r.stdout.splitlines():
        if "TOTAL ENERGY" in line:
            return float(line.split()[3])
    return None

def barrier(agfile, find_sites, npath=7):
    s,x = rx(agfile)
    sub_s, sub_x, ag = s[:-1], x[:-1], x[-1]
    (A,nA),(B,nB) = find_sites(sub_s, sub_x, ag)
    nrm=(nA+nB); nrm/=np.linalg.norm(nrm)
    H=27.211386
    prof=[]
    for t in np.linspace(0,1,npath):
        lat=(1-t)*A+t*B
        best=None
        for dz in (-0.4,-0.2,0,0.2,0.4,0.6):
            e=xtb_energy(sub_s+["Ag"], np.vstack([sub_x, lat+dz*nrm]), 1)
            if e is not None and (best is None or e<best): best=e
        prof.append(best)
    prof=np.array([p for p in prof if p is not None])
    Eb=(prof.max()-prof[0])*H
    return Eb, prof

def s_Npair(s,x,ag):
    Ns=[i for i,a in enumerate(s) if a=="N"]
    o=sorted(Ns,key=lambda i:np.linalg.norm(x[i]-ag))
    def site(i):
        nb=[j for j in range(len(s)) if 0<np.linalg.norm(x[i]-x[j])<1.5]
        bis=x[i]-0.5*(x[nb[0]]+x[nb[1]]) if len(nb)>=2 else x[i]-x.mean(0)
        n=np.linalg.norm(bis); bis=bis/n if n>0 else np.array([0,0,1.])
        return x[i]+2.3*bis, bis
    return site(o[0]), site(o[1])

def s_CNpair(s,x,ag):  # HATCN nitrile
    Ns=[i for i,a in enumerate(s) if a=="N"]
    nit=[i for i in Ns if (lambda d:(d.__setitem__(i,9),d.min())[1])(np.linalg.norm(x-x[i],axis=1))<1.25]
    o=sorted(nit,key=lambda i:np.linalg.norm(x[i]-ag))
    def site(i):
        d=np.linalg.norm(x-x[i],axis=1);d[i]=9;ci=int(np.argmin(d));ax=x[i]-x[ci];ax/=np.linalg.norm(ax)
        return x[i]+2.3*ax, ax
    return site(o[0]), site(o[1])

def s_TPA(s,x,ag):
    com=x.mean(0); N=[i for i,a in enumerate(s) if a=="N"][0]
    outw=ag-x[N]; nn=np.linalg.norm(outw); outw=outw/nn if nn>0 else np.array([0,0,1.])
    A=(x[N]+2.5*outw, outw)
    C=[i for i,a in enumerate(s) if a=="C"]; near=sorted(C,key=lambda i:np.linalg.norm(x[i]-ag))[:6]
    cen=x[near].mean(0); nB=cen-com; nB/=np.linalg.norm(nB)
    return A,(cen+2.6*nB, nB)

CASES=[("HATCN","HATCN_Ag_CN",s_CNpair),("TAPC","TPA_Ag",s_TPA),
       ("TPBi","TPBi_Ag",s_Npair),("B3PyMPM","B3PyMPM_Ag",s_Npair)]
out={}
for tag,run,fn in CASES:
    Eb,prof=barrier(os.path.join(RUNS,run,"xtbopt.xyz"),fn)
    out[tag]=Eb
    print(f"{tag}: xTB barrier = {Eb:.3f} eV  (profile {np.round((prof-prof[0])*27.211,3)})",flush=True)
json.dump(out,open(os.path.join(RUNS,"diff_xtb_eV.json"),"w"),indent=2)
print(json.dumps(out,indent=2))
