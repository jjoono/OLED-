"""Ag adatom surface diffusion barrier (item #2) on 3 surfaces:
  clean Al (Al13 cluster), AlOx (Al4O6), HATCN (nitrile site -> adjacent site).
Method: frozen substrate, lateral drag of Ag from site A toward B; at each lateral
position relax Ag height by a small z-scan; DFT PBE-D3BJ/def2-SVP single points.
Barrier = max(E_path) - E(site).  Screening-level (frozen substrate, cluster).
"""
import numpy as np, os, json, psi4

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")
psi4.set_memory("5 GB"); psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_diff.out"), False)
H2EV = 27.211386
BASEOPTS = {"basis":"def2-svp","scf_type":"df","reference":"uks","maxiter":300,"guess":"sad"}
ROB = dict(BASEOPTS, level_shift=1.0, level_shift_cutoff=1e-3, damping_percentage=15.0, maxiter=400)

def read_xyz(p):
    L=open(p).read().strip().splitlines(); n=int(L[0]); s,x=[],[]
    for l in L[2:2+n]:
        q=l.split(); s.append(q[0]); x.append(list(map(float,q[1:4])))
    return s, np.array(x)

def gstr(syms, xyz, mult, extra_ghost=None):
    st=f"0 {mult}\n"
    for i,(s,c) in enumerate(zip(syms,xyz)):
        st+=f"{s} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n"
    return st+"symmetry c1\nno_reorient\nno_com\n"

def sp(syms, xyz, mult, opts):
    psi4.set_options(opts)
    e=psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms,xyz,mult)))
    psi4.core.clean(); return e

def relaxed_point(syms, xyz_sub, ag_xy, z0, normal, mult, opts):
    """place Ag at lateral ag_xy, scan height along `normal` around z0, return min E."""
    best=None
    for dz in (-0.4,-0.2,0.0,0.2,0.4,0.7):
        ag = ag_xy + dz*normal
        xyz=np.vstack([xyz_sub, ag])
        try:
            e=sp(syms+["Ag"], xyz, mult, opts)
        except Exception:
            continue
        if best is None or e<best: best=e
    return best

def barrier(tag, sub_xyz_file, site_a, site_b, normal, mult, opts, npath=5):
    syms, xyz = read_xyz(sub_xyz_file)
    # substrate = all but last (Ag) if file is the *_Ag geometry; detect
    if syms[-1]=="Ag":
        syms=syms[:-1]; xyz=xyz[:-1]
    energies=[]
    for t in np.linspace(0,1,npath):
        lateral = (1-t)*site_a + t*site_b
        e = relaxed_point(syms, xyz, lateral, 0.0, normal, mult, opts)
        energies.append(e)
        es = "FAIL" if e is None else f"{e:.6f}"
        print(f"  {tag} t={t:.2f}: E={es}", flush=True)
    vals=[(i,e) for i,e in enumerate(energies) if e is not None]
    if len(vals)<2:
        print(f"{tag}: barrier FAILED (too few points)", flush=True); return None
    e0 = energies[0] if energies[0] is not None else vals[0][1]
    emax = max(e for _,e in vals)
    Eb=(emax-e0)*H2EV
    print(f"{tag}: diffusion barrier ~ {Eb:.3f} eV (site->max along path)", flush=True)
    return Eb

res={}

# ---- HATCN: nitrile N site -> adjacent nitrile N ----
syms,xyz = read_xyz(os.path.join(RUNS,"HATCN_Ag_CN","xtbopt.xyz"))
sub_s, sub_x = syms[:-1], xyz[:-1]
# find the two nitrile N nearest to the Ag in optimized geom
agpos = xyz[-1]
Nidx=[i for i,s in enumerate(sub_s) if s=="N"]
# nitrile N = N with single close C neighbor (~1.16A)
nitrile=[]
for i in Nidx:
    d=np.linalg.norm(sub_x-sub_x[i],axis=1); d[i]=9
    if d.min()<1.25: nitrile.append(i)
# site A = nearest nitrile to Ag; site B = next nearest nitrile
order=sorted(nitrile, key=lambda i: np.linalg.norm(sub_x[i]-agpos))
A_i, B_i = order[0], order[1]
# place sites 2.3 A out along each C#N axis (approx Ag on-site positions)
def cn_site(i):
    d=np.linalg.norm(sub_x-sub_x[i],axis=1); d[i]=9; ci=int(np.argmin(d))
    ax=sub_x[i]-sub_x[ci]; ax/=np.linalg.norm(ax)
    return sub_x[i]+2.3*ax, ax
siteA,axA = cn_site(A_i); siteB,axB = cn_site(B_i)
# normal ~ average CN axis (out-of-plane relaxation direction)
normal=(axA+axB); normal/=np.linalg.norm(normal)
res["HATCN"]=barrier("HATCN", os.path.join(RUNS,"HATCN_Ag_CN","xtbopt.xyz"),
                     siteA, siteB, normal, 2, ROB, npath=5)
json.dump(res,open(os.path.join(RUNS,"diffusion_barrier_eV.json"),"w"),indent=2)

# ---- AlOx (Al4O6): O site -> adjacent O site ----
syms,xyz=read_xyz(os.path.join(RUNS,"Al4O6_Ag","xtbopt.xyz"))
sub_s,sub_x=syms[:-1],xyz[:-1]; agpos=xyz[-1]
Oidx=[i for i,s in enumerate(sub_s) if s=="O"]
com=sub_x.mean(axis=0)
order=sorted(Oidx,key=lambda i: np.linalg.norm(sub_x[i]-agpos))
A_i,B_i=order[0],order[1]
def o_site(i):
    outw=sub_x[i]-com; outw/=np.linalg.norm(outw); return sub_x[i]+2.2*outw, outw
siteA,nA=o_site(A_i); siteB,nB=o_site(B_i)
normal=(nA+nB); normal/=np.linalg.norm(normal)
res["AlOx"]=barrier("AlOx", os.path.join(RUNS,"Al4O6_Ag","xtbopt.xyz"),
                    siteA, siteB, normal, 2, ROB, npath=5)
json.dump(res,open(os.path.join(RUNS,"diffusion_barrier_eV.json"),"w"),indent=2)

# ---- clean Al (Al13): top site -> adjacent surface atom (hollow/bridge) ----
syms,xyz=read_xyz(os.path.join(RUNS,"Al13_Ag","xtbopt.xyz"))
sub_s,sub_x=syms[:-1],xyz[:-1]; agpos=xyz[-1]
com=sub_x.mean(axis=0)
# two outermost Al atoms nearest Ag = path across facet
order=sorted(range(len(sub_s)),key=lambda i: np.linalg.norm(sub_x[i]-agpos))
A_i,B_i=order[0],order[1]
def al_site(i):
    outw=sub_x[i]-com; outw/=np.linalg.norm(outw); return sub_x[i]+2.6*outw, outw
siteA,nA=al_site(A_i); siteB,nB=al_site(B_i)
normal=(nA+nB); normal/=np.linalg.norm(normal)
res["cleanAl"]=barrier("cleanAl", os.path.join(RUNS,"Al13_Ag","xtbopt.xyz"),
                       siteA, siteB, normal, 1, ROB, npath=5)
json.dump(res,open(os.path.join(RUNS,"diffusion_barrier_eV.json"),"w"),indent=2)
print(json.dumps(res,indent=2))
