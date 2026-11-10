"""Diagnostic A: root-cause of 3-angle inconsistency in #2 (perovskite).
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

Window 650-1050nm (transparent perovskite, opaque Si, low depol).
Variants: 0 baseline / conv Delta-sign / 1 per-angle independent /
          2 per-angle thickness / 3 per-angle angle-offset / 4 both."""
import numpy as np, ellipsometry_fit as ef
from scipy.optimize import minimize
base=str(DATA_DIR)
angles=[65.0,70.0,75.0]
Ng=ef.load_grahil_nk(str(OUT_DIR / 'GraHIL_genosc_nk.csv'))
wl,psi_m,del_m=ef.load_data_xlsx(base+r'\#2.xlsx',angles,650,1050)
Wd=ef.depol_weights(wl,angles,base+r'\#2_depol.xlsx',depol_cut=3.0,depol_soft=1.5)
print('window 650-1050: N=%d, mean depol-weight=%.3f'%(len(wl),Wd.mean()))

def model_one(ang, dr, f, db):
    n_air=np.ones(len(wl),dtype=complex)
    n=None  # set by caller closure
def psidel(ang, Np, dr, f, db):
    n_air=np.ones(len(wl),dtype=complex)
    layers=[n_air, ef.bruggeman_ema(Np,n_air,f), Np, ef.n_SiO2(wl), ef.n_Si(wl)]
    d=[dr,db,2.0]
    rp,rs=ef._tmm(wl,layers,d,ang); rho=rp/rs
    return np.degrees(np.arctan(np.abs(rho))), np.degrees(np.angle(rho))

def wrmse(sel_angles, drs, fs, dbs, A, B, dths=None, flipdel=False):
    n=A+B/wl**2; Np=(n+0j).astype(complex)
    num=0.0; den=0.0
    for idx,(ai,ang) in enumerate(sel_angles):
        th=ang+(dths[idx] if dths is not None else 0.0)
        pc,dc=psidel(th,Np,drs[idx],fs[idx],dbs[idx])
        dm=(360-del_m[:,ai]) if flipdel else del_m[:,ai]
        dd=(dc-dm+180)%360-180
        w=Wd[:,ai]
        num+=np.sum(w*((pc-psi_m[:,ai])**2+dd**2)); den+=2*np.sum(w)
    return np.sqrt(num/den)

def fit(fun,x0s):
    best=None
    for x0 in x0s:
        r=minimize(fun,x0,method='Nelder-Mead',options={'maxiter':20000,'xatol':1e-6,'fatol':1e-9,'adaptive':True})
        if best is None or r.fun<best.fun: best=r
    return best

C=lambda v,lo,hi:min(max(v,lo),hi)
seeds=[[98,0.23,259,2.04,20000],[60,0.40,240,2.10,15000],[120,0.15,270,2.00,25000]]

# --- variant 0: joint baseline (shared everything)
def f0(x):
    dr,f,db,A,B=C(x[0],20,160),C(x[1],0.05,0.7),C(x[2],180,320),C(x[3],1.6,2.6),C(x[4],0,8e4)
    all3=[(i,a) for i,a in enumerate(angles)]
    return wrmse(all3,[dr]*3,[f]*3,[db]*3,A,B)
r0=fit(f0,seeds)
dr,f,db,A,B=C(r0.x[0],20,160),C(r0.x[1],0.05,0.7),C(r0.x[2],180,320),C(r0.x[3],1.6,2.6),C(r0.x[4],0,8e4)
print('\n[V0 baseline joint]  RMSE=%.3f  dr=%.1f f=%.2f db=%.1f n800=%.3f'%(r0.fun,dr,f,db,A+B/800**2))

# --- Delta-sign convention check
def f0f(x):
    dr,f,db,A,B=C(x[0],20,160),C(x[1],0.05,0.7),C(x[2],180,320),C(x[3],1.6,2.6),C(x[4],0,8e4)
    all3=[(i,a) for i,a in enumerate(angles)]
    return wrmse(all3,[dr]*3,[f]*3,[db]*3,A,B,flipdel=True)
r0f=fit(f0f,seeds)
print('[V0-conv Delta→360-Delta] RMSE=%.3f  (should be much worse if convention OK)'%r0f.fun)

# --- variant 1: per-angle independent
print('\n[V1 per-angle independent]')
for i,a in enumerate(angles):
    def f1(x,i=i,a=a):
        dr,f,db,A,B=C(x[0],20,160),C(x[1],0.05,0.7),C(x[2],180,320),C(x[3],1.6,2.6),C(x[4],0,8e4)
        return wrmse([(i,a)],[dr],[f],[db],A,B)
    r1=fit(f1,seeds)
    dr,f,db,A,B=C(r1.x[0],20,160),C(r1.x[1],0.05,0.7),C(r1.x[2],180,320),C(r1.x[3],1.6,2.6),C(r1.x[4],0,8e4)
    print('  %.0f°: RMSE=%.3f  dr=%.1f f=%.2f db=%.1f  n800=%.3f n1000=%.3f'%(a,r1.fun,dr,f,db,A+B/800**2,A+B/1000**2))

# --- variant 2: per-angle thickness (shared optics)
def f2(x):
    dr,f=C(x[0],20,160),C(x[1],0.05,0.7)
    d1,d2,d3=C(x[2],180,320),C(x[3],180,320),C(x[4],180,320)
    A,B=C(x[5],1.6,2.6),C(x[6],0,8e4)
    all3=[(i,a) for i,a in enumerate(angles)]
    return wrmse(all3,[dr]*3,[f]*3,[d1,d2,d3],A,B)
s2=[[98,0.23,259,259,259,2.04,20000],[60,0.4,240,250,260,2.1,15000]]
r2=fit(f2,s2)
d1,d2,d3=C(r2.x[2],180,320),C(r2.x[3],180,320),C(r2.x[4],180,320)
print('\n[V2 per-angle d]  RMSE=%.3f  d65=%.1f d70=%.1f d75=%.1f  (spread %.1fnm)  n800=%.3f'%(
      r2.fun,d1,d2,d3,max(d1,d2,d3)-min(d1,d2,d3),C(r2.x[5],1.6,2.6)+C(r2.x[6],0,8e4)/800**2))

# --- variant 3: per-angle angle offset (shared d)
def f3(x):
    dr,f,db,A,B=C(x[0],20,160),C(x[1],0.05,0.7),C(x[2],180,320),C(x[3],1.6,2.6),C(x[4],0,8e4)
    dt=[C(x[5],-1.5,1.5),C(x[6],-1.5,1.5),C(x[7],-1.5,1.5)]
    all3=[(i,a) for i,a in enumerate(angles)]
    return wrmse(all3,[dr]*3,[f]*3,[db]*3,A,B,dths=dt)
s3=[s+[0,0,0] for s in seeds]
r3=fit(f3,s3)
dt=[C(r3.x[5],-1.5,1.5),C(r3.x[6],-1.5,1.5),C(r3.x[7],-1.5,1.5)]
print('[V3 angle offsets]  RMSE=%.3f  δθ=%.2f/%.2f/%.2f°'%(r3.fun,dt[0],dt[1],dt[2]))

# --- variant 4: both
def f4(x):
    dr,f=C(x[0],20,160),C(x[1],0.05,0.7)
    d1,d2,d3=C(x[2],180,320),C(x[3],180,320),C(x[4],180,320)
    A,B=C(x[5],1.6,2.6),C(x[6],0,8e4)
    dt=[C(x[7],-1.5,1.5),C(x[8],-1.5,1.5),C(x[9],-1.5,1.5)]
    all3=[(i,a) for i,a in enumerate(angles)]
    return wrmse(all3,[dr]*3,[f]*3,[d1,d2,d3],A,B,dths=dt)
s4=[[98,0.23,259,259,259,2.04,20000,0,0,0]]
r4=fit(f4,s4)
d1,d2,d3=C(r4.x[2],180,320),C(r4.x[3],180,320),C(r4.x[4],180,320)
dt=[C(r4.x[7],-1.5,1.5),C(r4.x[8],-1.5,1.5),C(r4.x[9],-1.5,1.5)]
print('[V4 d+offsets]  RMSE=%.3f  d=%.1f/%.1f/%.1f  δθ=%.2f/%.2f/%.2f'%(
      r4.fun,d1,d2,d3,dt[0],dt[1],dt[2]))

# GraHIL #1 depol sanity in its fit band
wg,dg=ef.load_depol(base+r'\#1_depol;.xlsx',angles)
m=(wg>=400)&(wg<=1000)
print('\n[GraHIL #1 depol in 400-1000nm band] mean|depol|=%.2f%% max=%.2f%% -> Psi/Delta reliable'%(
      np.abs(dg[m]).mean(),np.abs(dg[m]).max()))
