"""Diagnostic C: branch discrimination on the reliable angles (65, 70).
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

Does adding a free low-index buried layer (partially-damaged GraHIL)
recover the physical branch (n800~2.0-2.1, d~250-270) at equal fit quality?
Window 650-1050nm."""
import numpy as np, ellipsometry_fit as ef
from scipy.optimize import minimize
base=str(DATA_DIR)
angles=[65.0,70.0,75.0]
wl,psi_m,del_m=ef.load_data_xlsx(base+r'\#2.xlsx',angles,650,1050)
Wd=ef.depol_weights(wl,angles,base+r'\#2_depol.xlsx',depol_cut=3.0,depol_soft=1.5)
C=lambda v,lo,hi:min(max(v,lo),hi)

def psidel(ang, Np, db, dlow, nlow, dr, f):
    n_air=np.ones(len(wl),dtype=complex)
    layers=[n_air]
    dl=[]
    if dr>0.5:
        layers.append(ef.bruggeman_ema(Np,n_air,f)); dl.append(dr)
    layers.append(Np); dl.append(db)
    if dlow>0.5:
        layers.append((nlow+0j)*np.ones(len(wl))); dl.append(dlow)
    layers+= [ef.n_SiO2(wl), ef.n_Si(wl)]
    dl.append(2.0)
    rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
    return np.degrees(np.arctan(np.abs(rho))), np.degrees(np.angle(rho))

def wrmse(sel,db,dlow,nlow,A,B,dr=0.0,f=0.3):
    n=A+B/wl**2; Np=(n+0j).astype(complex)
    num=0.0;den=0.0
    for ai,ang in sel:
        pc,dc=psidel(ang,Np,db,dlow,nlow,dr,f)
        dd=(dc-del_m[:,ai]+180)%360-180
        w=Wd[:,ai]
        num+=np.sum(w*((pc-psi_m[:,ai])**2+dd**2)); den+=2*np.sum(w)
    return np.sqrt(num/den)

def fit(fun,x0s):
    best=None
    for x0 in x0s:
        r=minimize(fun,x0,method='Nelder-Mead',options={'maxiter':25000,'xatol':1e-6,'fatol':1e-9,'adaptive':True})
        if best is None or r.fun<best.fun: best=r
    return best

sels={'65':[(0,65.0)],'70':[(1,70.0)],'65+70':[(0,65.0),(1,70.0)],'all3':[(0,65.0),(1,70.0),(2,75.0)]}
for name in ['65','70','65+70','all3']:
    sel=sels[name]
    # base: no low layer
    fb=lambda x: wrmse(sel,C(x[0],180,320),0,0,C(x[1],1.6,2.6),C(x[2],0,8e4))
    rb=fit(fb,[[226,2.35,5000],[270,2.05,15000],[240,2.2,10000]])
    db,A,B=C(rb.x[0],180,320),C(rb.x[1],1.6,2.6),C(rb.x[2],0,8e4)
    print('[%s base   ] RMSE=%6.3f d=%.1f n800=%.3f'%(name,rb.fun,db,A+B/800**2))
    # +low free
    fl=lambda x: wrmse(sel,C(x[0],180,320),C(x[1],0,90),C(x[2],1.35,1.95),C(x[3],1.6,2.6),C(x[4],0,8e4))
    rl=fit(fl,[[226,5,1.45,2.35,5000],[260,40,1.5,2.10,12000],[270,50,1.45,2.05,15000],[240,25,1.7,2.2,9000]])
    db,dlow,nlow,A,B=C(rl.x[0],180,320),C(rl.x[1],0,90),C(rl.x[2],1.35,1.95),C(rl.x[3],1.6,2.6),C(rl.x[4],0,8e4)
    print('[%s +low   ] RMSE=%6.3f d=%.1f dlow=%.1f nlow=%.2f n800=%.3f'%(name,rl.fun,db,dlow,nlow,A+B/800**2))
    # +low, d fixed 270
    ff=lambda x: wrmse(sel,270.0,C(x[0],0,90),C(x[1],1.35,1.95),C(x[2],1.6,2.6),C(x[3],0,8e4))
    rf=fit(ff,[[40,1.5,2.05,12000],[10,1.45,2.1,9000],[70,1.6,2.0,15000]])
    dlow,nlow,A,B=C(rf.x[0],0,90),C(rf.x[1],1.35,1.95),C(rf.x[2],1.6,2.6),C(rf.x[3],0,8e4)
    print('[%s d=270  ] RMSE=%6.3f dlow=%.1f nlow=%.2f n800=%.3f'%(name,rf.fun,dlow,nlow,A+B/800**2))
    # +low +rough (only for joint, else overfit)
    if '+' in name or name=='all3':
        fr=lambda x: wrmse(sel,C(x[0],180,320),C(x[1],0,90),C(x[2],1.35,1.95),C(x[3],1.6,2.6),C(x[4],0,8e4),
                           dr=C(x[5],0,140),f=C(x[6],0.05,0.7))
        rr=fit(fr,[[226,5,1.45,2.35,5000,20,0.3],[255,40,1.5,2.1,12000,80,0.25],[270,50,1.45,2.05,15000,100,0.2]])
        db,dlow,nlow,A,B=C(rr.x[0],180,320),C(rr.x[1],0,90),C(rr.x[2],1.35,1.95),C(rr.x[3],1.6,2.6),C(rr.x[4],0,8e4)
        dr_,f_=C(rr.x[5],0,140),C(rr.x[6],0.05,0.7)
        print('[%s +low+rg] RMSE=%6.3f d=%.1f dlow=%.1f nlow=%.2f rough=%.0f(f=%.2f) n800=%.3f'%(
              name,rr.fun,db,dlow,nlow,dr_,f_,A+B/800**2))
    print()
