"""Final: k-anchored (UV-Vis SI Fig.6a) + total thickness FREE within SEM
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

window [260,280]nm. Uncertainty re-estimated: optimizer restarts + angle-pair
+ blue-k sensitivity, all within the window."""
import numpy as np, ellipsometry_fit as ef, csv
from scipy.optimize import minimize
base=str(DATA_DIR)
scr=str(SCRATCH_DIR)
angles=[65.0,70.0,75.0]
wl,psi_m,del_m=ef.load_data_xlsx(base+r'\#2.xlsx',angles,380,1688)
Wd=ef.depol_weights(wl,angles,base+r'\#2_depol.xlsx',depol_cut=3.0,depol_soft=1.5)
CL=lambda v,lo,hi:min(max(v,lo),hi)
wl_si,k_si=[],[]
with open(scr+r'\k_SI_digitized.csv') as f:
    r=csv.reader(f); next(r)
    for row in r: wl_si.append(float(row[0])); k_si.append(float(row[1]))
wl_si=np.array(wl_si); k_si=np.array(k_si)

def stack_pd(Np,dr,f,db):
    n_air=np.ones(len(wl),dtype=complex)
    layers=[n_air, ef.bruggeman_ema(Np,n_air,f), Np, ef.n_SiO2(wl), ef.n_Si(wl)]
    dl=[dr,db,2.0]
    P=[];D=[]
    for ang in angles:
        rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
        P.append(np.degrees(np.arctan(np.abs(rho)))); D.append(np.degrees(np.angle(rho)))
    return np.array(P).T, np.array(D).T

def unpack(x):
    dr=CL(x[0],40,140); f=CL(x[1],0.08,0.60); db=CL(x[2],170,300)
    epsinf=CL(x[3],1.0,3.0)
    tlEg=CL(x[7],2.20,2.42); uvEg=CL(x[14],2.4,3.6)
    osc=[{'type':'TL','A':CL(x[4],0,300),'E0':CL(x[5],tlEg+0.05,3.3),'C':CL(x[6],0.05,1.2),'Eg':tlEg},
         {'type':'Gaussian','A':CL(x[8],0,5),'Ecen':CL(x[9],2.30,2.42),'Br':CL(x[10],0.03,0.10)},
         {'type':'TL','A':CL(x[11],0,400),'E0':CL(x[12],uvEg+0.05,5.5),'C':CL(x[13],0.2,2.5),'Eg':uvEg}]
    return dr,f,db,epsinf,osc

def huber(r,d=3.0):
    a=np.abs(r); return np.where(a<=d,r*r,2*d*a-d*d)

def make_loss(angq, k450_target=None):
    Wt=Wd*np.array(angq)[None,:]*np.where(wl<560,0.6,1.0)[:,None]
    sW=Wt.sum()
    def loss(x):
        dr,f,db,epsinf,osc=unpack(x)
        total=db+f*dr
        try: n,k=ef.genosc_nk(wl,osc,eps_inf=epsinf)
        except Exception: return 1e7
        if np.any(n<1.4) or np.any(n>2.9): return 1e7
        k_at=np.interp(wl_si,wl,k)
        pen=2e5*np.mean((k_at-k_si)**2)+3000.0*np.mean(k[wl>572]**2)
        pen+=50.0*(max(0.0,total-280.0)**2+max(0.0,260.0-total)**2)  # SEM window
        if k450_target is not None:
            pen+=1e5*(np.interp(450,wl,k)-k450_target)**2
        Np=(n+1j*k).astype(complex)
        pc,dc=stack_pd(Np,dr,f,db)
        dd=(dc-del_m+180)%360-180
        return (np.sum(Wt*huber(pc-psi_m))+np.sum(Wt*huber(dd)))/(2*sW)+pen
    return loss

def run(loss,x0s,it):
    best=None
    for x0 in x0s:
        r=minimize(loss,x0,method='Nelder-Mead',options={'maxiter':it,'xatol':1e-6,'fatol':1e-9,'adaptive':True})
        if best is None or r.fun<best.fun: best=r
    return best

seeds=[[90,0.25,255,1.9, 30,2.50,0.30,2.30, 0.6,2.365,0.06, 60,3.6,1.2,2.6],
       [70,0.35,258,2.0, 15,2.45,0.20,2.33, 0.7,2.36,0.055, 90,3.8,1.0,2.7],
       [110,0.18,262,1.8, 50,2.55,0.40,2.25, 0.5,2.37,0.07, 40,3.4,1.5,2.5]]

lM=make_loss([1,1,0.35])
print('=== MAIN (window 260-280) ===')
rM=run(lM,seeds,25000); rM=run(lM,[rM.x],20000)
dr,f,db,epsinf,osc=unpack(rM.x)
print('MAIN loss=%.2f dr=%.1f f=%.2f db=%.1f TOTAL=%.1f epsinf=%.2f'%(rM.fun,dr,f,db,db+f*dr,epsinf))
for o in osc: print('  ',{kk:(round(v,3) if isinstance(v,float) else v) for kk,v in o.items()})
x_main=rM.x

variants={'MAIN':x_main}
# optimizer scatter: perturbed restart
rng=np.random.default_rng(1)
xp=np.array(x_main)*(1+0.05*rng.standard_normal(len(x_main)))
rR=run(lM,[xp],18000); variants['RESTART']=rR.x
# angle pair
rP=run(make_loss([1,1,0]),[x_main],14000); variants['P6570']=rP.x
# blue-k sensitivity within window
n0,k0=None,None
dr,f,db,epsinf,osc=unpack(x_main)
n0,k0=ef.genosc_nk(wl,osc,eps_inf=epsinf)
k450=float(np.interp(450,wl,k0))
rL=run(make_loss([1,1,0.35],k450*0.6),[x_main],14000); variants['K450LO']=rL.x
rH=run(make_loss([1,1,0.35],k450*1.4),[x_main],14000); variants['K450HI']=rH.x

curves={}
for nm,x in variants.items():
    dr,f,db,epsinf,osc=unpack(x)
    n,k=ef.genosc_nk(wl,osc,eps_inf=epsinf)
    curves[nm]=(n,k)
    print('%-8s total=%.1f  n600=%.3f n535=%.3f n450=%.3f'%(
          nm,db+f*dr,np.interp(600,wl,n),np.interp(535,wl,n),np.interp(450,wl,n)))
nM,kM=curves['MAIN']
alln=np.array([c[0] for c in curves.values()])
n_lo,n_hi=alln.min(0),alln.max(0)

print('\n wl  |  n    ( lo .. hi )   half-width')
for tw in [400,420,450,480,500,510,524,535,549,560,575,600]:
    i=int(np.argmin(np.abs(wl-tw)))
    print(' %4d | %.3f (%.3f..%.3f)  ±%.3f'%(wl[i],nM[i],n_lo[i],n_hi[i],(n_hi[i]-n_lo[i])/2))

# CE MSE of MAIN
dr,f,db,epsinf,osc=unpack(x_main)
Np=(nM+1j*kM).astype(complex)
pc,dc=stack_pd(Np,dr,f,db)
def ncs(p,d):
    p=np.deg2rad(p); d=np.deg2rad(d)
    return np.cos(2*p),np.sin(2*p)*np.cos(d),np.sin(2*p)*np.sin(d)
Nm,Cm,Sm=ncs(psi_m,del_m); Nc,Cc,Sc=ncs(pc,dc)
r2=(Nm-Nc)**2+(Cm-Cc)**2+(Sm-Sc)**2
msk=Wd>0.5
print('\nMSE_CE(reliable)=%.1f'%(1000*np.sqrt(r2[msk].sum()/(3*msk.sum()-15))))

with open(str(OUT_DIR / 'Perov_final_window_nk.csv'),'w',newline='') as fo:
    w=csv.writer(fo); w.writerow(['wl_nm','n','k','n_lo','n_hi'])
    for i in range(len(wl)):
        w.writerow(['%.3f'%wl[i],'%.5f'%nM[i],'%.5f'%kM[i],'%.5f'%n_lo[i],'%.5f'%n_hi[i]])
dd=(dc-del_m+180)%360-180
np.savez(str(OUT_DIR / 'Perov_final_window_fit.npz'),wl=wl,psi_m=psi_m,del_m=del_m,
         psi_c=pc,del_c=dc,n=nM,k=kM,n_lo=n_lo,n_hi=n_hi,W=Wd,
         d_bulk=db,d_rough=dr,f_rough=f,wl_si=wl_si,k_si=k_si,x_main=x_main)
print('saved Perov_final_window_nk.csv / .npz')
