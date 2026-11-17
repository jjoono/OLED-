"""Error budget for n in 400-600nm (k-anchored fit).
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

Variants: thickness 260/270/280nm (SEM +-10), blue-k(450nm) forced low/high
(k<510nm is NOT UV-Vis-anchored -> quantify its leverage on n)."""
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

Wt=Wd*np.array([1,1,0.35])[None,:]*np.where(wl<560,0.6,1.0)[:,None]
sW=Wt.sum()

def make_loss(d_target, k450_target=None):
    def loss(x):
        dr,f,db,epsinf,osc=unpack(x)
        try: n,k=ef.genosc_nk(wl,osc,eps_inf=epsinf)
        except Exception: return 1e7
        if np.any(n<1.4) or np.any(n>2.9): return 1e7
        k_at=np.interp(wl_si,wl,k)
        pen=2e5*np.mean((k_at-k_si)**2)+3000.0*np.mean(k[wl>572]**2)
        pen+=50.0*(db+f*dr-d_target)**2
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

s1=[90,0.25,248,1.9, 30,2.50,0.30,2.30, 0.6,2.365,0.06, 60,3.6,1.2,2.6]
s2=[70,0.35,246,2.0, 15,2.45,0.20,2.33, 0.7,2.36,0.055, 90,3.8,1.0,2.7]

print('=== reference (270nm) ===')
rr=run(make_loss(270.0),[s1,s2],22000); rr=run(make_loss(270.0),[rr.x],15000)
x_ref=rr.x
dr,f,db,epsinf,osc=unpack(x_ref)
n_ref,k_ref=ef.genosc_nk(wl,osc,eps_inf=epsinf)
print('ref loss=%.2f total=%.1f k450=%.3f'%(rr.fun,db+f*dr,np.interp(450,wl,k_ref)))

k450_ref=float(np.interp(450,wl,k_ref))
cases=[('D260',260.0,None),('D280',280.0,None),
       ('K450LO',270.0,k450_ref*0.6),('K450HI',270.0,k450_ref*1.4)]
curves={'REF':(n_ref,k_ref)}
for name,dt,k4 in cases:
    r=run(make_loss(dt,k4),[x_ref],14000)
    dr,f,db,epsinf,osc=unpack(r.x)
    n,k=ef.genosc_nk(wl,osc,eps_inf=epsinf)
    curves[name]=(n,k)
    print('%-6s loss=%.2f total=%.1f k450=%.3f'%(name,r.fun,db+f*dr,np.interp(450,wl,k)))

print('\nΔn vs REF at key wavelengths:')
print(' wl  | n_ref |  D260    D280   | K450LO  K450HI')
for tw in [400,420,450,480,510,535,550,560,580,600]:
    i=int(np.argmin(np.abs(wl-tw)))
    row=' %4d | %.3f |'%(wl[i],n_ref[i])
    for nm in ['D260','D280']:
        row+=' %+.3f '%(curves[nm][0][i]-n_ref[i])
    row+=' |'
    for nm in ['K450LO','K450HI']:
        row+=' %+.3f '%(curves[nm][0][i]-n_ref[i])
    print(row)

np.savez(str(OUT_DIR / 'Perov_errbudget.npz'),wl=wl,
         **{('n_'+nm):curves[nm][0] for nm in curves},
         **{('k_'+nm):curves[nm][1] for nm in curves})
print('saved Perov_errbudget.npz')
