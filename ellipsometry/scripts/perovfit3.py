import numpy as np, ellipsometry_fit as ef, csv
from scipy.optimize import minimize
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

base=str(DATA_DIR)
fp=base+r'\#2.xlsx'; depf=base+r'\#2_depol.xlsx'
angles=[65,70,75]
Ng=ef.load_grahil_nk(str(OUT_DIR / 'GraHIL_genosc_nk.csv'))
wl,psi_m,del_m=ef.load_data_xlsx(fp,angles,380,1688)
Wd=ef.depol_weights(wl,angles,depf,depol_cut=3.0,depol_soft=1.5)

D_GRAHIL=15.0   # thin, DMSO-thinned/damaged GraHIL (SEM shows very thin red Buf-HIL line)

def stack(wlq, Np, d_rough, f_rough, d_bulk):
    n_air=np.ones(len(wlq),dtype=complex)
    Ndam=ef.bruggeman_ema(Ng(wlq), Np, 0.0)   # damaged GraHIL = index-matched (100% perov)
    Nrough=ef.bruggeman_ema(Np, n_air, f_rough)  # f_rough = perovskite vol-fraction in rough layer
    layers=[n_air, Nrough, Np, Ndam, ef.n_SiO2(wlq), ef.n_Si(wlq)]
    d_list=[d_rough, d_bulk, D_GRAHIL, 2.0]
    po,do=[],[]
    for ang in angles:
        rp,rs=ef._tmm(wlq,layers,d_list,ang); rho=rp/rs
        po.append(np.degrees(np.arctan(np.abs(rho)))); do.append(np.degrees(np.angle(rho)))
    return np.array(po).T, np.array(do).T

def oscs(p):
    epsinf=max(p[0],1.0)
    Eg=min(max(p[4],2.05),2.45); uEg=min(max(p[11],2.4),3.6)
    return epsinf,[
       {'type':'TL','A':min(max(p[1],0),300),'E0':min(max(p[2],Eg+0.05),3.3),'C':min(max(p[3],0.05),1.5),'Eg':Eg},
       {'type':'Gaussian','A':min(max(p[5],0),5),'Ecen':min(max(p[6],2.20),2.45),'Br':min(max(p[7],0.03),0.09)},
       {'type':'TL','A':min(max(p[8],0),400),'E0':min(max(p[9],uEg+0.05),5.5),'C':min(max(p[10],0.2),2.5),'Eg':uEg}]

def unpack(x):
    d_rough=min(max(x[0],30),150); f_rough=min(max(x[1],0.2),0.65); d_bulk=min(max(x[2],120),260)
    return d_rough,f_rough,d_bulk,x[3:]

def resid(x):
    d_rough,f_rough,d_bulk,op=unpack(x)
    epsinf,L=oscs(op)
    try: n,k=ef.genosc_nk(wl,L,eps_inf=epsinf)
    except Exception: return 1e6
    if np.any(n<1.5) or np.any(n>2.8): return 1e6
    Np=(n+1j*k).astype(complex)
    pc,dc=stack(wl,Np,d_rough,f_rough,d_bulk)
    dd=(dc-del_m+180)%360-180
    mse=(np.sum(Wd*(pc-psi_m)**2)+np.sum(Wd*dd**2))/(2*np.sum(Wd))
    ktr=80.0*np.mean((k[wl>620])**2)
    prior=0.02*(d_bulk+f_rough*d_rough-270)**2   # perovskite mass ~270nm-equiv
    return mse+ktr+prior

osc0=[1.0,50,2.42,0.42,2.05, 1.47,2.24,0.09, 38,4.9,0.2,2.4]
best=None
for dr0 in [60,90,120]:
  for fr0 in [0.35,0.5]:
    for db0 in [180,210]:
      r=minimize(resid,[dr0,fr0,db0]+osc0,method='Nelder-Mead',
                 options={'maxiter':40000,'xatol':1e-6,'fatol':1e-9,'adaptive':True})
      if best is None or r.fun<best[0]: best=(r.fun,r.x)
r=minimize(resid,best[1],method='Nelder-Mead',options={'maxiter':60000,'xatol':1e-7,'fatol':1e-10,'adaptive':True})
d_rough,f_rough,d_bulk,op=unpack(r.x)
epsinf,L=oscs(op)
n,k=ef.genosc_nk(wl,L,eps_inf=epsinf); Np=(n+1j*k).astype(complex)
pc,dc=stack(wl,Np,d_rough,f_rough,d_bulk)
dd=(dc-del_m+180)%360-180
wmse=np.sqrt((np.sum(Wd*(pc-psi_m)**2)+np.sum(Wd*dd**2))/(2*np.sum(Wd)))
rel=Wd>0.5
rmse_rel=np.sqrt(np.mean(((pc-psi_m)[rel])**2+(dd[rel])**2))
bm=(wl[:,None]<600)&(Wd>0.5)
blue_psi=np.sqrt(np.mean((pc-psi_m)[bm]**2))
print('=== thin damaged GraHIL (%.0fnm) + free large roughness ==='%D_GRAHIL)
print('rough=%.1f nm (perov frac f=%.2f, => %.0f%% void)  bulk=%.1f nm'%(d_rough,f_rough,100*(1-f_rough),d_bulk))
print('perov mass-equiv = bulk+f*rough = %.0f nm'%(d_bulk+f_rough*d_rough))
for tw in [450,500,535,575,600,700,900]:
    print('  %5d n=%.3f k=%.3f'%(tw,np.interp(tw,wl,n),np.interp(tw,wl,k)))
print('max k(620-1600)=%.4f'%max(k[(wl>620)&(wl<1600)]))
print('weighted RMSE=%.2f  reliable(w>0.5)=%.2f  BLUE(<600)Psi=%.2f'%(wmse,rmse_rel,blue_psi))
print('(prev: weighted 30, reliable 43, blue Psi 46)')
with open(str(OUT_DIR / 'Perov_v3_nk.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['wl','n','k'])
    for a,b,c in zip(wl,n,k): w.writerow(['%.3f'%a,'%.5f'%b,'%.5f'%c])
np.savez(str(OUT_DIR / 'Perov_v3_fit.npz'),wl=wl,psi_m=psi_m,del_m=del_m,psi_c=pc,del_c=dc,n=n,k=k,W=Wd,
         d_bulk=d_bulk,d_rough=d_rough,f_rough=f_rough,d_grahil=D_GRAHIL)
print('saved Perov_v3_nk.csv, Perov_v3_fit.npz')
