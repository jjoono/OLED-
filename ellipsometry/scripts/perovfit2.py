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

D_GRAHIL=75.0   # SEM/paper physical GraHIL thickness (fixed)
DR=80.0         # roughness EMA (fixed, from consistency scan)

def stack_psidelta(wlq, Np, d_bulk, f_grahil_frac):
    # Air / rough-EMA(perov+void 50/50) / perov bulk / damaged-GraHIL(EMA GraHIL+perov, 75nm) / SiO2 / Si
    n_air=np.ones(len(wlq),dtype=complex)
    Ngr=Ng(wlq)
    Ndam=ef.bruggeman_ema(Ngr, Np, f_grahil_frac)   # f_grahil_frac = GraHIL vol-fraction; low=heavily damaged
    layers=[n_air, ef.bruggeman_ema(Np,n_air,0.5), Np, Ndam, ef.n_SiO2(wlq), ef.n_Si(wlq)]
    d_list=[DR, d_bulk, D_GRAHIL, 2.0]
    po,do=[],[]
    for ang in angles:
        rp,rs=ef._tmm(wlq,layers,d_list,ang); rho=rp/rs
        po.append(np.degrees(np.arctan(np.abs(rho)))); do.append(np.degrees(np.angle(rho)))
    return np.array(po).T, np.array(do).T

def oscs(p):
    epsinf=max(p[0],1.0)
    Eg=min(max(p[4],2.05),2.45); uEg=min(max(p[11],2.4),3.6)
    L=[{'type':'TL','A':min(max(p[1],0),300),'E0':min(max(p[2],Eg+0.05),3.3),'C':min(max(p[3],0.05),1.5),'Eg':Eg},
       {'type':'Gaussian','A':min(max(p[5],0),5),'Ecen':min(max(p[6],2.20),2.45),'Br':min(max(p[7],0.03),0.09)},
       {'type':'TL','A':min(max(p[8],0),400),'E0':min(max(p[9],uEg+0.05),5.5),'C':min(max(p[10],0.2),2.5),'Eg':uEg}]
    return epsinf,L

def resid(x):
    d_bulk=min(max(x[0],120),260); fgr=min(max(x[1],0.0),1.0)
    epsinf,L=oscs(x[2:])
    try: n,k=ef.genosc_nk(wl,L,eps_inf=epsinf)
    except Exception: return 1e6
    if np.any(n<1.5) or np.any(n>2.8): return 1e6
    Np=(n+1j*k).astype(complex)
    pc,dc=stack_psidelta(wl,Np,d_bulk,fgr)
    dd=(dc-del_m+180)%360-180
    mse=(np.sum(Wd*(pc-psi_m)**2)+np.sum(Wd*dd**2))/(2*np.sum(Wd))
    ktr=80.0*np.mean((k[wl>620])**2)
    prior=0.02*(d_bulk+0.5*DR+0.3*D_GRAHIL-270)**2  # soft: perov(bulk+½rough)≈270 above GraHIL
    return mse+ktr+prior

# seed oscillators from prior fit; bulk seed ~190 (so perov total ~230), fgr seed 0.3 (damaged)
osc0=[1.0,50,2.42,0.42,2.05, 1.47,2.24,0.09, 38,4.9,0.2,2.4]
best=None
for db0 in [170,195,220]:
  for fgr0 in [0.15,0.4,0.7]:
    r=minimize(resid,[db0,fgr0]+osc0,method='Nelder-Mead',options={'maxiter':40000,'xatol':1e-6,'fatol':1e-9,'adaptive':True})
    if best is None or r.fun<best[0]: best=(r.fun,r.x)
r=minimize(resid,best[1],method='Nelder-Mead',options={'maxiter':60000,'xatol':1e-7,'fatol':1e-10,'adaptive':True})
x=r.x
d_bulk=min(max(x[0],120),260); fgr=min(max(x[1],0.0),1.0)
epsinf,L=oscs(x[2:])
n,k=ef.genosc_nk(wl,L,eps_inf=epsinf); Np=(n+1j*k).astype(complex)
pc,dc=stack_psidelta(wl,Np,d_bulk,fgr)
dd=(dc-del_m+180)%360-180
rmse=np.sqrt(np.mean((pc-psi_m)**2+dd**2))
wmse=np.sqrt((np.sum(Wd*(pc-psi_m)**2)+np.sum(Wd*dd**2))/(2*np.sum(Wd)))
rel=Wd>0.5
rmse_rel=np.sqrt(np.mean(((pc-psi_m)[rel])**2+(dd[rel])**2))
Nd=ef.bruggeman_ema(Ng(np.array([633.0])),np.array([np.interp(633,wl,n)+1j*np.interp(633,wl,k)]),fgr)
print('=== GraHIL PRESENT (fixed %.0fnm) as damaged intermix EMA ==='%D_GRAHIL)
print('perov bulk=%.1f nm  rough=%.0f nm  GraHIL=%.0f nm (fixed)  perov total=%.0f'%(d_bulk,DR,D_GRAHIL,d_bulk+0.5*DR))
print('GraHIL vol-fraction in damaged layer f_GraHIL=%.2f (=> %.0f%% perovskite intermixed)'%(fgr,100*(1-fgr)))
print('damaged-GraHIL index @633nm: n=%.2f (pristine GraHIL 1.42, perovskite ~2.15)'%Nd[0].real)
for tw in [450,500,535,575,600,700,900]:
    print('  %5d n=%.3f k=%.3f'%(tw,np.interp(tw,wl,n),np.interp(tw,wl,k)))
print('max k(620-1600)=%.4f'%max(k[(wl>620)&(wl<1600)]))
print('RMSE all=%.2f  weighted=%.2f  reliable(w>0.5)=%.2f'%(rmse,wmse,rmse_rel))
print('(compare: GraHIL=0 model gave weighted 31, reliable 43)')
with open(str(OUT_DIR / 'Perov_grahil75_nk.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['wl','n','k'])
    for a,b,c in zip(wl,n,k): w.writerow(['%.3f'%a,'%.5f'%b,'%.5f'%c])
np.savez(str(OUT_DIR / 'Perov_grahil75_fit.npz'),wl=wl,psi_m=psi_m,del_m=del_m,psi_c=pc,del_c=dc,n=n,k=k,W=Wd,
         d_bulk=d_bulk,d_rough=DR,d_grahil=D_GRAHIL,fgr=fgr)
print('saved Perov_grahil75_nk.csv, Perov_grahil75_fit.npz')
