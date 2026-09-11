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
D,DR,DI=230.0,80.0,0.0

def oscs(p):
    epsinf=max(p[0],1.0)
    Eg=min(max(p[4],2.05),2.45); uEg=min(max(p[11],2.4),3.6)
    L=[{'type':'TL','A':min(max(p[1],0),300),'E0':min(max(p[2],Eg+0.05),3.3),'C':min(max(p[3],0.05),1.5),'Eg':Eg},
       {'type':'Gaussian','A':min(max(p[5],0),5),'Ecen':min(max(p[6],2.20),2.45),'Br':min(max(p[7],0.03),0.09)},
       {'type':'TL','A':min(max(p[8],0),400),'E0':min(max(p[9],uEg+0.05),5.5),'C':min(max(p[10],0.2),2.5),'Eg':uEg}]
    return epsinf,L

def resid(p):
    epsinf,L=oscs(p)
    try: n,k=ef.genosc_nk(wl,L,eps_inf=epsinf)
    except Exception: return 1e6
    if np.any(n<1.5) or np.any(n>2.8): return 1e6
    Np=(n+1j*k).astype(complex)
    pc,dc=ef.calc_psi_delta_perov(wl,Np,Ng,DR,D,DI,0.0,2,angles,f_damage=0.5)
    dd=(dc-del_m+180)%360-180
    mse=(np.sum(Wd*(pc-psi_m)**2)+np.sum(Wd*dd**2))/(2*np.sum(Wd))
    ktr=80.0*np.mean((k[wl>620])**2)
    return mse+ktr

# p = [epsinf,TLA,TLE0,TLC,TLEg, excA,excEcen,excBr, uvA,uvE0,uvC,uvEg]
seeds=[[2.0,40,2.45,0.25,2.28, 0.4,2.33,0.05, 60,3.6,1.2,2.5],
       [2.3,60,2.40,0.20,2.25, 0.6,2.30,0.05, 80,3.8,1.0,2.6],
       [1.8,30,2.50,0.30,2.30, 0.3,2.35,0.06, 40,3.5,1.5,2.5]]
best=None
for s in seeds:
    r=minimize(resid,s,method='Nelder-Mead',options={'maxiter':40000,'xatol':1e-6,'fatol':1e-9,'adaptive':True})
    r=minimize(resid,r.x,method='Nelder-Mead',options={'maxiter':40000,'xatol':1e-7,'fatol':1e-10,'adaptive':True})
    if best is None or r.fun<best[0]: best=(r.fun,r.x)
epsinf,L=oscs(best[1])
n,k=ef.genosc_nk(wl,L,eps_inf=epsinf)
Np=(n+1j*k).astype(complex)
pc,dc=ef.calc_psi_delta_perov(wl,Np,Ng,DR,D,DI,0.0,2,angles,f_damage=0.5)
dd=(dc-del_m+180)%360-180
rmse=np.sqrt(np.mean((pc-psi_m)**2+dd**2))
wmse=np.sqrt((np.sum(Wd*(pc-psi_m)**2)+np.sum(Wd*dd**2))/(2*np.sum(Wd)))
rel=wl<1000
rmse_rel=np.sqrt(np.mean((pc[rel]-psi_m[rel])**2+dd[rel]**2))
print('thickness fixed: bulk=%.0f rough=%.0f (total=%.0f)'%(D,DR,D+0.5*DR))
print('eps_inf=%.3f'%epsinf)
for o in L: print('  ',{kk:(round(vv,3) if isinstance(vv,float) else vv) for kk,vv in o.items()})
for tw in [400,450,500,535,550,580,600,650,700,900,1200]:
    print('  %5d n=%.3f k=%.3f'%(tw,np.interp(tw,wl,n),np.interp(tw,wl,k)))
print('max k(620-1600)=%.4f'%max(k[(wl>620)&(wl<1600)]))
print('RMSE all=%.2f  weighted=%.2f  reliable(<1000)=%.2f'%(rmse,wmse,rmse_rel))
with open(str(OUT_DIR / 'Perov_final_nk.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['wl','n','k'])
    for a,b,c in zip(wl,n,k): w.writerow(['%.3f'%a,'%.5f'%b,'%.5f'%c])
np.savez(str(OUT_DIR / 'Perov_final_fit.npz'),wl=wl,psi_m=psi_m,del_m=del_m,psi_c=pc,del_c=dc,n=n,k=k,W=Wd)
print('saved Perov_final_nk.csv, Perov_final_fit.npz')
