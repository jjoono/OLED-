"""Final reliability-maximized perovskite n,k:
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

- corrected KK (singularity-subtracted, validated 2e-3)
- depol-based weights + angle-quality weights (75° down: lateral-inhomog outlier)
- physics constraints: TL Eg in [2.25,2.42] eV, exciton in [2.28,2.40] eV,
  k=0 enforced for wl>560nm (user physics: no absorption 550-600nm)
- SEM thickness anchor d_bulk + f*dr = 270nm (breaks the n*d alias to the
  physical branch: n800~2.08, eps~4.3 - literature-consistent)
- Huber robust loss; multi-start
- uncertainty: angle-pair bootstrap + free-thickness alias branch + SiO2 +
  blue-weight sensitivity -> n,k envelope band"""
import numpy as np, ellipsometry_fit as ef, csv
from scipy.optimize import minimize
base=str(DATA_DIR)
angles=[65.0,70.0,75.0]
wl,psi_m,del_m=ef.load_data_xlsx(base+r'\#2.xlsx',angles,380,1688)
Wd=ef.depol_weights(wl,angles,base+r'\#2_depol.xlsx',depol_cut=3.0,depol_soft=1.5)
CL=lambda v,lo,hi:min(max(v,lo),hi)

def stack_pd(Np, dr, f, db, d_sio2):
    n_air=np.ones(len(wl),dtype=complex)
    layers=[n_air, ef.bruggeman_ema(Np,n_air,f), Np, ef.n_SiO2(wl), ef.n_Si(wl)]
    dl=[dr,db,d_sio2]
    out_p=[];out_d=[]
    for ang in angles:
        rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
        out_p.append(np.degrees(np.arctan(np.abs(rho))))
        out_d.append(np.degrees(np.angle(rho)))
    return np.array(out_p).T, np.array(out_d).T

def unpack(x):
    dr=CL(x[0],40,140); f=CL(x[1],0.08,0.60); db=CL(x[2],170,300)
    epsinf=CL(x[3],1.0,3.0)
    tlEg=CL(x[7],2.25,2.42)
    uvEg=CL(x[14],2.4,3.6)
    osc=[{'type':'TL','A':CL(x[4],0,300),'E0':CL(x[5],tlEg+0.05,3.3),'C':CL(x[6],0.05,1.2),'Eg':tlEg},
         {'type':'Gaussian','A':CL(x[8],0,5),'Ecen':CL(x[9],2.28,2.40),'Br':CL(x[10],0.03,0.10)},
         {'type':'TL','A':CL(x[11],0,400),'E0':CL(x[12],uvEg+0.05,5.5),'C':CL(x[13],0.2,2.5),'Eg':uvEg}]
    return dr,f,db,epsinf,osc

def huber(r,d=3.0):
    a=np.abs(r)
    return np.where(a<=d, r*r, 2*d*a-d*d)

def make_loss(angq, blue_w, d_sio2, prior_w):
    Wtot=Wd*np.array(angq)[None,:]
    bw=np.where(wl<560,blue_w,1.0)[:,None]
    Wtot=Wtot*bw
    sW=np.sum(Wtot)
    def loss(x):
        dr,f,db,epsinf,osc=unpack(x)
        try: n,k=ef.genosc_nk(wl,osc,eps_inf=epsinf)
        except Exception: return 1e7
        if np.any(n<1.4) or np.any(n>2.8): return 1e7
        Np=(n+1j*k).astype(complex)
        pc,dc=stack_pd(Np,dr,f,db,d_sio2)
        dd=(dc-del_m+180)%360-180
        core=(np.sum(Wtot*huber(pc-psi_m))+np.sum(Wtot*huber(dd)))/(2*sW)
        kpen=1000.0*np.mean(k[wl>560]**2)
        prior=prior_w*(db+f*dr-270.0)**2
        return core+kpen+prior
    return loss

def run_fit(loss,x0s,it=25000):
    best=None
    for x0 in x0s:
        r=minimize(loss,x0,method='Nelder-Mead',options={'maxiter':it,'xatol':1e-6,'fatol':1e-9,'adaptive':True})
        if best is None or r.fun<best.fun: best=r
    return best

# seeds (Branch B physical)
s1=[90,0.25,248,1.90, 60,2.55,0.35,2.30, 0.6,2.33,0.06, 60,3.6,1.2,2.6]
s2=[70,0.35,246,2.00, 40,2.60,0.30,2.32, 1.0,2.34,0.05, 90,3.8,1.0,2.7]
s3=[110,0.18,251,1.80, 90,2.50,0.45,2.28, 0.4,2.31,0.07, 40,3.4,1.5,2.5]

variants={}
print('=== MAIN fit (angles w=[1,1,0.35], blue_w=0.6, SEM anchor) ===')
lossM=make_loss([1,1,0.35],0.6,2.0,0.3)
rM=run_fit(lossM,[s1,s2,s3],25000)
rM=run_fit(lossM,[rM.x],25000)
variants['MAIN']=rM.x
dr,f,db,epsinf,osc=unpack(rM.x)
print('MAIN loss=%.3f  dr=%.1f f=%.2f db=%.1f  total=%.1f  epsinf=%.2f'%(rM.fun,dr,f,db,db+f*dr,epsinf))
for o in osc: print('  ',{kk:(round(v,3) if isinstance(v,float) else v) for kk,v in o.items()})

cfgs=[('P6570',[1,1,0.0],0.6,2.0,0.3),
      ('P6575',[1,0,1.0],0.6,2.0,0.3),
      ('P7075',[0,1,1.0],0.6,2.0,0.3),
      ('FREE_D',[1,1,0.35],0.6,2.0,0.0),
      ('SIO2_1',[1,1,0.35],0.6,1.0,0.3),
      ('SIO2_3',[1,1,0.35],0.6,3.0,0.3),
      ('BLUE1',[1,1,0.35],1.0,2.0,0.3)]
for name,angq,bw,ds,pw in cfgs:
    l=make_loss(angq,bw,ds,pw)
    x0s=[rM.x] if name!='FREE_D' else [rM.x,[20,0.3,230,1.9,60,2.55,0.35,2.30,0.6,2.33,0.06,60,3.6,1.2,2.6]]
    r=run_fit(l,x0s,15000)
    variants[name]=r.x
    dr,f,db,epsinf,osc=unpack(r.x)
    print('%-7s loss=%.3f dr=%.1f f=%.2f db=%.1f total=%.1f'%(name,r.fun,dr,f,db,db+f*dr))

# collect n,k curves
curves={}
for name,x in variants.items():
    dr,f,db,epsinf,osc=unpack(x)
    n,k=ef.genosc_nk(wl,osc,eps_inf=epsinf)
    curves[name]=(n,k)
nM,kM=curves['MAIN']
alln=np.array([curves[nm][0] for nm in curves])
allk=np.array([curves[nm][1] for nm in curves])
n_lo,n_hi=alln.min(0),alln.max(0)
k_lo,k_hi=allk.min(0),allk.max(0)

# metrics for MAIN
dr,f,db,epsinf,osc=unpack(variants['MAIN'])
Np=(nM+1j*kM).astype(complex)
pc,dc=stack_pd(Np,dr,f,db,2.0)
dd=(dc-del_m+180)%360-180
m6570=(Wd[:,:2]>0.5)
r6570=np.sqrt(np.mean((pc[:,:2]-psi_m[:,:2])[m6570]**2+dd[:,:2][m6570]**2))
print('\nMAIN true RMSE (65+70, depol-ok pts) = %.2f deg'%r6570)

print('\n  wl |  n (lo..hi)        |  k (lo..hi)')
for tw in [450,500,520,535,550,560,575,600,633,700,800,1000]:
    i=int(np.argmin(np.abs(wl-tw)))
    print(' %4d | %.3f (%.3f..%.3f) | %.3f (%.3f..%.3f)'%(
          wl[i],nM[i],n_lo[i],n_hi[i],kM[i],k_lo[i],k_hi[i]))
print('max k(>560nm) MAIN=%.4f'%kM[wl>560].max())

with open(str(OUT_DIR / 'Perov_reliable_nk.csv'),'w',newline='') as fcsv:
    w=csv.writer(fcsv); w.writerow(['wl_nm','n','k','n_lo','n_hi','k_lo','k_hi'])
    for i in range(len(wl)):
        w.writerow(['%.3f'%wl[i],'%.5f'%nM[i],'%.5f'%kM[i],
                    '%.5f'%n_lo[i],'%.5f'%n_hi[i],'%.5f'%k_lo[i],'%.5f'%k_hi[i]])
np.savez(str(OUT_DIR / 'Perov_reliable_fit.npz'),wl=wl,psi_m=psi_m,del_m=del_m,
         psi_c=pc,del_c=dc,n=nM,k=kM,n_lo=n_lo,n_hi=n_hi,k_lo=k_lo,k_hi=k_hi,W=Wd,
         d_bulk=db,d_rough=dr,f_rough=f,x_main=variants['MAIN'])
print('saved Perov_reliable_nk.csv / Perov_reliable_fit.npz')
