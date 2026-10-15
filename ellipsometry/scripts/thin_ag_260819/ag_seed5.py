"""QCM says both seeds are 5 nm.  Thickness and index are degenerate for a layer
this thin - only (n^2-1)*d is really measured - so ask the decisive question:
if d is FORCED to 5 nm, what index does the data then demand, and is that index
physically possible for HATCN / MoOx?
  HATCN (organic, literature n550 ~ 1.75-1.80)
  MoOx  (literature n550 ~ 2.0-2.2)
"""
import numpy as np, json
from scipy.optimize import least_squares
import ag_load as L, ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef

D_OX = 2.0
def seedN(p, wl):
    E=1239.841984/wl; e1,e2 = osc.tauc_lorentz(E,p[1],p[2],p[3],p[4])
    N=np.sqrt((p[0]+e1+1j*e2).astype(complex)); return np.where(N.imag<0,-N,N)
LO=np.array([1.0,  1.0,0.10,3.5,2.5]); HI=np.array([6.0, 400.,6.00,9.0,5.5])

def fit(sheet, d_fix, nstart=16, seed=13):
    wl,P,D,_ = L.load(sheet); m=(wl>=260)&(wl<=1080); wl,P,D=wl[m],P[m],D[m]
    ox,si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    M=[cf.ncs(P[:,i],D[:,i]) for i in range(5)]; amb=np.ones(len(wl))
    def r(x):
        Nf=seedN(x,wl); o=[]
        for a,an in enumerate(L.ANG):
            rp,rs=ef._tmm(wl,[amb,Nf,ox,si],[d_fix,D_OX],an)
            rr=rp/rs; q,dl=np.arctan(np.abs(rr)),np.angle(rr)
            o+=[np.cos(2*q)-M[a][0],np.sin(2*q)*np.cos(dl)-M[a][1],np.sin(2*q)*np.sin(dl)-M[a][2]]
        return np.concatenate(o)
    p0=np.array([2.2,60.,1.5,5.5,3.5]); rng=np.random.default_rng(seed); best=None
    for k in range(nstart):
        s=p0 if k==0 else (LO+rng.random(len(LO))*(HI-LO) if k%4==3
                           else np.clip(p0*(1+0.3*rng.standard_normal(len(p0))),LO+1e-9,HI-1e-9))
        try: b=least_squares(r,np.clip(s,LO+1e-9,HI-1e-9),bounds=(LO,HI),x_scale='jac',max_nfev=1200)
        except Exception: continue
        if best is None or b.cost<best.cost: best=b
    n=len(best.fun)//3
    N=seedN(best.x,np.array([550.,633.]))
    return 1000*np.sqrt(np.sum(best.fun**2)/(3*n-len(LO))), N.real[0], N.real[1], best.x

for sheet,name,lit in [('1-5','HATCN','1.75-1.80 (organic)'),('1-13','MoOx','2.0-2.2')]:
    print('=== %s  (QCM 5 nm; literature n550 %s) ==='%(name,lit))
    print('  d fixed   MSE    n550    n633   (n^2-1)*d')
    for d in [4.5, 5.0, 5.5, 6.0, 6.42, 7.0, 7.5, 8.0]:
        M,n5,n6,x = fit(sheet,d)
        print('  %5.2f  %6.2f  %6.3f  %6.3f  %8.2f'%(d,M,n5,n6,(n5**2-1)*d), flush=True)
    print()
