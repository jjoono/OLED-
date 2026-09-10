"""How much of the Ag thickness offset comes from my fixed seed thickness?
Refit the Ag layer with the seed held at several values and watch d_Ag move."""
import numpy as np, json
from scipy.optimize import least_squares
import ag_load as L, ce_fit as cf, ellipsometry_fit as ef, ag_fit as AF, ag_final as AG

R = json.load(open('ag_final_result.json'))

def refit(sh, d_seed, nstart=8, seed=31):
    name,_,_ = L.SPEC[sh]; _, ps = AF.SEED[name]
    wl,P,D,_ = L.load(sh); m=(wl>=260)&(wl<=1080); wl,P,D=wl[m],P[m],D[m]
    ox,si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    Ns = AF.seedN(ps,wl); amb=np.ones(len(wl)); M=[cf.ncs(P[:,i],D[:,i]) for i in range(5)]
    def r(p):
        Na=AG.agN(p,wl); Nr=ef.bruggeman_ema50(Na,amb); o=[]
        for a,an in enumerate(L.ANG):
            rp,rs=ef._tmm(wl,[amb,Nr,Na,Ns,ox,si],[p[1],p[0],d_seed,AF.D_OX],an)
            rr=rp/rs; q,dl=np.arctan(np.abs(rr)),np.angle(rr)
            o+=[np.cos(2*q)-M[a][0],np.sin(2*q)*np.cos(dl)-M[a][1],np.sin(2*q)*np.sin(dl)-M[a][2]]
        return np.concatenate(o)
    p0=np.clip(np.array(R[sh]['p']),AG.LO+1e-9,AG.HI-1e-9)
    rng=np.random.default_rng(seed); best=None
    for k in range(nstart):
        s=p0 if k==0 else np.clip(p0*(1+0.2*rng.standard_normal(len(p0))),AG.LO+1e-9,AG.HI-1e-9)
        try: b=least_squares(r,s,bounds=(AG.LO,AG.HI),x_scale='jac',max_nfev=1000)
        except Exception: continue
        if best is None or b.cost<best.cost: best=b
    n=len(best.fun)//3
    return 1000*np.sqrt(np.sum(best.fun**2)/(3*n-len(AG.LO))), best.x

HB=6.582119569e-16; EPS0=8.8541878128e-12
wpf=lambda p: np.sqrt(HB**2/(EPS0*(p[3]/100.0)*(p[4]*1e-15)))
REF=max(wpf(np.array(R[s]['p'])) for s in R if R[s]['seed']=='HATCN')
print('HATCN series, Ag thickness vs the ASSUMED seed thickness')
print('  seed=4.5 nm        seed=6.42 (used)     seed=8.5 nm')
print('%-5s %3s | %6s %6s | %6s %6s | %6s %6s'%('id','nom','d_Ag','MSE','d_Ag','MSE','d_Ag','MSE'))
tot={}
for sh in ['1-6','1-7','1-8','2-5','2-6','2-8']:
    row=[]; 
    for ds in (4.5, 6.42, 8.5):
        M,x = refit(sh, ds); row += [x[0], M]
        tot.setdefault(ds,[]).append(((x[0]+x[1]/2)*(wpf(x)/REF)**2)/L.SPEC[sh][2])
    print('%-5s %3d | %6.2f %6.2f | %6.2f %6.2f | %6.2f %6.2f'%(sh,L.SPEC[sh][2],*row), flush=True)
print()
for ds,v in tot.items():
    print('  seed=%.2f nm -> mass-equivalent / nominal = %.2f +- %.2f'%(ds,np.mean(v),np.std(v)))
