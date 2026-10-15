"""Refine the underlayer thickness, then check the +-10% question there.

Fixing every Ag thickness at its deposited value and scanning the seed gave
mean MSE 3.14 (seed 6.42) -> 2.53 (7.50) -> 3.43 (8.50): a real minimum near
7.5 nm, where nominal thicknesses fit as well as free ones did before (2.46).
Bracket it more finely, then free d_Ag at the optimum and see where it lands.
"""
import numpy as np, json
from scipy.optimize import least_squares
import ag_load as L, ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef, ag_final as AG

S=json.load(open('ag_seed_result.json')); R=json.load(open('ag_final_result.json'))
SP=[r for r in S['HATCN'] if abs(r['d_ox']-2.0)<1e-9][0]['p'][1:]
HB=6.582119569e-16; EPS0=8.8541878128e-12
wpf=lambda p: np.sqrt(HB**2/(EPS0*(p[3]/100.0)*(p[4]*1e-15)))
def seedN(p,wl):
    E=1239.841984/wl; e1,e2=osc.tauc_lorentz(E,p[1],p[2],p[3],p[4])
    N=np.sqrt((p[0]+e1+1j*e2).astype(complex)); return np.where(N.imag<0,-N,N)

def go(sh,d_seed,fix_d,nstart=6,seed=61):
    wl,P,D,_=L.load(sh); m=(wl>=260)&(wl<=1080); wl,P,D=wl[m],P[m],D[m]
    ox,si=cf._mat('NTVE_JAW',wl),cf._mat('SI_JAW',wl)
    Ns=seedN(SP,wl); amb=np.ones(len(wl)); M=[cf.ncs(P[:,i],D[:,i]) for i in range(5)]
    free=list(range(1,14)) if fix_d else list(range(0,14))
    p_all=np.array(R[sh]['p']).copy()
    if fix_d: p_all[0]=float(L.SPEC[sh][2])
    def r(x):
        p=p_all.copy(); p[free]=x
        Na=AG.agN(p,wl); Nr=ef.bruggeman_ema50(Na,amb); o=[]
        for a,an in enumerate(L.ANG):
            rp,rs=ef._tmm(wl,[amb,Nr,Na,Ns,ox,si],[p[1],p[0],d_seed,2.0],an)
            rr=rp/rs; q,dl=np.arctan(np.abs(rr)),np.angle(rr)
            o+=[np.cos(2*q)-M[a][0],np.sin(2*q)*np.cos(dl)-M[a][1],np.sin(2*q)*np.sin(dl)-M[a][2]]
        return np.concatenate(o)
    lo,hi=AG.LO[free],AG.HI[free]; x0=np.clip(p_all[free],lo+1e-9,hi-1e-9)
    rng=np.random.default_rng(seed); best=None
    for k in range(nstart):
        s=x0 if k==0 else np.clip(x0*(1+0.2*rng.standard_normal(len(x0))),lo+1e-9,hi-1e-9)
        try: b=least_squares(r,s,bounds=(lo,hi),x_scale='jac',max_nfev=900)
        except Exception: continue
        if best is None or b.cost<best.cost: best=b
    p=p_all.copy(); p[free]=best.x; n=len(best.fun)//3
    return 1000*np.sqrt(np.sum(best.fun**2)/(3*n-len(free)-1)), p

sheets=[s for s in L.ORDER if L.SPEC[s][0]=='HATCN' and L.SPEC[s][2]>0]
print('Ag FIXED at deposited value, finer seed scan')
print('%-8s %-8s | %s | %7s'%('d_seed','under.',' '.join('%5d'%L.SPEC[s][2] for s in sheets),'meanMSE'))
best=None
for dsd in [7.0, 7.25, 7.5, 7.75, 8.0]:
    ms=[go(sh,dsd,True)[0] for sh in sheets]
    mm=np.mean(ms)
    if best is None or mm<best[0]: best=(mm,dsd)
    print('%-8.2f %-8.2f | %s | %7.2f'%(dsd,dsd+2.0,' '.join('%5.2f'%v for v in ms),mm),flush=True)
print('\nbest seed = %.2f nm (underlayer %.2f nm), mean MSE %.2f'%(best[1],best[1]+2.0,best[0]))
print('\nnow FREE the Ag thickness at that seed:')
print('%-5s %4s | %6s %6s %7s %6s'%('id','nom','d_Ag','d/nom','MSE','f_Ag'))
rat=[]
for sh in sheets:
    M,p = go(sh,best[1],False)
    f=(wpf(p)/9.29)**2; rat.append(p[0]/L.SPEC[sh][2])
    print('%-5s %4d | %6.2f %6.2f %7.2f %6.2f'%(sh,L.SPEC[sh][2],p[0],rat[-1],M,f),flush=True)
rat=np.array(rat)
print('  d/nominal = %.3f +- %.3f ; within +-10%%: %d/%d'%(rat.mean(),rat.std(),
      int(np.sum(np.abs(rat-1)<=0.10)),len(rat)))
