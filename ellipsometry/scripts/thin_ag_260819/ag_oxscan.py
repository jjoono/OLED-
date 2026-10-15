"""Is the native oxide the culprit?

d_ox and d_seed trade off almost exactly in the bare-seed fit (sum ~8.3 nm for
HATCN).  A thinner oxide means a thicker seed, which pushes the fitted Ag
thickness DOWN.  So scan the oxide, take the seed thickness AND dispersion that
the bare sample gives for that oxide, refit every Ag sample, and ask: is there
an oxide value where all Ag thicknesses land within +-10% of the deposited
(nominal) value at the same time?
"""
import numpy as np, json
from scipy.optimize import least_squares
import ag_load as L, ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef, ag_final as AG

S = json.load(open('ag_seed_result.json'))
R = json.load(open('ag_final_result.json'))
HB=6.582119569e-16; EPS0=8.8541878128e-12
wpf = lambda p: np.sqrt(HB**2/(EPS0*(p[3]/100.0)*(p[4]*1e-15)))

def seedN(p, wl):
    E=1239.841984/wl
    e1,e2 = osc.tauc_lorentz(E, p[1], p[2], p[3], p[4])
    N=np.sqrt((p[0]+e1+1j*e2).astype(complex)); return np.where(N.imag<0,-N,N)

def refit(sh, d_ox, d_seed, sp, nstart=6, seed=41):
    wl,P,D,_ = L.load(sh); m=(wl>=260)&(wl<=1080); wl,P,D=wl[m],P[m],D[m]
    ox,si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    Ns = seedN(sp, wl); amb=np.ones(len(wl)); M=[cf.ncs(P[:,i],D[:,i]) for i in range(5)]
    def r(p):
        Na=AG.agN(p,wl); Nr=ef.bruggeman_ema50(Na,amb); o=[]
        for a,an in enumerate(L.ANG):
            rp,rs=ef._tmm(wl,[amb,Nr,Na,Ns,ox,si],[p[1],p[0],d_seed,d_ox],an)
            rr=rp/rs; q,dl=np.arctan(np.abs(rr)),np.angle(rr)
            o+=[np.cos(2*q)-M[a][0],np.sin(2*q)*np.cos(dl)-M[a][1],np.sin(2*q)*np.sin(dl)-M[a][2]]
        return np.concatenate(o)
    p0=np.clip(np.array(R[sh]['p']),AG.LO+1e-9,AG.HI-1e-9)
    rng=np.random.default_rng(seed); best=None
    for k in range(nstart):
        s=p0 if k==0 else np.clip(p0*(1+0.18*rng.standard_normal(len(p0))),AG.LO+1e-9,AG.HI-1e-9)
        try: b=least_squares(r,s,bounds=(AG.LO,AG.HI),x_scale='jac',max_nfev=900)
        except Exception: continue
        if best is None or b.cost<best.cost: best=b
    n=len(best.fun)//3
    return 1000*np.sqrt(np.sum(best.fun**2)/(3*n-len(AG.LO))), best.x

GRP='HATCN'
sheets=[s for s in L.ORDER if L.SPEC[s][0]==GRP and L.SPEC[s][2]>0]
print('\n%-9s %-9s | %s | %6s %6s %s'%('d_ox','d_seed',
      ' '.join('%5d'%L.SPEC[s][2] for s in sheets),'meanMSE','within','mass/nom'))
for rec in S[GRP]:
    dox=rec['d_ox']; dsd=rec['p'][0]; sp=rec['p'][1:]
    ds=[]; ms=[]; mass=[]
    for sh in sheets:
        M,x = refit(sh,dox,dsd,sp); ds.append(x[0]); ms.append(M)
        mass.append((x[0]+x[1]/2)*(wpf(x)/9.29)**2/L.SPEC[sh][2])
    rat=np.array(ds)/np.array([L.SPEC[s][2] for s in sheets])
    within=int(np.sum(np.abs(rat-1)<=0.10))
    print('%-9.1f %-9.2f | %s | %6.2f %4d/%d  %.2f+-%.2f'%(dox,dsd,
          ' '.join('%5.2f'%v for v in ds), np.mean(ms), within, len(sheets),
          np.mean(mass), np.std(mass)), flush=True)
    print('%-9s %-9s | %s | (d/nominal)'%('','',' '.join('%5.2f'%v for v in rat)), flush=True)
