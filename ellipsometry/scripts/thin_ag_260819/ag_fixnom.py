"""The decisive form of the question: FIX every Ag thickness at its deposited
(nominal) value and scan the total underlayer thickness (native oxide + seed).
Is there an underlayer where nominal thicknesses fit as well as free ones?

Only the SUM matters - swapping 1 nm of oxide for 1 nm of seed left the Ag
thickness unchanged to 0.02 nm, because the two are optically interchangeable.
So d_ox is held at 2.0 and the seed carries the scan.
"""
import numpy as np, json
from scipy.optimize import least_squares
import ag_load as L, ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef, ag_final as AG

S=json.load(open('ag_seed_result.json')); R=json.load(open('ag_final_result.json'))
SP = [r for r in S['HATCN'] if abs(r['d_ox']-2.0)<1e-9][0]['p'][1:]   # seed dispersion
def seedN(p,wl):
    E=1239.841984/wl; e1,e2=osc.tauc_lorentz(E,p[1],p[2],p[3],p[4])
    N=np.sqrt((p[0]+e1+1j*e2).astype(complex)); return np.where(N.imag<0,-N,N)

def run(sh, d_seed, d_ag, nstart=6, seed=53):
    wl,P,D,_=L.load(sh); m=(wl>=260)&(wl<=1080); wl,P,D=wl[m],P[m],D[m]
    ox,si=cf._mat('NTVE_JAW',wl),cf._mat('SI_JAW',wl)
    Ns=seedN(SP,wl); amb=np.ones(len(wl)); M=[cf.ncs(P[:,i],D[:,i]) for i in range(5)]
    free=list(range(1,14))                       # everything except d_Ag
    p_all=np.array(R[sh]['p']).copy(); p_all[0]=d_ag
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
        try: b=least_squares(r,s,bounds=(lo,hi),x_scale='jac',max_nfev=800)
        except Exception: continue
        if best is None or b.cost<best.cost: best=b
    n=len(best.fun)//3
    return 1000*np.sqrt(np.sum(best.fun**2)/(3*n-len(free)-1))

sheets=[s for s in L.ORDER if L.SPEC[s][0]=='HATCN' and L.SPEC[s][2]>0]
print('Ag thickness FIXED at the deposited value; scanning the seed thickness')
print('(oxide fixed 2.0 nm; bare-seed sample says seed = 6.42 nm)')
print('%-8s %-8s | %s | %7s'%('d_seed','under.',' '.join('%5d'%L.SPEC[s][2] for s in sheets),'meanMSE'))
for dsd in [6.42, 7.5, 8.5, 9.5, 10.5]:
    ms=[run(sh,dsd,float(L.SPEC[sh][2])) for sh in sheets]
    print('%-8.2f %-8.2f | %s | %7.2f'%(dsd,dsd+2.0,' '.join('%5.2f'%v for v in ms),np.mean(ms)),flush=True)
print()
print('reference: same samples with d_Ag FREE at seed 6.42 -> mean MSE %.2f'
      % np.mean([R[s]['mse'] for s in sheets]))
