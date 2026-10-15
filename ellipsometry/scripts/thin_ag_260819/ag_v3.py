"""v3: find the underlayer thickness each seed series actually wants, then refit.

The bare-seed samples put the underlayer (native oxide + seed) at 8.4 nm for
HATCN, but with that value the fitted Ag thicknesses run 12-22% over the
deposited values.  Pinning Ag at nominal and scanning the seed gives a clear
minimum at 7.5 nm (underlayer 9.5), and freeing Ag there lands all 7 samples at
1.010 +- 0.052 of nominal with equal or better MSE.  Do the same for MoOx, then
refit everything at each series' own optimum.
"""
import numpy as np, json
from scipy.optimize import least_squares
import ag_load as L, ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef, ag_final as AG

S=json.load(open('ag_seed_result.json')); R=json.load(open('ag_final_result.json'))
SP={k:[r for r in S[k] if abs(r['d_ox']-2.0)<1e-9][0]['p'][1:] for k in ('HATCN','MoOx')}
HB=6.582119569e-16; EPS0=8.8541878128e-12
wpf=lambda p: np.sqrt(HB**2/(EPS0*(p[3]/100.0)*(p[4]*1e-15)))

def seedN(p,wl):
    E=1239.841984/wl; e1,e2=osc.tauc_lorentz(E,p[1],p[2],p[3],p[4])
    N=np.sqrt((p[0]+e1+1j*e2).astype(complex)); return np.where(N.imag<0,-N,N)

def go(sh,d_seed,fix_d,nstart=8,seed=71):
    name=L.SPEC[sh][0]
    wl,P,D,_=L.load(sh); m=(wl>=260)&(wl<=1080); wl,P,D=wl[m],P[m],D[m]
    ox,si=cf._mat('NTVE_JAW',wl),cf._mat('SI_JAW',wl)
    Ns=seedN(SP[name],wl); amb=np.ones(len(wl)); M=[cf.ncs(P[:,i],D[:,i]) for i in range(5)]
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

OPT={}
for grp,grid in [('HATCN',[7.50]), ('MoOx',[6.93,7.5,8.0,8.5,9.0,9.5])]:
    sheets=[s for s in L.ORDER if L.SPEC[s][0]==grp and L.SPEC[s][2]>0]
    print('%s: Ag fixed at deposited value, seed scan'%grp)
    best=None
    for dsd in grid:
        ms=[go(sh,dsd,True)[0] for sh in sheets]
        mm=np.mean(ms)
        if best is None or mm<best[0]: best=(mm,dsd)
        print('  d_seed=%.2f (underlayer %.2f)  mean MSE %.2f   [%s]'%(
            dsd,dsd+2.0,mm,' '.join('%.2f'%v for v in ms)),flush=True)
    OPT[grp]=best[1]
    print('  -> %s optimum seed = %.2f nm\n'%(grp,best[1]),flush=True)
    json.dump(OPT,open('ag_v3_seed.json','w'))

out={}
print('final refit with Ag FREE at each optimum')
wp=np.array([450.,550.,633.,800.,1000.])
for sh in [s for s in L.ORDER if L.SPEC[s][2]>0]:
    grp=L.SPEC[sh][0]
    M,p=go(sh,OPT[grp],False,nstart=12)
    N=AG.agN(p,wp); w=wpf(p)
    out[sh]=dict(sheet=sh,seed=grp,ag_nom=L.SPEC[sh][2],mse=float(M),d_seed=OPT[grp],
                 p=[float(v) for v in p],n=[float(v) for v in N.real],
                 k=[float(v) for v in N.imag],wp=float(w))
    print('  %-5s %-6s nom=%2d  d=%5.2f (%.2f x nom)  MSE=%5.2f  n633=%.3f  hw_p=%.2f'
          %(sh,grp,L.SPEC[sh][2],p[0],p[0]/L.SPEC[sh][2],M,N.real[2],w),flush=True)
    json.dump(out,open('ag_v3_result.json','w'),indent=1)
print('done')
