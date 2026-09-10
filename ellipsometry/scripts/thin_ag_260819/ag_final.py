"""Final Ag model: Einf + Drude + 3 Gaussians over 260-1080 nm.

Residuals were concentrated in 260-350 nm, where Ag's interband structure is,
and two Gaussians could not follow it.  Cutting that band drops MSE much further
(1-16: 9.0 -> 2.1) but lets the thickness drift by +2.5 nm and away from nominal,
so the UV is carrying real thickness information - better to describe it than to
discard it.  A third Gaussian keeps the thickness and improves the fit.
"""
import numpy as np, json
from scipy.optimize import least_squares
import ag_load as L, ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef, ag_fit as AF

R0 = json.load(open('ag_polish_result.json'))
NG = 3
def agN(p, wl, ng=NG):
    E = 1239.841984/wl
    e1, e2 = osc.drude_rt(E, p[3], p[4])
    for i in range(ng):
        a,b = osc.gaussian(E, p[5+3*i], p[6+3*i], p[7+3*i]); e1 += a; e2 += b
    N = np.sqrt((p[2]+e1+1j*e2).astype(complex)); return np.where(N.imag<0,-N,N)

LO = np.array([0.5,0.0,0.5,3e-7,0.4,  0.0,0.15,3.5,  0.0,0.15,1.3,  0.0,0.15,4.2])
HI = np.array([30.,8.0,6.0,1e-1,50.,  30.,4.00,6.5,  30.,3.00,3.5,  30.,4.00,7.5])
NM = ['d_Ag','rough','Einf','rho','tau','G1A','G1Br','G1En','G2A','G2Br','G2En','G3A','G3Br','G3En']

def fit(sh, nstart=30, seed=23):
    name,_,_ = L.SPEC[sh]; d_seed, ps = AF.SEED[name]
    wl,P,D,_ = L.load(sh); m=(wl>=260)&(wl<=1080); wl,P,D=wl[m],P[m],D[m]
    ox,si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    Ns = AF.seedN(ps,wl); amb=np.ones(len(wl)); M=[cf.ncs(P[:,i],D[:,i]) for i in range(5)]
    def r(p):
        Na=agN(p,wl); Nr=ef.bruggeman_ema50(Na,amb); o=[]
        for a,an in enumerate(L.ANG):
            rp,rs=ef._tmm(wl,[amb,Nr,Na,Ns,ox,si],[p[1],p[0],d_seed,AF.D_OX],an)
            rr=rp/rs; q,dl=np.arctan(np.abs(rr)),np.angle(rr)
            o+=[np.cos(2*q)-M[a][0],np.sin(2*q)*np.cos(dl)-M[a][1],np.sin(2*q)*np.sin(dl)-M[a][2]]
        return np.concatenate(o)
    p0 = np.clip(np.array(R0[sh]['p'][:11] + [2.0,1.0,5.5]), LO+1e-9, HI-1e-9)
    rng=np.random.default_rng(seed); best=None
    for k in range(nstart):
        s = p0 if k==0 else (LO+rng.random(len(LO))*(HI-LO) if k%4==3
                             else np.clip(p0*(1+0.3*rng.standard_normal(len(p0))),LO+1e-9,HI-1e-9))
        try: b=least_squares(r,np.clip(s,LO+1e-9,HI-1e-9),bounds=(LO,HI),x_scale='jac',max_nfev=1600)
        except Exception: continue
        if best is None or b.cost<best.cost: best=b
    res=best.fun; n=len(res)//3
    M_=1000*np.sqrt(np.sum(res**2)/(3*n-len(LO)))
    s2=np.sum(res**2)/(len(res)-len(LO)); J=np.zeros((len(res),len(LO)))
    for j in range(len(LO)):
        h=max(abs(best.x[j])*1e-5,1e-9); a=best.x.copy(); a[j]+=h; c=best.x.copy(); c[j]-=h
        J[:,j]=(r(a)-r(c))/(2*h)
    sig=np.sqrt(np.diag(np.linalg.pinv(J.T@J)*s2))*1.645
    return M_, best.x, sig

if __name__=='__main__':
    wp=np.array([450.,550.,633.,800.,1000.])
    out={}
    HB=6.582119569e-16; EPS0=8.8541878128e-12
    for sh in [s for s in L.ORDER if L.SPEC[s][2]>0]:
        M,x,sig = fit(sh); N=agN(x,wp)
        w=np.sqrt(HB**2/(EPS0*(x[3]/100.0)*(x[4]*1e-15)))
        out[sh]=dict(sheet=sh,seed=L.SPEC[sh][0],ag_nom=L.SPEC[sh][2],mse=float(M),
                     p=[float(v) for v in x],err=[float(v) for v in sig],
                     n=[float(v) for v in N.real],k=[float(v) for v in N.imag],wp=float(w))
        print('%-5s %-6s nom=%2d  MSE %5.2f (was %5.2f)  d=%5.2f (was %5.2f)  n633=%.3f  hw_p=%.2f'
              %(sh,L.SPEC[sh][0],L.SPEC[sh][2],M,R0[sh]['mse'],x[0],R0[sh]['p'][0],N.real[2],w),flush=True)
        json.dump(out,open('ag_final_result.json','w'),indent=1)
    print('done')
