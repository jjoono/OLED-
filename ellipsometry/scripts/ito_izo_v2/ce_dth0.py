"""If the angle offset is forced to 0, what does it cost and how much do n,k move?"""
import numpy as np, json
from scipy.optimize import least_squares
import ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef, tmm_fix, ce_v11 as v11

V = json.load(open('ce_v13_result.json'))
LO = np.array([0.0, 1.0, 5e-5, 1.0,  10., 0.10, 3.3, 2.60, 0.0, 0.20, 1.20])
HI = np.array([12., 5.0, 5e-1, 30.0, 1500., 6.0, 6.5, 4.10, 10., 3.00, 3.60])
def NfB(p, wl):
    E = 1239.841984/wl
    e1,e2 = osc.drude_rt(E,p[2],p[3])
    a,b = osc.tauc_lorentz(E,p[4],p[5],p[6],p[7]); e1+=a; e2+=b
    a,b = osc.gaussian(E,p[8],p[9],p[10]);         e1+=a; e2+=b
    N=np.sqrt((p[1]+e1+1j*e2).astype(complex)); return np.where(N.imag<0,-N,N)

def go(sh, dth, dscan):
    wl,Pm,Dm = cf.load(sh); m=(wl>=340)&(wl<=1080); wl,Pm,Dm=wl[m],Pm[m],Dm[m]
    ox,si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    # re-optimise thickness for this angle offset (they are aliased)
    bd = min(((np.sqrt(np.mean(v11.chain(wl,Pm,Dm,ox,si,d,0.0,dth,stride=14)[:,3]**2)), d) for d in dscan))
    d = bd[1]
    M=[cf.ncs(Pm[:,i],Dm[:,i]) for i in range(3)]; amb=np.ones(len(wl))
    def r(p):
        Nf=NfB(p,wl); Nr=ef.bruggeman_ema50(Nf,amb); o=[]
        for a,an in enumerate(cf.ANG0+dth):
            rp,rs=tmm_fix.tmm(wl,[amb,Nr,Nf,ox,si],[p[0],d,cf.D_OX],an)
            rr=np.conj(rp/rs); ps,dl=np.arctan(np.abs(rr)),np.angle(rr)
            o+=[np.cos(2*ps)-M[a][0],np.sin(2*ps)*np.cos(dl)-M[a][1],np.sin(2*ps)*np.sin(dl)-M[a][2]]
        return np.concatenate(o)
    p0=np.clip(np.array(V[sh]['p']),LO+1e-9,HI-1e-9)
    rng=np.random.default_rng(7); best=None
    for k in range(16):
        s=p0 if k==0 else np.clip(p0*(1+0.25*rng.standard_normal(len(p0))),LO+1e-9,HI-1e-9)
        try: b=least_squares(r,s,bounds=(LO,HI),x_scale='jac',max_nfev=1400)
        except Exception: continue
        if best is None or b.cost<best.cost: best=b
    return cf.mse(best.fun,13), d, best.x, bd[0]

wp=np.array([450.,550.,633.,800.,1000.])
print('%-5s %-22s %7s %7s %8s   %s'%('','angle offset','d(nm)','MSE','chainRMS','n / k at 450,550,633,800,1000 nm'))
for sh in ['#1','#2','#5']:
    ref=None
    for dth,tag in [(V[sh]['dth'],'fitted %+.2f deg'%V[sh]['dth']),(0.0,'forced 0.00 deg')]:
        M,d,p,cr = go(sh,dth,[float(x) for x in range(50,62)])
        N=NfB(p,wp)
        if ref is None: ref=N
        print('%-5s %-22s %7.1f %7.2f %8.5f   n '%(sh,tag,d,M,cr)+' '.join('%.3f'%v for v in N.real))
        print('%-5s %-22s %7s %7s %8s   k '%('','','','','')+' '.join('%.4f'%v for v in N.imag)
              + ('   <- dn %.3f dk %.4f'%(np.abs(N.real-ref.real).max(),np.abs(N.imag-ref.imag).max()) if dth==0.0 else ''))
