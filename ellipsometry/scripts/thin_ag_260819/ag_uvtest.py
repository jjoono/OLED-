"""The residual sits in 260-350 nm, where Ag's interband structure lives.
Two remedies: (a) cut the deep UV, (b) add a third oscillator to describe it.
Check what each costs and whether the visible n,k moves."""
import numpy as np, json
from scipy.optimize import least_squares
import ag_load as L, ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef, ag_fit as AF

R = json.load(open('ag_polish_result.json'))

def agN3(p, wl, ng):
    E = 1239.841984/wl
    e1, e2 = osc.drude_rt(E, p[3], p[4])
    for i in range(ng):
        a,b = osc.gaussian(E, p[5+3*i], p[6+3*i], p[7+3*i]); e1 += a; e2 += b
    N = np.sqrt((p[2]+e1+1j*e2).astype(complex)); return np.where(N.imag<0,-N,N)

def bounds(ng):
    lo=[0.5,0.0,0.5,3e-7,0.4]; hi=[30.,8.0,6.0,1e-1,50.]
    for i in range(ng):
        lo += [0.0,0.15,1.3 if i==1 else (3.5 if i==0 else 4.2)]
        hi += [30.,4.00,3.5 if i==1 else (6.5 if i==0 else 7.5)]
    return np.array(lo), np.array(hi)

def run(sh, lo_w, ng, nstart=22, seed=17):
    name,_,_ = L.SPEC[sh]; d_seed, ps = AF.SEED[name]
    wl,P,D,_ = L.load(sh); m=(wl>=lo_w)&(wl<=1080); wl,P,D=wl[m],P[m],D[m]
    ox,si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    Ns = AF.seedN(ps,wl); amb=np.ones(len(wl)); M=[cf.ncs(P[:,i],D[:,i]) for i in range(5)]
    def r(p):
        Na=agN3(p,wl,ng); Nr=ef.bruggeman_ema50(Na,amb); o=[]
        for a,an in enumerate(L.ANG):
            rp,rs=ef._tmm(wl,[amb,Nr,Na,Ns,ox,si],[p[1],p[0],d_seed,AF.D_OX],an)
            rr=rp/rs; q,dl=np.arctan(np.abs(rr)),np.angle(rr)
            o+=[np.cos(2*q)-M[a][0],np.sin(2*q)*np.cos(dl)-M[a][1],np.sin(2*q)*np.sin(dl)-M[a][2]]
        return np.concatenate(o)
    LO,HI = bounds(ng)
    p0 = np.array(R[sh]['p'][:11].copy() + ([2.0,1.0,5.5] if ng==3 else []))
    p0 = np.clip(p0, LO+1e-9, HI-1e-9)
    rng=np.random.default_rng(seed); best=None
    for k in range(nstart):
        s = p0 if k==0 else (LO+rng.random(len(LO))*(HI-LO) if k%4==3
                             else np.clip(p0*(1+0.3*rng.standard_normal(len(p0))),LO+1e-9,HI-1e-9))
        try: b=least_squares(r,np.clip(s,LO+1e-9,HI-1e-9),bounds=(LO,HI),x_scale='jac',max_nfev=1400)
        except Exception: continue
        if best is None or b.cost<best.cost: best=b
    n=len(best.fun)//3
    return 1000*np.sqrt(np.sum(best.fun**2)/(3*n-len(LO))), best.x

wp=np.array([450.,633.,800.])
print('%-5s %-6s %3s | %-22s %-22s %-22s'%('id','seed','nom','260-1080, 2 gauss','350-1080, 2 gauss','260-1080, 3 gauss'))
for sh in ['1-6','1-16','2-16','2-8']:
    base=R[sh]['mse']; out=['%6.2f  n633=%.3f'%(base, R[sh]['n'][2])]
    for lo_w,ng in [(350.,2),(260.,3)]:
        M,x = run(sh, lo_w, ng)
        N = agN3(x, wp, ng)
        out.append('%6.2f  n633=%.3f d=%.2f'%(M, N.real[1], x[0]))
    print('%-5s %-6s %3d | %-22s %-22s %-22s'%(sh,L.SPEC[sh][0],L.SPEC[sh][2],*out), flush=True)
