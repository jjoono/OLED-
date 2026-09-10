"""Final model per sample at the corrected geometry.

The 3-oscillator fits pushed the Gaussian down to En=0.6 eV, where it merely
duplicates the Drude, so try Drude+TL alone as well and keep whichever wins.
Then lock whatever the covariance says is unidentifiable, so the confidence
intervals CompleteEASE reports mean something.
"""
import numpy as np, json
from scipy.optimize import least_squares
import ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef, tmm_fix

V = json.load(open('ce_v12_result.json'))

def NfB(p, wl, ng):                 # p=[rough,Einf,rho,tau,TLa,TLbr,TLeo,TLeg,(Ga,Gbr,Gen)]
    E = 1239.841984/wl
    e1, e2 = osc.drude_rt(E, p[2], p[3])
    a, b = osc.tauc_lorentz(E, p[4], p[5], p[6], p[7]); e1 += a; e2 += b
    if ng: 
        a, b = osc.gaussian(E, p[8], p[9], p[10]); e1 += a; e2 += b
    N = np.sqrt((p[1] + e1 + 1j*e2).astype(complex))
    return np.where(N.imag < 0, -N, N)

LO2 = np.array([0.0, 1.0, 5e-5, 1.0,  10., 0.10, 3.3, 2.60])
HI2 = np.array([12., 5.0, 5e-1, 30.0, 1500., 8.0, 6.5, 4.10])
LO3 = np.r_[LO2, [0.0, 0.10, 0.30]]
HI3 = np.r_[HI2, [10., 8.00, 5.00]]
NM  = ['rough','Einf','rho','tau','TL_Amp','TL_Br','TL_Eo','TL_Eg','G_Amp','G_Br','G_En']

def build(sh, ng, nstart=26, seed=17):
    d, dth = V[sh]['d'], V[sh]['dth']
    wl, Pm, Dm = cf.load(sh)
    m = (wl>=340)&(wl<=1080); wl,Pm,Dm = wl[m],Pm[m],Dm[m]
    ox, si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    M = [cf.ncs(Pm[:,i],Dm[:,i]) for i in range(3)]
    amb = np.ones(len(wl))
    def r(p):
        Nf = NfB(p, wl, ng); Nr = ef.bruggeman_ema50(Nf, amb); out=[]
        for a, an in enumerate(cf.ANG0+dth):
            rp, rs = tmm_fix.tmm(wl, [amb,Nr,Nf,ox,si], [p[0], d, cf.D_OX], an)
            rr = np.conj(rp/rs); ps, dl = np.arctan(np.abs(rr)), np.angle(rr)
            out += [np.cos(2*ps)-M[a][0], np.sin(2*ps)*np.cos(dl)-M[a][1],
                    np.sin(2*ps)*np.sin(dl)-M[a][2]]
        return np.concatenate(out)
    lo, hi = (LO3, HI3) if ng else (LO2, HI2)
    p0 = np.array(V[sh]['p'])[:len(lo)]
    rng = np.random.default_rng(seed); best=None
    for k in range(nstart):
        s = p0 if k==0 else (lo+rng.random(len(lo))*(hi-lo) if k%5==4
                             else np.clip(p0*(1+0.3*rng.standard_normal(len(p0))), lo+1e-9, hi-1e-9))
        try: b = least_squares(r, np.clip(s,lo+1e-9,hi-1e-9), bounds=(lo,hi), x_scale='jac', max_nfev=1600)
        except Exception: continue
        if best is None or b.cost<best.cost: best=b
    return best, r, d, dth, lo, hi

if __name__ == '__main__':
    out={}
    for sh in ['#1','#2','#3','#4','#5']:
        cand=[]
        for ng in (0,1):
            b, r, d, dth, lo, hi = build(sh, ng)
            cand.append((cf.mse(b.fun, len(lo)+2), ng, b, r, d, dth, lo, hi))
        cand.sort(key=lambda z: z[0])
        M, ng, b, r, d, dth, lo, hi = cand[0]
        p = b.x
        res = r(p); s2 = np.sum(res**2)/(len(res)-len(p))
        J = np.zeros((len(res), len(p)))
        for j in range(len(p)):
            h=max(abs(p[j])*1e-5,1e-9); a=p.copy(); a[j]+=h; c=p.copy(); c[j]-=h
            J[:,j]=(r(a)-r(c))/(2*h)
        sig = np.sqrt(np.diag(np.linalg.pinv(J.T@J)*s2))*1.645
        rel = 100*sig/np.maximum(np.abs(p),1e-12)
        out[sh]=dict(sheet=sh, d=d, dth=dth, ng=ng, mse=float(M), p=[float(v) for v in p],
                     err=[float(v) for v in sig], rel=[float(v) for v in rel],
                     mse_2osc=float(cand[1][0] if cand[1][1]!=ng else M),
                     unident=[NM[j] for j in range(len(p)) if rel[j]>50])
        print('%s  d=%.1f dth=%+.2f  %d-osc  MSE=%5.2f  (other option %5.2f)'
              %(sh,d,dth,2+ng,M,cand[1][0]), flush=True)
        print('    '+'  '.join('%s=%.4g(%.0f%%)'%(NM[j],p[j],rel[j]) for j in range(len(p))), flush=True)
        json.dump(out, open('ce_final2_result.json','w'), indent=1)
    print('done')
