"""v13: keep the third oscillator INSIDE the measured window.

With G_Br allowed up to 9 eV the optimiser escapes to En3 = 0.34 eV (3600 nm),
far outside 340-1080 nm, where only the tail is seen: CompleteEASE then reports
Amp3 +-458 and En3 +-107 even though MSE is 1.83.  Capping G_Br at 3 eV pulls it
back to En3 ~ 2.3 eV (530 nm, a plausible ITO sub-gap defect band) and every
parameter becomes determined to a few percent, for dMSE ~ +0.2.
"""
import numpy as np, json
from scipy.optimize import least_squares
import ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef, tmm_fix

V = json.load(open('ce_final2_result.json'))
NM = ['rough','Einf','rho','tau','TL_Amp','TL_Br','TL_Eo','TL_Eg','G_Amp','G_Br','G_En']
LO = np.array([0.0, 1.0, 5e-5, 1.0,  10., 0.10, 3.3, 2.60, 0.0, 0.20, 1.20])
HI = np.array([12., 5.0, 5e-1, 30.0, 1500., 6.0, 6.5, 4.10, 10., 3.00, 3.60])

def NfB(p, wl):
    E = 1239.841984/wl
    e1, e2 = osc.drude_rt(E, p[2], p[3])
    a, b = osc.tauc_lorentz(E, p[4], p[5], p[6], p[7]); e1 += a; e2 += b
    a, b = osc.gaussian(E, p[8], p[9], p[10]);          e1 += a; e2 += b
    N = np.sqrt((p[1] + e1 + 1j*e2).astype(complex))
    return np.where(N.imag < 0, -N, N)

out = {}
for sh in ['#1','#2','#3','#4','#5']:
    d, dth = V[sh]['d'], V[sh]['dth']
    wl, Pm, Dm = cf.load(sh)
    m = (wl>=340)&(wl<=1080); wl,Pm,Dm = wl[m],Pm[m],Dm[m]
    ox, si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    M = [cf.ncs(Pm[:,i],Dm[:,i]) for i in range(3)]; amb = np.ones(len(wl))
    def r(p):
        Nf = NfB(p, wl); Nr = ef.bruggeman_ema50(Nf, amb); o=[]
        for a, an in enumerate(cf.ANG0+dth):
            rp, rs = tmm_fix.tmm(wl, [amb,Nr,Nf,ox,si], [p[0], d, cf.D_OX], an)
            rr = np.conj(rp/rs); ps, dl = np.arctan(np.abs(rr)), np.angle(rr)
            o += [np.cos(2*ps)-M[a][0], np.sin(2*ps)*np.cos(dl)-M[a][1],
                  np.sin(2*ps)*np.sin(dl)-M[a][2]]
        return np.concatenate(o)
    p0 = np.clip(np.array(V[sh]['p']), LO+1e-9, HI-1e-9)
    p0[10] = 2.3; p0[9] = min(p0[9], 2.0)
    rng = np.random.default_rng(31); best=None
    for k in range(30):
        s = p0 if k==0 else (LO+rng.random(len(LO))*(HI-LO) if k%5==4
                             else np.clip(p0*(1+0.3*rng.standard_normal(len(p0))), LO+1e-9, HI-1e-9))
        try: b = least_squares(r, np.clip(s,LO+1e-9,HI-1e-9), bounds=(LO,HI), x_scale='jac', max_nfev=1600)
        except Exception: continue
        if best is None or b.cost<best.cost: best=b
    p = best.x; res = best.fun
    s2 = np.sum(res**2)/(len(res)-len(p))
    J = np.zeros((len(res), len(p)))
    for j in range(len(p)):
        h=max(abs(p[j])*1e-5,1e-9); a=p.copy(); a[j]+=h; c=p.copy(); c[j]-=h
        J[:,j]=(r(a)-r(c))/(2*h)
    sig = np.sqrt(np.diag(np.linalg.pinv(J.T@J)*s2))*1.645
    rel = 100*sig/np.maximum(np.abs(p),1e-12)
    out[sh]=dict(sheet=sh,d=d,dth=dth,ng=1,mse=float(cf.mse(res,13)),
                 p=[float(v) for v in p],err=[float(v) for v in sig],rel=[float(v) for v in rel],
                 atb=[NM[j] for j in range(len(p))
                      if abs(p[j]-LO[j])<1e-6*max(1,abs(LO[j])) or abs(p[j]-HI[j])<1e-6*max(1,abs(HI[j]))])
    print('%s  MSE=%5.2f  d=%.1f dth=%+.2f  worst uncertainty %.0f%%  at-bound:%s'
          %(sh,out[sh]['mse'],d,dth,rel[1:].max(),out[sh]['atb'] or '-'),flush=True)
    print('    '+'  '.join('%s=%.4g±%.2g(%.0f%%)'%(NM[j],p[j],sig[j],rel[j]) for j in range(1,len(p))),flush=True)
    json.dump(out, open('ce_v13_result.json','w'), indent=1)
print('done')
