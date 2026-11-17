"""Does starting at 370 nm instead of 340 nm matter?
Refit each sample on both windows and compare parameters, uncertainties and n,k."""
import numpy as np, json
from scipy.optimize import least_squares
import ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef, tmm_fix

V = json.load(open('ce_v13_result.json'))
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

def fit(sh, lo_w, lock_drude=False):
    d, dth = V[sh]['d'], V[sh]['dth']
    wl, Pm, Dm = cf.load(sh)
    m = (wl>=lo_w)&(wl<=1080); wl,Pm,Dm = wl[m],Pm[m],Dm[m]
    ox, si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    M = [cf.ncs(Pm[:,i],Dm[:,i]) for i in range(3)]; amb = np.ones(len(wl))
    p_all = np.array(V[sh]['p'])
    free = [j for j in range(1,11) if not (lock_drude and j in (2,3))]
    def full(x):
        q = p_all.copy(); q[free] = x; return q
    def r(x):
        p = full(x); Nf = NfB(p, wl); Nr = ef.bruggeman_ema50(Nf, amb); o=[]
        for a, an in enumerate(cf.ANG0+dth):
            rp, rs = tmm_fix.tmm(wl, [amb,Nr,Nf,ox,si], [p[0], d, cf.D_OX], an)
            rr = np.conj(rp/rs); ps, dl = np.arctan(np.abs(rr)), np.angle(rr)
            o += [np.cos(2*ps)-M[a][0], np.sin(2*ps)*np.cos(dl)-M[a][1],
                  np.sin(2*ps)*np.sin(dl)-M[a][2]]
        return np.concatenate(o)
    lo, hi = LO[free], HI[free]
    x0 = np.clip(p_all[free], lo+1e-9, hi-1e-9)
    rng = np.random.default_rng(41); best=None
    for k in range(14):
        s = x0 if k==0 else np.clip(x0*(1+0.25*rng.standard_normal(len(x0))), lo+1e-9, hi-1e-9)
        try: b = least_squares(r, s, bounds=(lo,hi), x_scale='jac', max_nfev=1300)
        except Exception: continue
        if best is None or b.cost<best.cost: best=b
    p = full(best.x); res = best.fun
    s2 = np.sum(res**2)/(len(res)-len(free))
    J = np.zeros((len(res), len(free)))
    for j in range(len(free)):
        h=max(abs(best.x[j])*1e-5,1e-9); a=best.x.copy(); a[j]+=h; c=best.x.copy(); c[j]-=h
        J[:,j]=(r(a)-r(c))/(2*h)
    sig=np.sqrt(np.diag(np.linalg.pinv(J.T@J)*s2))*1.645
    rel={NM[free[j]]:100*sig[j]/max(abs(best.x[j]),1e-12) for j in range(len(free))}
    return cf.mse(res,len(free)+2), p, rel

wp = np.array([400.,450.,550.,633.,800.,1000.])
print('window comparison (all other settings identical)')
for sh in ['#1','#2','#3','#4','#5']:
    ld = (sh=='#5')
    out={}
    for lo_w in (340., 370., 400.):
        M,p,rel = fit(sh, lo_w, ld)
        out[lo_w]=(M,p,rel)
    p340 = out[340.][1]; N340 = NfB(p340, wp)
    print('\n%s' % sh)
    print('  window     MSE   TL_Eg(+-)        TL_Amp(+-)      worst%   max|dn| max|dk| vs 340nm fit')
    for lo_w in (340., 370., 400.):
        M,p,rel = out[lo_w]; N = NfB(p, wp)
        print('  %3.0f-1080 %6.2f   %.3f(%.0f%%)      %.1f(%.0f%%)      %4.0f%%   %.4f  %.4f'
              % (lo_w, M, p[7], rel['TL_Eg'], p[4], rel['TL_Amp'],
                 max(rel.values()), np.abs(N.real-N340.real).max(), np.abs(N.imag-N340.imag).max()))
