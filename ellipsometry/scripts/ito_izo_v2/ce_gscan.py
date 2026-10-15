"""The third oscillator ran off to En=0.34 eV (3600 nm), far outside the
340-1080 nm data, where only its tail is seen and Amp/Br/En are perfectly
degenerate (CompleteEASE reports Amp3 +-458, En3 +-107).
Profile the fit with En3 FIXED at a series of values inside/near the data and
see what MSE it actually costs to keep the oscillator where it is measurable."""
import numpy as np, json, sys
from scipy.optimize import least_squares
import ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef, tmm_fix

V = json.load(open('ce_final2_result.json'))

def NfB(p, wl):
    E = 1239.841984/wl
    e1, e2 = osc.drude_rt(E, p[2], p[3])
    a, b = osc.tauc_lorentz(E, p[4], p[5], p[6], p[7]); e1 += a; e2 += b
    a, b = osc.gaussian(E, p[8], p[9], p[10]);          e1 += a; e2 += b
    N = np.sqrt((p[1] + e1 + 1j*e2).astype(complex))
    return np.where(N.imag < 0, -N, N)

NM = ['rough','Einf','rho','tau','TL_Amp','TL_Br','TL_Eo','TL_Eg','G_Amp','G_Br','G_En']
LO = np.array([0.0, 1.0, 5e-5, 1.0,  10., 0.10, 3.3, 2.60, 0.0, 0.10, 0.30])
HI = np.array([12., 5.0, 5e-1, 30.0, 1500., 8.0, 6.5, 4.10, 10., 3.00, 6.00])

def run(sh, en_fix, nstart=12, seed=23):
    d, dth = V[sh]['d'], V[sh]['dth']
    wl, Pm, Dm = cf.load(sh)
    m = (wl>=340)&(wl<=1080); wl,Pm,Dm = wl[m],Pm[m],Dm[m]
    ox, si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    M = [cf.ncs(Pm[:,i],Dm[:,i]) for i in range(3)]; amb = np.ones(len(wl))
    free = [1,2,3,4,5,6,7,8,9] if en_fix else list(range(1,11))
    p_full = np.array(V[sh]['p']).copy(); p_full[0] = max(p_full[0], 0.0)
    if en_fix: p_full[10] = en_fix
    def full(x):
        q = p_full.copy(); q[free] = x; return q
    def r(x):
        p = full(x); Nf = NfB(p, wl); Nr = ef.bruggeman_ema50(Nf, amb); out=[]
        for a, an in enumerate(cf.ANG0+dth):
            rp, rs = tmm_fix.tmm(wl, [amb,Nr,Nf,ox,si], [p[0], d, cf.D_OX], an)
            rr = np.conj(rp/rs); ps, dl = np.arctan(np.abs(rr)), np.angle(rr)
            out += [np.cos(2*ps)-M[a][0], np.sin(2*ps)*np.cos(dl)-M[a][1],
                    np.sin(2*ps)*np.sin(dl)-M[a][2]]
        return np.concatenate(out)
    lo, hi = LO[free], HI[free]
    x0 = np.clip(p_full[free], lo+1e-9, hi-1e-9)
    rng = np.random.default_rng(seed); best=None
    for k in range(nstart):
        s = x0 if k==0 else np.clip(x0*(1+0.3*rng.standard_normal(len(x0))), lo+1e-9, hi-1e-9)
        try: b = least_squares(r, s, bounds=(lo,hi), x_scale='jac', max_nfev=1200)
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
    return cf.mse(res, len(free)+2), p, rel

sh = sys.argv[1] if len(sys.argv)>1 else '#1'
print('%s : third oscillator centre fixed at En3 =' % sh)
print('  En3(eV)  lambda(nm)    MSE     worst param uncertainty')
for en in [None, 0.8, 1.0, 1.2, 1.5, 1.8, 2.2, 2.8]:
    M, p, rel = run(sh, en)
    worst = max(rel.items(), key=lambda z: z[1])
    tag = 'free -> %.3f' % p[10] if en is None else '%.2f' % en
    lam = 1239.841984/p[10]
    print('  %-12s %7.0f  %7.3f    %s %.0f%%' % (tag, lam, M, worst[0], worst[1]))
