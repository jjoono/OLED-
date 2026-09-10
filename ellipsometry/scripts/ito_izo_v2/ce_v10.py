"""Refit with the CORRECTED transfer matrix.

tmm_fix.tmm has the physically right propagation sign (a 200 nm Ag film now
returns the bulk Fresnel reflectance instead of R>1).  That flips the sign of
Delta, which is a pure convention difference - the instrument reports the
conjugate - so rho is conjugated to put it back in the measured convention.
Question this answers: does the bug change the extracted n,k?
"""
import numpy as np, json
from scipy.optimize import least_squares
import ce_fit as cf, ce_fit2 as f2, ellipsometry_fit as ef, tmm_fix

FREE = [1, 3, 6, 9, 10, 11, 12]
NM = {1:'Roughness',3:'Einf',6:'Amp2(TL)',9:'Eg2(TL)',10:'Amp3(G)',11:'Br3(G)',12:'En3(G)'}
LO = {1:0.0, 3:1.0, 6:20., 9:3.00, 10:0.0, 11:0.10, 12:0.60}
HI = {1:12., 3:4.5, 6:600., 9:3.95, 10:8.0, 11:5.00, 12:3.40}

def resid_fixed(wl, Pm, Dm, ox, si):
    M = [cf.ncs(Pm[:, i], Dm[:, i]) for i in range(3)]
    amb = np.ones(len(wl))
    def r(p):
        Nf = f2.film_N(p, wl); Nr = ef.bruggeman_ema50(Nf, amb)
        out = []
        for a, an in enumerate(cf.ANG0 + p[2]):
            rp, rs = tmm_fix.tmm(wl, [amb, Nr, Nf, ox, si], [p[1], p[0], cf.D_OX], an)
            rr = np.conj(rp/rs)                    # instrument convention
            ps, dl = np.arctan(np.abs(rr)), np.angle(rr)
            out += [np.cos(2*ps)-M[a][0], np.sin(2*ps)*np.cos(dl)-M[a][1],
                    np.sin(2*ps)*np.sin(dl)-M[a][2]]
        return np.concatenate(out)
    return r

V9 = json.load(open('ce_v9_result.json'))
out = {}
for sh in ['#1','#2','#3','#4','#5']:
    p = np.array(V9[sh]['p'])
    wl, Pm, Dm = cf.load(sh)
    m = (wl>=340)&(wl<=1080); wl,Pm,Dm = wl[m],Pm[m],Dm[m]
    ox, si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    r0 = resid_fixed(wl, Pm, Dm, ox, si)
    mse_before = cf.mse(r0(p), 9)
    lo = np.array([LO[i] for i in FREE]); hi = np.array([HI[i] for i in FREE])
    def r(x):
        q = p.copy(); q[FREE] = x; return r0(q)
    x0 = np.clip(p[FREE], lo+1e-9, hi-1e-9)
    best=None; rng=np.random.default_rng(5)
    for k in range(12):
        s = x0 if k==0 else np.clip(x0*(1+0.3*rng.standard_normal(len(x0))), lo+1e-9, hi-1e-9)
        try: b=least_squares(r, s, bounds=(lo,hi), x_scale='jac', max_nfev=1200)
        except Exception: continue
        if best is None or b.cost<best.cost: best=b
    q = p.copy(); q[FREE] = best.x
    res = r0(q); M2 = cf.mse(res, len(FREE)+2)
    s2 = np.sum(res**2)/(len(res)-len(FREE))
    J = np.zeros((len(res), len(FREE)))
    for j,idx in enumerate(FREE):
        h=max(abs(q[idx])*1e-5,1e-9); a=q.copy(); a[idx]+=h; bb=q.copy(); bb[idx]-=h
        J[:,j]=(r0(a)-r0(bb))/(2*h)
    sig=np.sqrt(np.diag(np.linalg.pinv(J.T@J)*s2))*1.645
    wp=np.linspace(400,1000,7)
    n9=f2.film_N(p,wp); n10=f2.film_N(q,wp)
    out[sh]=dict(sheet=sh,p=q.tolist(),mse=float(M2),mse_v9=V9[sh]['mse'],
                 err={NM[i]:float(sig[j]) for j,i in enumerate(FREE)})
    print('%s  MSE %6.2f  (old-TMM fit was %6.2f; that solution scored %7.1f under the fixed TMM)'
          %(sh,M2,V9[sh]['mse'],mse_before), flush=True)
    print('    dn max %.4f  dk max %.4f  over 400-1000 nm  |  n550 %.3f -> %.3f   k550 %.4f -> %.4f'
          %(np.abs(n10.real-n9.real).max(), np.abs(n10.imag-n9.imag).max(),
            n9.real[np.argmin(abs(wp-550))], n10.real[np.argmin(abs(wp-550))],
            n9.imag[np.argmin(abs(wp-550))], n10.imag[np.argmin(abs(wp-550))]), flush=True)
    print('    '+'  '.join('%s=%.4g+-%.2g'%(NM[i],q[i],sig[j]) for j,i in enumerate(FREE)), flush=True)
    json.dump(out, open('ce_v10_result.json','w'), indent=1)
print('done')
