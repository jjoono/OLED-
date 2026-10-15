"""Seed layers: HATCN (1-5) and MoOx (1-13) on Si + native oxide.

Both are only 4-5 nm, so n and d are strongly correlated; the native oxide adds
a third correlated term.  So scan the assumed oxide thickness and report the
whole (d_seed, n550) family rather than a single number.
Dispersion: Einf + Tauc-Lorentz (KK consistent, k >= 0 by construction).
"""
import numpy as np, json
from scipy.optimize import least_squares
import ag_load as L, ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef

LOW, HIGH = 260.0, 1080.0

def seedN(p, wl):          # p = [Einf, TL_Amp, TL_Br, TL_Eo, TL_Eg]
    E = 1239.841984/wl
    e1, e2 = osc.tauc_lorentz(E, p[1], p[2], p[3], p[4])
    N = np.sqrt((p[0] + e1 + 1j*e2).astype(complex))
    return np.where(N.imag < 0, -N, N)

LO = np.array([1.0,  1.0, 0.10, 3.5, 2.5])
HI = np.array([4.0, 400., 6.00, 9.0, 5.5])

def fit(sheet, d_ox, nstart=20, seed=5, dlo=1.0, dhi=25.0):
    wl, P, D, _ = L.load(sheet)
    m = (wl>=LOW)&(wl<=HIGH); wl,P,D = wl[m],P[m],D[m]
    ox, si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    M = [cf.ncs(P[:,i],D[:,i]) for i in range(5)]; amb = np.ones(len(wl))
    def r(x):
        Nf = seedN(x[1:], wl); o=[]
        for a,an in enumerate(L.ANG):
            rp,rs = ef._tmm(wl,[amb,Nf,ox,si],[x[0],d_ox],an)
            rr=rp/rs; ps,dl=np.arctan(np.abs(rr)),np.angle(rr)
            o += [np.cos(2*ps)-M[a][0], np.sin(2*ps)*np.cos(dl)-M[a][1],
                  np.sin(2*ps)*np.sin(dl)-M[a][2]]
        return np.concatenate(o)
    lo = np.r_[dlo, LO]; hi = np.r_[dhi, HI]
    p0 = np.array([6.0, 2.2, 60., 1.5, 5.5, 3.5])
    rng = np.random.default_rng(seed); best=None
    for k in range(nstart):
        s = p0 if k==0 else (lo+rng.random(len(lo))*(hi-lo) if k%4==3
                             else np.clip(p0*(1+0.3*rng.standard_normal(len(p0))), lo+1e-9, hi-1e-9))
        try: b = least_squares(r, np.clip(s,lo+1e-9,hi-1e-9), bounds=(lo,hi), x_scale='jac', max_nfev=1500)
        except Exception: continue
        if best is None or b.cost<best.cost: best=b
    n=len(best.fun)//3
    return 1000*np.sqrt(np.sum(best.fun**2)/(3*n-len(lo))), best.x

if __name__ == '__main__':
    out={}
    for sheet,name in [('1-5','HATCN'),('1-13','MoOx')]:
        print('=== %s (%s, nominal %d nm) ==='%(sheet,name,L.SPEC[sheet][1]))
        print('  d_ox(nm)  MSE    d_seed(nm)   n550    n633    k350   Einf   TL(Amp,Br,Eo,Eg)')
        rec=[]
        for dox in [1.0, 1.5, 2.0, 2.5, 3.0]:
            M,x = fit(sheet, dox)
            N = seedN(x[1:], np.array([350.,550.,633.]))
            print('  %5.1f   %6.2f  %8.2f   %6.3f  %6.3f  %6.3f  %5.2f  %.1f %.2f %.2f %.2f'
                  %(dox,M,x[0],N.real[1],N.real[2],N.imag[0],x[1],x[2],x[3],x[4],x[5]))
            rec.append(dict(d_ox=dox,mse=float(M),p=[float(v) for v in x],
                            n550=float(N.real[1]),n633=float(N.real[2])))
        out[name]=rec
        json.dump(out, open('ag_seed_result.json','w'), indent=1)
