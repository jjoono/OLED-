"""Ag on HATCN / MoOx, 14 samples.

Stack: air / roughness(Bruggeman Ag+void) / Ag / seed / native oxide 2 nm / Si
The seed is held at the dispersion and thickness fitted from the bare-seed
samples 1-5 and 1-13 (oxide = 2.0 nm branch).
Ag dispersion: Einf + Drude + Gaussian(interband, UV) + Gaussian(visible).
The second Gaussian is what lets a DISCONTINUOUS island film be described - a
localized plasmon sits in the visible, which is the percolation signature.
"""
import numpy as np, json
from scipy.optimize import least_squares
import ag_load as L, ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef

LOW, HIGH = 260.0, 1080.0
D_OX = 2.0
SEED = {'HATCN': (6.42, [2.65, 19.5, 0.54, 3.95, 2.54]),
        'MoOx':  (6.93, [2.58, 44.6, 2.37, 4.42, 2.99])}

def seedN(p, wl):
    E = 1239.841984/wl
    e1, e2 = osc.tauc_lorentz(E, p[1], p[2], p[3], p[4])
    N = np.sqrt((p[0] + e1 + 1j*e2).astype(complex))
    return np.where(N.imag < 0, -N, N)

#              d_Ag  rough  Einf  rho     tau   G1A  G1Br G1En  G2A  G2Br G2En
LO = np.array([0.5,  0.0,   0.5, 3e-7,   0.4,  0.0, 0.15, 3.5,  0.0, 0.15, 1.3])
HI = np.array([30.,  8.0,   6.0, 1e-1,  50.0, 30.0, 4.00, 6.5, 30.0, 3.00, 3.5])

def agN(p, wl):
    E = 1239.841984/wl
    e1, e2 = osc.drude_rt(E, p[3], p[4])
    a, b = osc.gaussian(E, p[5], p[6], p[7]); e1 += a; e2 += b
    a, b = osc.gaussian(E, p[8], p[9], p[10]); e1 += a; e2 += b
    N = np.sqrt((p[2] + e1 + 1j*e2).astype(complex))
    return np.where(N.imag < 0, -N, N)

def fit(sheet, nstart=34, seed=9):
    name, dseed_nom, dag_nom = L.SPEC[sheet]
    d_seed, ps = SEED[name]
    wl, P, D, _ = L.load(sheet)
    m = (wl>=LOW)&(wl<=HIGH); wl,P,D = wl[m],P[m],D[m]
    ox, si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    Ns = seedN(ps, wl)
    M = [cf.ncs(P[:,i],D[:,i]) for i in range(5)]; amb = np.ones(len(wl))
    def r(p):
        Na = agN(p, wl); Nr = ef.bruggeman_ema50(Na, amb); o=[]
        for a,an in enumerate(L.ANG):
            rp,rs = ef._tmm(wl,[amb,Nr,Na,Ns,ox,si],[p[1],p[0],d_seed,D_OX],an)
            rr=rp/rs; ps_,dl=np.arctan(np.abs(rr)),np.angle(rr)
            o += [np.cos(2*ps_)-M[a][0], np.sin(2*ps_)*np.cos(dl)-M[a][1],
                  np.sin(2*ps_)*np.sin(dl)-M[a][2]]
        return np.concatenate(o)
    p0 = np.array([float(dag_nom), 1.0, 3.0, 5e-6, 12., 3.0, 1.5, 4.6, 1.0, 0.8, 2.5])
    p0 = np.clip(p0, LO+1e-9, HI-1e-9)
    rng = np.random.default_rng(seed); best=None
    for k in range(nstart):
        s = p0 if k==0 else (LO+rng.random(len(LO))*(HI-LO) if k%4==3
                             else np.clip(p0*(1+0.35*rng.standard_normal(len(p0))), LO+1e-9, HI-1e-9))
        s[0] = np.clip(dag_nom*(1+0.4*rng.standard_normal()), 0.5, 30.) if k else s[0]
        try: b = least_squares(r, np.clip(s,LO+1e-9,HI-1e-9), bounds=(LO,HI), x_scale='jac', max_nfev=1800)
        except Exception: continue
        if best is None or b.cost<best.cost: best=b
    n = len(best.fun)//3
    return 1000*np.sqrt(np.sum(best.fun**2)/(3*n-len(LO))), best.x, r

if __name__ == '__main__':
    import sys
    sheets = sys.argv[1:] if len(sys.argv)>1 else [s for s in L.ORDER if L.SPEC[s][2]>0]
    out = {}
    try: out = json.load(open('ag_fit_result.json'))
    except Exception: pass
    wp = np.array([450.,550.,633.,800.,1000.])
    for sh in sheets:
        M, x, _ = fit(sh)
        N = agN(x, wp)
        out[sh] = dict(sheet=sh, seed=L.SPEC[sh][0], ag_nom=L.SPEC[sh][2],
                       mse=float(M), p=[float(v) for v in x],
                       n=[float(v) for v in N.real], k=[float(v) for v in N.imag])
        print('%-5s %-6s Ag_nom=%2d  d_Ag=%5.2f rough=%4.2f  MSE=%6.2f | n633=%.3f k633=%.3f | rho=%.3g tau=%.1f'
              %(sh, L.SPEC[sh][0], L.SPEC[sh][2], x[0], x[1], M, N.real[2], N.imag[2], x[3], x[4]), flush=True)
        json.dump(out, open('ag_fit_result.json','w'), indent=1)
    print('done')
