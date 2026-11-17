"""Full redo with the CORRECTED transfer matrix (tmm_fix), 340-1080 nm.
  1  geometry scan (d, roughness, angle offset) scored by model-free chain residual
  2  chain n,k at the best geometry
  3  Einf + Drude + Tauc-Lorentz + Gaussian fitted to the raw NCS at locked geometry
"""
import numpy as np, json, time
from scipy.optimize import least_squares
import ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef, tmm_fix

LOW, HIGH = 340.0, 1080.0
D_OX = cf.D_OX

def rho_model(wl, Nf, rough, d, ox, si, ang):
    amb = np.ones(len(wl)); Nr = ef.bruggeman_ema50(Nf, amb)
    rp, rs = tmm_fix.tmm(wl, [amb, Nr, Nf, ox, si], [rough, d, D_OX], ang)
    return np.conj(rp/rs)

def chain(wl, Pm, Dm, ox, si, d, rough, dth, stride=1):
    M = [cf.ncs(Pm[:, i], Dm[:, i]) for i in range(3)]
    ang = cf.ANG0 + dth
    idx = np.arange(len(wl)-1, -1, -stride)
    out = np.empty((len(idx), 4)); g = np.array([2.0, 0.02])
    for t, i in enumerate(idx):
        w1 = wl[i:i+1]; o1, s1 = ox[i:i+1], si[i:i+1]
        meas = np.array([[M[a][c][i] for c in range(3)] for a in range(3)]).ravel()
        def res(v):
            Nf = np.array([complex(v[0], max(v[1], 0.0))]); o = []
            for an in ang:
                r = rho_model(w1, Nf, rough, d, o1, s1, an)
                ps, dl = np.arctan(np.abs(r)), np.angle(r)
                o += [np.cos(2*ps)[0], (np.sin(2*ps)*np.cos(dl))[0], (np.sin(2*ps)*np.sin(dl))[0]]
            return np.array(o) - meas
        b = least_squares(res, g, bounds=([0.05, 0.0], [6.0, 8.0]), xtol=1e-13, ftol=1e-13)
        g = b.x; out[t] = (wl[i], b.x[0], b.x[1], np.sqrt(2*b.cost/9))
    return out[::-1]

def film_N(p, wl):
    """p = [rough, Einf, rho, tau, TL_A, TL_Br, TL_Eo, TL_Eg, G_A, G_Br, G_En]"""
    E = 1239.841984/wl
    e1, e2 = osc.drude_rt(E, p[2], p[3])
    a, b = osc.tauc_lorentz(E, p[4], p[5], p[6], p[7]); e1 += a; e2 += b
    a, b = osc.gaussian(E, p[8], p[9], p[10]);          e1 += a; e2 += b
    N = np.sqrt((p[1] + e1 + 1j*e2).astype(complex))
    return np.where(N.imag < 0, -N, N)
#            rough Einf  rho    tau   TLa   TLbr  TLeo  TLeg   Ga   Gbr   Gen
LO = np.array([0.0, 1.0, 1e-4,  2.0,  20., 0.10, 3.4, 2.80,  0.0, 0.10, 0.60])
HI = np.array([12., 4.5, 5e-1, 15.0, 1500., 6.0, 6.5, 4.10,  10., 6.00, 4.50])

def fit_osc(sh, wl, Pm, Dm, ox, si, d, dth, nstart=22, seed=13):
    M = [cf.ncs(Pm[:, i], Dm[:, i]) for i in range(3)]
    def r(p):
        Nf = film_N(p, wl); out = []
        for a, an in enumerate(cf.ANG0 + dth):
            rr = rho_model(wl, Nf, p[0], d, ox, si, an)
            ps, dl = np.arctan(np.abs(rr)), np.angle(rr)
            out += [np.cos(2*ps)-M[a][0], np.sin(2*ps)*np.cos(dl)-M[a][1],
                    np.sin(2*ps)*np.sin(dl)-M[a][2]]
        return np.concatenate(out)
    p0 = np.array([3.0, 2.5, 2e-3, 6.582, 300., 1.0, 4.3, 3.60, 0.3, 1.0, 2.5])
    rng = np.random.default_rng(seed); best = None
    for k in range(nstart):
        s = p0 if k == 0 else (LO + rng.random(len(LO))*(HI-LO) if k % 5 == 4
                               else np.clip(p0*(1+0.35*rng.standard_normal(len(p0))), LO+1e-9, HI-1e-9))
        try:
            b = least_squares(r, np.clip(s, LO+1e-9, HI-1e-9), bounds=(LO, HI),
                              x_scale='jac', max_nfev=1500)
        except Exception as e:
            print('   start %d failed: %r' % (k, e), flush=True); continue
        if best is None or b.cost < best.cost: best = b
    return best, r

if __name__ == '__main__':
    R = {}
    for sh in ['#1', '#2', '#3', '#4', '#5']:
        t0 = time.time()
        wl, Pm, Dm = cf.load(sh)
        m = (wl >= LOW) & (wl <= HIGH); wl, Pm, Dm = wl[m], Pm[m], Dm[m]
        ox, si = cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl)
        bg = None
        for d in [44., 46., 48., 50., 52., 54.]:
            for rg in [0.0, 2.5, 5.0]:
                for dth in [-0.2, 0.0, 0.2, 0.4, 0.6]:
                    c = chain(wl, Pm, Dm, ox, si, d, rg, dth, stride=14)
                    s = np.sqrt(np.mean(c[:, 3]**2))
                    if bg is None or s < bg[0]: bg = (s, d, rg, dth)
        _, d, rg0, dth = bg
        b, r = fit_osc(sh, wl, Pm, Dm, ox, si, d, dth)
        M = cf.mse(b.fun, len(LO)+2)
        ch = chain(wl, Pm, Dm, ox, si, d, b.x[0], dth, stride=4)
        R[sh] = dict(sheet=sh, d=d, dth=dth, chain_rms=float(bg[0]), mse=float(M),
                     p=[float(v) for v in b.x],
                     chain=[[float(v) for v in row] for row in ch])
        print('%s  scan rms=%.5f d=%.0f rough=%.1f dth=%+.1f | MSE=%6.2f rough_fit=%.2f  %.0fs'
              % (sh, bg[0], d, rg0, dth, M, b.x[0], time.time()-t0), flush=True)
        print('    Einf=%.3f rho=%.4g tau=%.2f | TL A=%.1f Br=%.3f Eo=%.3f Eg=%.3f | G A=%.4f Br=%.3f En=%.3f'
              % tuple(b.x[1:]), flush=True)
        json.dump(R, open('ce_v11_result.json', 'w'), indent=1)
    print('done')
