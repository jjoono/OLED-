"""Direct Gen-Osc fit (Einf + Drude + Tauc-Lorentz + Gaussian) to the raw data
over the range where a homogeneous single layer is actually valid."""
import numpy as np, json, sys
from scipy.optimize import least_squares
import ce_fit as cf, ellipsometry_fit as ef, ce_osc as osc

#      d   rough  dth   Einf   rho     tau   TL_A  TL_C  TL_E0 TL_Eg  G_A   G_Br  G_En
LO = np.array([30., 0.0, -1.0, 0.8, 8e-5,  1.5,  20., 0.05, 3.3, 3.00, 0.0, 0.15, 1.0])
HI = np.array([90., 12.,  1.0, 4.5, 2e-2, 20.0, 500., 4.0,  5.2, 4.00, 3.0, 4.0,  3.3])

def film_N(p, wl):
    E = 1239.841984 / wl
    e1, e2 = osc.drude_rt(E, p[4], p[5])
    a, b = osc.tauc_lorentz(E, p[6], p[7], p[8], p[9]); e1 += a; e2 += b
    a, b = osc.gaussian(E, p[10], p[11], p[12]);        e1 += a; e2 += b
    N = np.sqrt((p[3] + e1 + 1j*e2).astype(complex))
    return np.where(N.imag < 0, -N, N)

def resid_factory(wl, Pm, Dm, ox, si):
    M = [cf.ncs(Pm[:, i], Dm[:, i]) for i in range(3)]
    def r(p):
        Nf = film_N(p, wl)
        Nr = ef.bruggeman_ema50(Nf, np.ones_like(Nf)); amb = np.ones_like(Nf)
        out = []
        for a, an in enumerate(cf.ANG0 + p[2]):
            rp, rs = ef._tmm(wl, [amb, Nr, Nf, ox, si], [p[1], p[0], cf.D_OX], an)
            rr = rp/rs; ps, dl = np.arctan(np.abs(rr)), np.angle(rr)
            out += [np.cos(2*ps)-M[a][0], np.sin(2*ps)*np.cos(dl)-M[a][1],
                    np.sin(2*ps)*np.sin(dl)-M[a][2]]
        return np.concatenate(out)
    return r

def fit(sheet, lo, hi, nstart=14, seed=1):
    wl, Pm, Dm = cf.load(sheet)
    m = (wl >= lo) & (wl <= hi)
    wl, Pm, Dm = wl[m], Pm[m], Dm[m]
    ox, si = cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl)
    r = resid_factory(wl, Pm, Dm, ox, si)
    rng = np.random.default_rng(seed)
    P0 = np.array([51., 2., 0.3, 2.0, 6e-4, 6.0, 150., 1.0, 4.1, 3.45, 0.05, 1.0, 1.9])
    best = None
    for k in range(nstart):
        s = P0.copy()
        if k:
            s = LO + rng.random(len(LO))*(HI-LO) if k % 4 == 0 else \
                np.clip(P0*(1+0.35*rng.standard_normal(len(P0))), LO+1e-9, HI-1e-9)
            s[0] = rng.uniform(44, 60); s[2] = rng.uniform(-0.6, 1.0)
        s = np.clip(s, LO+1e-9, HI-1e-9)
        try:
            b = least_squares(r, s, bounds=(LO, HI), method='trf', x_scale='jac', max_nfev=700)
        except Exception:
            continue
        if best is None or b.cost < best.cost: best = b
    return wl, best, cf.mse(best.fun, len(LO)), r

NAMES = ['d','rough','dth','Einf','rho','tau','TL_Amp','TL_Br','TL_Eo','TL_Eg','G_Amp','G_Br','G_En']
if __name__ == '__main__':
    for lo, hi in [(400, 1689)]:
        wl, b, M, _ = fit('#1', lo, hi, nstart=14)
        at = [NAMES[i] for i, v in enumerate(b.x)
              if abs(v-LO[i]) < 1e-6*max(1,abs(LO[i])) or abs(v-HI[i]) < 1e-6*max(1,abs(HI[i]))]
        print('#1  %4d-%4d nm   MSE = %7.3f   at-bound: %s' % (lo, hi, M, at or 'none'))
        print('    ' + '  '.join('%s=%.4g' % (n, v) for n, v in zip(NAMES, b.x)))
        np.save('ce_p1.npy', b.x)
