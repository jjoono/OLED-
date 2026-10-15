"""Oscillator fit at LOCKED, transmittance-validated geometry.

Letting d/roughness/angle-offset float lets the optimiser buy MSE by sliding
onto the n*d alias branch (d=46, dth=+0.73) that measured transmittance already
rejected.  So the geometry is fixed and only the dielectric function is fitted;
the remaining question is how many oscillators it takes to represent it.
"""
import numpy as np, json, math, sys
from scipy.optimize import least_squares
import ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef

LO_W, HI_W = 340.0, 1080.0
HBAR_EVS, EPS0, HB_FS = 6.582119569e-16, 8.8541878128e-12, 0.6582119569
GEO = {'#1': (50.0, 5.0, 0.38), '#2': (50.5, 3.0, 0.37),
       '#3': (42.0, 2.0, 0.18), '#4': (42.0, 2.0, -0.14), '#5': (42.0, 2.0, 0.75)}

def unpack(p, ng):
    """p = [rough, Einf, rho, tau, TLa, TLbr, TLeo, TLeg, (Ga,Gbr,Gen)*ng]"""
    return p[1], p[2:4], p[4:8], [p[8+3*i:11+3*i] for i in range(ng)]

def film_N(p, wl, ng):
    E = 1239.841984 / wl
    einf, dr, tl, gs = unpack(p, ng)
    e1, e2 = osc.drude_rt(E, dr[0], dr[1])
    a, b = osc.tauc_lorentz(E, *tl); e1 += a; e2 += b
    for g in gs:
        a, b = osc.gaussian(E, *g); e1 += a; e2 += b
    N = np.sqrt((einf + e1 + 1j*e2).astype(complex))
    return np.where(N.imag < 0, -N, N)

def bounds(ng):
    lo = [0.0, 1.0, 8e-5, 3.0,  20., 0.30, 3.5, 3.00]
    hi = [12., 4.5, 2e-2, 10.0, 400., 4.0, 4.8, 3.90]
    for i in range(ng):
        lo += [0.0, 0.10, 0.60]
        hi += [8.0, 5.00, 3.40]
    return np.array(lo), np.array(hi)

def start(sh, ng):
    q = json.load(open('genosc_params.json'))[{'#1':'1','#2':'2','#3':'3','#4':'4','#5':'5'}[sh]]['p']
    tau = HB_FS / 0.10
    rho = HBAR_EVS**2 / (q[1]*EPS0*(tau*1e-15)) * 100.0
    s = [GEO[sh][1], q[0], rho, tau, q[2], q[4], q[3], q[5],
         q[6], q[8]*2*math.sqrt(math.log(2.0)), q[7]]
    for i in range(1, ng):
        s += [0.3, 1.5, 1.2 + 0.9*i]
    return np.array(s)

def fit(sh, ng, nstart=14, seed=11):
    d, rg, dth = GEO[sh]
    wl, Pm, Dm = cf.load(sh)
    m = (wl >= LO_W) & (wl <= HI_W); wl, Pm, Dm = wl[m], Pm[m], Dm[m]
    ox, si = cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl)
    M = [cf.ncs(Pm[:, i], Dm[:, i]) for i in range(3)]
    amb = np.ones(len(wl))
    def r(p):
        Nf = film_N(p, wl, ng)
        Nr = ef.bruggeman_ema50(Nf, amb)
        out = []
        for a, an in enumerate(cf.ANG0 + dth):
            rp, rs = ef._tmm(wl, [amb, Nr, Nf, ox, si], [p[0], d, cf.D_OX], an)
            rr = rp/rs; ps, dl = np.arctan(np.abs(rr)), np.angle(rr)
            out += [np.cos(2*ps)-M[a][0], np.sin(2*ps)*np.cos(dl)-M[a][1],
                    np.sin(2*ps)*np.sin(dl)-M[a][2]]
        return np.concatenate(out)
    lo, hi = bounds(ng); p0 = np.clip(start(sh, ng), lo+1e-9, hi-1e-9)
    rng = np.random.default_rng(seed); best = None
    for k in range(nstart):
        s = p0 if k == 0 else np.clip(p0*(1+0.25*rng.standard_normal(len(p0))), lo+1e-9, hi-1e-9)
        try: b = least_squares(r, s, bounds=(lo, hi), x_scale='jac', max_nfev=1500)
        except Exception: continue
        if best is None or b.cost < best.cost: best = b
    return best, cf.mse(best.fun, len(p0)+3), (d, rg, dth)

NAMES = ['rough','Einf','rho','tau','TL_Amp','TL_Br','TL_Eo','TL_Eg','G_Amp','G_Br','G_En']
if __name__ == '__main__':
    R = {}
    for sh in ['#1','#2','#3','#4','#5']:
        b, M, g = fit(sh, 1, nstart=16)
        lo, hi = bounds(1)
        at=[NAMES[i] for i,v in enumerate(b.x)
            if abs(v-lo[i])<1e-6*max(1,abs(lo[i])) or abs(v-hi[i])<1e-6*max(1,abs(hi[i]))]
        R[sh]=dict(sheet=sh, mse=float(M), d=g[0], dth=g[2],
                   p=[g[0], float(b.x[0]), g[2]]+[float(v) for v in
                      [b.x[1],b.x[2],b.x[3],b.x[4],b.x[5],b.x[6],b.x[7],b.x[8],b.x[9],b.x[10]]],
                   at_bound=at)
        print('%s  MSE=%7.3f  d=%.1f(fix) rough=%.2f dth=%+.2f(fix)  bound:%s'
              %(sh,M,g[0],b.x[0],g[2],at or '-'), flush=True)
        print('    '+'  '.join('%s=%.5g'%(n,v) for n,v in zip(NAMES,b.x)), flush=True)
        json.dump(R, open('ce_final_result.json','w'), indent=1)
    print('done')
