"""layer2: perovskite nanocrystal film on GraHIL on Si, 205-1000 nm.

What the data ask for (gp_explore notes in the README): the 108 nm GraHIL of
layer1 is still there, and the perovskite is not one homogeneous layer.

  air
  EMA roughness      d_r     dense perovskite / void 50:50
  dense layer        d_t     eps_p(E), a KK-consistent Gaussian basis
  lower layer        d_b     Bruggeman(eps_p, void, f_b) + extinction k_sc (633 / lambda)^p
                             (p comes out ~4: scattering by grains and voids, not absorption)
  GraHIL             d_g     n,k of layer1 (gp_grahil 'iso'), 108.3 nm
  native oxide 1.5 nm / Si

Data: N, C, S at 65 / 70 / 75 deg, and the unpolarised reflectance R at the
same angles, R2 = (I2 / I1) R1 with R1 the layer1 model (the two sheets were
measured with the same lamp; one common scale s_R is fitted).

Stages
  A  structure from the transparent range (570-1000 nm), Sellmeier perovskite
  B  n,k of the opaque range (205-430 nm) point by point under the roughness of A
  C  the KK basis started from A + B (+ the absorbance-derived k of the paper's
     SI Fig. 6a, 510-565 nm, as a start only), then everything refitted together

Variants (VARIANTS) repeat C with one assumption changed each; their spread is
the uncertainty of the n,k that is reported."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import gp_load as G, gp_model as P, gp_grahil as GH
from tmm_fix import tmm
from scipy.optimize import least_squares, lsq_linear

W = _os.path.join(ELLIPS_OUT, '')
LO, HI, STEP = 205.0, 1000.0, 2
D_G, D_OX = 108.3, 1.5
WR = 0.5                                   # weight of a reflectance residual against an N,C,S residual
SMOOTH = 0.004
E_UV = 11.0
# Gaussian centres (eV): nothing below 2.24 eV (the film is transparent there), 0.03 eV steps through
# the band edge and exciton, coarser into the UV
_c = np.r_[np.arange(2.24, 2.801, 0.03), np.arange(2.90, 4.61, 0.12), np.arange(4.85, 7.31, 0.25)]
CEN, WID = _c, 1.6 * np.gradient(_c)

STRUCT = ['d_t', 'd_b', 'f_b', 'd_r', 'k_sc', 'p', 's_R']
S_LO = dict(d_t=30.0, d_b=30.0, f_b=0.15, d_r=0.0, k_sc=0.0, p=0.0, s_R=0.94, d_g=40.0, sig_s=0.0, n_m=1.0, k_top=0.0, f_b2=0.1, d_x=0.0)
S_HI = dict(d_t=260.0, d_b=330.0, f_b=1.0, d_r=25.0, k_sc=0.6, p=8.0, s_R=1.08, d_g=160.0, sig_s=40.0, n_m=1.65, k_top=0.2, f_b2=1.0, d_x=120.0)
S_0 = dict(d_t=120.0, d_b=160.0, f_b=0.5, d_r=4.0, k_sc=0.25, p=4.0, s_R=1.015, d_g=D_G, sig_s=10.0, n_m=1.3, k_top=0.01, f_b2=0.5, d_x=10.0)
M_SUB = 5
VARIANTS = {                                # name: extra structure parameters, options
    'base':            ([], {}),
    'GraHIL free':     (['d_g'], {}),
    'surface spread':  (['sig_s'], {}),
    'organic matrix':  (['n_m'], {}),
    'top scatters':    (['k_top'], {}),
    'no reflectance':  ([], dict(wr=0.0)),
    '65+70 only':      ([], dict(angles=[0, 1])),
    '70+75 only':      ([], dict(angles=[1, 2])),
    'from 230 nm':     ([], dict(lo=230.0)),
    # refinements of the lower layer: is the residual of the transparent range a too-sharp profile?
    'graded lower':    (['f_b2'], {}),                       # perovskite fraction f_b at its top -> f_b2 at its bottom
    'soft boundary':   (['d_x'], {}),                        # linear transition of thickness d_x between dense and lower
    'graded+GraHIL':   (['f_b2', 'd_g'], {}),
    'graded+spread':   (['f_b2', 'sig_s'], {}),
    'all':             (['f_b2', 'd_x', 'd_g', 'sig_s'], {}),
}


# ------------------------------------------------------------------ data
def grahil_N(wl):
    r = [v for v in json.load(open(W + 'grahil_fits.json'))['variants'] if v['variant'] == 'iso'][0]
    kb = P.KKBasis(wl, GH.CEN, GH.WID, GH.E_UV); v = GH.unpack(np.array(r['x']), 'iso', kb.n)
    return kb.N(v['a'], v['a_uv'], v['einf']), v['d']


def setup(lo=LO, hi=HI, step=STEP):
    d2 = G.load('layer2', lo, hi, step); d1 = G.load('layer1', lo, hi, step); wl = d2['wl']; ang = d2['ang']
    air, ox, si = P.n_air(wl), P.n_ox(wl), P.n_si(wl); Ng, dg1 = grahil_N(wl)
    R1 = []
    for th in ang:
        rp, rs = tmm(wl, [air, Ng, ox, si], [dg1, D_OX], th); R1.append((np.abs(rp)**2 + np.abs(rs)**2) / 2)
    Rm = d2['I'] / d1['I'] * np.array(R1).T
    return dict(wl=wl, ang=ang, meas=P.data_ncs(d2), Rm=Rm, air=air, ox=ox, si=si, Ng=Ng, kb=P.KKBasis(wl, CEN, WID, E_UV), E=P.HC / wl)


# ------------------------------------------------------------------ model
def perov_layers(wl, v, Np):
    """The perovskite as a list of (index, thickness, perovskite fraction), top (air / ETL side) first:
    dense layer, optional linear transition, lower layer (optionally graded f_b -> f_b2)."""
    one = np.ones(len(wl), complex)
    sc = (633.0 / wl)**v.get('p', 4.0)
    nm = v.get('n_m', 1.0) * one
    Nt = Np + 1j * v.get('k_top', 0.0) * sc
    low = lambda f, w=1.0: P.bruggeman(Np, nm, f) + 1j * w * v['k_sc'] * sc         # lower-layer medium at fraction f
    dx = min(v.get('d_x', 0.0), 1.8 * v['d_t'], 1.8 * v['d_b'])                     # transition, taken half from each side
    dx = dx if dx > 0.5 else 0.0
    out = [(Nt, v['d_t'] - dx / 2, 1.0)]
    for i in range(M_SUB if dx else 0):                                            # dense (f = 1, no scatter) -> lower (f_b)
        u = (i + 0.5) / M_SUB; f = 1 + (v['f_b'] - 1) * u; out.append((low(f, u), dx / M_SUB, f))
    db = v['d_b'] - dx / 2
    if 'f_b2' in v:
        for i in range(M_SUB):
            u = (i + 0.5) / M_SUB; f = v['f_b'] + (v['f_b2'] - v['f_b']) * u; out.append((low(f), db / M_SUB, f))
    else:
        out.append((low(v['f_b']), db, v['f_b']))
    return out


def stack_r(S, v, Np, angles=None):
    """r_p, r_s [n_wl, n_ang] for dense index Np and structure v (dict)."""
    wl = S['wl']; ang = S['ang'] if angles is None else S['ang'][angles]
    pl = perov_layers(wl, v, Np)
    L = [S['air'], P.bruggeman(pl[0][0], S['air'], 0.5)] + [x[0] for x in pl] + [S['Ng'], S['ox'], S['si']]
    D = [v['d_r']] + [x[1] for x in pl] + [v.get('d_g', D_G), D_OX]
    if v.get('sig_s', 0.0) > 0:
        return P.rr_coh2(wl, ang, L, D, 1, 2, v['sig_s'], 0.0, nodes=7)
    RP, RS = [], []
    for th in ang:
        rp, rs = tmm(wl, L, D, th); RP.append(rp); RS.append(rs)
    return np.stack(RP, 1), np.stack(RS, 1)


def residuals(S, v, Np, wr=WR, angles=None):
    rp, rs = stack_r(S, v, Np, angles)
    m = P.ncs_from_r(rp, rs); R = (np.abs(rp)**2 + np.abs(rs)**2) / 2
    a = slice(None) if angles is None else angles
    meas = tuple(x[:, a] for x in S['meas'])
    return P.resid(m, meas), (R * v.get('s_R', 1.0) - S['Rm'][:, a]).ravel() * wr


# ------------------------------------------------------------------ stage A: structure, transparent range
def stage_a(n_start=60, seed=0):
    S = setup(570.0, 1000.0, 3); wl = S['wl']; names = STRUCT[:-1] + ['einf', 'Ssel', 'lam0']
    lo = np.array([S_LO[k] for k in STRUCT[:-1]] + [1.0, 0.0, 380.0]); hi = np.array([S_HI[k] for k in STRUCT[:-1]] + [4.5, 3.5, 535.0])
    base = np.array([S_0[k] for k in STRUCT[:-1]] + [2.6, 1.3, 430.0]); rng = np.random.default_rng(seed); best = None
    def res(x):
        v = dict(zip(names, x)); v['s_R'] = 1.0 / 1.015
        Np = np.sqrt(v['einf'] + v['Ssel'] * wl**2 / (wl**2 - v['lam0']**2) + 0j)
        a, b = residuals(S, v, Np); return np.r_[a, b]
    for k in range(n_start):
        x0 = base * (1 + 0.6 * (rng.random(len(lo)) - 0.5)) if k % 2 else lo + rng.random(len(lo)) * (hi - lo)
        try: r = least_squares(res, np.clip(x0, lo + 1e-6, hi - 1e-6), bounds=(lo, hi), x_scale='jac', max_nfev=300)
        except Exception: continue
        if best is None or r.cost < best.cost: best = r
    v = dict(zip(names, [float(t) for t in best.x])); n3 = 9 * len(wl)
    v['mse_ncs'] = float(P.mse_ncs(best.fun[:n3], len(lo))); v['rms_R'] = float(1e3 * np.sqrt(np.mean((best.fun[n3:] / WR)**2)))
    return v


# ------------------------------------------------------------------ stage B: opaque range, point by point
def stage_b(d_r):
    from scipy.sparse import lil_matrix
    S = setup(LO, 430.0, 3); wl = S['wl']; nw = len(wl); ang = S['ang']
    rho = np.tan(np.deg2rad(G.load('layer2', LO, 430.0, 3)['psi'])) * np.exp(-1j * np.deg2rad(G.load('layer2', LO, 430.0, 3)['del']))
    th = np.deg2rad(ang); e0 = np.mean([np.sin(t)**2 * (1 + np.tan(t)**2 * ((1 - rho[:, a]) / (1 + rho[:, a]))**2) for a, t in enumerate(th)], 0)
    n0 = np.sqrt(e0); n0 = np.where(n0.imag < 0, -n0, n0)
    Jp = lil_matrix((9 * nw, 2 * nw), dtype=int)
    for i in range(nw):
        for c in range(3):
            for a_ in range(3): Jp[c * nw * 3 + i * 3 + a_, i] = 1; Jp[c * nw * 3 + i * 3 + a_, nw + i] = 1
    def res(x):
        n = x[:nw] + 1j * x[nw:]
        return P.resid(P.ncs(wl, ang, [S['air'], P.bruggeman(n, S['air'], 0.5), n], [d_r]), S['meas'])
    r = least_squares(res, np.r_[n0.real * 1.05, n0.imag], jac_sparsity=Jp, bounds=(np.r_[np.full(nw, 0.8), np.zeros(nw)], np.r_[np.full(nw, 4.0), np.full(nw, 3.0)]), max_nfev=80)
    return wl, r.x[:nw] + 1j * r.x[nw:]


# ------------------------------------------------------------------ stage C: KK basis, everything together
def init_basis(va, wl_uv, N_uv):
    """amplitudes that reproduce stage A (transparent), stage B (UV) and, in between, a smooth join
    through the paper's absorbance-derived k (start values only)."""
    wt = np.arange(570.0, 1001.0, 10.0); Nt = np.sqrt(va['einf'] + va['Ssel'] * wt**2 / (wt**2 - va['lam0']**2) + 0j)
    wm = np.arange(440.0, 566.0, 5.0)
    kp = np.interp(wm, [440, 480, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565], [0.30, 0.18, 0.115, 0.119, 0.125, 0.132, 0.125, 0.082, 0.035, 0.014, 0.006, 0.003, 0.0015, 0.0004])
    nmid = np.interp(wm, [430, 500, 530, 545, 570], [N_uv[-1].real, 2.15, 2.35, 2.40, Nt[0].real])
    wl = np.r_[wl_uv, wm, wt]; N = np.r_[N_uv, nmid + 1j * kp, Nt]; eps = N**2
    kb = P.KKBasis(wl, CEN, WID, E_UV)
    A = np.c_[np.ones(len(wl)), kb.pole, kb.G]                   # unknowns: einf, a_uv, a
    wgt = np.r_[np.ones(len(wl_uv)), np.full(len(wm), 2.0), np.full(len(wt), 3.0)]
    M = np.r_[(A.real * wgt[:, None]), (A.imag * wgt[:, None])]; y = np.r_[eps.real * wgt, eps.imag * wgt]
    D2 = np.c_[np.zeros((len(CEN) - 2, 2)), np.diff(np.eye(len(CEN)), 2, axis=0)] * 0.05
    r = lsq_linear(np.r_[M, D2], np.r_[y, np.zeros(len(CEN) - 2)], bounds=(np.r_[1.0, 0.0, np.zeros(len(CEN))], np.r_[4.0, 400.0, np.full(len(CEN), 30.0)]))
    return r.x


def stage_c(variant, va, x_basis):
    extra, opt = VARIANTS[variant]
    S = setup(opt.get('lo', LO), HI, STEP); nb = S['kb'].n; names = STRUCT + extra; ns = len(names)
    wr = opt.get('wr', WR); angles = opt.get('angles')
    lo = np.r_[[S_LO[k] for k in names], 1.0, 0.0, np.zeros(nb)]; hi = np.r_[[S_HI[k] for k in names], 4.0, 400.0, np.full(nb, 30.0)]
    x0 = np.r_[[va.get(k, S_0[k]) for k in names], x_basis]
    x0 = np.clip(x0, lo + 1e-9, hi - 1e-9)
    D2 = np.diff(np.eye(nb), 2, axis=0)
    def split(x):
        v = dict(zip(names, x[:ns])); return v, S['kb'].N(x[ns + 2:], x[ns + 1], x[ns])
    def res(x):
        v, Np = split(x); a, b = residuals(S, v, Np, wr, angles)
        return np.r_[a, b, SMOOTH * (D2 @ x[ns + 2:])]
    r = least_squares(res, x0, bounds=(lo, hi), x_scale='jac', max_nfev=120)
    x = r.x; v, Np = split(x); a, b = residuals(S, v, Np, wr if wr > 0 else WR, angles)
    wl = S['wl']; na = len(S['ang'] if angles is None else angles); rr = a.reshape(3, len(wl), na)
    band = lambda p, q: float(1e3 * np.sqrt(np.mean(rr[:, (wl >= p) & (wl < q), :]**2)))
    J = r.jac; s2 = np.sum(r.fun**2) / max(r.fun.size - x.size, 1)
    ci = 1.645 * np.sqrt(np.clip(np.diag(np.linalg.pinv(J.T @ J) * s2), 0, None))
    return dict(variant=variant, names=names, x=[float(t) for t in x], struct={k: float(t) for k, t in v.items()},
                struct_ci={k: float(c) for k, c in zip(names, ci[:ns])},
                mse_ncs=float(P.mse_ncs(a, len(x))), rms_R=float(1e3 * np.sqrt(np.mean((b / (wr if wr > 0 else WR))**2))),
                bands={'205-300': band(205, 300), '300-450': band(300, 450), '450-570': band(450, 570), '570-1000': band(570, 1001)},
                nfev=int(r.nfev), lo=float(opt.get('lo', LO)))


def nk_of(r, wl):
    ns = len(r['names']); x = np.array(r['x']); kb = P.KKBasis(wl, CEN, WID, E_UV)
    return kb.N(x[ns + 2:], x[ns + 1], x[ns])


def run_variant(a):
    return stage_c(*a)


def main():
    from multiprocessing import Pool
    only = sys.argv[1:]                                     # variant names: rerun just these, starting from the saved base fit
    old = json.load(open(W + 'perov_fits.json')) if _os.path.exists(W + 'perov_fits.json') else None
    if only and old:
        va = old['stage_a']; b = [r for r in old['variants'] if r['variant'] == 'base'][0]
        va = dict(va, **b['struct']); xb = np.array(b['x'][len(b['names']):])
    else:
        va = stage_a(); print('stage A (570-1000 nm, Sellmeier):', {k: round(v, 3) for k, v in va.items()})
        wl_uv, N_uv = stage_b(va['d_r']); xb = init_basis(va, wl_uv, N_uv)
    todo = only if only else list(VARIANTS)
    with Pool(4) as pool:
        res = pool.map(run_variant, [(k, va, xb) for k in todo])
    if only and old:
        keep = [r for r in old['variants'] if r['variant'] not in todo]; res = keep + res; va = old['stage_a']
    json.dump(dict(stage_a=va, variants=res), open(W + 'perov_fits.json', 'w'), indent=1)
    q = np.array([300.0, 400.0, 450.0, 500.0, 520.0, 530.0, 540.0, 550.0, 600.0, 633.0, 800.0, 1000.0])
    print('%-16s %6s %6s | %s | d_t    d_b    f_b   d_r   k_sc  p    s_R   extra' % ('variant', 'MSE', 'rms R', '205-300 300-450 450-570 570-1000'))
    for r in res:
        s = r['struct']
        s = dict(s, f_b=s['f_b'])
        print('%-16s %6.2f %6.1f | %7.1f %7.1f %7.1f %8.1f | %5.1f %6.1f %5.2f %5.1f %5.3f %4.1f %5.3f  %s'
              % (r['variant'], r['mse_ncs'], r['rms_R'], *[r['bands'][k] for k in ('205-300', '300-450', '450-570', '570-1000')], s['d_t'], s['d_b'], s['f_b'], s['d_r'], s['k_sc'], s['p'], s['s_R'],
                 '  '.join('%s=%.3g' % (k, s[k]) for k in r['names'][len(STRUCT):])))
    print('\ndense perovskite n + ik:')
    print('%-16s ' % 'wavelength' + ' '.join('%11d' % w for w in q))
    for r in res:
        N = nk_of(r, q); print('%-16s ' % r['variant'] + ' '.join('%5.3f+%.3fi' % (z.real, z.imag) for z in N))


if __name__ == '__main__':
    main()
