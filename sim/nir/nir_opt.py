"""NIR bottom-emitting OLED with a low-loss multi-pass extractor: global optimisation of the stack.

  Ag 100 nm | ETL d_E | EML 20 nm (dipole at centre, Theta 2/3) | HTL d_H | transparent electrode (TE) [| seed 5 nm]
  | chirped DBR (L/H x 4, every thickness free, 0 allowed) | substrate n_sub (1.5 or 1.8) | hemispherical MLA (n_MLA = n_sub)

Organics: isotropic, lossless, n = 1.75.  Emitter: Gaussian, lambda0 760 nm, FWHM 80 nm (photon-number weights).
Low PLQY q: per exciton EQE = q * sum_l w F eta_sub eta_ext / (1 - q + q * sum_l w F)   (F = P_tot of the CPS model,
i.e. the Purcell factor relative to the dipole in bulk EML; eta_sub = air + sub of cps2 (the tunnelled near field is part of 'sub'); eta_ext = recycling series with the
LightTools BSDF of the hemispherical MLA and the stack reflectance R(theta) from the substrate).
TE options (600-900 nm data):
  ITO  ITO_NIR_ccho (measured, highly doped)          d >= 44 nm  (rho 2.2e-4 Ohm cm, Rs <= 50 Ohm/sq)
  IZO  IZO_NIR (measured)                              d >= 80 nm  (rho 4e-4)
  IOH  hydrogenated In2O3, Drude model (eps_inf 4.0, N 1.5e20 cm-3, mu 100 cm2/Vs, m* 0.3) d >= 84 nm (rho 4.2e-4)
  Ag   thin Ag on HATCN (measured joint n,k of 8/10/12 nm films, interpolated in thickness), 8-12 nm, + seed 5 nm
  Au   thin Au, bulk literature n,k (Au_liter_NIR; optimistic for 8-15 nm films), + seed 5 nm
DBR: TiO2 (TiO2_NIR) / SiO2 (l_SiO2_SJ_RTP), deposited on the substrate before the OLED.
Ag reflector: McPeak to 800 nm, Drude fit beyond.
Usage: python3 nir_opt.py <TE> <n_sub> <ndbr> [q]   -> results/<TE>_n<n_sub>_dbr<ndbr>_q<q>.json"""
import sys, json, os, numpy as np, scipy.io as sio
from scipy.optimize import differential_evolution, minimize
sys.path[:0] = ['/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/mla']
import cps2, series

L = np.arange(600, 901)                      # nm, 1-nm grid
M = sio.loadmat('/home/user/OLED-/sim/table1_author/nk_JH_total.mat', squeeze_me=True, struct_as_record=False)['material']
g = lambda n: np.asarray(getattr(M, n))
def vis401(n):                               # 400-800 nm arrays -> 600-900, flat beyond 800
    a = g(n)[200:]; return np.concatenate([a, np.full(100, a[-1])])
# Ag reflector: McPeak <= 800 nm, Drude beyond (fit 650-800 nm)
agv = g('l_Ag_McPeak')[250:401]; lv = np.arange(650, 801); w_ = 1239.84 / lv; eps = agv ** 2
A_ = np.c_[np.ones_like(w_), -1 / w_ ** 2]; einf, wp2 = np.linalg.lstsq(A_, eps.real, rcond=None)[0]
gam = np.median(eps.imag * w_ ** 3 / wp2)
def drude(lam, einf, wp2, gam):
    w = 1239.84 / lam; return np.sqrt(einf - wp2 / (w * (w + 1j * gam)))
AG = np.concatenate([g('l_Ag_McPeak')[200:401], drude(np.arange(801, 901), einf, wp2, gam)])
TIO2 = g('TiO2_NIR').astype(complex); SIO2 = vis401('l_SiO2_SJ_RTP')
def ioh(lam, einf=4.0, N=1.5e26, mu=0.01, ms=0.3):
    e, m0, e0, c = 1.602e-19, 9.109e-31, 8.854e-12, 2.998e8
    w = 2 * np.pi * c / (lam * 1e-9); wp2 = N * e ** 2 / (e0 * ms * m0); gm = e / (ms * m0 * mu)
    return np.sqrt(einf - wp2 / (w * (w + 1j * gm)))
TCO = {'ITO': (g('ITO_NIR_ccho'), 44.0), 'IZO': (g('IZO_NIR'), 80.0), 'IOH': (ioh(L.astype(float)), 84.0)}
def joint(f):
    a = np.genfromtxt('/home/user/OLED-/sim/teoled_opt/nk_joint/' + f, delimiter=',', names=True)
    return np.interp(L, a['wavelength_nm'], a['n_joint']) + 1j * np.interp(L, a['wavelength_nm'], a['k_joint'])
AGT = {8.0: joint('Ag_HATCN_2-6_8nm_nk.csv'), 10.0: joint('Ag_HATCN_2-7_10nm_nk.csv'), 12.0: joint('Ag_HATCN_2-8_12nm_nk.csv')}
def agthin(d, j):
    d = min(max(d, 8.0), 12.0); lo = 8.0 if d <= 10 else 10.0; t = (d - lo) / 2.0
    return (1 - t) * AGT[lo][j] + t * AGT[lo + 2][j]
AU = g('Au_liter_NIR')
NORG, THETA, DEML = 1.75, 2 / 3, 20.0
TH = np.arange(90) + 0.5
B = sio.loadmat('/home/user/OLED-/sim/mla/lt_hemisphere_bsdf.mat', squeeze_me=True)['BSDF_MLA']
def ext(nsub, kind='mla'):
    if kind == 'mla':
        b = B[:, :, int(round((nsub - 1.30) / 0.05))]; return b[:90].sum(0), b[90:180][::-1]
    c = np.cos(np.radians(TH)); s2 = (nsub * np.sin(np.radians(TH))) ** 2; ct = np.sqrt(np.clip(1 - s2, 0, 1))
    Rf = np.where(s2 < 1, 0.5 * (((nsub * c - ct) / (nsub * c + ct)) ** 2 + ((nsub * ct - c) / (nsub * ct + c)) ** 2), 1.0)
    return 1 - Rf, np.diag(Rf)
def spectrum(step):
    l = np.arange(680, 841, step); s = 80 / 2.3548; w = np.exp(-0.5 * ((l - 760) / s) ** 2); return l, w / w.sum()

def stack(x, te, ndbr, lam):
    j = int(lam - 600); dE, dH, dT = x[:3]; dbr = x[3:3 + 2 * ndbr]
    below = [(NORG, NORG, dH)]
    if te in TCO: n = TCO[te][0][j]; below.append((n, n, dT))
    else:
        n = agthin(dT, j) if te == 'Ag' else AU[j]; below += [(n, n, dT), (NORG, NORG, 5.0)]
    for k, d in enumerate(dbr):
        if d > 0.5: n = SIO2[j] if k % 2 == 0 else TIO2[j]; below.append((n, n, d))
    return below

def point(args):
    x, te, ndbr, nsub, lam, npts, kind = args
    S = cps2.Stack(float(lam), (NORG, NORG), DEML, DEML / 2, above=[(NORG, NORG, x[0]), (AG[int(lam - 600)],) * 2 + (100.0,)],
                   below=stack(x, te, ndbr, lam), n_sub=nsub, h=THETA)
    r = cps2.solve_pol(S, npts=npts); F = float(np.real(r['P_tot']))
    es = float(np.real(r['air'] + r['sub']))          # 'sub' already contains the tunnelled near field
    P = np.real(cps2.sub_angular(S, TH)); R = np.real(cps2.stack_reflectance(S, TH))
    BT, BR = ext(nsub, kind); ee = series.eta_ext(BT, BR, R, P, 150)[0]
    wP = P / max(P.sum(), 1e-30)
    return F, es, ee, float(wP @ (1 - R)), float(np.real(r['spp'])), float(np.real(r['abs']))

def evaluate(x, te, ndbr, nsub, step=10, npts=1500, kind='mla', mapper=map, qs=(0.3,)):
    l, w = spectrum(step)
    res = np.array(list(mapper(point, [(tuple(x), te, ndbr, nsub, lam, npts, kind) for lam in l])))
    F, es, ee = res[:, 0], res[:, 1], res[:, 2]
    num = float(w @ (F * es * ee)); Fb = float(w @ F)
    out = {f'EQE_q{q}': q * num / (1 - q + q * Fb) for q in qs}
    out.update(F=Fb, eta_sub=float(w @ (F * es)) / Fb, eta_ext=num / float(w @ (F * es)),
               A_round=float(w @ res[:, 3]), spp=float(w @ res[:, 4]), absorb=float(w @ res[:, 5]))
    return out

def bounds(te, ndbr):
    b = [(20.0, 400.0), (20.0, 400.0)]
    b.append((TCO[te][1], 300.0) if te in TCO else ((8.0, 12.0) if te == 'Ag' else (8.0, 15.0)))
    for k in range(2 * ndbr): b.append((0.0, 300.0) if k % 2 == 0 else (0.0, 200.0))
    return b

_CASE = {}
def obj(x):
    c = _CASE; return -evaluate(x, c['te'], c['ndbr'], c['nsub'], step=10, npts=1500, kind=c['kind'], qs=(c['q'],))[f"EQE_q{c['q']}"]

if __name__ == '__main__':
    from multiprocessing import Pool
    te, nsub, ndbr = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]); q = float(sys.argv[4]) if len(sys.argv) > 4 else 0.3
    kind = sys.argv[5] if len(sys.argv) > 5 else 'mla'
    _CASE.update(te=te, nsub=nsub, ndbr=ndbr, q=q, kind=kind)
    bd = bounds(te, ndbr); dim = len(bd)
    with Pool(4) as pool:
        r = differential_evolution(obj, bd, seed=1, popsize=max(6, 40 // dim), maxiter=35 if dim > 5 else 30, tol=1e-4,
                                   polish=False, workers=pool.map, updating='deferred', init='latinhypercube')
        f5 = lambda x: -evaluate(np.clip(x, [b[0] for b in bd], [b[1] for b in bd]), te, ndbr, nsub, step=5, npts=1500, kind=kind,
                                 mapper=pool.map, qs=(q,))[f'EQE_q{q}']
        nm = minimize(f5, r.x, method='Nelder-Mead', options=dict(maxiter=60 * dim, xatol=0.5, fatol=1e-5))
        xb = np.clip(nm.x, [b[0] for b in bd], [b[1] for b in bd])
        fin = evaluate(xb, te, ndbr, nsub, step=2, npts=3000, kind=kind, mapper=pool.map, qs=(0.1, 0.3, 0.5, 1.0))
    names = ['d_ETL', 'd_HTL', 'd_TE'] + [('SiO2_' if k % 2 == 0 else 'TiO2_') + str(k // 2 + 1) for k in range(2 * ndbr)]
    out = dict(TE=te, n_sub=nsub, ndbr=ndbr, q_opt=q, extractor=kind, de=-r.fun, nm5=-nm.fun, **fin,
               structure={k: round(float(v), 1) for k, v in zip(names, xb)})
    os.makedirs('/home/user/OLED-/sim/nir/results', exist_ok=True)
    fn = '/home/user/OLED-/sim/nir/results/%s_n%.1f_dbr%d_q%.1f_%s.json' % (te, nsub, ndbr, q, kind)
    json.dump(out, open(fn, 'w'), indent=1); print(json.dumps(out), flush=True)
