"""On-axis luminance of top-emitting stacks with a front DBR angular filter, without (A) and with (B) a back scattering
layer + Ag reflector that redistributes the light the DBR sends back, relative to a conventional bottom-emitting glass OLED.

Emitter: Ir(ppy)2acac spectrum, Theta 0.76, PLQY 1, dipoles at the centre of a 25 nm EML; organics birefringent (measured).
  REF : glass n 1.52 | IZO 50 | HATCN 5 | TAPC d_H | EML | B3PyMPM d_E | Ag 100      (on-axis in air via the flat glass face)
  A   : air | DBR (ZnS/LiF, N pairs, ZnS next to the IZO, quarter-wave at lam_c) | IZO 50 | HATCN 5 | TAPC d_H | EML |
        B3PyMPM d_E | IZO 50 | Ag 100
  B   : as A, but IZO 50 | scattering layer (n_s 1.8, lossless, full angular randomisation) | Ag 100 at the back.
        Direct part: coherent dipole emission into air (scattering medium as a semi-infinite half-space).
        Recycled part: power into the scattering medium P_s -> Ag (R_b, Lambertian redistribution) -> front stack
        (T_f, R_f seen from n_s, coherent, cos-sin weighted) -> ... ;  Phi = P_s R_b / (1 - R_b <R_f>),
        on-axis radiance in air = Phi T_f(0) / (pi n_s^2).
Luminance = sum over lambda of I(lam) V(lam) x (on-axis radiance per emitted photon).  Ratios to REF are reported.
Usage: python3 front_lum.py <case REF|A|B> <N pairs>   -> results/<case>_N<N>.json"""
import sys, os, json, numpy as np, scipy.io as sio
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize
sys.path[:0] = ['/home/user/OLED-/sim/diag', '/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/be_cone', '/home/user/OLED-/sim/mla']
import cps2, dev3
D = dev3.D; IZO = dev3.IZO
M = sio.loadmat('/home/user/OLED-/sim/table1_author/nk_JH_total.mat', squeeze_me=True, struct_as_record=False)['material']
ZNS = np.asarray(M.l_ZnS); LIF = np.asarray(M.l_LiF)              # 400-800 nm, 1 nm
NS, THETA, NG = 1.8, 0.76, 1.52
I = np.clip(np.asarray(D.spec.l_I_Irppy2acac, float)[:301], 0, None)
V = lambda l: 1.019 * np.exp(-285.4 * (l / 1000 - 0.559) ** 2)
TH = np.arange(90) + 0.5; CW = np.cos(np.radians(TH)) * np.sin(np.radians(TH)); CW /= CW.sum()

def dbr(N, lamc, j):
    if N == 0: return []
    jc = int(round(lamc - 400)); dH = lamc / 4 / ZNS[jc].real; dL = lamc / 4 / LIF[jc].real
    return [(ZNS[j], ZNS[j], dH), (LIF[j], LIF[j], dL)] * N
def front_layers(j, dH, N, lamc):
    return [(D.TPo[j], D.TPe[j], dH), (D.HAT[j], D.HAT[j], 5.0), (IZO[j], IZO[j], 50.0)] + dbr(N, lamc, j)
def tmm(ns, ds, n0, nexit, lam, kx, pol):
    k0 = 2 * np.pi / lam
    def q(n):
        v = np.sqrt(complex(n) ** 2 * k0 ** 2 - kx ** 2 + 0j); return v if v.imag > 0 or (v.imag == 0 and v.real >= 0) else -v
    eta = (lambda n: q(n) / k0) if pol == 's' else (lambda n: complex(n) ** 2 * k0 / q(n))
    Mm = np.eye(2, dtype=complex)
    for n, d in zip(ns, ds):
        de = q(n) * d; e = eta(n); Mm = Mm @ np.array([[np.cos(de), -1j * np.sin(de) / e], [-1j * e * np.sin(de), np.cos(de)]])
    e0, es = eta(n0), eta(nexit); B, C = Mm @ np.array([1, es]); den = e0 * B + C
    return abs((e0 * B - C) / den) ** 2, 4 * e0.real * es.real / abs(den) ** 2
def per_sr(S, th=0.25):
    return float(np.real(cps2.sub_angular(S, [th]))[0]) * 180 / np.pi / (2 * np.pi * np.sin(np.radians(th)))

def point(args):
    """per wavelength: on-axis radiance in air per unit emitted power, plus bookkeeping"""
    case, j, dE, dH, N, lamc, npts = args
    lam = float(400 + j); EML = (D.EMo[j], D.EMe[j])
    if case == 'REF':
        S = cps2.Stack(lam, EML, 25.0, 12.5, above=[(D.B3o[j], D.B3e[j], dE), (D.AG[j], D.AG[j], 100.0)],
                       below=[(D.TPo[j], D.TPe[j], dH), (D.HAT[j], D.HAT[j], 5.0), (IZO[j], IZO[j], 50.0)], n_sub=NG, h=THETA)
        r = cps2.solve_pol(S, npts=npts); P = float(np.real(r['P_tot']))
        Rg = ((NG - 1) / (NG + 1)) ** 2; Rl = float(np.real(cps2.stack_reflectance(S, [0.25]))[0])
        L0 = per_sr(S) / P * (1 - Rg) / NG ** 2 / (1 - Rg * Rl)
        return L0, float(np.real(r['air'])), 0.0
    back = [(D.B3o[j], D.B3e[j], dE), (IZO[j], IZO[j], 50.0)]
    fr = front_layers(j, dH, N, lamc)
    if case == 'A':
        S = cps2.Stack(lam, EML, 25.0, 12.5, above=back + [(D.AG[j], D.AG[j], 100.0)], below=fr, n_sub=1.0, h=THETA)
        r = cps2.solve_pol(S, npts=npts); P = float(np.real(r['P_tot']))
        return per_sr(S) / P, float(np.real(r['air'] + r['sub'])), 0.0
    # case B
    S2 = cps2.Stack(lam, EML, 25.0, 12.5, above=fr, below=back, n_sub=NS, n_top=1.0, h=THETA)     # viewpoint: scattering medium
    r2 = cps2.solve_pol(S2, npts=npts); P = float(np.real(r2['P_tot'])); Ps = float(np.real(r2['air'] + r2['sub']))   # 'sub' excludes the air-cone part
    S1 = cps2.Stack(lam, EML, 25.0, 12.5, above=back, below=fr, n_sub=1.0, n_top=NS, h=THETA)      # viewpoint: air
    Ld = per_sr(S1) / P
    Pair = float(np.sum(np.real(cps2.sub_angular(S1, TH)))) / P                                     # direct power into air
    ns_ = [l[0] for l in back[::-1]] + [D.EMo[j]] + [l[0] for l in fr]; ds_ = [l[2] for l in back[::-1]] + [25.0] + [l[2] for l in fr]
    Tf = np.zeros(90); Rf = np.zeros(90); Rb = np.zeros(90)
    for i, th in enumerate(TH):
        kx = 2 * np.pi / lam * NS * np.sin(np.radians(th))
        for pol in 's', 'p':
            R_, T_ = tmm(ns_, ds_, NS, 1.0, lam, kx, pol); Tf[i] += T_ / 2; Rf[i] += R_ / 2
            Rb[i] += tmm([D.AG[j]], [100.0], NS, 1.0, lam, kx, pol)[0] / 2
    return Ld, Pair, (Ps, float(CW @ Tf), float(CW @ Rf), float(CW @ Rb), float(Tf[0]))

def evaluate(case, x, N, step=10, npts=1500, mapper=map, Rsc=0.0, full=False):
    dE, dH, lamc = (list(x) + [550.0])[:3]
    js = list(range(0, 301, step)); res = list(mapper(point, [(case, j, dE, dH, N, lamc, npts) for j in js]))
    lam = np.array(js) + 400.0; w = I[js] * V(lam); w = w / w.sum(); wp = I[js] * lam; wp = wp / wp.sum()
    if case != 'B':
        L = np.array([r[0] for r in res]); E = np.array([r[1] for r in res])
        return dict(L0=float(w @ L), EQE=float(wp @ E))
    Ld = np.array([r[0] for r in res]); Pa = np.array([r[1] for r in res]); X = np.array([r[2] for r in res])
    Ps, Tm, Rm, Rbm, T0 = X.T; Rb = Rbm * (1 - Rsc)
    Phi = Ps * Rb / (1 - Rb * Rm); Lr = Phi * T0 / (np.pi * NS ** 2)
    out = dict(L0=float(w @ (Ld + Lr)), L0_direct=float(w @ Ld), L0_recycled=float(w @ Lr), EQE=float(wp @ (Pa + Phi * Tm)),
               EQE_direct=float(wp @ Pa), P_s=float(wp @ Ps), Rf_mean=float(wp @ Rm), Tf_mean=float(wp @ Tm), Rb=float(wp @ Rb))
    return out

class Obj:
    def __init__(self, case, N, Rsc=0.0): self.a = (case, N, Rsc)
    def __call__(self, x): return -evaluate(self.a[0], x, self.a[1], Rsc=self.a[2])['L0']

if __name__ == '__main__':
    case, N = sys.argv[1], int(sys.argv[2]); os.makedirs('results', exist_ok=True)
    bd = [(20.0, 300.0), (20.0, 300.0)] + ([(440.0, 800.0)] if case != 'REF' and N > 0 else [])
    with Pool(4) as pool:
        r = differential_evolution(Obj(case, N), bd, seed=1, popsize=12, maxiter=25, tol=1e-4, polish=False, workers=pool.map, updating='deferred')
        lo = np.array([b[0] for b in bd]); hi = np.array([b[1] for b in bd])
        m = minimize(lambda x: -evaluate(case, np.clip(x, lo, hi), N, step=5, mapper=pool.map)['L0'], r.x, method='Nelder-Mead',
                     options=dict(maxiter=80 * len(bd), xatol=0.5, fatol=1e-6))
        xb = np.clip(m.x, lo, hi); fin = evaluate(case, xb, N, step=5, npts=3000, mapper=pool.map)
        if case == 'B':
            for rs in (0.02, 0.05): fin['L0_Rsc%.2f' % rs] = evaluate(case, xb, N, step=5, npts=3000, mapper=pool.map, Rsc=rs)['L0']
    out = dict(case=case, N=N, d_ETL=round(float(xb[0]), 1), d_HTL=round(float(xb[1]), 1), lam_c=round(float(xb[2]), 1) if len(xb) > 2 else None, **fin)
    json.dump(out, open('results/%s_N%d.json' % (case, N), 'w'), indent=1); print(json.dumps(out), flush=True)
