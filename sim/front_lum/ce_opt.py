"""Maximum on-axis current efficiency (cd/A): top-emitting OLED with a front ZnS/LiF DBR and a back diffuser
(scattering layer + Ag), stack and scattering layer optimised together.
Emitter: Ir(ppy)2acac spectrum, PLQY q 0.98, Theta 0.76, centre of a 25 nm EML; birefringent organics; Purcell effect
included (photons per electron = q sum p F eta / (1 - q + q sum p F)); charge balance and spin factor 1.
  air | DBR (N pairs, quarter-wave at lam_c, ZnS next to the IZO) | IZO d_If | HATCN 5 | TAPC d_H | EML 25 | B3PyMPM d_E |
  IZO d_Ib | scattering layer (index n_s; HG phase function g, optical thickness S, absorption mu_a d) | Ag 100
Back diffuser: Monte Carlo redistribution matrices B[o, i] (scat_ag_mc.py, 2-deg bins in the layer); recycling
  U = (I - B Rf)^-1 B P_s ; on-axis radiance (air) = sum_{0-6 deg} U Tf / (Omega n_s^2) + direct coherent part.
CE = 683 sum_l V(l) (1239.84 / l) x photons per electron per sr on axis.   Reference: glass n 1.52 | IZO 50 | HATCN 5 |
TAPC | EML | B3PyMPM | Ag 100 (bottom emission, flat glass), same emitter, thicknesses optimised for on-axis CE.
Usage: python3 ce_opt.py REF  |  python3 ce_opt.py B <N> <mua>      -> results_ce/*.json"""
import sys, os, json, numpy as np
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize
import front_lum as FL, cps2
D, IZO, ZNS, LIF = FL.D, FL.IZO, FL.ZNS, FL.LIF
Q, THETA, NG = 0.98, 0.76, 1.52
NB = 45; TB = np.arange(NB) * 2 + 1.0; OM = 2 * np.pi * (1 - np.cos(np.radians(6.0)))
BZ = np.load('/home/user/OLED-/sim/front_lum/scat_ag_B.npz')
SG = [(S, g) for S in (0.0, 0.3, 1.0, 3.0, 10.0, 30.0) for g in (0.0, 0.6, 0.9) if not (S == 0 and g > 0)]
I = FL.I; V = FL.V
def dbr(N, lamc, j):
    if N == 0: return []
    jc = int(round(lamc - 400)); return [(ZNS[j], ZNS[j], lamc / 4 / ZNS[jc].real), (LIF[j], LIF[j], lamc / 4 / LIF[jc].real)] * N
def point(args):
    case, j, x, N, npts = args
    lam = float(400 + j); EML = (D.EMo[j], D.EMe[j])
    if case == 'REF':
        dE, dH = x
        S = cps2.Stack(lam, EML, 25.0, 12.5, above=[(D.B3o[j], D.B3e[j], dE), (D.AG[j], D.AG[j], 100.0)],
                       below=[(D.TPo[j], D.TPe[j], dH), (D.HAT[j], D.HAT[j], 5.0), (IZO[j], IZO[j], 50.0)], n_sub=NG, h=THETA)
        r = cps2.solve_pol(S, npts=npts); P = float(np.real(r['P_tot']))
        Rg = ((NG - 1) / (NG + 1)) ** 2; Rl = float(np.real(cps2.stack_reflectance(S, [0.25]))[0])
        return P, FL.per_sr(S) / P * (1 - Rg) / NG ** 2 / (1 - Rg * Rl), float(np.real(r['air'])), None
    dE, dH, lamc, dIf, dIb, ns = x
    fr = [(D.TPo[j], D.TPe[j], dH), (D.HAT[j], D.HAT[j], 5.0), (IZO[j], IZO[j], dIf)] + dbr(N, lamc, j)
    back = [(D.B3o[j], D.B3e[j], dE), (IZO[j], IZO[j], dIb)]
    S2 = cps2.Stack(lam, EML, 25.0, 12.5, above=fr, below=back, n_sub=ns, n_top=1.0, h=THETA)
    r2 = cps2.solve_pol(S2, npts=npts); P = float(np.real(r2['P_tot']))
    Ps = np.real(cps2.sub_angular(S2, TB)) * 2.0 / P
    S1 = cps2.Stack(lam, EML, 25.0, 12.5, above=back, below=fr, n_sub=1.0, n_top=ns, h=THETA)
    Ld = FL.per_sr(S1) / P; Pair = float(np.sum(np.real(cps2.sub_angular(S1, np.arange(90) + 0.5)))) / P
    nsl = [l[0] for l in back[::-1]] + [D.EMo[j]] + [l[0] for l in fr]; dsl = [l[2] for l in back[::-1]] + [25.0] + [l[2] for l in fr]
    Tf = np.zeros(NB); Rf = np.zeros(NB)
    for i, th in enumerate(TB):
        kx = 2 * np.pi / lam * ns * np.sin(np.radians(th))
        for pol in 's', 'p':
            R_, T_ = FL.tmm(nsl, dsl, ns, 1.0, lam, kx, pol); Tf[i] += T_ / 2; Rf[i] += R_ / 2
    return P, Ld, Pair, (Ps, Tf, Rf, ns)
def combine(case, js, res, mua=None):
    lam = np.array(js) + 400.0; p = I[js] * lam; p = p / p.sum()
    F = np.array([r[0] for r in res]); den = 1 - Q + Q * (p @ F); k = Q * p * F / den          # photons per electron per lambda
    conv = 683 * V(lam) * 1239.84 / lam
    if case == 'REF':
        L = np.array([r[1] for r in res]); E = np.array([r[2] for r in res])
        return dict(CE=float(conv * k @ L), EQE=float(k @ E), F=float(p @ F))
    Ld = np.array([r[1] for r in res]); Pa = np.array([r[2] for r in res]); best = None
    for S, g in SG:
        B = BZ['S%g_g%g_a%g' % (S, g, mua)]; Lr = np.zeros(len(js)); Er = np.zeros(len(js))
        for t, r in enumerate(res):
            Ps, Tf, Rf, ns = r[3]
            U = np.linalg.solve(np.eye(NB) - B * Rf[None, :], B @ Ps)
            Lr[t] = (U[:3] * Tf[:3]).sum() / OM / ns ** 2; Er[t] = U @ Tf
        CE = float(conv * k @ (Ld + Lr))
        if best is None or CE > best['CE']:
            best = dict(CE=CE, CE_direct=float(conv * k @ Ld), EQE=float(k @ (Pa + Er)), S=S, g=g, F=float(p @ F))
    return best
def evaluate(case, x, N=0, mua=0.01, step=20, npts=1200, mapper=map, muas=None):
    js = list(range(0, 301, step)); res = list(mapper(point, [(case, j, tuple(x), N, npts) for j in js]))
    if case == 'REF': return combine(case, js, res)
    if muas: return {m: combine(case, js, res, m) for m in muas}
    return combine(case, js, res, mua)
class Obj:
    def __init__(self, *a): self.a = a
    def __call__(self, x):
        case, N, mua = self.a; xx = list(x)
        if case == 'B' and N == 0: xx = xx[:2] + [550.0] + xx[2:]
        return -evaluate(case, xx, N, mua)['CE']
if __name__ == '__main__':
    case = sys.argv[1]; N = int(sys.argv[2]) if len(sys.argv) > 2 else 0; mua = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01
    os.makedirs('results_ce', exist_ok=True)
    if case == 'REF': bd = [(20.0, 400.0), (20.0, 400.0)]
    else: bd = [(20.0, 400.0), (20.0, 400.0)] + ([(440.0, 800.0)] if N else []) + [(30.0, 150.0), (30.0, 150.0), (1.5, 2.1)]
    lo = np.array([b[0] for b in bd]); hi = np.array([b[1] for b in bd])
    full = lambda x: (list(x[:2]) + [550.0] + list(x[2:])) if case == 'B' and N == 0 else list(x)
    with Pool(4) as pool:
        r = differential_evolution(Obj(case, N, mua), bd, seed=1, popsize=8, maxiter=30, tol=1e-4, polish=False, workers=pool.map, updating='deferred')
        m = minimize(lambda x: -evaluate(case, full(np.clip(x, lo, hi)), N, mua, step=10, mapper=pool.map)['CE'], r.x, method='Nelder-Mead',
                     options=dict(maxiter=60 * len(bd), xatol=0.5, fatol=1e-4))
        xb = full(np.clip(m.x, lo, hi))
        fin = evaluate(case, xb, N, step=5, npts=3000, mapper=pool.map, muas=(0.0, 0.01, 0.03) if case == 'B' else None)
    names = ['d_ETL', 'd_HTL'] if case == 'REF' else ['d_ETL', 'd_HTL', 'lam_c', 'd_IZO_front', 'd_IZO_back', 'n_s']
    out = dict(case=case, N=N, mua_opt=mua, structure={k: round(float(v), 2) for k, v in zip(names, xb)},
               results=({str(k): v for k, v in fin.items()} if case == 'B' else fin))
    fn = 'results_ce/%s.json' % ('REF' if case == 'REF' else 'B_N%d_mua%g' % (N, mua))
    json.dump(out, open(fn, 'w'), indent=1); print(json.dumps(out), flush=True)
