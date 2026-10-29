"""Fabricable variants of the angle-selective recycling OLED (top emission), co-optimised for on-axis current efficiency.
Common: substrate | Ag 100 | back diffuser (S <= 10, g, mu_a d 0.01) + high-index planarizer (n_s 1.6-1.9) | bottom electrode
IZO d_Ib (30-80) | organics | top electrode | ZnS/LiF DBR (N pairs, quarter-wave at lam_c) | air.  Emitter as ce_opt
(Ir(ppy)2acac, PLQY 0.98, Theta 0.76, EML 25 nm, Purcell).  Transport layers 30-250 nm (doped).
  F1 (conventional top emission): IZO | HATCN 5 | TAPC d_B | EML | B3PyMPM d_F | HATCN 5 (seed) | Ag d_Ag 8-12 (measured
      HATCN-seeded thin-Ag n,k) | DBR
  F2 (inverted, sputtered top IZO): IZO | B3PyMPM d_B | EML | TAPC d_F | MoO3 15 | IZO d_If 30-80 | DBR
  TREF (conventional top-emitting AMOLED reference, no diffuser/DBR): Ag 100 | IZO 10 | HATCN 5 | TAPC d_B | EML | B3PyMPM d_F |
      HATCN 5 | Ag d_Ag | capping layer (n 1.9) d_C 20-120 | air
Usage: python3 ce_fab.py TREF | python3 ce_fab.py F1|F2 <N>   -> results_fab/*.json"""
import sys, os, json, numpy as np, scipy.io as sio
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize
import ce_opt as C, front_lum as FL, cps2
D, IZO, ZNS, LIF = C.D, C.IZO, C.ZNS, C.LIF
Mt = sio.loadmat('/home/user/OLED-/sim/table1_author/nk_JH_total.mat', squeeze_me=True, struct_as_record=False)['material']
MOO3 = np.asarray(Mt.l_MoO3_JH)
L = np.arange(400, 701)
def joint(f):
    a = np.genfromtxt('/home/user/OLED-/sim/teoled_opt/nk_joint/' + f, delimiter=',', names=True)
    return np.interp(L, a['wavelength_nm'], a['n_joint']) + 1j * np.interp(L, a['wavelength_nm'], a['k_joint'])
AGT = {6.0: joint('Ag_HATCN_1-8_6nm_nk.csv'), 7.0: joint('Ag_HATCN_2-5_7nm_nk.csv'), 8.0: joint('Ag_HATCN_2-6_8nm_nk.csv'),
       10.0: joint('Ag_HATCN_2-7_10nm_nk.csv'), 12.0: joint('Ag_HATCN_2-8_12nm_nk.csv')}
AGD = sorted(AGT)
def agthin(d, j):                       # measured HATCN-seeded thin-Ag n,k, interpolated in thickness (6-12 nm)
    d = min(max(d, AGD[0]), AGD[-1]); i = max(k for k in range(len(AGD) - 1) if AGD[k] <= d) if d < AGD[-1] else len(AGD) - 2
    lo, hi = AGD[i], AGD[i + 1]; t = (d - lo) / (hi - lo)
    return (1 - t) * AGT[lo][j] + t * AGT[hi][j]
C.SG = [(S, g) for S in (0.3, 1.0, 3.0, 10.0) for g in (0.0, 0.6, 0.9)]          # realistic diffusers (S <= 10)
def layers(case, j, x):
    if case == 'F1':
        dF, dB, lamc, dAg, dIb, ns = x
        fr = [(D.B3o[j], D.B3e[j], dF), (D.HAT[j], D.HAT[j], 5.0), (agthin(dAg, j),) * 2 + (dAg,)]
        back = [(D.TPo[j], D.TPe[j], dB), (D.HAT[j], D.HAT[j], 5.0), (IZO[j], IZO[j], dIb)]
    else:
        dF, dB, lamc, dIf, dIb, ns = x
        fr = [(D.TPo[j], D.TPe[j], dF), (MOO3[j], MOO3[j], 15.0), (IZO[j], IZO[j], dIf)]
        back = [(D.B3o[j], D.B3e[j], dB), (IZO[j], IZO[j], dIb)]
    return fr, back, lamc, ns
def point(args):
    case, j, x, N, npts = args
    lam = float(400 + j); EML = (D.EMo[j], D.EMe[j])
    if case == 'TREF':
        dF, dB, dAg, dC = x
        S = cps2.Stack(lam, EML, 25.0, 12.5,
                       above=[(D.TPo[j], D.TPe[j], dB), (D.HAT[j], D.HAT[j], 5.0), (IZO[j], IZO[j], 10.0), (D.AG[j], D.AG[j], 100.0)],
                       below=[(D.B3o[j], D.B3e[j], dF), (D.HAT[j], D.HAT[j], 5.0), (agthin(dAg, j),) * 2 + (dAg,), (1.9, 1.9, dC)], n_sub=1.0, h=C.THETA)
        r = cps2.solve_pol(S, npts=npts); P = float(np.real(r['P_tot']))
        return P, FL.per_sr(S) / P, float(np.real(r['air'] + r['sub'])), None
    fr, back, lamc, ns = layers(case, j, x); fr = fr + C.dbr(N, lamc, j)
    S2 = cps2.Stack(lam, EML, 25.0, 12.5, above=fr, below=back, n_sub=ns, n_top=1.0, h=C.THETA)
    r2 = cps2.solve_pol(S2, npts=npts); P = float(np.real(r2['P_tot']))
    Ps = np.real(cps2.sub_angular(S2, C.TB)) * 2.0 / P
    S1 = cps2.Stack(lam, EML, 25.0, 12.5, above=back, below=fr, n_sub=1.0, n_top=ns, h=C.THETA)
    Ld = FL.per_sr(S1) / P; Pair = float(np.sum(np.real(cps2.sub_angular(S1, np.arange(90) + 0.5)))) / P
    nsl = [l[0] for l in back[::-1]] + [D.EMo[j]] + [l[0] for l in fr]; dsl = [l[2] for l in back[::-1]] + [25.0] + [l[2] for l in fr]
    Tf = np.zeros(C.NB); Rf = np.zeros(C.NB)
    for i, th in enumerate(C.TB):
        kx = 2 * np.pi / lam * ns * np.sin(np.radians(th))
        for pol in 's', 'p':
            R_, T_ = FL.tmm(nsl, dsl, ns, 1.0, lam, kx, pol); Tf[i] += T_ / 2; Rf[i] += R_ / 2
    return P, Ld, Pair, (Ps, Tf, Rf, ns)
def evaluate(case, x, N=0, step=20, npts=1200, mapper=map, mua=0.01):
    js = list(range(0, 301, step)); res = list(mapper(point, [(case, j, tuple(x), N, npts) for j in js]))
    return C.combine('REF' if case == 'TREF' else 'B', js, res, mua)
class Obj:
    def __init__(self, *a): self.a = a
    def __call__(self, x): return -evaluate(self.a[0], list(x), self.a[1])['CE']
if __name__ == '__main__':
    case = sys.argv[1]; N = int(sys.argv[2]) if len(sys.argv) > 2 else 0; os.makedirs('results_fab', exist_ok=True)
    if case == 'TREF': bd = [(30.0, 250.0), (30.0, 250.0), (8.0, 12.0), (20.0, 120.0)]; names = ['d_ETL', 'd_HTL', 'd_Ag', 'd_CPL']
    elif case == 'F1': bd = [(30.0, 250.0), (30.0, 250.0), (440.0, 800.0), (6.0, 12.0), (30.0, 80.0), (1.6, 1.9)]; names = ['d_ETL', 'd_HTL', 'lam_c', 'd_Ag', 'd_IZO_b', 'n_s']
    else: bd = [(30.0, 250.0), (30.0, 250.0), (440.0, 800.0), (30.0, 80.0), (30.0, 80.0), (1.6, 1.9)]; names = ['d_HTL', 'd_ETL', 'lam_c', 'd_IZO_f', 'd_IZO_b', 'n_s']
    lo = np.array([b[0] for b in bd]); hi = np.array([b[1] for b in bd])
    with Pool(4) as pool:
        r = differential_evolution(Obj(case, N), bd, seed=1, popsize=8, maxiter=30, tol=1e-4, polish=False, workers=pool.map, updating='deferred')
        m = minimize(lambda x: -evaluate(case, list(np.clip(x, lo, hi)), N, step=10, mapper=pool.map)['CE'], r.x, method='Nelder-Mead',
                     options=dict(maxiter=60 * len(bd), xatol=0.3, fatol=1e-4))
        xb = list(np.clip(m.x, lo, hi))
        fin = {m_: evaluate(case, xb, N, step=5, npts=3000, mapper=pool.map, mua=m_) for m_ in ((0.0, 0.01, 0.03) if case != 'TREF' else (0.01,))}
    out = dict(case=case, N=N, structure={k: round(float(v), 2) for k, v in zip(names, xb)}, results={str(k): v for k, v in fin.items()})
    fn = 'results_fab/%s.json' % (case if case == 'TREF' else '%s_N%d' % (case, N))
    json.dump(out, open(fn, 'w'), indent=1); print(json.dumps(out), flush=True)
