"""Inverted bottom-emitting green OLED vs the conventional stack: optimum thicknesses for three outcoupling cases.
  INV : substrate | ITO d_ITO | B3PyMPM:Cs2CO3 d_n | B3PyMPM 20 | TCTA:B3PyMPM:Ir(ppy)2acac 25 | TCTA 10 | TAPC d_H | HATCN 10 | Ag 100
  NORM: substrate | ITO d_ITO | HATCN 5 | TAPC d_H | TCTA 10 | EML 25 | B3PyMPM d_E (undoped 20 + doped) | Ag 100
Cases: 'planar' (glass n 1.52, flat), 'mla15' (glass + hemispherical MLA film, n 1.5, BSDF lt_hemisphere),
       'mla177' (HNPS n 1.77 + close-packed hemispherical MLA, full back reflector = OLED stack reflectance).
Emitter Ir(ppy)2acac (photon spectrum), Theta 0.76, PLQY 1, dipoles at the EML centre; measured birefringent organics,
ITO l_ITO_SNU_temp, Ag McPeak.  EQE = sum_l w eta_sub eta_ext (matrix series).
Usage: python3 inv_opt.py <INV|NORM> <planar|mla15|mla177> [fixITO]   -> results/*.json"""
import sys, os, json, numpy as np, scipy.io as sio
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize
sys.path[:0] = ['/home/user/OLED-/sim/be_cone', '/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/mla']
import device_eqe as D, cps2, series
TH = np.arange(90) + 0.5; THETA = 0.76
def flat(ns):
    c = np.cos(np.radians(TH)); s2 = (ns * np.sin(np.radians(TH))) ** 2; ct = np.sqrt(np.clip(1 - s2, 0, 1))
    Rf = np.where(s2 < 1, 0.5 * (((ns * c - ct) / (ns * c + ct)) ** 2 + ((ns * ct - c) / (ns * ct + c)) ** 2), 1.0); return 1 - Rf, np.diag(Rf)
Bh = sio.loadmat('/home/user/OLED-/sim/mla/lt_hemisphere_bsdf.mat', squeeze_me=True)['BSDF_MLA'][:, :, 4]     # n 1.50
EXT = {'planar': (1.52, flat(1.52)), 'mla15': (1.50, (Bh[:90].sum(0), Bh[90:180][::-1])), 'mla177': (1.77, D.EXT['hemisphere (close-packed)'])}
def stack(case, j, x, ns):
    dn, dH, dITO = x
    EML = (D.EMo[j], D.EMe[j]); ito = (D.ITO[j], D.ITO[j], dITO)
    if case == 'INV':
        below = [(D.B3o[j], D.B3e[j], 20.0), (D.B3o[j], D.B3e[j], dn), ito]
        above = [(D.TCo[j], D.TCe[j], 10.0), (D.TPo[j], D.TPe[j], dH), (D.HAT[j], D.HAT[j], 10.0), (D.AG[j], D.AG[j], 100.0)]
    else:
        below = [(D.TCo[j], D.TCe[j], 10.0), (D.TPo[j], D.TPe[j], dH), (D.HAT[j], D.HAT[j], 5.0), ito]
        above = [(D.B3o[j], D.B3e[j], 20.0 + dn), (D.AG[j], D.AG[j], 100.0)]
    return cps2.Stack(float(D.LAM[j]), EML, 25.0, 12.5, above=above, below=below, n_sub=ns, h=THETA)
def point(a):
    case, ext, j, x, npts = a
    ns, (BT, BR) = EXT[ext]; S = stack(case, j, x, ns); r = cps2.solve_pol(S, npts=npts)
    es = float(np.real(r['air'] + r['sub']))
    P = np.real(cps2.sub_angular(S, TH)); R = np.real(cps2.stack_reflectance(S, TH))
    return es, series.eta_ext(BT, BR, R, P, 300)[0], float(np.real(r['spp'])), float(np.real(r['abs']))
def evaluate(case, ext, x, step=10, npts=1500, mapper=map):
    js = list(range(0, 301, step)); res = np.array(list(mapper(point, [(case, ext, j, tuple(x), npts) for j in js])))
    w = D.WP[js] / D.WP[js].sum(); es, ee, spp, ab = res.T
    return dict(EQE=float(w @ (es * ee)), eta_sub=float(w @ es), eta_ext=float(w @ (es * ee)) / float(w @ es), spp=float(w @ spp), absorb=float(w @ ab))
class Obj:
    def __init__(self, *a): self.a = a
    def __call__(self, x):
        case, ext, fix = self.a; xx = list(x) + ([150.0] if fix else [])
        return -evaluate(case, ext, xx)['EQE']
if __name__ == '__main__':
    case, ext = sys.argv[1], sys.argv[2]; fix = len(sys.argv) > 3
    bd = [(10.0, 250.0), (20.0, 400.0)] + ([] if fix else [(30.0, 150.0)])
    with Pool(4) as pool:
        r = differential_evolution(Obj(case, ext, fix), bd, seed=1, popsize=12, maxiter=30, tol=1e-5, polish=False, workers=pool.map, updating='deferred')
        lo = np.array([b[0] for b in bd]); hi = np.array([b[1] for b in bd])
        f = lambda x: list(np.clip(x, lo, hi)) + ([150.0] if fix else [])
        m = minimize(lambda x: -evaluate(case, ext, f(x), step=5, mapper=pool.map)['EQE'], r.x, method='Nelder-Mead', options=dict(maxiter=150, xatol=0.5, fatol=1e-5))
        xb = f(m.x); fin = evaluate(case, ext, xb, step=2, npts=3000, mapper=pool.map)
    names = (['d_ETL_doped', 'd_HTL_TAPC', 'd_ITO'] if case == 'INV' else ['d_ETL_doped', 'd_HTL_TAPC', 'd_ITO'])
    out = dict(case=case, ext=ext, fixITO=fix, structure={k: round(float(v), 1) for k, v in zip(names, xb)}, **{k: round(v, 4) for k, v in fin.items()})
    os.makedirs('results', exist_ok=True); json.dump(out, open('results/%s_%s%s.json' % (case, ext, '_ITO150' if fix else ''), 'w'), indent=1)
    print(json.dumps(out), flush=True)
