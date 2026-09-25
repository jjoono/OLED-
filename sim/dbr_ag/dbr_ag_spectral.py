"""Reflector comparison weighted by what the device really delivers: the Ir(ppy)2acac spectrum
(photon number) and the angle- and wavelength-resolved substrate-delivered power P_sub(theta, lam)
of a dipole in the Fig. 2 generic stack, then the matrix series with the hemispherical MLA
(LightTools BSDF, n = 1.5).  The series is split pass by pass into what leaves through the lens,
what the stack absorbs and what it transmits out of the back.
Stack: glass 1.5 / ITO 150 / HTL 200 / EML 20 (dipole at centre) / ETL 200 (all n = 1.8) / reflector."""
import numpy as np, scipy.io as sio, sys, os, csv
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.join(H, '..', 'design_rule4'), os.path.join(H, '..', 'fig3c')]
import fig2d as F, cps2
LAM = np.arange(430.0, 701.0, 5.0); IX = [int(np.argmin(abs(F.LAM - l))) for l in LAM]
TH = np.arange(90) + 0.5
B = sio.loadmat(os.path.join(H, '..', 'mla', 'lt_hemisphere_bsdf.mat'))['BSDF_MLA'][:, :, 4]   # n_MLA = 1.5
BT = B[:90].sum(0); BR = B[90:180][::-1]
SPEC = np.clip(np.interp(LAM, F.LAM, F.GREEN), 0, None) * LAM; SPEC /= SPEC.sum()   # photon number
def refl(j):
    ito, ag, al, zns, lif = F.ITO[j], F.AG[j], F.AL[j], F.ZNS[j], F.LIF[j]
    ch = lambda n, dz, dl, f: [(zns if k % 2 == 0 else lif, d) for k, (_, d) in enumerate(F.dbr_chirp(n, dz, dl, f)) if d > 0]
    return {'Al': [(al, 100.0)], 'Ag': [(ag, 100.0)],
            'ITO50 + chirped DBR (air)': [(ito, 50.0)] + ch(10, 55.0, 91.0, 1.40),
            'ITO50 + chirped DBR + Ag': [(ito, 50.0)] + ch(10, 55.0, 91.0, 1.40) + [(ag, 100.0)],
            'ITO50 + 2 pairs + Ag (re-opt)': [(ito, 50.0)] + ch(2, 50.0, 100.0, 1.0) + [(ag, 100.0)]}
NAMES = list(refl(150).keys())
def one(a):
    name, i = a; j = IX[i]; lam = LAM[i]; ref = refl(j)[name]
    S = cps2.Stack(lam, (1.8, 1.8), 20.0, 10.0, above=[(1.8, 1.8, 200.0)] + [(n, n, d) for n, d in ref],
                   below=[(1.8, 1.8, 200.0), (F.ITO[j], F.ITO[j], 150.0)], n_sub=1.5, n_top=1.0)
    r = cps2.solve_pol(S, npts=6000); es = float(np.real(r['air'] + r['sub']))
    P = np.real(cps2.sub_angular(S, TH)); R = np.real(cps2.stack_reflectance(S, TH))
    lay = [(np.array([F.ITO[j]]), 150.0), (np.array([1.8+0j]), 420.0)] + [(np.array([n]), d) for n, d in ref]
    _, T = F.RT(lay, np.array([1.5+0j]), np.array([1.0+0j]), np.radians(TH), np.array([lam])); T = T[0]
    v = P / P.sum(); out = ab = tr = 0.0
    for _ in range(200):
        out += BT @ v; w = BR @ v; tr += T @ w; ab += (1 - R - T) @ w; v = R * w
    Pn = P / P.sum()
    return name, i, es, out, ab, tr, Pn @ (1 - R - T), Pn @ T
if __name__ == '__main__':
    jobs = [(n, i) for n in NAMES for i in range(len(LAM))]
    with Pool(4) as p: res = p.map(one, jobs)
    rows = [['reflector', 'eta_sub', 'Aprime_first_pass_absorbed_%', 'Aprime_first_pass_transmitted_%', 'Aprime_first_pass_%',
             'final_absorbed_%', 'final_transmitted_%', 'eta_ext', 'EQE']]
    for n in NAMES:
        R_ = sorted([r for r in res if r[0] == n], key=lambda r: r[1]); a = np.array([r[2:] for r in R_])
        es, out, ab, tr, a1, t1 = a.T; wsub = SPEC * es; ws = wsub / wsub.sum()
        rows.append([n, SPEC @ es, 100*ws @ a1, 100*ws @ t1, 100*ws @ (a1+t1), 100*ws @ ab, 100*ws @ tr, ws @ out, SPEC @ (es*out)])
        print('%-30s eta_sub %.3f | first pass A: abs %.2f tr %.2f tot %.2f | final: abs %.2f tr %.2f | eta_ext %.3f EQE %.3f' % tuple([n] + rows[-1][1:]))
    with open(os.path.join(H, 'dbr_ag_spectral.csv'), 'w', newline='') as f: csv.writer(f).writerows(rows)
