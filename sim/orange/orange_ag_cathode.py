"""Orange preprint device with the Al cathode replaced by Ag 100 nm, ETL (B3PyMPM) thickness swept.
Stack and models as orange_reflectors.py (IZO 50 | HATCN/TAPC x2 | TCTA 10 | EML 25 | B3PyMPM d_ETL | Liq 1 | cathode 100);
the organics (incl. the ETL) also cover the area outside the pixel under the auxiliary reflector (DBR 4.5 pairs ZnS 70 / LiF 115,
unchanged, or Ag 100). Finite geometry: 3 x 3 mm pixel (9 mm^2, as fabricated) on a 24 x 24 mm substrate, thickness T."""
import numpy as np, sys, csv
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/orange']
import orange_reflectors as O
import orange_bounds as OB

ETLS = (55, 80, 110, 150, 200)


def organics(j, d):
    L = O.organics(j); L[7] = (O.B3o[j], float(d)); return L


def back_R(j, kind, d):
    lam = float(O.LAM[j]); L = organics(j, d)
    if kind == 'DBR':
        for k in range(9): L = L + [(O.ZNS[j], 70.0) if k % 2 == 0 else (O.LIF[j], 115.0)]
    elif kind == 'Ag':
        L = L + [(O.AG[j], 100.0)]
    return O.rt_unpol(L, lam)


def device(args):
    cat, d, j = args
    m = O.AL[j] if cat == 'Al' else O.AG[j]
    S = O.cps2.Stack(float(O.LAM[j]), (O.EMo[j], O.EMe[j]), 25.0, 12.5,
                     above=[(O.B3o[j], O.B3e[j], float(d)), (1.70, 1.70, 1.0), (m, m, 100.0)],
                     below=[(O.TCo[j], O.TCe[j], 10.0), (O.TPo[j], O.TPe[j], 95.0), (O.HAT[j], O.HAT[j], 5.0),
                            (O.TPo[j], O.TPe[j], 95.0), (O.HAT[j], O.HAT[j], 5.0), (O.IZO[j], O.IZO[j], 50.0)],
                     n_sub=O.NS, h=O.THETA)
    r = O.cps2.solve_pol(S, npts=3000)
    return float(np.real(r['air'] + r['sub'])), np.real(O.cps2.sub_angular(S, O.TH)), np.real(O.cps2.stack_reflectance(S, O.TH))


if __name__ == '__main__':
    wp = (O.LAM * O.SPEC)[O.LAMJ]; wp = wp / wp.sum()
    T = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
    rows = []
    for cat, d in [('Al', 55)] + [('Ag', d) for d in ETLS]:
        with Pool(4) as p: res = p.map(device, [(cat, d, j) for j in O.LAMJ])
        ES = np.array([r[0] for r in res]); PS = np.array([r[1] for r in res]); RS = np.array([r[2] for r in res])
        es = float(wp @ ES); WL = wp * ES / (wp @ ES)
        Pn = PS / PS.sum(1, keepdims=True); Ap = float(WL @ (1 - (Pn * RS).sum(1)))
        for kind in ('DBR', 'Ag'):
            RB = np.array([back_R(j, kind, d) for j in O.LAMJ]); AB = float(WL @ (1 - (Pn * RB).sum(1)))
            with Pool(4) as p:
                e = float(np.mean(p.map(OB.mcg, [(10 * k + 1, 120000, 1.5, T, PS, RS, RB, WL) for k in range(8)])))
            rows.append(dict(cathode=cat, d_ETL=d, aux=kind, T_mm=T, eta_sub=es, A_pixel=Ap, A_aux=AB, eta_ext=e,
                             EQE=O.PLQY * es * e, EQE_PLQY1=es * e))
            print('%s ETL %3d aux %-3s | eta_sub %.3f  A\'pix %.3f  A\'aux %.3f | eta_ext %.3f  EQE %.3f'
                  % (cat, d, kind, es, Ap, AB, e, O.PLQY * es * e), flush=True)
    with open('/home/user/OLED-/sim/orange/orange_ag_cathode_T%.1f.csv' % T, 'w', newline='') as fh:
        wr = csv.DictWriter(fh, list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
