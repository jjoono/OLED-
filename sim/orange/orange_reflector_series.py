"""Orange preprint device (Al cathode, ETL 55 nm): auxiliary-reflector series measured as
none 51% | Al mirror behind an air gap 61% | DBR 2.5 pairs ~74% | DBR 4.5 pairs 77%.
Back reflectance outside the pixel (substrate side, unpolarised, 400-700 nm):
  none : glass | IZO | organics | air
  Al   : the same stack, then an air gap (incoherent) and an Al mirror (Al_JO, Fresnel from air)
  DBRn : the same stack + ZnS 70 / LiF 115 alternating, starting and ending with ZnS (2.5 pairs = 5 layers, 4.5 = 9)
A'_aux = 1 - <R_back> weighted by the substrate-delivered P_sub(lambda, theta); EQE from the finite-pixel Monte Carlo
(3 x 3 mm pixel, 24 x 24 mm substrate, thickness T) of orange_bounds.mcg."""
import numpy as np, sys, csv
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/orange']
import orange_reflectors as O
import orange_bounds as OB
import orange_ag_cathode as A

MEAS = {'none': 51.0, 'Al': 61.0, 'DBR2.5': 74.0, 'DBR4.5': 77.0}


def rt_pol(layers, lam, pol, n_in=O.NS):
    """R, T (power) vs internal angle O.TH from the substrate side, exit medium air"""
    neff = n_in * np.sin(np.radians(O.TH))
    def eta(n):
        c = np.sqrt(1 - (neff / n) ** 2 + 0j); c = np.where(np.imag(c) > 0, np.conj(c), c)
        return (n * c if pol == 's' else n / c)
    Mx = np.broadcast_to(np.eye(2, dtype=complex), neff.shape + (2, 2)).copy()
    for n, d in layers:
        n = np.conj(complex(n)); e = eta(n)
        c = np.sqrt(1 - (neff / n) ** 2 + 0j); c = np.where(np.imag(c) > 0, np.conj(c), c)
        dl = 2 * np.pi * n * c * d / lam
        L = np.empty_like(Mx); L[..., 0, 0] = np.cos(dl); L[..., 0, 1] = 1j * np.sin(dl) / e
        L[..., 1, 0] = 1j * e * np.sin(dl); L[..., 1, 1] = np.cos(dl); Mx = Mx @ L
    eo, ei = eta(1.0 + 0j), eta(complex(n_in))
    Bq = Mx[..., 0, 0] + Mx[..., 0, 1] * eo; Cq = Mx[..., 1, 0] + Mx[..., 1, 1] * eo
    R = np.abs((ei * Bq - Cq) / (ei * Bq + Cq)) ** 2
    T = 4 * np.real(ei) * np.real(eo) / np.abs(ei * Bq + Cq) ** 2
    return np.clip(R, 0, 1), np.clip(np.nan_to_num(T), 0, 1)


def r_al_air(lam, j, pol):
    sa = np.clip(O.NS * np.sin(np.radians(O.TH)), 0, 0.999999); ca = np.sqrt(1 - sa ** 2)
    m = np.conj(complex(O.AL[j])); cm = np.sqrt(1 - (sa / m) ** 2 + 0j)
    r = (ca - m * cm) / (ca + m * cm) if pol == 's' else (m * ca - cm) / (m * ca + cm)
    return np.abs(r) ** 2


def back_R(j, kind, d=55):
    lam = float(O.LAM[j]); L = A.organics(j, d)
    if kind.startswith('DBR'):
        nl = {'DBR2.5': 5, 'DBR4.5': 9}[kind]
        L = L + [(O.ZNS[j], 70.0) if k % 2 == 0 else (O.LIF[j], 115.0) for k in range(nl)]
    out = 0
    for pol in ('s', 'p'):
        R, T = rt_pol(L, lam, pol)
        if kind == 'Al':                       # incoherent air gap + Al mirror; stack reflectance from air taken as R
            Ra = r_al_air(lam, j, pol); R = R + T * T * Ra / (1 - R * Ra)
        out = out + 0.5 * R
    return np.clip(out, 0, 1)


if __name__ == '__main__':
    T = float(sys.argv[1]) if len(sys.argv) > 1 else 0.7
    wp = (O.LAM * O.SPEC)[O.LAMJ]; wp = wp / wp.sum()
    with Pool(4) as p: res = p.map(A.device, [('Al', 55, j) for j in O.LAMJ])
    ES = np.array([r[0] for r in res]); PS = np.array([r[1] for r in res]); RS = np.array([r[2] for r in res])
    es = float(wp @ ES); WL = wp * ES / (wp @ ES); Pn = PS / PS.sum(1, keepdims=True)
    rows = []
    for kind in MEAS:
        RB = np.array([back_R(j, kind) for j in O.LAMJ]); AB = float(WL @ (1 - (Pn * RB).sum(1)))
        with Pool(4) as p:
            e = float(np.mean(p.map(OB.mcg, [(10 * k + 1, 120000, 1.5, T, PS, RS, RB, WL) for k in range(8)])))
        rows.append(dict(aux=kind, T_mm=T, A_aux=AB, eta_sub=es, eta_ext=e, EQE_model=100 * O.PLQY * es * e, EQE_meas=MEAS[kind]))
        print('%-6s A\'aux %.3f | model EQE %.1f%% | measured %.1f%%' % (kind, AB, 100 * O.PLQY * es * e, MEAS[kind]), flush=True)
    with open('/home/user/OLED-/sim/orange/orange_reflector_series_T%.1f.csv' % T, 'w', newline='') as fh:
        wr = csv.DictWriter(fh, list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
