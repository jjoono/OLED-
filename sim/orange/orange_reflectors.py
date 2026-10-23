"""Preprint orange device on the high-index MLA substrate, full-back auxiliary reflector: DBR (4.5 pairs ZnS 70 / LiF 115)
vs Ag 100 nm, with the Al cathode (as fabricated) and with an Ag cathode.
  substrate n 1.77 + close-packed hemisphere MLA (n 1.77; own ray tracer - the HNPS01 BSDF is not in the repo)
  | IZO 50 (l_IZO_SJ_RTP) | HATCN 5 | TAPC 95 | HATCN 5 | TAPC 95 | TCTA 10 | TCTA:B3PyMPM:Ir(dmppy-ph)2tmd 25 | B3PyMPM 55
  | Liq 1 | Al (or Ag) 100 | air;   Theta 0.865, PLQY 0.98 (results also for PLQY 1), l_I_Irdmppyph2tmd, 400-700 nm.
Outside the pixel the back is: glass | IZO 50 | the organic stack | (nothing / DBR / Ag 100) | air.
(1) infinite pixel (matrix series; auxiliary reflector irrelevant because the cathode is opaque);
(2) finite: 24 x 24 x 0.5 mm substrate, 2 x 2 mm pixel at the centre, absorbing side walls (Monte Carlo, as be_cone)."""
import numpy as np, scipy.io as sio, sys, csv
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/mla', '/home/user/OLED-/sim/be_cone']
import cps2, series

M = sio.loadmat('/home/user/OLED-/sim/table1_author/nk_JH_total.mat', squeeze_me=True, struct_as_record=False)
mat, spec = M['material'], M['spectrum']
LAM = np.arange(400, 701)
def lib(f): return np.asarray(getattr(mat, f))[:301]
IZO, HAT, AG, AL, ZNS, LIF = lib('l_IZO_SJ_RTP'), lib('l_HATCN'), lib('l_Ag_McPeak'), lib('l_Al_JO'), lib('l_ZnS'), lib('l_LiF')
TPo, TPe, TCo, TCe = lib('l_TAPC_o_JO'), lib('l_TAPC_e_JO'), lib('l_TCTA_o_JO'), lib('l_TCTA_e_JO')
EMo, EMe, B3o, B3e = lib('l_TCTA_B3_o_JO'), lib('l_TCTA_B3_e_JO'), lib('l_B3_o_JO'), lib('l_B3_e_JO')
SPEC = np.clip(np.asarray(spec.l_I_Irdmppyph2tmd, float)[:301], 0, None)
NS, THETA, PLQY = 1.77, 0.865, 0.98
TH = np.arange(90) + 0.5
C = np.load('/home/user/OLED-/sim/mla/cone_bsdf.npz'); B = C['hemi_cp_n1.77']; BT, BR = B[:90].sum(0), B[90:180][::-1]
BRn = BR / np.maximum(BR.sum(0), 1e-12); BRc = np.cumsum(BRn, 0)
LAMJ = list(range(0, 301, 10))


def rt_unpol(layers, lam, n_in=NS):
    """unpolarised R vs internal angle from the substrate side; layers [(n, d)] from the substrate outwards, air behind"""
    neff = n_in * np.sin(np.radians(TH)); out = 0
    for pol in ('s', 'p'):
        def eta(n):
            c = np.sqrt(1 - (neff / n) ** 2 + 0j); c = np.where(np.imag(c) > 0, np.conj(c), c)
            return (n * c if pol == 's' else n / c), c
        Mx = np.broadcast_to(np.eye(2, dtype=complex), neff.shape + (2, 2)).copy()
        for n, d in layers:
            n = np.conj(complex(n)); e, c = eta(n); dl = 2 * np.pi * n * c * d / lam
            L = np.empty_like(Mx); L[..., 0, 0] = np.cos(dl); L[..., 0, 1] = 1j * np.sin(dl) / e
            L[..., 1, 0] = 1j * e * np.sin(dl); L[..., 1, 1] = np.cos(dl); Mx = Mx @ L
        eo, _ = eta(1.0 + 0j); ei, _ = eta(complex(n_in))
        Bq = Mx[..., 0, 0] + Mx[..., 0, 1] * eo; Cq = Mx[..., 1, 0] + Mx[..., 1, 1] * eo
        out = out + 0.5 * np.abs((ei * Bq - Cq) / (ei * Bq + Cq)) ** 2
    return np.clip(out, 0, 1)


def organics(j):   # from the IZO outwards (n_o for the planar reflectance)
    return [(IZO[j], 50), (HAT[j], 5), (TPo[j], 95), (HAT[j], 5), (TPo[j], 95), (TCo[j], 10), (EMo[j], 25), (B3o[j], 55), (1.70, 1)]


def back_R(j, kind):
    lam = float(LAM[j]); L = organics(j)
    if kind == 'DBR':
        for k in range(9): L = L + [(ZNS[j], 70.0) if k % 2 == 0 else (LIF[j], 115.0)]
    elif kind == 'Ag':
        L = L + [(AG[j], 100.0)]
    return rt_unpol(L, lam)


def device(args):
    cat, j = args
    m = AL[j] if cat == 'Al' else AG[j]
    S = cps2.Stack(float(LAM[j]), (EMo[j], EMe[j]), 25.0, 12.5, above=[(B3o[j], B3e[j], 55.0), (1.70, 1.70, 1.0), (m, m, 100.0)],
                   below=[(TCo[j], TCe[j], 10.0), (TPo[j], TPe[j], 95.0), (HAT[j], HAT[j], 5.0), (TPo[j], TPe[j], 95.0),
                          (HAT[j], HAT[j], 5.0), (IZO[j], IZO[j], 50.0)], n_sub=NS, h=THETA)
    r = cps2.solve_pol(S, npts=3000)
    return float(np.real(r['air'] + r['sub'])), np.real(cps2.sub_angular(S, TH)), np.real(cps2.stack_reflectance(S, TH))


def mc(args):
    seed, N, ES, PS, RS, RB, WL = args
    rng = np.random.default_rng(seed)
    j = rng.choice(len(WL), N, p=WL); Pc = np.cumsum(PS / PS.sum(1, keepdims=True), 1)
    ib = (Pc[j] < rng.random(N)[:, None]).sum(1).clip(0, 89)
    x = (rng.random(N) - 0.5) * 2; y = (rng.random(N) - 0.5) * 2; ph = rng.random(N) * 2 * np.pi; w = np.ones(N)
    top = back = edge = 0.0; alive = np.ones(N, bool)
    for _ in range(800):
        a = np.flatnonzero(alive)
        if a.size == 0: break
        for half in (0, 1):
            th = np.radians(TH[ib[a]] + rng.random(a.size) - 0.5); dd = 0.5 * np.tan(th)
            x[a] += dd * np.cos(ph[a]); y[a] += dd * np.sin(ph[a])
            out = (np.abs(x[a]) > 12) | (np.abs(y[a]) > 12); edge += w[a[out]].sum(); w[a[out]] = 0
            if half == 0:
                bt = BT[ib[a]]; top += (w[a] * bt).sum(); w[a] *= 1 - bt; u = rng.random(a.size)
                ib[a] = (BRc[:, ib[a]].T < u[:, None]).sum(1).clip(0, 89); ph[a] = rng.random(a.size) * 2 * np.pi
        inpix = (np.abs(x[a]) <= 1) & (np.abs(y[a]) <= 1)
        Rb = np.where(inpix, RS[j[a], ib[a]], RB[j[a], ib[a]])
        back += (w[a] * np.where(inpix, 0, 1 - RB[j[a], ib[a]])).sum()
        w[a] *= Rb; alive &= w > 1e-5
    return top / N, back / N, edge / N


if __name__ == '__main__':
    wp = (LAM * SPEC)[LAMJ]; wp = wp / wp.sum()
    rows = []
    for cat in ('Al', 'Ag'):
        with Pool(4) as p: res = p.map(device, [(cat, j) for j in LAMJ])
        ES = np.array([r[0] for r in res]); PS = np.array([r[1] for r in res]); RS = np.array([r[2] for r in res])
        es = float(wp @ ES); WL = wp * ES / (wp @ ES)
        ext_inf = float(sum(WL[k] * series.eta_ext(BT, BR, RS[k], PS[k], 400)[0] for k in range(len(LAMJ))))
        Pn = PS / PS.sum(1, keepdims=True); Ap = float(WL @ (1 - (Pn * RS).sum(1)))
        print('%s cathode: eta_sub %.3f | cathode-stack A\' %.3f | infinite pixel: eta_ext %.3f, EQE %.3f (PLQY 1: %.3f)'
              % (cat, es, Ap, ext_inf, PLQY * es * ext_inf, es * ext_inf), flush=True)
        rows.append(dict(cathode=cat, back='infinite pixel', eta_sub=es, eta_ext=ext_inf, EQE=PLQY * es * ext_inf, EQE_PLQY1=es * ext_inf,
                         back_loss=0, edge_loss=0))
        for kind in ('none', 'DBR', 'Ag'):
            RB = np.array([back_R(j, kind) for j in LAMJ])
            AB = float(WL @ (1 - (Pn * RB).sum(1)))
            with Pool(4) as p:
                out = np.mean(p.map(mc, [(100 * k + 1, 150000, ES, PS, RS, RB, WL) for k in range(8)]), 0)
            rows.append(dict(cathode=cat, back=kind, eta_sub=es, eta_ext=out[0], EQE=PLQY * es * out[0], EQE_PLQY1=es * out[0],
                             back_loss=out[1], edge_loss=out[2]))
            print('   back outside the pixel: %-4s (loss per hit %.3f) | 2x2 mm pixel on 24x24x0.5 mm: eta_ext %.3f  EQE %.3f (PLQY 1: %.3f) | back %.3f edge %.3f'
                  % (kind, AB, out[0], PLQY * es * out[0], es * out[0], out[1], out[2]), flush=True)
    with open('orange_reflectors.csv', 'w', newline='') as fh:
        wr = csv.DictWriter(fh, list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
