"""Orange preprint device with the electrode layout of the sketch (dimensions read from the sketch, approximate):
25 x 25 mm substrate (n 1.77, MLA on the whole front, thickness T, absorbing edges); IZO anode over the upper part
(y > 0.6 mm, origin at the substrate centre); five Al cathode stripes 3 mm wide at x = -9.4, -4.7, 0, 4.7, 9.4 mm running
from y = 3.6 down to y = -10.3 mm; lit pixel = centre stripe x IZO overlap (3 x 3 mm, y 0.6..3.6); organics everywhere;
auxiliary reflector (none / Al mirror behind an air gap / DBR 2.5 / DBR 4.5) over a centred 23 x 21 mm area.
Back regions: Al with IZO (pixel stack, cps2), Al without IZO, aux with / without IZO, bare with / without IZO."""
import numpy as np, sys, csv
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/orange']
import orange_reflectors as O
import orange_ag_cathode as A
import orange_reflector_series as S

HALF, DX = 12.5, 0.05
GX = np.arange(-HALF, HALF, DX) + DX / 2
XX, YY = np.meshgrid(GX, GX, indexing='ij')
IZO = YY > 0.6
AL = np.zeros_like(IZO)
for xc in (-9.4, -4.7, 0.0, 4.7, 9.4): AL |= (np.abs(XX - xc) <= 1.5) & (YY <= 3.6) & (YY >= -10.3)
AUX = (np.abs(XX) <= 11.5) & (np.abs(YY) <= 10.5)
# region index: 0 Al+IZO, 1 Al no IZO, 2 aux+IZO, 3 aux no IZO, 4 bare+IZO, 5 bare no IZO
REG = np.where(AL, np.where(IZO, 0, 1), np.where(AUX, np.where(IZO, 2, 3), np.where(IZO, 4, 5))).astype(np.int8)
PIX = (np.abs(XX) <= 1.5) & (YY > 0.6) & (YY <= 3.6)


def stack(j, izo, back):
    L = A.organics(j, 55)
    if not izo: L = L[1:]
    if back == 'Al': L = L + [(O.AL[j], 100.0)]
    elif back.startswith('DBR'):
        L = L + [(O.ZNS[j], 70.0) if k % 2 == 0 else (O.LIF[j], 115.0) for k in range({'DBR2.5': 5, 'DBR4.5': 9}[back])]
    return L


def back_R(j, izo, kind):
    """kind: 'cath' (Al cathode), 'none', 'AlMirror', 'DBR2.5', 'DBR4.5'"""
    lam = float(O.LAM[j]); out = 0
    L = stack(j, izo, 'Al' if kind == 'cath' else (kind if kind.startswith('DBR') else 'none'))
    for pol in ('s', 'p'):
        R, T = S.rt_pol(L, lam, pol)
        if kind == 'AlMirror':
            Ra = S.r_al_air(lam, j, pol); R = R + T * T * Ra / (1 - R * Ra)
        out = out + 0.5 * R
    return np.clip(out, 0, 1)


def mc(args):
    seed, N, T, PS, RTAB, WL = args          # RTAB[region, lambda, theta]
    rng = np.random.default_rng(seed)
    j = rng.choice(len(WL), N, p=WL); Pc = np.cumsum(PS / PS.sum(1, keepdims=True), 1)
    ib = (Pc[j] < rng.random(N)[:, None]).sum(1).clip(0, 89)
    x = (rng.random(N) - 0.5) * 3.0; y = 0.6 + rng.random(N) * 3.0
    ph = rng.random(N) * 2 * np.pi; w = np.ones(N); top = 0.0; alive = np.ones(N, bool)
    for _ in range(800):
        a = np.flatnonzero(alive)
        if a.size == 0: break
        for h in (0, 1):
            th = np.radians(O.TH[ib[a]] + rng.random(a.size) - 0.5); dd = T * np.tan(th)
            x[a] += dd * np.cos(ph[a]); y[a] += dd * np.sin(ph[a])
            out = (np.abs(x[a]) > HALF) | (np.abs(y[a]) > HALF); w[a[out]] = 0
            if h == 0:
                bt = O.BT[ib[a]]; top += (w[a] * bt).sum(); w[a] *= 1 - bt; u = rng.random(a.size)
                ib[a] = (O.BRc[:, ib[a]].T < u[:, None]).sum(1).clip(0, 89); ph[a] = rng.random(a.size) * 2 * np.pi
        ix = np.clip(((x[a] + HALF) / DX).astype(int), 0, len(GX) - 1); iy = np.clip(((y[a] + HALF) / DX).astype(int), 0, len(GX) - 1)
        w[a] *= RTAB[REG[ix, iy], j[a], ib[a]]; alive &= w > 1e-5
    return top / N


if __name__ == '__main__':
    T = float(sys.argv[1]) if len(sys.argv) > 1 else 0.7
    wp = (O.LAM * O.SPEC)[O.LAMJ]; wp = wp / wp.sum()
    with Pool(4) as p: res = p.map(A.device, [('Al', 55, j) for j in O.LAMJ])
    ES = np.array([r[0] for r in res]); PS = np.array([r[1] for r in res]); RS = np.array([r[2] for r in res])
    es = float(wp @ ES); WL = wp * ES / (wp @ ES); Pn = PS / PS.sum(1, keepdims=True)
    frac = {k: float(np.mean(REG == i)) for i, k in enumerate(('Al+IZO', 'Al', 'aux+IZO', 'aux', 'bare+IZO', 'bare'))}
    print('area fractions', {k: round(v, 3) for k, v in frac.items()})
    Rcath0 = RS; Rcath1 = np.array([back_R(j, False, 'cath') for j in O.LAMJ])
    Rn = {iz: np.array([back_R(j, iz, 'none') for j in O.LAMJ]) for iz in (True, False)}
    rows = []
    for kind in S.MEAS:
        kk = {'none': 'none', 'Al': 'AlMirror', 'DBR2.5': 'DBR2.5', 'DBR4.5': 'DBR4.5'}[kind]
        Ra = {iz: np.array([back_R(j, iz, kk) for j in O.LAMJ]) for iz in (True, False)}
        RTAB = np.stack([Rcath0, Rcath1, Ra[True], Ra[False], Rn[True], Rn[False]])
        Aa = {iz: float(WL @ (1 - (Pn * Ra[iz]).sum(1))) for iz in (True, False)}
        with Pool(4) as p:
            e = float(np.mean(p.map(mc, [(10 * k + 1, 120000, T, PS, RTAB, WL) for k in range(8)])))
        rows.append(dict(aux=kind, T_mm=T, A_aux_IZO=Aa[True], A_aux_noIZO=Aa[False], eta_ext=e,
                         EQE_model=100 * O.PLQY * es * e, EQE_meas=S.MEAS[kind]))
        print('%-6s A\'aux IZO %.3f / no IZO %.3f | model EQE %.1f%% | measured %.1f%%'
              % (kind, Aa[True], Aa[False], 100 * O.PLQY * es * e, S.MEAS[kind]), flush=True)
    Ac = {'Al+IZO': float(WL @ (1 - (Pn * Rcath0).sum(1))), 'Al': float(WL @ (1 - (Pn * Rcath1).sum(1)))}
    print('cathode regions A\'', {k: round(v, 3) for k, v in Ac.items()})
    with open('/home/user/OLED-/sim/orange/orange_layout_T%.1f.csv' % T, 'w', newline='') as fh:
        wr = csv.DictWriter(fh, list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
