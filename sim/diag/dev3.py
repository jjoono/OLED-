"""User's green BE device, three geometries, finite substrate (24 x 24 x 0.5 mm, 2 x 2 mm pixel at (-2.48, 1.5) mm,
absorbing side walls):
  stack: substrate | IZO 50 | HATCN 5 | TAPC 70 | HATCN 5 | TAPC 70 | HATCN 5 | TAPC 70 | TCTA 10 | TCTA:B3PyMPM:Ir(ppy)2acac 25
         | B3PyMPM 120*f | B3PyMPM:Cs2CO3 10 wt% 60*f (k_d) | B3PyMPM:Cs2CO3 50 wt% 5 (n 1.7, k_cs) | Ag 100
  (1) reference: flat glass n 1.52, no MLA (front-side emission; matrix series with the flat-interface BSDF)
  (2) MLA substrate n 1.77 (cone BSDF), back reflection only where the cathode pattern is (no auxiliary reflector)
  (3) as (2) with a full back reflector (OLED-stack reflectance everywhere)
Theta 0.76, PLQY 1, Ir(ppy)2acac."""
import numpy as np, sys
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/be_cone', '/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/mla']
import device_eqe as D, cps2, series
import cathode_as_reflector as CR

TH = D.TH; IZO = D.lib('l_IZO_SJ_RTP'); BTc, BRc_ = D.EXT['cone']
BRn = BRc_ / np.maximum(BRc_.sum(0), 1e-12); BRC = np.cumsum(BRn, 0)
LAMJ = list(range(0, 301, 10)); THETA = 0.76
I = np.clip(np.asarray(D.spec.l_I_Irppy2acac, float)[:301], 0, None); WP = (D.LAM * I)[LAMJ]; WP = WP / WP.sum()
def flat(ns):
    c = np.cos(np.radians(TH)); s2 = (ns * np.sin(np.radians(TH))) ** 2; ct = np.sqrt(np.clip(1 - s2, 0, 1))
    Rf = np.where(s2 < 1, 0.5 * (((ns * c - ct) / (ns * c + ct)) ** 2 + ((ns * ct - c) / (ns * ct + c)) ** 2), 1.0)
    return 1 - Rf, np.diag(Rf)


def stack(args):
    j, ns, f, kd, kcs, theta = args
    B3o, B3e, EMo, EMe, TCo, TCe, TPo, TPe, HAT, AG = (D.B3o[j], D.B3e[j], D.EMo[j], D.EMe[j], D.TCo[j], D.TCe[j], D.TPo[j], D.TPe[j], D.HAT[j], D.AG[j])
    cs = 1.7 + 1j * kcs
    S = cps2.Stack(float(D.LAM[j]), (EMo, EMe), 25.0, 12.5,
                   above=[(B3o, B3e, 120.0 * f), (B3o + 1j * kd, B3e + 1j * kd, 60.0 * f), (cs, cs, 5.0), (AG, AG, 100.0)],
                   below=[(TCo, TCe, 10.0), (TPo, TPe, 70.0), (HAT, HAT, 5.0), (TPo, TPe, 70.0), (HAT, HAT, 5.0), (TPo, TPe, 70.0),
                          (HAT, HAT, 5.0), (IZO[j], IZO[j], 50.0)], n_sub=ns, h=theta)
    r = cps2.solve_pol(S, npts=3000)
    return float(np.real(r['air'] + r['sub'])), np.real(cps2.sub_angular(S, TH)), np.real(cps2.stack_reflectance(S, TH))


def mc(args):
    seed, N, PS, RS, WL, mode = args
    rng = np.random.default_rng(seed)
    j = rng.choice(len(WL), N, p=WL); Pc = np.cumsum(PS / PS.sum(1, keepdims=True), 1)
    ib = (Pc[j] < rng.random(N)[:, None]).sum(1).clip(0, 89)
    x = -2.48 + (rng.random(N) - 0.5) * 2.0; y = 1.5 + (rng.random(N) - 0.5) * 1.96
    ph = rng.random(N) * 2 * np.pi; w = np.ones(N); top = 0.0; alive = np.ones(N, bool)
    RF = flat(1.77)[1].diagonal()
    for _ in range(800):
        a = np.flatnonzero(alive)
        if a.size == 0: break
        for h in (0, 1):
            th = np.radians(TH[ib[a]] + rng.random(a.size) - 0.5); dd = 0.5 * np.tan(th)
            x[a] += dd * np.cos(ph[a]); y[a] += dd * np.sin(ph[a])
            out = (np.abs(x[a]) > 12) | (np.abs(y[a]) > 12); w[a[out]] = 0
            if h == 0:
                bt = BTc[ib[a]]; top += (w[a] * bt).sum(); w[a] *= 1 - bt; u = rng.random(a.size)
                ib[a] = (BRC[:, ib[a]].T < u[:, None]).sum(1).clip(0, 89); ph[a] = rng.random(a.size) * 2 * np.pi
        Rd = RS[j[a], ib[a]]
        if mode == 'cathode':
            ix = np.clip(((x[a] + 12) / CR.DX).astype(int), 0, len(CR.GX) - 1); iy = np.clip(((y[a] + 12) / CR.DX).astype(int), 0, len(CR.GX) - 1)
            Rd = np.where(CR.MET[ix, iy], Rd, RF[ib[a]])
        w[a] *= Rd; alive &= w > 1e-5
    return top / N


def three(pool, f=1.0, kd=0.0, kcs=0.0, theta=THETA, nref=1.52, N=60000, nrun=4):
    out = {}
    res = pool.map(stack, [(j, nref, f, kd, kcs, theta) for j in LAMJ])
    ES = np.array([r[0] for r in res]); PS = np.array([r[1] for r in res]); RS = np.array([r[2] for r in res])
    es = WP @ ES; WL = WP * ES / es; BT, BR = flat(nref)
    out['ref'] = float(sum(WP[k] * ES[k] * series.eta_ext(BT, BR, RS[k], PS[k], 300)[0] for k in range(len(LAMJ))))
    out['es_ref'] = float(es)
    res = pool.map(stack, [(j, 1.77, f, kd, kcs, theta) for j in LAMJ])
    ES = np.array([r[0] for r in res]); PS = np.array([r[1] for r in res]); RS = np.array([r[2] for r in res])
    es = WP @ ES; WL = WP * ES / es; out['es'] = float(es)
    for mode in ('cathode', 'full'):
        e = np.mean(pool.map(mc, [(17 * k + 1, N, PS, RS, WL, mode) for k in range(nrun)]))
        out[mode] = float(es * e)
    return out


if __name__ == '__main__':
    meas = np.array([0.16, 0.51, 0.59])
    with Pool(4) as p:
        o = three(p, N=150000, nrun=8)
        m = np.array([o['ref'], o['cathode'], o['full']])
        print('nominal: eta_sub(glass) %.3f eta_sub(1.77) %.3f | EQE ref %.3f  MLA %.3f  MLA+reflector %.3f' % (o['es_ref'], o['es'], *m))
        print('measured / model:', np.round(meas / m, 3), '| ratios measured MLA/ref %.2f, refl/MLA %.2f ; model %.2f, %.2f'
              % (meas[1] / meas[0], meas[2] / meas[1], m[1] / m[0], m[2] / m[1]))
        o2 = three(p, nref=1.77, N=150000, nrun=8)
        print('(reference on n 1.77 flat instead of n 1.52 glass: %.3f)' % o2['ref'])
