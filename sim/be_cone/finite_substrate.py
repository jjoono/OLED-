"""Finite substrate: 24 x 24 mm, 0.5 mm thick, n 1.77, cone array on the whole top face, 2 x 2 mm pixel at the centre of the
bottom face, side walls perfectly absorbing.  Monte-Carlo ray tracing in the slab (weights), spectrally resolved:
  source  : uniform over the pixel, lambda and theta sampled from the substrate-delivered power P_sub(lambda, theta)
            (Theta 0.76, PLQY 1, l_I_Irppy2acac), azimuth uniform
  top     : cone-array BSDF (escape B_T(theta); return polar angle from B_R, azimuth re-randomised)
  bottom  : inside the pixel -> OLED stack reflectance R(lambda, theta), specular
            outside the pixel -> scenario A: bare glass/air (Fresnel, TIR above 34.4 deg; transmitted light = back loss)
                                 scenario B: same reflectance as the OLED stack (wide reflector, upper bound)
  walls   : |x| or |y| > 12 mm -> absorbed
Check: scenario B with walls removed reproduces the matrix-series eta_ext."""
import numpy as np, sys, csv
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/be_cone', '/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/mla']
import device_eqe as D, cps2, series

THETA, NS, T, HALF, PIX = 0.76, 1.77, 0.5, 12.0, 1.0          # mm
LAMJ = list(range(0, 301, 10))
BT, BR = D.EXT['cone']; TH = D.TH
BRn = BR / np.maximum(BR.sum(0), 1e-12); BRc = np.cumsum(BRn, 0)
c = np.cos(np.radians(TH)); s2 = (NS * np.sin(np.radians(TH))) ** 2; ct = np.sqrt(np.clip(1 - s2, 0, 1))
RF = np.where(s2 < 1, 0.5 * (((NS * c - ct) / (NS * c + ct)) ** 2 + ((NS * ct - c) / (NS * ct + c)) ** 2), 1.0)


def stack(j):
    S = cps2.Stack(float(D.LAM[j]), (D.EMo[j], D.EMe[j]), 25.0, 12.5,
                   above=[(D.B3o[j], D.B3e[j], 10.0), (D.B3o[j], D.B3e[j], 170.0), (D.AG[j], D.AG[j], 100.0)],
                   below=[(D.TCo[j], D.TCe[j], 10.0), (D.TPo[j], D.TPe[j], 200.0), (D.HAT[j], D.HAT[j], 5.0), (D.ITO[j], D.ITO[j], 50.0)],
                   n_sub=NS, h=THETA)
    r = cps2.solve_pol(S, npts=3000)
    return float(np.real(r['air'] + r['sub'])), np.real(cps2.sub_angular(S, TH)), np.real(cps2.stack_reflectance(S, TH))


def trace(args):
    seed, N, scen, walls, ES, PS, RS, wl = args
    rng = np.random.default_rng(seed)
    j = rng.choice(len(wl), N, p=wl)
    Pc = np.cumsum(PS / PS.sum(1, keepdims=True), 1)
    ib = (Pc[j] < rng.random(N)[:, None]).sum(1).clip(0, 89)
    x = (rng.random(N) - 0.5) * 2 * PIX; y = (rng.random(N) - 0.5) * 2 * PIX
    ph = rng.random(N) * 2 * np.pi; w = np.ones(N)
    top = np.zeros(N); back = np.zeros(N); edge = np.zeros(N); absb = np.zeros(N)
    alive = np.ones(N, bool)
    for _ in range(600):
        a = np.flatnonzero(alive)
        if a.size == 0: break
        th = np.radians(TH[ib[a]] + rng.random(a.size) - 0.5)
        d = T * np.tan(th)
        x[a] += d * np.cos(ph[a]); y[a] += d * np.sin(ph[a])                       # bottom -> top
        if walls:
            out = (np.abs(x[a]) > HALF) | (np.abs(y[a]) > HALF)
            edge[a[out]] += w[a[out]]; w[a[out]] = 0
        bt = BT[ib[a]]; top[a] += w[a] * bt; w[a] *= 1 - bt                          # top: escape / return
        u = rng.random(a.size)
        ib[a] = (BRc[:, ib[a]].T < u[:, None]).sum(1).clip(0, 89)
        ph[a] = rng.random(a.size) * 2 * np.pi
        th = np.radians(TH[ib[a]] + rng.random(a.size) - 0.5)
        d = T * np.tan(th)
        x[a] += d * np.cos(ph[a]); y[a] += d * np.sin(ph[a])                       # top -> bottom
        if walls:
            out = (np.abs(x[a]) > HALF) | (np.abs(y[a]) > HALF)
            edge[a[out]] += w[a[out]]; w[a[out]] = 0
        inpix = (np.abs(x[a]) <= PIX) & (np.abs(y[a]) <= PIX)
        Rd = RS[j[a], ib[a]]
        if scen == 'A':
            Rb = np.where(inpix, Rd, RF[ib[a]])
            back[a] += np.where(inpix, 0, w[a] * (1 - RF[ib[a]]))
            absb[a] += np.where(inpix, w[a] * (1 - Rd), 0)
        else:
            Rb = Rd; absb[a] += w[a] * (1 - Rd)
        w[a] *= Rb                                                                    # specular: theta, phi kept
        alive &= w > 1e-5
    return top.sum(), back.sum(), edge.sum(), absb.sum(), w.sum(), N


if __name__ == '__main__':
    with Pool(4) as p:
        res = p.map(stack, LAMJ)
    ES = np.array([r[0] for r in res]); PS = np.array([r[1] for r in res]); RS = np.array([r[2] for r in res])
    I = np.clip(np.asarray(D.spec.l_I_Irppy2acac, float)[:301], 0, None)
    wp = (D.LAM * I)[LAMJ]; wp = wp / wp.sum()
    es = float(wp @ ES); wl = wp * ES / (wp @ ES)                                     # spectrum of substrate-delivered photons
    eta_series = float(sum(wl[k] * series.eta_ext(BT, BR, RS[k], PS[k], 400)[0] for k in range(len(LAMJ))))
    print('eta_sub %.4f | infinite-substrate series eta_ext %.4f -> EQE %.4f' % (es, eta_series, es * eta_series))
    rows = []
    for scen, walls, lab in (('B', False, 'check: wide reflector, no walls (= series)'), ('B', True, 'B: wide reflector (= OLED R everywhere), walls absorbing'),
                             ('A', True, 'A: Ag only in the 2x2 mm pixel, bare glass elsewhere, walls absorbing'), ('A', False, 'A without walls')):
        with Pool(4) as p:
            out = np.array(p.map(trace, [(100 * k + 7, 250000, scen, walls, ES, PS, RS, wl) for k in range(8)]))
        tot = out.sum(0); n = tot[5]
        f = dict(top=tot[0] / n, back=tot[1] / n, edge=tot[2] / n, absorbed=tot[3] / n, residual=tot[4] / n)
        rows.append(dict(case=lab, eta_ext=f['top'], EQE=es * f['top'], back_loss=f['back'], edge_loss=f['edge'], absorbed=f['absorbed'], residual=f['residual']))
        print('%-72s eta_ext %.4f  EQE %.4f | back %.4f edge %.4f absorbed %.4f residual %.5f' % (lab, f['top'], es * f['top'], f['back'], f['edge'], f['absorbed'], f['residual']), flush=True)
    with open('finite_substrate.csv', 'w', newline='') as fh:
        wr = csv.DictWriter(fh, list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
