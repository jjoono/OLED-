"""Goniometer measurement of the cone-array device (24 x 24 x 0.5 mm substrate, 2 x 2 mm pixel, absorbing edges):
 1. far-field anisotropy of the hexagonal cone array (ray tracer, isotropic substrate radiance in, exit direction out):
    I(theta) along the nearest-neighbour azimuth (0 deg) and the 30 deg azimuth vs the azimuthal mean
 2. where on the front face the light leaves (Monte Carlo of finite_substrate.py, scenarios A and B)
 3. a goniometer with a 75 mm^2 detector at L = 50 / 80 mm rotating the sample about the pixel centre: measured I(theta)
    and the hemisphere-integrated EQE relative to the true value (also with the rotation axis 1 mm off the pixel centre)."""
import numpy as np, sys, csv
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/be_cone', '/home/user/OLED-/sim/mla']
import raytrace as rt
import finite_substrate as F

R_UM, H_UM, GAP = 1.4, 1.8, 0.2
TB = np.linspace(0, 90, 46); TC = 0.5 * (TB[1:] + TB[:-1])          # 2-deg air-angle bins


def exitdirs(seed, N=400000):
    rng = np.random.default_rng(seed)
    esc, _, Dd = rt._trace_chunk(H_UM / R_UM, 1.77, N, 80, rng, None, 'cone', (2 * R_UM + GAP) / R_UM, return_dir=True)
    d = Dd[esc]; th = np.degrees(np.arccos(np.clip(d[:, 2], -1, 1))); ph = np.degrees(np.arctan2(d[:, 1], d[:, 0])) % 60
    ph = np.where(ph > 30, 60 - ph, ph)                                 # fold onto 0..30 (C6v)
    return th, ph


def farfield():
    with Pool(4) as p: out = p.map(exitdirs, range(10))
    th = np.concatenate([o[0] for o in out]); ph = np.concatenate([o[1] for o in out])
    sa = 2 * np.pi * (np.cos(np.radians(TB[:-1])) - np.cos(np.radians(TB[1:])))          # solid angle per ring
    def prof(sel, frac):                                                                    # radiant intensity per sr
        h, _ = np.histogram(th[sel], TB); return h / (sa * frac)
    Iavg = prof(np.ones(th.size, bool), 1.0)
    I0 = prof(ph < 3, 3 / 30); I30 = prof(ph > 27, 3 / 30)
    tot = (Iavg * sa).sum()
    return Iavg / tot, I0 / tot, I30 / tot, sa


G = np.arange(-12, 12.001, 0.25)


def _spatial_run(args):
    scen, seed = args
    rng = np.random.default_rng(seed); N = 150000
    j = rng.choice(len(WL), N, p=WL); Pc = np.cumsum(PS / PS.sum(1, keepdims=True), 1)
    ib = (Pc[j] < rng.random(N)[:, None]).sum(1).clip(0, 89)
    x = (rng.random(N) - 0.5) * 2 * F.PIX; y = (rng.random(N) - 0.5) * 2 * F.PIX; ph = rng.random(N) * 2 * np.pi; w = np.ones(N)
    H = np.zeros((len(G) - 1, len(G) - 1)); alive = np.ones(N, bool)
    for _ in range(600):
        a = np.flatnonzero(alive)
        if a.size == 0: break
        for half in (0, 1):
            th = np.radians(F.TH[ib[a]] + rng.random(a.size) - 0.5); dd = F.T * np.tan(th)
            x[a] += dd * np.cos(ph[a]); y[a] += dd * np.sin(ph[a])
            out = (np.abs(x[a]) > F.HALF) | (np.abs(y[a]) > F.HALF); w[a[out]] = 0
            if half == 0:
                bt = F.BT[ib[a]]; h2, _, _ = np.histogram2d(x[a], y[a], [G, G], weights=w[a] * bt); H += h2
                w[a] *= 1 - bt; u = rng.random(a.size)
                ib[a] = (F.BRc[:, ib[a]].T < u[:, None]).sum(1).clip(0, 89); ph[a] = rng.random(a.size) * 2 * np.pi
        inpix = (np.abs(x[a]) <= F.PIX) & (np.abs(y[a]) <= F.PIX); Rd = RS[j[a], ib[a]]
        w[a] *= np.where(inpix, Rd, F.RF[ib[a]]) if scen == 'A' else Rd
        alive &= w > 1e-5
    return H / N


def spatial(scen, n=8):
    """front-face exit positions (weights) on a 0.25 mm grid"""
    with Pool(4) as p: Hs = p.map(_spatial_run, [(scen, 1000 * k + 3) for k in range(n)])
    return G, np.mean(Hs, 0)


def gonio(g, H, Iavg, L, offset=0.0, rpd=np.sqrt(75 / np.pi)):
    xc = 0.5 * (g[1:] + g[:-1]); X, Y = np.meshgrid(xc, xc, indexing='ij'); m = H > H.max() * 1e-4
    sx, sy, sw = X[m] - offset, Y[m], H[m] / H[m].sum()
    r = rpd * np.sqrt((np.arange(60) + 0.5) / 60); a = np.arange(60) * 2.399963          # sunflower points on the PD
    pu, pv = r * np.cos(a), r * np.sin(a); dA = np.pi * rpd ** 2 / 60
    Ith = lambda t: np.interp(t, TC, Iavg / np.interp(0, TC, Iavg))                        # relative I(theta) per point
    I0 = Iavg / ((Iavg * SA).sum())
    tg = np.arange(1, 89, 2.0); S = []
    for t in tg:
        tr = np.radians(t); C = L * np.array([np.sin(tr), 0, np.cos(tr)])
        e1 = np.array([np.cos(tr), 0, -np.sin(tr)]); e2 = np.array([0, 1, 0]); nrm = -C / L
        P = C[None, :] + pu[:, None] * e1[None, :] + pv[:, None] * e2[None, :]               # (60, 3)
        v = P[None, :, :] - np.stack([sx, sy, np.zeros_like(sx)], 1)[:, None, :]
        d = np.linalg.norm(v, axis=2); th = np.degrees(np.arccos(np.clip(v[..., 2] / d, -1, 1)))
        cosd = np.clip(-(v @ nrm) / d, 0, None)
        S.append((sw[:, None] * np.interp(th, TC, I0) * cosd / d ** 2 * dA).sum())
    S = np.array(S); Imeas = S / (np.pi * rpd ** 2 / L ** 2)
    est = (Imeas * 2 * np.pi * np.sin(np.radians(tg)) * np.radians(2.0)).sum()
    return tg, Imeas, est


if __name__ == '__main__':
    with Pool(4) as p: res = p.map(F.stack, F.LAMJ)
    ES = np.array([r[0] for r in res]); PS = np.array([r[1] for r in res]); RS = np.array([r[2] for r in res])
    I = np.clip(np.asarray(F.D.spec.l_I_Irppy2acac, float)[:301], 0, None); wp = (F.D.LAM * I)[F.LAMJ]; wp /= wp.sum()
    WL = wp * ES / (wp @ ES)
    Iavg, I0, I30, SA = farfield()
    tot = lambda Ip: (Ip * SA).sum()
    print('far field (isotropic substrate radiance): I(theta)/I(0) avg:', ' '.join('%.2f' % (Iavg[k] / Iavg[0]) for k in (0, 7, 14, 22, 29, 37)), '(theta 1,15,29,45,59,75)')
    print('Lambertian cos               :', ' '.join('%.2f' % np.cos(np.radians(TC[k])) for k in (0, 7, 14, 22, 29, 37)))
    for lab, Ip in (('0 deg (nearest neighbour)', I0), ('30 deg', I30)):
        print('azimuth %-26s I/I_avg at theta 1/15/29/45/59/75: %s | hemisphere total if this line is assumed for all azimuths: %.3f'
              % (lab, ' '.join('%.2f' % (Ip[k] / Iavg[k]) for k in (0, 7, 14, 22, 29, 37)), tot(Ip)))
    rows = []
    for scen in ('A', 'B'):
        g, H = spatial(scen)
        xc = 0.5 * (g[1:] + g[:-1]); X, Y = np.meshgrid(xc, xc, indexing='ij'); rr = np.hypot(X, Y); Hn = H / H.sum()
        rs = np.sort(rr.ravel()); cum = np.cumsum(Hn.ravel()[np.argsort(rr.ravel())])
        r50, r90, r99 = (rs[np.searchsorted(cum, q)] for q in (0.5, 0.9, 0.99))
        inpix = Hn[(np.abs(X) <= 1) & (np.abs(Y) <= 1)].sum()
        print('scenario %s: front emission inside the 2x2 mm pixel %.2f | radius containing 50/90/99%%: %.1f / %.1f / %.1f mm' % (scen, inpix, r50, r90, r99))
        for L in (50.0, 80.0):
            pt = gonio(np.array([-0.125, 0.125]), np.ones((1, 1)), Iavg, L)[2]               # point source = quadrature reference
            for off in (0.0, 0.5, 1.0):
                tp, Ip, ep = gonio(g, H, Iavg, L, off); tm, Imn, em = gonio(g, H, Iavg, L, -off)   # +theta and -theta scans
                ref = np.interp(tp, TC, Iavg / tot(Iavg)) * (Ip[0] / np.interp(tp[0], TC, Iavg / tot(Iavg)))
                k45 = np.argmin(abs(tp - 45)); k75 = np.argmin(abs(tp - 75))
                print('   L %2.0f mm, axis offset %.1f mm: one-side scan %.4f | both-side average %.4f | point-source ref %.4f -> '
                      'extended/point one-side %.4f, averaged %.4f'
                      % (L, off, ep, 0.5 * (ep + em), pt, ep / pt, 0.5 * (ep + em) / pt))
                rows.append(dict(scenario=scen, L_mm=L, axis_offset_mm=off, one_side=ep / pt, both_sides=0.5 * (ep + em) / pt, r50=r50, r90=r90, in_pixel=inpix))
    with open('gonio_check.csv', 'w', newline='') as fh:
        wr = csv.DictWriter(fh, list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
