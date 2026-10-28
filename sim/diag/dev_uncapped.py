"""MLA devices (n 1.77 cone MLA, finite 24 x 24 x 0.5 mm substrate, centre-left pixel) with a realistic back-face map:
  cathode metal                      -> OLED-stack reflectance (dev3.stack, fitted stack, no extra absorption in the pixel)
  organic area not under the cathode -> glass | organics (+ absorber) | TFE 30 | parylene-C (incoherent) | back
  outside the organic area           -> glass | TFE 30 | parylene-C (incoherent) | back
  back = air (no auxiliary reflector; TIR above the parylene/air critical angle)  or  MoOx 5 (k_MoOx) | Ag 100 (n + 0.099)
Absorber scenarios: none / HATCN k in the uncapped organics / MoOx k / both.   Measured: MLA 0.51, MLA + reflector 0.579.
Organic area (measured mask): pattern x 6.10-14.80, y 5.15-17.85 mm (21 x 23 pattern centred on the substrate)."""
import sys, numpy as np
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/diag', '/home/user/OLED-/sim/be_cone', '/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/mla']
import dev3, cps2
CR = dev3.CR; D = dev3.D; TH = dev3.TH; LAMJ = dev3.LAMJ; WP = dev3.WP; BTc = dev3.BTc; BRC = dev3.BRC
NS = 1.77; F_, G_, THETA, DIZO = 172.5 / 185, 160.0 / 210, 0.714, 12.7
ORG = (CR.XX >= 6.10 - 10.5) & (CR.XX <= 14.80 - 10.5) & (CR.YY >= 5.15 - 11.5) & (CR.YY <= 17.85 - 11.5)
def tmm(ns, ds, n0, nexit, lam, kx, pol):
    k0 = 2 * np.pi / lam
    def q(n):
        v = np.sqrt(complex(n) ** 2 * k0 ** 2 - kx ** 2 + 0j); return v if v.imag > 0 or (v.imag == 0 and v.real >= 0) else -v
    eta = (lambda n: q(n) / k0) if pol == 's' else (lambda n: complex(n) ** 2 * k0 / q(n))
    M = np.eye(2, dtype=complex)
    for n, d in zip(ns, ds):
        de = q(n) * d; e = eta(n); M = M @ np.array([[np.cos(de), -1j * np.sin(de) / e], [-1j * e * np.sin(de), np.cos(de)]])
    e0, es = eta(n0), eta(nexit); B, C = M @ np.array([1, es]); den = e0 * B + C
    return abs((e0 * B - C) / den) ** 2, 4 * e0.real * es.real / abs(den) ** 2
def back_map(args):
    """R(theta) seen from the substrate for the uncapped stacks, s/p averaged; parylene treated incoherently"""
    j, kh, kmo, refl, org = args
    lam = float(D.LAM[j]); npar = 1.62 + 0.0065 / (lam / 1000) ** 2; hat = D.HAT[j] + 1j * kh
    ns = ([hat, D.TPo[j], hat, D.TPo[j], hat, D.TPo[j], D.TCo[j], D.EMo[j], D.B3o[j], D.B3o[j], 1.7] if org else []) + [1.75]
    ds = ([5, 70 * G_, 5, 70 * G_, 5, 70 * G_, 10, 25, 120 * F_, 60 * F_, 5] if org else []) + [30]
    out = []
    for th in TH:
        kx = 2 * np.pi / lam * NS * np.sin(np.radians(th)); Rt = 0
        for pol in 's', 'p':
            R1, T1 = tmm(ns, ds, NS, npar, lam, kx, pol); R1b, _ = tmm(ns[::-1], ds[::-1], npar, NS, lam, kx, pol)
            if refl: R2, _ = tmm([2.1 + 1j * kmo, D.AG[j] + 0.099], [5, 100], npar, 1.0, lam, kx, pol)
            else: R2, _ = tmm([], [], npar, 1.0, lam, kx, pol)
            Rt += 0.5 * (R1 + T1 * T1 * R2 / (1 - R1b * R2))
        out.append(Rt)
    return np.array(out)
def mc(args):
    seed, N, PS, RS, RO, RN, WL = args
    rng = np.random.default_rng(seed)
    j = rng.choice(len(WL), N, p=WL); Pc = np.cumsum(PS / PS.sum(1, keepdims=True), 1)
    ib = (Pc[j] < rng.random(N)[:, None]).sum(1).clip(0, 89)
    x = -2.48 + (rng.random(N) - 0.5) * 2.0; y = 1.5 + (rng.random(N) - 0.5) * 1.96
    ph = rng.random(N) * 2 * np.pi; w = np.ones(N); top = 0.0; alive = np.ones(N, bool)
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
        ix = np.clip(((x[a] + 12) / CR.DX).astype(int), 0, len(CR.GX) - 1); iy = np.clip(((y[a] + 12) / CR.DX).astype(int), 0, len(CR.GX) - 1)
        Rd = np.where(CR.MET[ix, iy], RS[j[a], ib[a]], np.where(ORG[ix, iy], RO[j[a], ib[a]], RN[j[a], ib[a]]))
        w[a] *= Rd; alive &= w > 1e-5
    return top / N
if __name__ == '__main__':
    with Pool(4) as p:
        res = p.map(dev3.stack, [(j, NS, F_, 0.0, 0.0, THETA, G_, 1.0, DIZO, 1.0) for j in LAMJ])
        ES = np.array([r[0] for r in res]); PS = np.array([r[1] for r in res]); RS = np.array([r[2] for r in res])
        es = WP @ ES; WL = WP * ES / es
        print('eta_sub(1.77) %.3f' % es, flush=True)
        for name, kh, kmo in (('none', 0, 0), ('HATCN k 0.1', 0.1, 0), ('HATCN k 0.2', 0.2, 0), ('MoOx k 1.0', 0, 1.0), ('HATCN 0.2 + MoOx 1.0', 0.2, 1.0)):
            row = []
            for refl in (False, True):
                RO = np.array(p.map(back_map, [(j, kh, kmo, refl, True) for j in LAMJ]))
                RN = np.array(p.map(back_map, [(j, 0.0, kmo, refl, False) for j in LAMJ]))
                e = np.mean(p.map(mc, [(17 * k + 1, 100000, PS, RS, RO, RN, WL) for k in range(8)]))
                row.append(es * e)
            print('%-22s MLA %.3f  MLA+reflector %.3f   (measured 0.51 / 0.579)' % (name, *row), flush=True)
