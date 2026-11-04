"""Forward-directed NIR OLED: multi-pass extraction with an angle-selective top DBR on a scattering slab.

Geometry (bottom emission):  device stack | substrate n_sub | scattering slab (index-matched, HG g, S) | DBR (TiO2/SiO2,
free thicknesses) | air.  Light returned by the DBR is re-scattered by the slab, goes back to the device (stack
reflectance R(theta)) and is recycled.  Composite outcoupler (slab + DBR), all in 1-deg polar bins of the slab medium:
   u = (I - R1 Rd)^-1 T1 v,   out = Td u,   back = R1 v + T1 Rd u
Recycling series as in the paper: U = (I - BR_eff R_LED)^-1 v0, air output = BT_eff U.
Air distribution D[a] (1-deg air bins, photons per electron) with PLQY q and Purcell factor F:
   D = q sum_l w_l F_l es_l D_l / (1 - q + q sum_l w_l F_l).
Metrics: EQE = sum D;  C(t) = photons per electron inside the air cone of half angle t;  Leq(t) = C(t)/sin^2 t
(= EQE of a Lambertian source with the same mean radiance inside t);  Erel(z, a) = on-axis irradiance at distance z
from a uniform disc source of radius a, in photons per electron (contact with a large source -> EQE)."""
import numpy as np, glob, os, sys
sys.path.insert(0, '/home/user/OLED-/sim/nir')
import nir_opt as N
TH = np.arange(90) + 0.5
SLABS = sorted(glob.glob('/home/user/OLED-/sim/nir/collim/slabs/*.npz'))
def load_slab(f): d = np.load(f); return d['T'], d['R'], float(d['S']), float(d['g'])
def load_dev(f):
    d = np.load(f); return {k: d[k] for k in d.files}
def spectrum(lams, lam0, fwhm):
    w = np.exp(-0.5 * ((lams - lam0) / (fwhm / 2.3548)) ** 2); return w / w.sum()
# ---------------------------------------------------------------- thin-film R/T (unpolarised), incidence from medium n_in at polar angles th (deg)
def tmm_RT(n_in, layers, lam, th, n_out=1.0):
    s = n_in * np.sin(np.radians(th)) + 0j
    def cosang(n): return np.sqrt(1 - (s / n) ** 2 + 0j)
    ci, co = cosang(n_in), cosang(n_out); co = np.where(co.imag < 0, -co, co)
    out = []
    for pol in ('s', 'p'):
        eta = lambda n, c: n * c if pol == 's' else n / c
        ei, eo = eta(n_in, ci), eta(n_out, co)
        m11 = np.ones_like(s); m12 = np.zeros_like(s); m21 = np.zeros_like(s); m22 = np.ones_like(s)
        for n, d in layers:
            if d <= 0.5: continue
            c = cosang(n); c = np.where(c.imag < 0, -c, c); e = eta(n, c); dl = 2 * np.pi * n * d * c / lam
            a11, a12, a21, a22 = np.cos(dl), 1j * np.sin(dl) / e, 1j * e * np.sin(dl), np.cos(dl)
            m11, m12, m21, m22 = m11 * a11 + m12 * a21, m11 * a12 + m12 * a22, m21 * a11 + m22 * a21, m21 * a12 + m22 * a22
        Bv = m11 + m12 * eo; Cv = m21 + m22 * eo
        r = (ei * Bv - Cv) / (ei * Bv + Cv); t = 2 * ei / (ei * Bv + Cv)
        out.append((np.abs(r) ** 2, np.real(eo) / np.real(ei) * np.abs(t) ** 2))
    R = 0.5 * (out[0][0] + out[1][0]); T = 0.5 * (out[0][1] + out[1][1]); return R, T
# ---------------------------------------------------------------- slab-angle bins -> air-angle bins (power uniform in sin^2 within a bin)
def to_air(n):
    C = np.zeros((90, 90)); ea = np.sin(np.radians(np.arange(91))) ** 2
    for k in range(90):
        lo, hi = n ** 2 * np.sin(np.radians(k)) ** 2, n ** 2 * np.sin(np.radians(k + 1)) ** 2
        if lo >= 1: continue
        hi_c = min(hi, 1.0)
        for a in range(90):
            ov = max(0.0, min(hi_c, ea[a + 1]) - max(lo, ea[a]))
            if ov > 0: C[a, k] += ov / (hi - lo)
    return C       # columns sum to the fraction of the bin inside the escape cone (rest cannot reach air; Td is ~0 there)
def dbr_layers(x, j):
    return [((N.TIO2[j] if k % 2 == 0 else N.SIO2[j]), float(d)) for k, d in enumerate(x)]
# ---------------------------------------------------------------- air distribution per substrate photon for one wavelength
def air_dist(P, R, BT, BR):
    v0 = P / P.sum(); M = BR * R[None, :] if False else BR * R[:, None]
    # M[i, j] = BR[i, j] R[i]: incident j -> returned in i -> device reflectance at i
    U = np.linalg.solve(np.eye(90) - M, v0); return BT @ U
def composite(T1, R1, Rd, Td, C):
    A = np.linalg.solve(np.eye(90) - R1 * Rd[None, :], T1)           # (I - R1 diag(Rd))^-1 T1
    BT = C @ (Td[:, None] * A); BR = R1 + T1 @ (Rd[:, None] * A); return BT, BR
def evaluate(dev, mode, lam0, fwhm, q, step=2, x=None, slab=None):
    lams = dev['lam']; sel = np.arange(0, len(lams), step // 2); w = spectrum(lams[sel], lam0, fwhm)
    n = float(dev['n_sub']); Cn = to_air(n); D = np.zeros(90); Fb = 0.0
    if mode == 'mla':
        Bm = N.B[:, :, int(round((n - 1.30) / 0.05))]; BT, BR = Bm[:90], Bm[90:180][::-1]
    if mode == 'slabdbr': T1, R1 = slab[0], slab[1]
    for wi, i in zip(w, sel):
        lam = float(lams[i]); j = int(lam - 600); P, R = dev['P'][i], dev['R'][i]
        if mode == 'flat':
            Rf, Tf = tmm_RT(n, [], lam, TH); BT, BR = Cn * Tf[None, :], np.diag(Rf)
        elif mode == 'flatdbr':
            Rf, Tf = tmm_RT(n, dbr_layers(x, j), lam, TH); BT, BR = Cn * Tf[None, :], np.diag(Rf)
        elif mode == 'slabdbr':
            Rd, Td = tmm_RT(n, dbr_layers(x, j), lam, TH) if x is not None and len(x) else tmm_RT(n, [], lam, TH)
            BT, BR = composite(T1, R1, Rd, Td, Cn)
        D += wi * dev['F'][i] * dev['es'][i] * air_dist(P, R, BT, BR); Fb += wi * dev['F'][i]
    return q * D / (1 - q + q * Fb)
def cone(D, t):
    s2 = np.sin(np.radians(np.arange(91))) ** 2; cum = np.r_[0, np.cumsum(D)]; return float(np.interp(np.sin(np.radians(t)) ** 2, s2, cum))
def irradiance(D, z, a, nr=400):
    sr = 2 * np.pi * (np.cos(np.radians(np.arange(90))) - np.cos(np.radians(np.arange(1, 91)))); I = D / sr     # per sr, per unit area
    rho = (np.arange(nr) + 0.5) * a / nr; th = np.degrees(np.arctan(rho / z)); Ith = np.interp(th, TH, I)
    return float(np.sum(2 * np.pi * rho * (a / nr) * Ith * np.cos(np.radians(th)) / (z ** 2 + rho ** 2)))
def metrics(D):
    m = dict(EQE=float(D.sum()))
    for t in (5, 10, 20, 30): m['C%d' % t] = cone(D, t); m['Leq%d' % t] = cone(D, t) / np.sin(np.radians(t)) ** 2
    for z in (1, 3, 10): m['E_z%d_a1' % z] = irradiance(D, z, 1.0)
    return m
