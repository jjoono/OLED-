"""Goniometer with a bare photodiode (no imaging optics): square PD of side s at D = 59 mm from the rotation axis
point C = (0, 0, z_c); the sample rotates about an axis through C. In the observation plane of each photon
(axisymmetric pixel), a photon leaving the front face at P with polar angle theta_a hits the PD for goniometer angles
theta_g in [theta_a - asin((s/2 - b)/D), theta_a + asin((s/2 + b)/D)], b = (P - C).u (parallax), and only if its
offset perpendicular to the observation plane is <= s/2. The measured EQE is the integral of the measured intensity,
EQE_m = sum_photons e * |window within [0, theta_max]| / (2 asin(s / 2D)); the measured angular profile is binned at the
window centre. PE = EQE_m / (pi I_m(0-5 deg))."""
import json, numpy as np
import simulate as SM, gonio as G
from simulate import TAN
D = 59.0; SIZES = (3.0, 5.0, 10.0); AXES = G.AXES; TMAX = (90, 80, 70)
dv, L, wq = G.dv, G.L, G.wq


def run(S, PL, SA, PA, Rin, Rout, P, t_mm, Wp_mm, N=60000, seed=0, maxit=5000, a_mm=1.0, W_mm=12.5):
    rng = np.random.default_rng(seed); a, W, Wp = a_mm / t_mm, W_mm / t_mm, Wp_mm / t_mm
    th = np.minimum(np.searchsorted(np.cumsum(P) / P.sum(), rng.random(N)), 89)
    r0 = a * np.sqrt(rng.random(N)); f0 = 2 * np.pi * rng.random(N)
    x = r0 * np.cos(f0); y = r0 * np.sin(f0); phi = 2 * np.pi * rng.random(N); g = np.ones(N)
    x += TAN[th] * np.cos(phi); y += TAN[th] * np.sin(phi)
    zc = np.array(list(AXES.values())) * t_mm
    E = np.zeros((len(AXES), len(SIZES), len(TMAX))); H = np.zeros((len(AXES), len(SIZES), 90)); tot = 0.0
    for it in range(maxit):
        out = (np.abs(x) > W) | (np.abs(y) > W); k = ~out; th, x, y, phi, g = th[k], x[k], y[k], phi[k], g[k]
        if g.size == 0: break
        ins = (np.abs(x) <= Wp) & (np.abs(y) <= Wp)
        pp = np.where(ins, S['p'][th], PL['p'][th]); e = g * pp; tot += e.sum()
        u = rng.random(g.size); j = np.empty_like(th)
        j[ins] = np.minimum((SA['cA'][th[ins]] < u[ins, None]).sum(1), 89); j[~ins] = np.minimum((PA['cA'][th[~ins]] < u[~ins, None]).sum(1), 89)
        tha = np.radians(j + rng.random(g.size)); pha = phi.copy()
        if ins.any() and 'cT' in SA:
            c = SA['cT'][th[ins], j[ins]]; kk = np.minimum((c < rng.random(ins.sum())[:, None]).sum(1), SA['NA'] - 1)
            pha[ins] = phi[ins] + (kk + rng.random(ins.sum())) * 2 * np.pi / SA['NA']
        X = x * t_mm; Y = y * t_mm
        xp = X * np.cos(pha) + Y * np.sin(pha); yp = -X * np.sin(pha) + Y * np.cos(pha)
        for ia, z in enumerate(zc):
            b = xp * np.cos(tha) - (t_mm - z) * np.sin(tha)
            for isz, s in enumerate(SIZES):
                okY = np.abs(yp) <= s / 2; w0 = 2 * np.arcsin(s / (2 * D))
                lo = tha - np.arcsin(np.clip((s / 2 - b) / D, -1, 1)); hi = tha + np.arcsin(np.clip((s / 2 + b) / D, -1, 1))
                for itm, tm in enumerate(TMAX):
                    ov = np.clip(np.minimum(hi, np.radians(tm)) - np.maximum(lo, 0.0), 0, None)
                    E[ia, isz, itm] += (e * okY * ov).sum() / w0
                cen = np.degrees(tha + b / D)
                H[ia, isz] += np.histogram(cen, np.arange(91), weights=e * okY)[0]
        g = g * (1 - pp); u = rng.random(g.size); new = th.copy()
        if ins.any():
            new[ins] = np.minimum((S['cdf'][th[ins]] < u[ins, None]).sum(1), 89)
            c = S['cond'][th[ins], new[ins]]; kk = np.minimum((c < rng.random(ins.sum())[:, None]).sum(1), S['NA'] - 1)
            phi[ins] = phi[ins] + (kk + rng.random(ins.sum())) * 2 * np.pi / S['NA']
        th = new; x += TAN[th] * np.cos(phi); y += TAN[th] * np.sin(phi)
        out = (np.abs(x) > W) | (np.abs(y) > W); g[out] = 0.0
        g = g * np.where(np.hypot(x, y) <= a, Rin[th], Rout[th])
        x += TAN[th] * np.cos(phi); y += TAN[th] * np.sin(phi)
        k = g > 1e-7; th, x, y, phi, g = th[k], x[k], y[k], phi[k], g[k]
    return E / N, H / N, tot / N


import os
PATCH_MM = float(os.environ.get('PATCH_MM', '10'))          # full width of the square lens patch
SFX = '' if PATCH_MM == 10 else f'_patch{PATCH_MM:g}'

if __name__ == '__main__':
    PL = SM.planar_kernel(); PA = G.kernel_air(PL); KER = {'bare': (PL, PA)}
    for nm in ('hemisphere', 'freeform'):
        K = SM.lens_kernel(nm); K['KT'] = np.load(f'kernel_{nm}.npz')['KT']; KER[nm] = (K, G.kernel_air(K))
    T = {'bare': 0.7, 'hemisphere': 1.4, 'freeform': 1.4}
    res = dict(D_mm=D, sizes_mm=list(SIZES), axes=list(AXES), theta_max=list(TMAX))
    dOm = 2 * np.pi * (1 - np.cos(np.radians(5)))
    for nm, (S, SA) in KER.items():
        for o in ('Al', 'noAl'):
            Ea = np.zeros((len(AXES), len(SIZES), len(TMAX))); Ha = np.zeros((len(AXES), len(SIZES), 90)); tot = 0.0
            for k in range(0, len(L), 2):
                Rin = dv['R_pix'][k]; Rout = Rin if o == 'Al' else dv['R_noAl'][k]
                Ek, Hk, tk = run(S, PL, SA, PA, Rin, Rout, dv['P'][k], T[nm], 0.0 if nm == 'bare' else PATCH_MM / 2, seed=200 + k)
                Ea += wq[k] * Ek; Ha += wq[k] * Hk; tot += wq[k] * tk
            sc = wq.sum() / wq[::2].sum(); Ea *= sc; Ha *= sc; tot *= sc
            PE = Ea[:, :, 0] / (np.pi * Ha[:, :, :5].sum(2) / dOm)
            res[f'{nm}|{o}'] = dict(total=tot, EQE_meas=Ea.tolist(), PE=PE.tolist(), hist=Ha.tolist())
            print(nm, o, f'total {tot*100:.1f}%', ' | '.join(f"{ax}: " + ', '.join(f"s{s:g} {Ea[i, j, 0]*100:.1f}/PE {PE[i, j]:.2f}" for j, s in enumerate(SIZES)) for i, ax in enumerate(AXES)), flush=True)
    json.dump(res, open(f'gonio_pd{SFX}.json', 'w'))
