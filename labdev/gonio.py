"""Goniometer collection geometry for the lab devices.
The detector sees, at each observation direction d, a circular measuring field of radius r0 perpendicular to d,
centred on the rotation axis point C = (0, 0, z_c) above the pixel centre (z = 0: emitting plane).
A photon leaving the front face (z = t) at P with direction d is counted if the distance from C to the ray line
is <= r0, i.e. |(P - C) x d| <= r0. Integrating the measured intensity over all directions then gives the
'measured EQE' as the sum of counted photons. Rotation axis: device plane (z_c = 0), middle (t/2), front/lens surface (t).
Air azimuth: planar face keeps the substrate azimuth; lenslet arrays use the traced transmitted-azimuth kernel KT.
Axisymmetry of the pixel lets every photon be evaluated in its own observation plane."""
import json, numpy as np
import simulate as SM
from simulate import TAN
dv = np.load('device.npz'); L = dv['lam']; s = dv['s']; F = dv['F']; fsub = dv['f_sub']
Q = SM.Q; wq = Q * s * F * fsub / (Q * np.sum(s * F) + 1 - Q)
R0 = np.array([0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0])
AXES = {'device plane': 0.0, 'middle': 0.5, 'front/lens surface': 1.0}       # fraction of t
TMAX = (90, 80, 70)                                                           # goniometer angle limits (deg)


def kernel_air(K):
    BTn = K['BT'] / np.maximum(K['p'][:, None], 1e-15); cA = np.cumsum(BTn, 1); cA[:, -1] = 1.0
    out = dict(cA=cA)
    if 'KT' in K:
        KT = K['KT']; c = np.cumsum(KT, 2) / np.maximum(KT.sum(2, keepdims=True), 1e-15); c[..., -1] = 1.0; out['cT'] = c; out['NA'] = KT.shape[2]
    return out


def run(S, PL, SA, PA, Rin, Rout, P, t_mm, Wp_mm, N=60000, seed=0, maxit=5000, a_mm=1.0, W_mm=12.5):
    rng = np.random.default_rng(seed); a, W, Wp = a_mm / t_mm, W_mm / t_mm, Wp_mm / t_mm
    th = np.minimum(np.searchsorted(np.cumsum(P) / P.sum(), rng.random(N)), 89)
    r0 = a * np.sqrt(rng.random(N)); f0 = 2 * np.pi * rng.random(N)
    x = r0 * np.cos(f0); y = r0 * np.sin(f0); phi = 2 * np.pi * rng.random(N); g = np.ones(N)
    x += TAN[th] * np.cos(phi); y += TAN[th] * np.sin(phi)
    acc = np.zeros((len(AXES), len(R0), len(TMAX))); tot = 0.0
    zc = np.array(list(AXES.values())) * t_mm
    for it in range(maxit):
        out = (np.abs(x) > W) | (np.abs(y) > W); k = ~out; th, x, y, phi, g = th[k], x[k], y[k], phi[k], g[k]
        if g.size == 0: break
        ins = (np.abs(x) <= Wp) & (np.abs(y) <= Wp)
        pp = np.where(ins, S['p'][th], PL['p'][th]); e = g * pp; tot += e.sum()
        # direction of the escaping photon
        u = rng.random(g.size); j = np.empty_like(th)
        j[ins] = np.minimum((SA['cA'][th[ins]] < u[ins, None]).sum(1), 89); j[~ins] = np.minimum((PA['cA'][th[~ins]] < u[~ins, None]).sum(1), 89)
        tha = np.radians(j + rng.random(g.size)); pha = phi.copy()
        if ins.any() and 'cT' in SA:
            c = SA['cT'][th[ins], j[ins]]; kk = np.minimum((c < rng.random(ins.sum())[:, None]).sum(1), SA['NA'] - 1)
            pha[ins] = phi[ins] + (kk + rng.random(ins.sum())) * 2 * np.pi / SA['NA']
        X = x * t_mm; Y = y * t_mm
        xp = X * np.cos(pha) + Y * np.sin(pha); yp = -X * np.sin(pha) + Y * np.cos(pha)    # observation plane = x'z
        st, ct = np.sin(tha), np.cos(tha)
        for ia, z in enumerate(zc):
            vz = t_mm - z
            # |v x d| with v = (xp, yp, vz), d = (st, 0, ct)
            cx = yp * ct; cy = vz * st - xp * ct; cz = -yp * st
            dist = np.sqrt(cx ** 2 + cy ** 2 + cz ** 2)
            for it2, tm in enumerate(TMAX):
                okang = tha <= np.radians(tm)
                acc[ia, :, it2] += np.array([(e * ((dist <= r) & okang)).sum() for r in R0])
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
    return acc / N, tot / N


if __name__ == '__main__':
    PL = SM.planar_kernel(); PA = kernel_air(PL)
    KER = {'bare': (PL, PA)}
    for nm in ('hemisphere', 'freeform'):
        K = SM.lens_kernel(nm); K['KT'] = np.load(f'kernel_{nm}.npz')['KT']; KER[nm] = (K, kernel_air(K))
    T = {'bare': 0.7, 'hemisphere': 1.4, 'freeform': 1.4}; res = {'r0_mm': R0.tolist(), 'axes': list(AXES), 'theta_max': list(TMAX)}
    for nm, (S, SA) in KER.items():
        for o in ('Al', 'noAl'):
            A = np.zeros((len(AXES), len(R0), len(TMAX))); tot = 0.0; ks = range(0, len(L), 2)
            for k in ks:
                Rin = dv['R_pix'][k]; Rout = Rin if o == 'Al' else dv['R_noAl'][k]
                a, t = run(S, PL, SA, PA, Rin, Rout, dv['P'][k], T[nm], 0.0 if nm == 'bare' else 5.0, seed=100 + k)
                A += wq[k] * a; tot += wq[k] * t
            sc = wq.sum() / wq[::2].sum(); A *= sc; tot *= sc
            res[f'{nm}|{o}'] = dict(total=tot, measured=A.tolist())
            print(nm, o, f'total {tot*100:.1f}%', ' | '.join(f"{ax}: r0=1.5 {A[i, 4, 0]*100:.1f} r0=2 {A[i, 6, 0]*100:.1f}" for i, ax in enumerate(AXES)), flush=True)
    json.dump(res, open('gonio.json', 'w'))
