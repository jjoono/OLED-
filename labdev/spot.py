"""Where the light leaves the front face: encircled fraction of the front-face output (all angles, and 0-10 deg in air)
versus radius from the pixel centre, for the six lab cases (spectrally weighted as in simulate.py)."""
import json, numpy as np
import simulate as SM
from simulate import TAN, mc as _unused
dv = np.load('device.npz'); L = dv['lam']; s = dv['s']; F = dv['F']; fsub = dv['f_sub']
Q = SM.Q; wq = Q * s * F * fsub / (Q * np.sum(s * F) + 1 - Q)
RB = np.linspace(0, 12.5 * np.sqrt(2), 361)                 # radius bins (mm)


def mc_pos(S, PL, Rin, Rout, P, t_mm, a_mm=1.0, W_mm=12.5, Wp_mm=5.0, N=60000, seed=0, maxit=5000):
    rng = np.random.default_rng(seed); a, W, Wp = a_mm / t_mm, W_mm / t_mm, Wp_mm / t_mm
    th = np.minimum(np.searchsorted(np.cumsum(P) / P.sum(), rng.random(N)), 89)
    r0 = a * np.sqrt(rng.random(N)); f0 = 2 * np.pi * rng.random(N)
    x = r0 * np.cos(f0); y = r0 * np.sin(f0); phi = 2 * np.pi * rng.random(N); g = np.ones(N)
    x += TAN[th] * np.cos(phi); y += TAN[th] * np.sin(phi)
    hall = np.zeros(RB.size - 1); hfwd = np.zeros(RB.size - 1)
    fS = S['BT'][:, :10].sum(1) / np.maximum(S['p'], 1e-15); fP = PL['BT'][:, :10].sum(1) / np.maximum(PL['p'], 1e-15)
    for it in range(maxit):
        out = (np.abs(x) > W) | (np.abs(y) > W); k = ~out; th, x, y, phi, g = th[k], x[k], y[k], phi[k], g[k]
        if g.size == 0: break
        ins = (np.abs(x) <= Wp) & (np.abs(y) <= Wp)
        pp = np.where(ins, S['p'][th], PL['p'][th]); e = g * pp; rr = np.hypot(x, y) * t_mm
        hall += np.histogram(rr, RB, weights=e)[0]; hfwd += np.histogram(rr, RB, weights=e * np.where(ins, fS[th], fP[th]))[0]
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
    return hall / N, hfwd / N


PL = SM.planar_kernel(); KER = {'bare': PL, 'hemisphere': SM.lens_kernel('hemisphere'), 'freeform': SM.lens_kernel('freeform')}
T = {'bare': 0.7, 'hemisphere': 1.4, 'freeform': 1.4}; out = {'r_mm': RB[1:].tolist()}
for nm, S in KER.items():
    for o in ('Al', 'noAl'):
        Ha = np.zeros(RB.size - 1); Hf = np.zeros(RB.size - 1)
        for k in range(0, len(L), 2):                           # every 20 nm
            Rin = dv['R_pix'][k]; Rout = Rin if o == 'Al' else dv['R_noAl'][k]
            a, f = mc_pos(S, PL, Rin, Rout, dv['P'][k], T[nm], Wp_mm=(0.0 if nm == 'bare' else 5.0), seed=k)
            Ha += wq[k] * a; Hf += wq[k] * f
        Ha *= wq.sum() / wq[::2].sum(); Hf *= wq.sum() / wq[::2].sum()
        out[f'{nm}|{o}'] = dict(EQE_within=np.cumsum(Ha).tolist(), EQE0_10_within=np.cumsum(Hf).tolist())
        ce = np.cumsum(Ha)
        print(f'{nm:10s} {o:5s} EQE total {ce[-1]*100:.1f}% | within r = 1, 1.5, 2, 3, 5 mm: ' + ', '.join(f'{np.interp(r, RB[1:], ce)*100:.1f}' for r in (1, 1.5, 2, 3, 5)), flush=True)
json.dump(out, open('spot.json', 'w'))
