"""Lab-device simulation: bare substrate vs hemispherical vs freeform MLA (finite geometry, weighted Monte Carlo).

Geometry: emitting disc r = 1 mm; 25 x 25 mm glass, blackened sides (absorbing); bare = 0.7 mm glass, planar front;
MLA cases = 0.7 mm device glass + 0.7 mm lens glass with index-matching oil (1.4 mm), lenses on the central 10 x 10 mm,
planar glass/air outside. Back face: full stack under the pixel; outside the pixel either the Al reflector (same stack,
non-emitting) or ITO/organics/air. Front-face encounters use the traced 1-degree kernels (escape into air bins, return
polar angle and change of lateral azimuth); planar glass/air uses Fresnel. Output counts only light leaving the front
face (goniometer). Spectral and Purcell weighting: EQE = q sum_l s F f_sub eta_out / (q sum_l s F + 1 - q), q = 0.98.
Writes results.json and results.npz (air angular distributions)."""
import sys, json, numpy as np
sys.path.insert(0, '../selectivity'); sys.path.insert(0, '.')
NS = 1.51; TAN = np.tan(np.radians(np.arange(90) + 0.5)); W_ET = np.diff(np.sin(np.radians(np.arange(91))) ** 2); W_ET /= W_ET.sum()
Q = 0.98; E = np.arange(91)


def planar_kernel(n=NS, sub=400):
    BT = np.zeros((90, 90))
    for i in range(90):
        t = np.radians(np.linspace(i, i + 1, sub + 2)[1:-1]); wt = np.sin(2 * t); s = n * np.sin(t); ok = s < 1
        ci = np.cos(t); ct = np.sqrt(np.clip(1 - s ** 2, 0, 1))
        rs = (n * ci - ct) / (n * ci + ct); rp = (n * ct - ci) / (n * ct + ci); T = np.where(ok, 1 - 0.5 * (rs ** 2 + rp ** 2), 0)
        a = np.degrees(np.arcsin(np.clip(s, 0, 1))).astype(int).clip(0, 89)
        np.add.at(BT[i], a[ok], (T * wt)[ok] / wt.sum())
    return dict(BT=BT, p=BT.sum(1), cdf=np.tile(np.eye(90).cumsum(1), (1, 1))[:, :] if False else np.cumsum(np.eye(90), 1), cond=None)


def lens_kernel(name):
    z = np.load(f'kernel_{name}.npz'); BT, K = z['BT'], z['K']; NA = K.shape[2]
    M = K.sum(2); cdf = np.cumsum(M / np.maximum(M.sum(1, keepdims=True), 1e-15), 1); cdf[:, -1] = 1.0
    cond = np.cumsum(K, 2) / np.maximum(K.sum(2, keepdims=True), 1e-15); cond[..., -1] = 1.0
    return dict(BT=BT, p=BT.sum(1), cdf=cdf, cond=cond, NA=NA)


def mc(S, PL, Rin, Rout, P, t_mm, a_mm=1.0, W_mm=12.5, Wp_mm=5.0, N=100000, seed=0, maxit=5000):
    """Returns escape weight per incident bin for the structure (inside patch) and the planar face, and side loss."""
    rng = np.random.default_rng(seed); a, W, Wp = a_mm / t_mm, W_mm / t_mm, Wp_mm / t_mm
    th = np.minimum(np.searchsorted(np.cumsum(P) / P.sum(), rng.random(N)), 89)
    r0 = a * np.sqrt(rng.random(N)); f0 = 2 * np.pi * rng.random(N)
    x = r0 * np.cos(f0); y = r0 * np.sin(f0); phi = 2 * np.pi * rng.random(N); g = np.ones(N)
    x += TAN[th] * np.cos(phi); y += TAN[th] * np.sin(phi)
    eS = np.zeros(90); eP = np.zeros(90); side = 0.0
    for it in range(maxit):
        out = (np.abs(x) > W) | (np.abs(y) > W); side += g[out].sum()
        k = ~out; th, x, y, phi, g = th[k], x[k], y[k], phi[k], g[k]
        if g.size == 0: break
        ins = (np.abs(x) <= Wp) & (np.abs(y) <= Wp)
        pp = np.where(ins, S['p'][th], PL['p'][th])
        np.add.at(eS, th[ins], (g * pp)[ins]); np.add.at(eP, th[~ins], (g * pp)[~ins]); g = g * (1 - pp)
        u = rng.random(g.size); new = th.copy()
        if ins.any():
            new[ins] = np.minimum((S['cdf'][th[ins]] < u[ins, None]).sum(1), 89)
            c = S['cond'][th[ins], new[ins]]; kk = np.minimum((c < rng.random(ins.sum())[:, None]).sum(1), S['NA'] - 1)
            phi[ins] = phi[ins] + (kk + rng.random(ins.sum())) * 2 * np.pi / S['NA']
        th = new                                                           # planar: specular, azimuth kept
        x += TAN[th] * np.cos(phi); y += TAN[th] * np.sin(phi)
        out = (np.abs(x) > W) | (np.abs(y) > W); side += g[out].sum(); g[out] = 0.0
        g = g * np.where(np.hypot(x, y) <= a, Rin[th], Rout[th])
        x += TAN[th] * np.cos(phi); y += TAN[th] * np.sin(phi)
        k = g > 1e-7; th, x, y, phi, g = th[k], x[k], y[k], phi[k], g[k]
    return eS / N, eP / N, side / N


def infinite(S, R, P):
    M = S['cdf'].copy(); M[:, 1:] -= S['cdf'][:, :-1]                      # M[i, r]
    BR = (M * (1 - S['p'])[:, None]).T
    return float(np.linalg.solve(np.eye(90) - R[:, None] * BR, P) @ S['p'])


def ybar(l):
    g = lambda x, m, s1, s2: np.exp(-0.5 * ((x - m) / np.where(x < m, s1, s2)) ** 2)
    return 0.821 * g(l, 568.8, 46.9, 40.5) + 0.286 * g(l, 530.9, 16.3, 31.1)


if __name__ == '__main__':
    dv = np.load('device.npz'); L = dv['lam']; s = dv['s']; F = dv['F']; fsub = dv['f_sub']
    PL = planar_kernel(); KER = {'bare': PL, 'hemisphere': lens_kernel('hemisphere'), 'freeform': lens_kernel('freeform')}
    T_MM = {'bare': 0.7, 'hemisphere': 1.4, 'freeform': 1.4}
    dOm = 2 * np.pi * (np.cos(np.radians(E[:-1])) - np.cos(np.radians(E[1:])))
    norm = Q * np.sum(s * F) + 1 - Q
    res = {}; ang = {}
    for nm, S in KER.items():
        for outside in ('Al', 'noAl'):
            air = np.zeros((len(L), 90)); eta = np.zeros(len(L)); side = np.zeros(len(L)); einf = np.zeros(len(L))
            for k in range(len(L)):
                P, Rin = dv['P'][k], dv['R_pix'][k]; Rout = Rin if outside == 'Al' else dv['R_noAl'][k]
                eS, eP, sd = mc(S, PL, Rin, Rout, P, T_MM[nm], Wp_mm=(0.0 if nm == 'bare' else 5.0), seed=k)
                BTn = lambda K: K['BT'] / np.maximum(K['p'][:, None], 1e-15)
                air[k] = eS @ BTn(S) + eP @ BTn(PL); eta[k] = air[k].sum(); side[k] = sd; einf[k] = infinite(S, Rin, P)
            wq = Q * s * F * fsub / norm                                   # photons into the substrate per electron, per lambda
            A = wq[:, None] * air                                          # photons per electron into each air bin
            I0 = (A[:, :5].sum(1) / dOm[:5].sum())                         # photons / sr / electron near 0 deg (0-5 deg average)
            cdA = float(np.sum(I0 * 683 * ybar(L) * 1239.84 / L))
            band = {b: float(A[:, lo:hi].sum()) for b, (lo, hi) in {'0-20': (0, 20), '20-40': (20, 40), '40-60': (40, 60), '60-80': (60, 80), '80-90': (80, 90)}.items()}
            key = f'{nm}|{outside}'
            res[key] = dict(EQE=float(A.sum()), cdA=cdA, EQE_lambertian=float(np.pi * I0.sum()), bands=band,
                            eta_out_weighted=float(np.sum(wq * eta) / wq.sum()), side_loss_weighted=float(np.sum(wq * side) / wq.sum()),
                            eta_out_infinite_weighted=float(np.sum(wq * einf) / wq.sum()), substrate_fraction_weighted=float(wq.sum()))
            ang[key] = A.sum(0) / dOm                                        # photons / sr / electron vs air angle
            print(key, json.dumps({k: (round(v, 4) if isinstance(v, float) else {b: round(x, 4) for b, x in v.items()}) for k, v in res[key].items()}), flush=True)
    res['_meta'] = dict(PLQY=Q, purcell_weighted=float(np.sum(s * F)), q_eff=float(Q * np.sum(s * F) / norm))
    json.dump(res, open('results.json', 'w'), indent=1); np.savez('results.npz', **{k.replace('|', '_'): v for k, v in ang.items()})
