"""Finite lateral size: the Markov chain on angle extended to (angle, lateral position) by weighted Monte Carlo.

Geometry (lengths in units of the substrate thickness t): emitting disc of radius a on the back face, square substrate
of half-width W, structure on a centred square patch of half-width Wp (planar glass/air outside), device reflectance
R(theta) of the OLED stack inside radius rc on the back face and R_out(theta) outside it. Light that reaches the
substrate edge is counted as lost. Each encounter with the front face uses the 1-degree return kernel of the structure
(escape p_i, return-angle distribution M[:, i]); between encounters the light moves laterally by 2 t tan(theta).
The azimuth of the lateral motion after an encounter with the structure is kept, randomized, or (hemispherical MLA)
drawn from the traced azimuth-resolved kernel (../multipass_reset/hftrace_az.py). With W, Wp, rc -> infinity the
estimator reproduces eta_sta of the main text.
Benchmark geometry: a = 1 mm, t = 1.295 mm, 25 x 25 mm patch = substrate, stack reflectance on the whole back face."""
import json, numpy as np, mixing
src = open('oe_figs.py').read(); exec(src[:src.index('# ------------------------------------------------------------------ Fig. 1')])
TAN = np.tan(np.radians(np.arange(90) + 0.5)); w = W1 / W1.sum(); I = np.eye(90)
T_MM = 1.295


def kernel(BT, BR, variant='actual'):
    p = BT.sum(1)
    if variant == 'actual':
        M = BR / np.maximum(BR.sum(0, keepdims=True), 1e-15)
    else:
        q = w * (1 - p) if variant == 'exit_mix' else w.copy()                 # exit_mix: reciprocal; rand: fully randomized
        M = np.repeat((q / q.sum())[:, None], 90, 1)
    cdf = np.cumsum(M, 0).T; cdf[:, -1] = 1.0                                 # cdf[i, r]
    return p, cdf


AZ = None
def load_az():
    global AZ
    z = np.load('../multipass_reset/azkernel_hemiMLA.npz'); K = z['K']; NA = K.shape[2]
    cond = np.cumsum(K, 2) / np.maximum(K.sum(2, keepdims=True), 1e-15); cond[..., -1] = 1.0
    AZ = dict(cond=cond, NA=NA, p=z['BT'], cdf=np.cumsum(K.sum(2) / np.maximum(K.sum((1, 2))[:, None], 1e-15), 1))
    AZ['cdf'][:, -1] = 1.0
    return AZ


def _fresnel_glass_air(c):
    """Unpolarized reflectance of glass (1.51) -> air for incidence cosine c on a side facet."""
    n = 1.51; s2 = n * n * (1 - c * c); tir = s2 >= 1; ct = np.sqrt(np.clip(1 - s2, 0, 1))
    rs = (n * c - ct) / (n * c + ct + 1e-30); rp = (n * ct - c) / (n * ct + c + 1e-30)
    return np.where(tir, 1.0, 0.5 * (rs * rs + rp * rp))


def _edges(x, y, phi, g, th, W, mode):
    """Apply the substrate side faces to positions outside |x|, |y| <= W. Returns updated arrays and the power that
    leaves through (fresnel) or is absorbed at (absorb) the sides."""
    lost = 0.0
    if mode == 'absorb':
        out = (np.abs(x) > W) | (np.abs(y) > W); lost = g[out].sum(); g = np.where(out, 0.0, g)
        return x, y, phi, g, lost
    st = np.sin(np.radians(th + 0.5))
    for _ in range(200):
        ox = np.abs(x) > W; oy = np.abs(y) > W
        if not (ox.any() or oy.any()): break
        for o, ax in ((ox, 0), (oy, 1)):
            if not o.any(): continue
            c = st[o] * np.abs(np.cos(phi[o]) if ax == 0 else np.sin(phi[o]))
            Rf = 1.0 if mode == 'mirror' else _fresnel_glass_air(c)
            lost += (g[o] * (1 - Rf)).sum(); g[o] = g[o] * Rf
            if ax == 0:
                x[o] = np.sign(x[o]) * 2 * W - x[o]; phi[o] = np.pi - phi[o]
            else:
                y[o] = np.sign(y[o]) * 2 * W - y[o]; phi[o] = -phi[o]
    return x, y, phi, g, lost


def mc(struct, planar, Rin, Rout, P, a=1 / T_MM, W=12.5 / T_MM, Wp=None, rc=np.inf, az='random', N=100000, seed=0, maxit=4000, edge='absorb'):
    """struct, planar: (p, cdf). edge: 'absorb', 'fresnel' (glass/air side faces; transmitted power is returned as the
    edge value), or 'mirror'. Returns (eta through the front face, power lost at or leaving through the sides)."""
    rng = np.random.default_rng(seed); Wp = W if Wp is None else Wp
    pS, cS = struct; pP, cP = planar
    th = np.minimum(np.searchsorted(np.cumsum(P) / P.sum(), rng.random(N)), 89)
    r0 = a * np.sqrt(rng.random(N)); f0 = 2 * np.pi * rng.random(N)
    x = r0 * np.cos(f0); y = r0 * np.sin(f0); phi = 2 * np.pi * rng.random(N); g = np.ones(N)
    x += TAN[th] * np.cos(phi); y += TAN[th] * np.sin(phi)
    esc = 0.0; side = 0.0
    for it in range(maxit):
        x, y, phi, g, l = _edges(x, y, phi, g, th, W, edge); side += l
        keep = g > 1e-7; th, x, y, phi, g = th[keep], x[keep], y[keep], phi[keep], g[keep]
        if g.size == 0: break
        ins = (np.abs(x) <= Wp) & (np.abs(y) <= Wp)
        pp = np.where(ins, pS[th], pP[th]); esc += (g * pp).sum(); g = g * (1 - pp)
        u = rng.random(g.size); new = np.empty_like(th)
        for msk, c in ((ins, cS), (~ins, cP)):
            if msk.any(): new[msk] = np.minimum((c[th[msk]] < u[msk, None]).sum(1), 89)
        if az == 'random':
            phi = np.where(ins, 2 * np.pi * rng.random(g.size), phi)
        elif az == 'traced' and ins.any():
            c = AZ['cond'][th[ins], new[ins]]; k = (c < rng.random(ins.sum())[:, None]).sum(1)
            phi[ins] = phi[ins] + (np.minimum(k, AZ['NA'] - 1) + rng.random(ins.sum())) * 2 * np.pi / AZ['NA']
        th = new
        x += TAN[th] * np.cos(phi); y += TAN[th] * np.sin(phi)
        x, y, phi, g, l = _edges(x, y, phi, g, th, W, edge); side += l
        g = g * np.where(np.hypot(x, y) <= rc, Rin[th], Rout[th])
        x += TAN[th] * np.cos(phi); y += TAN[th] * np.sin(phi)
        alive = g > 1e-7
        th, x, y, phi, g = th[alive], x[alive], y[alive], phi[alive], g[alive]
    return esc / N, side / N


def dev(name, i):
    return DV[f'{name}_P'][i].reshape(90, 2).sum(1), DV[f'{name}_R'][i].reshape(90, 2).mean(1)


SEL = ('planar', 'MLA hemi_total', 'MLA free_total', 'A_shape0.7', 'W_slope0.2', 'W_slope1.0', 'scatter g0.8_t1.0')
AZDEF = {'planar': 'keep', 'MLA hemi_total': 'traced'}
A0, W0 = 1 / T_MM, 12.5 / T_MM


def ref_inf(BT, BR, P, R, variant):
    pT = BT.sum(1)
    if variant == 'actual': return float(np.linalg.solve(I - R[:, None] * BR, P) @ pT)
    if variant == 'exit_mix': return float(mixing.reduced(pT, W1, R, P, 1.0))
    return float(mixing.closed(pT, W1, R))


if __name__ == '__main__':
    import mixing
    load_az(); out = {'geometry': dict(t_mm=T_MM, a_mm=1.0, W_mm=12.5), 'benchmark': [], 'size': [], 'pixel': []}
    BTp, BRp = LIB['planar'][1:]; PL = kernel(BTp, BRp); p_pl = BTp.sum(1)
    for nm in ('MLA hemi_total', 'W_slope0.2'):           # check of the infinite limit
        BT, BR = LIB[nm][1:]; P, R = dev('Al', 75)
        print(f'check {nm}: MC infinite {mc(kernel(BT, BR), PL, R, R, P, W=np.inf)[0]:.4f}, mirror sides {mc(kernel(BT, BR), PL, R, R, P, edge="mirror")[0]:.4f} vs eta_sta {ref_inf(BT, BR, P, R, "actual"):.4f}', flush=True)
    # ---- 1. benchmark geometry
    for d in DEVS:
        P, R = dev(d, 75)
        for nm in SEL:
            BT, BR = LIB[nm][1:]
            for var in ('actual', 'exit_mix', 'rand'):
                azs = [AZDEF.get(nm, 'random')] if var != 'actual' or nm == 'planar' else [AZDEF.get(nm, 'random'), 'keep'] + (['random'] if nm == 'MLA hemi_total' else [])
                if var != 'actual': azs = ['random']
                for az in azs:
                    for em in ('absorb', 'fresnel', 'mirror'):
                        e, ed = mc(kernel(BT, BR, var), PL, R, R, P, az=az, seed=1, edge=em)
                        out['benchmark'].append(dict(device=d, structure=nm, variant=var, az=az, edge_mode=em, eta=e, edge=ed, eta_inf=ref_inf(BT, BR, P, R, var)))
        print('benchmark', d, 'done', flush=True)
    # ---- 2. size of substrate = patch (stack reflectance everywhere), Al 550 nm
    P, R = dev('Al', 75)
    for nm in ('planar', 'MLA hemi_total', 'W_slope0.2', 'scatter g0.8_t1.0'):
        BT, BR = LIB[nm][1:]
        for var in ('actual', 'exit_mix', 'rand'):
            az = AZDEF.get(nm, 'random') if var == 'actual' else 'random'
            for Wt in (1, 2, 4, 8, 9.65, 16, 32, 64, 128, 256):
                e, ed = mc(kernel(BT, BR, var), PL, R, R, P, W=Wt, az=az, seed=2, edge='fresnel')
                out['size'].append(dict(structure=nm, variant=var, az=az, W_over_t=Wt, eta=e, edge=ed, eta_inf=ref_inf(BT, BR, P, R, var)))
    print('size done', flush=True)
    # ---- 3. stack reflectance only under the emitting disc (radius a); outside black or glass/air
    for rout_name, Rout in (('black', np.zeros(90)), ('glass/air', 1 - p_pl)):
        for nm in ('planar', 'MLA hemi_total', 'W_slope0.2', 'scatter g0.8_t1.0'):
            BT, BR = LIB[nm][1:]
            for var in ('actual', 'exit_mix', 'rand'):
                az = AZDEF.get(nm, 'random') if var == 'actual' else 'random'
                for at in (0.05, 0.1, 0.2, 0.5, A0, 2, 5, 10, 20, 50):
                    e, ed = mc(kernel(BT, BR, var), PL, R, Rout, P, a=at, W=at + W0, rc=at, az=az, seed=3, edge='fresnel')
                    out['pixel'].append(dict(outside=rout_name, structure=nm, variant=var, az=az, a_over_t=at, eta=e, edge=ed, eta_inf=ref_inf(BT, BR, P, R, var)))
    print('pixel done', flush=True)
    json.dump(out, open('data/finite_size.json', 'w'), indent=1)
