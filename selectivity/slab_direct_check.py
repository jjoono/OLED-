"""Independent check of the adding-method scattering-layer matrices: direct Monte Carlo of an index-matched HG slab
(n = 1.51) with the flat Fresnel exit face to air included in the same trace.  Rays enter from the substrate side at
an angle drawn uniformly in sin^2 within each 1-deg bin, scatter, and are followed until they escape into air
(binned by air angle) or return into the substrate (binned by substrate angle).  Compared with the library matrices
(slab MC + adding with the Fresnel face, 0.5 deg -> 1 deg).  Output: data/slab_direct_check.json"""
import json, numpy as np
from numba import njit, prange
import mixing
src = open('oe_figs.py').read(); exec(src[:src.index('# ------------------------------------------------------------------ Fig. 1')])

@njit(cache=True)
def _hg(g, u):
    if abs(g) < 1e-6:
        return 2 * u - 1
    s = (1 - g * g) / (1 - g + 2 * g * u)
    return (1 + g * g - s * s) / (2 * g)

@njit(parallel=True)
def run(S, g, N, seed):
    n = 1.51; BT = np.zeros((90, 90)); BR = np.zeros((90, 90))
    for i in prange(90):
        np.random.seed(seed + i)
        a2 = np.sin(np.radians(i)) ** 2; b2 = np.sin(np.radians(i + 1)) ** 2
        for k in range(N):
            th = np.arcsin(np.sqrt(a2 + (b2 - a2) * np.random.random())); ph = 2 * np.pi * np.random.random()
            ux = np.sin(th) * np.cos(ph); uy = np.sin(th) * np.sin(ph); uz = np.cos(th); z = 0.0
            for ev in range(100000):
                d = -np.log(np.random.random() + 1e-300) / S
                zn = z + d * uz
                if zn >= 1.0:                                   # reaches the exit face
                    ci = uz; st = n * np.sqrt(max(0.0, 1 - ci * ci))
                    if st >= 1.0:
                        Rf = 1.0
                    else:
                        ct = np.sqrt(1 - st * st)
                        rs = (n * ci - ct) / (n * ci + ct); rp = (n * ct - ci) / (n * ct + ci); Rf = 0.5 * (rs * rs + rp * rp)
                    if np.random.random() >= Rf:
                        j = min(int(np.degrees(np.arcsin(min(1.0, st)))), 89); BT[i, j] += 1.0; break
                    z = 1.0; uz = -uz                            # specular reflection back into the slab
                    continue
                if zn <= 0.0:
                    j = min(int(np.degrees(np.arccos(min(1.0, -uz)))), 89); BR[i, j] += 1.0; break
                z = zn
                ct = _hg(g, np.random.random()); stt = np.sqrt(max(0.0, 1 - ct * ct)); p2 = 2 * np.pi * np.random.random()
                if abs(uz) > 0.99999:
                    ux, uy, uz = stt * np.cos(p2), stt * np.sin(p2), ct * np.sign(uz)
                else:
                    dd = np.sqrt(1 - uz * uz)
                    nx = stt * (ux * uz * np.cos(p2) - uy * np.sin(p2)) / dd + ux * ct
                    ny = stt * (uy * uz * np.cos(p2) + ux * np.sin(p2)) / dd + uy * ct
                    nz = -stt * np.cos(p2) * dd + uz * ct
                    ux, uy, uz = nx, ny, nz
    return BT / N, BR / N

out = {}
for key, g, tau in (('scatter g0.8_t1.0', 0.8, 1.0), ('scatter g0.0_t0.25', 0.0, 0.25)):
    S = tau / (1 - g)
    BTd, BRd = run(S, g, 20000, 11); BRd = BRd.T                   # BR[out, in]
    _, BTl, BRl = LIB[key]
    P = DV['Al_P'][75].reshape(90, 2).sum(1); R = DV['Al_R'][75].reshape(90, 2).mean(1)
    r = {}
    for lab, BT, BR in (('adding', BTl, BRl), ('direct', BTd, BRd)):
        pT = BT.sum(1); l2 = mixing.lam2_cut(BR); pi, pt = mixing.pinf_cut(BR)
        x = np.linalg.solve(np.eye(90) - R[:, None] * BR, P); air = x @ BT
        r[lab] = dict(p=float(pT @ (W1 / W1.sum())), lam2=l2, N_mix=1 / (1 - l2), xi=pi / pt, eta=float(air.sum()), S4060=float(air[40:60].sum() / air.sum()))
    r['max_abs_dp_theta'] = float(np.abs(BTl.sum(1) - BTd.sum(1))[:85].max())
    Ml = BRl / np.maximum(BRl.sum(0, keepdims=True), 1e-12); Md = BRd / np.maximum(BRd.sum(0, keepdims=True), 1e-12)
    r['M_rel_frobenius_diff'] = float(np.linalg.norm((Ml - Md)[:85, :85]) / np.linalg.norm(Ml[:85, :85]))
    out[key] = r
    print(key, {k: (round(v, 4) if isinstance(v, float) else {a: round(b, 4) for a, b in v.items()}) for k, v in r.items()}, flush=True)
json.dump(out, open('data/slab_direct_check.json', 'w'), indent=1)
