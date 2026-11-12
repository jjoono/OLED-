"""Where the mixing happens: complete reciprocal mixing at the escape surface versus at the device (back) side.

Exit-side mixing  : the structure returns light with the rank-one reciprocal kernel (1 - p_i) q_r, q ~ w (1 - p);
                    device reflects specularly with R(theta)            -> eta_mix (main text, Sec. 2.4)
Device-side mixing: the structure keeps its own return matrix B_R; the device reflects with the rank-one reciprocal
                    kernel D(i -> j) = R_i R_j w_j / sum(w R) (column sums R_i, w_i D(i->j) = w_j D(j->i))
                    -> mixed light ~ w R, escape p_D = p + Cov_w(R, p) / <R>_w
Fully randomized  : p / [p + (1 - p) A'].
Part 1: the 38 structures with the three OLED stacks (500/550/600 nm).
Part 2: generic slab of index n, isotropic internal source, planar exit face (Fresnel), angle-independent back
        reflectance R; no device physics.
Writes data/placement.json; panels (c) and (d) of main-text Fig. 2 are drawn by oe_figs.py."""
import os, json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from openpyxl import Workbook
import sel, mixing

src = open('oe_figs.py').read()
exec(src[:src.index('# ------------------------------------------------------------------ Fig. 1')])   # LIB, DV, W1, helpers

w = W1 / W1.sum(); I = np.eye(90)


def cases(BT, BR, dev, i):
    P = DV[f'{dev}_P'][i].reshape(90, 2).sum(1); R = DV[f'{dev}_R'][i].reshape(90, 2).mean(1); pT = BT.sum(1)
    eta = mixing.full(pT, BR, R, P) if False else float(np.linalg.solve(I - R[:, None] * BR, P) @ BT.sum(1))
    exit_mix = mixing.reduced(pT, W1, R, P, 1.0)
    D = np.outer(R * w, R) / (w * R).sum()                         # D[j, i]
    dev_mix = float(pT @ np.linalg.solve(I - D @ BR, P))
    rand = mixing.closed(pT, W1, R)
    p = float(pT @ w); cov = float((w * (R - R @ w) * (pT - p)).sum())
    q = w * (1 - pT); q /= q.sum()
    return dict(eta=eta, exit_mix=exit_mix, dev_mix=dev_mix, rand=rand, p=p, p_mix=float(q @ pT), p_D=p + cov / float(R @ w))


out = {'structures': {}}
for nm, (c, BT, BR) in LIB.items():
    for dev in DEVS:
        for i, lam in WLS:
            out['structures'][f'{nm}|{dev}|{lam:.0f}'] = dict(cls=c, **cases(BT, BR, dev, i))

# ---- Part 2: generic slab
def fresnel_p(n, edges):
    p = np.zeros(len(edges) - 1)
    for k in range(len(p)):
        t = np.linspace(np.radians(edges[k]), np.radians(edges[k + 1]), 201)[1:-1]; s = n * np.sin(t); ok = s < 1
        ci = np.cos(t); ct = np.sqrt(np.clip(1 - s ** 2, 0, 1))
        rs = (n * ci - ct) / (n * ci + ct); rp = (n * ct - ci) / (n * ct + ci)
        T = np.where(ok, 1 - 0.5 * (rs ** 2 + rp ** 2), 0); wt = np.sin(2 * t); p[k] = (T * wt).sum() / wt.sum()
    return p

E = np.linspace(0, 90, 361); ws = np.diff(np.sin(np.radians(E)) ** 2); ws /= ws.sum(); Is = np.eye(len(ws))
slab = []
for n in np.round(np.arange(1.4, 3.61, 0.1), 2):
    p = fresnel_p(n, E); v = ws.copy()
    for R in (0.90, 0.98):
        none = float(p @ np.linalg.solve(Is - np.diag(R * (1 - p)), v))
        q = ws * (1 - p); q /= q.sum(); pm = float(q @ p); pb = float(p @ ws)
        exit_mix = pb + (1 - pb) * R * pm / (1 - R * (1 - pm))     # isotropic source: first pass p, then mixed
        rand = pb / (pb + (1 - pb) * (1 - R))                      # = device-side complete mixing for constant R
        slab.append(dict(n=float(n), R=R, p=pb, p_mix=pm, none=none, exit_mix=exit_mix, dev_mix=rand, rand=rand))
out['slab'] = slab
json.dump(out, open('data/placement.json', 'w'), indent=1)

# ---- summary
S = out['structures']
def summ(dev, lam):
    rows = [v for k, v in S.items() if k.endswith(f'|{dev}|{lam}') and v['cls'] != 'planar']
    g_exit = np.array([v['exit_mix'] - v['eta'] for v in rows]) * 100; g_dev = np.array([v['dev_mix'] - v['eta'] for v in rows]) * 100
    return g_exit, g_dev
for dev in DEVS:
    for lam in ('500', '550', '600'):
        ge, gd = summ(dev, lam)
        print(f'{dev:5s} {lam}: exit-side gain {ge.min():+.1f}..{ge.max():+.1f} pp (median {np.median(ge):.1f}) | device-side gain {gd.min():+.1f}..{gd.max():+.1f} (median {np.median(gd):.1f})')
for nm in ('planar', 'MLA hemi_total', 'MLA free_total', 'A_shape1.0', 'W_slope0.2', 'W_slope1.0', 'scatter g0.8_t1.0'):
    v = S[f'{nm}|Al|550']; print(f"{nm:18s} Al550 eta {v['eta']:.3f} exit {v['exit_mix']:.3f} device {v['dev_mix']:.3f} rand {v['rand']:.3f} | p {v['p']:.3f} p_mix {v['p_mix']:.3f} p_D {v['p_D']:.3f}")
for r in slab:
    if r['n'] in (1.5, 2.4, 3.5): print(f"slab n {r['n']} R {r['R']}: p {r['p']:.3f} p_mix {r['p_mix']:.3f} | none {r['none']:.3f} exit {r['exit_mix']:.3f} device/rand {r['rand']:.3f}")

