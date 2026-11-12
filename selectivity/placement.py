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
Writes data/placement.json and ../figures_oe/fig7_placement.{png,pdf,xlsx}."""
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

# ---- Figure 7
INK2 = '#52514e'
CLS7 = {'planar': ('#7d7c76', 's', 'planar'), 'scatter': ('#1baf7a', 'o', 'scattering layer'), 'MLA hemi': ('#2a78d6', 'D', 'hemispherical MLA'),
        'MLA free': ('#eb6834', '^', 'freeform MLA'), 'random lens': ('#6250d6', 'v', 'random MLA'), 'wavy': ('#e87ba4', 'P', 'wavy surface')}
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw=dict(wspace=0.32)); rows_a = []
xs, ys = [], []
for k, v in S.items():
    nm, dev, lam = k.split('|')
    if dev != 'Al' or lam != '550': continue
    x = v['p'] - v['p_mix']; y = (v['dev_mix'] - v['exit_mix']) * 100; col, mk, _ = CLS7[v['cls']]
    ax[0].plot(x, y, mk, color=col, ms=4, mec='white', mew=0.3); xs.append(x); ys.append(y)
    rows_a.append([nm, v['cls'], v['p'], v['p_mix'], x, v['eta'], v['exit_mix'], v['dev_mix'], v['rand'], y])
r = np.corrcoef(xs, ys)[0, 1]
for c, (col, mk, lab) in CLS7.items(): ax[0].plot([], [], mk, color=col, label=lab)
ax[0].set_xlabel('p − p_mix = σ²/(1 − p)'); ax[0].set_ylabel('device-side − exit-side mixing (pp of η_sta)')
ax[0].set_title(f'(a) placement of the mixing (Al, 550 nm), r = {r:.2f}', loc='left', fontsize=6.5); ax[0].legend(fontsize=4.8, frameon=False, loc='upper left', bbox_to_anchor=(0.0, 0.92))
for R, ls in ((0.90, '-'), (0.98, '--')):
    sl = [s for s in slab if s['R'] == R]; n = [s['n'] for s in sl]
    ax[1].plot(n, [s['none'] for s in sl], ls, color='#7d7c76', lw=1, label='no mixing' if R == 0.90 else None)
    ax[1].plot(n, [s['exit_mix'] for s in sl], ls, color='#eb6834', lw=1, label='complete mixing at the exit face' if R == 0.90 else None)
    ax[1].plot(n, [s['dev_mix'] for s in sl], ls, color='#2a78d6', lw=1, label='complete mixing at the back reflector (= fully randomized model)' if R == 0.90 else None)
ax[1].plot([], [], 'k-', lw=1, label='R = 0.90'); ax[1].plot([], [], 'k--', lw=1, label='R = 0.98')
ax[1].set_xlabel('refractive index n of the slab'); ax[1].set_ylabel('extraction efficiency'); ax[1].set_ylim(0, 1); ax[1].set_xlim(1.4, 3.6)
ax[1].set_title('(b) slab with planar exit face and back reflector', loc='left', fontsize=6.5); ax[1].legend(fontsize=4.8, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
fig.savefig('../figures_oe/fig7_placement.png', dpi=300, bbox_inches='tight'); fig.savefig('../figures_oe/fig7_placement.pdf', bbox_inches='tight'); plt.close(fig)
wb = Workbook(); ws_ = wb.active; ws_.title = 'a_structures_Al550'
ws_.append(['structure', 'class', 'p', 'p_mix', 'p - p_mix', 'eta_actual', 'eta_exit_side_mixing', 'eta_device_side_mixing', 'eta_fully_randomized', 'device - exit (pp)'])
[ws_.append(rw) for rw in rows_a]
w2 = wb.create_sheet('b_slab'); w2.append(['n', 'R', 'p', 'p_mix', 'eta_no_mixing', 'eta_exit_side_mixing', 'eta_back_side_mixing (= fully randomized)'])
[w2.append([s['n'], s['R'], s['p'], s['p_mix'], s['none'], s['exit_mix'], s['dev_mix']]) for s in slab]
wb.save('../figures_oe/fig7_placement.xlsx'); print('corr', r)
