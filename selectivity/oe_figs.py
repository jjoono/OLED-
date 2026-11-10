"""Final Figs. 1-4 (+ Excel source data) for the angular-mixing manuscript.  Fig. 5 is the benchmark figure of the
freeform study (make_figures.py).  N_mix uses lambda2 on theta < 85 deg (mixing.lam2_cut)."""
import os, glob, json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from openpyxl import Workbook
import sel, mixing

OUT = '../figures_oe'; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 7, 'axes.titlesize': 7.5, 'axes.spines.top': False, 'axes.spines.right': False})
INK2 = '#52514e'; N_SUB = 1.51
CLS = {'planar': ('#7d7c76', 's', 'planar'), 'scatter': ('#1baf7a', 'o', 'scattering layer'), 'MLA hemi': ('#2a78d6', 'D', 'hemispherical MLA'),
       'MLA free': ('#eb6834', '^', 'freeform MLA'), 'random lens': ('#6250d6', 'v', 'random MLA'), 'wavy': ('#e87ba4', 'P', 'wavy surface')}
DV = np.load('data/devices.npz'); SL = np.load('data/slabs.npz'); MB = np.load('../multipass_reset/bsdfs.npz')
BANDS = [(0, 20), (20, 40), (40, 60), (60, 80)]
DEVS = ('Al', 'Ag', 'LL15'); WLS = ((50, 500.0), (75, 550.0), (100, 600.0))


def save_xlsx(name, sheets):
    wb = Workbook(); wb.remove(wb.active)
    for title, header, rows in sheets:
        ws = wb.create_sheet(title[:31]); ws.append(header)
        for r in rows:
            ws.append([float(x) if isinstance(x, (np.floating, np.integer)) else x for x in r])
    wb.save(f'{OUT}/{name}.xlsx')


# ------------------------------------------------------------------ structure library on a common 1-deg grid
def to_1deg_air(BTsub, n=N_SUB):
    """Transmission per (input 0.5-deg bin) and exit substrate bin -> air-angle 1-deg histogram per input."""
    ta = np.degrees(np.arcsin(np.clip(n * np.sin(np.radians(sel.TH)), 0, 1)))
    ok = n * np.sin(np.radians(sel.TH)) < 1
    out = np.zeros((BTsub.shape[1], 90))
    idx = np.minimum(ta[ok].astype(int), 89)
    for k, j in zip(np.where(ok)[0], idx):
        out[:, j] += BTsub[k]
    return out                                                     # [in, air]


def library():
    lib = {}; I = np.eye(sel.NB); Rf, Tf = sel.tmm_RT([], 550.0, sel.TH, N_SUB)
    def slab(Ts, Rs):
        U = np.linalg.solve(I - Rs * Rf[None, :], Ts)
        return to_1deg_air(Tf[:, None] * U), Rs + (Ts * Rf[None, :]) @ U          # BTair[in, air], BR[out, in]
    def coarse(BTair, BR):                                                           # 0.5 -> 1 deg, Lambertian within bin
        w = sel.WL; W2 = w.reshape(90, 2); f = (W2 / W2.sum(1, keepdims=True)).ravel()
        BTc = (BTair * f[:, None]).reshape(90, 2, 90).sum(1)
        BRc = (BR * f[None, :]).reshape(90, 2, 90, 2).sum(axis=(1, 3))
        return BTc, BRc
    lib['planar'] = ('planar',) + coarse(*slab(I, np.zeros_like(I)))
    for k in ('g0.0_t0.25_w1', 'g0.0_t1.0_w1', 'g0.8_t0.25_w1', 'g0.8_t1.0_w1', 'g0.9_t1.0_w1', 'g0.8_t4.0_w1'):
        lib['scatter ' + k.replace('_w1', '')] = ('scatter',) + coarse(*slab(*SL[k]))
    for d in sorted({k.rsplit('_', 1)[0] for k in MB.files if k.endswith('_BT')} - {'free_0-20'}):
        lib['MLA ' + d] = ('MLA hemi' if d.startswith('hemi') else 'MLA free', MB[d + '_BT'], MB[d + '_BR'].T)
    for f in sorted(glob.glob('../multipass_reset/random_sweep/*.npz')):
        nm = os.path.basename(f)[:-4]; z = np.load(f)
        lib[nm] = ('wavy' if nm.startswith('W_') else 'random lens', z['BT'], z['BR'].T)
    return lib


W1 = sel.WL.reshape(90, 2).sum(1)


def evaluate(BTair, BR, dev, i):
    P = DV[f'{dev}_P'][i].reshape(90, 2).sum(1); R = DV[f'{dev}_R'][i].reshape(90, 2).mean(1)
    pT = BTair.sum(1)
    x = np.linalg.solve(np.eye(90) - R[:, None] * BR, P); air = x @ BTair
    lam2 = mixing.lam2_cut(BR)
    return dict(eta=float(air.sum()), S=[float(air[a:b].sum() / air.sum()) for a, b in BANDS],
                closed=mixing.closed(pT, W1, R), m0=mixing.reduced(pT, W1, R, P, 0.0), m1=mixing.reduced(pT, W1, R, P, 1.0), red=mixing.reduced(pT, W1, R, P, 1 - lam2),
                red_w=mixing.reduced(pT, W1, R, P, 1 - lam2, recip=False), pinf=PINF[id(BR)][0], p85=PINF[id(BR)][1],
                p=float(pT @ (W1 / W1.sum())), lam2=lam2, first=P @ BTair)


LIB = library()
PINF = {id(v[2]): mixing.pinf_cut(v[2]) for v in LIB.values()}
EV = {(nm, dev, lam): evaluate(BT, BR, dev, i) for nm, (c, BT, BR) in LIB.items() for dev in DEVS for i, lam in WLS}

# ------------------------------------------------------------------ Fig. 1
fig = plt.figure(figsize=(7.2, 2.7)); gs = fig.add_gridspec(1, 4, width_ratios=[1.35, 1, 1, 1], wspace=0.5)
a = fig.add_subplot(gs[0]); a.axis('off'); a.set_xlim(0, 10); a.set_ylim(0, 10)
a.add_patch(Rectangle((0.5, 0.6), 9, 1.0, color='#d9d8d3')); a.text(5, 1.1, "OLED: R(θ), round-trip loss A'", ha='center', va='center', fontsize=6)
a.add_patch(Rectangle((0.5, 1.6), 9, 4.6, color='#eef4fb')); a.text(9.15, 3.9, 'substrate', rotation=90, ha='center', va='center', fontsize=6, color=INK2)
for x0 in np.linspace(0.5, 8.5, 9):
    t = np.linspace(0, np.pi, 30); a.fill_between(x0 + 0.5 + 0.5 * np.cos(t)[::-1], 6.2, 6.2 + 0.9 * np.sin(t), color='#bcd3ef', lw=0)
a.text(5, 7.7, 'passive structure: p(θ), B_R(θᵢ→θᵣ)', ha='center', fontsize=6)
for (x, y, x2, y2, c) in ((2.2, 1.7, 3.2, 6.1, '#2a78d6'), (3.2, 6.1, 3.9, 9.0, '#2a78d6'), (3.2, 6.1, 4.6, 1.7, '#eb6834'), (4.6, 1.7, 6.0, 6.1, '#eb6834'), (6.0, 6.1, 6.6, 9.0, '#eb6834')):
    a.add_patch(FancyArrowPatch((x, y), (x2, y2), arrowstyle='-|>', mutation_scale=6, color=c, lw=0.9))
a.set_title('(a) multi-pass extraction as a\nMarkov chain on angle', loc='left', fontsize=6.5)
sh = []
for k, (nm, lab) in enumerate((('planar', 'planar'), ('MLA hemi_total', 'hemispherical MLA'), ('scatter g0.8_t1.0', 'scattering layer'))):
    ax = fig.add_subplot(gs[k + 1]); BR = LIB[nm][2]; M = kernel = BR / np.maximum(BR.sum(0, keepdims=True), 1e-12)
    ax.imshow(np.log10(M + 1e-4), origin='lower', extent=[0, 90, 0, 90], cmap='Blues', vmin=-3, vmax=0, aspect='auto')
    l2 = mixing.lam2_cut(BR)
    ax.set_title(f"({'bcd'[k]}) {lab}\n" + (f"λ₂ = {l2:.2f},  N_mix = {1 / (1 - l2):.1f}" if l2 < 0.999 else "λ₂ = 1,  N_mix = ∞"), loc='left', fontsize=6.5)
    ax.set_xlabel('incident angle θᵢ (deg)'); ax.set_ylabel('returned angle θᵣ (deg)' if k == 0 else '')
    sh.append((f'{lab} M', ['row = returned bin, col = incident bin (1 deg)'], M.tolist()))
fig.savefig(f'{OUT}/fig1_mixing_concept.png', dpi=300, bbox_inches='tight'); fig.savefig(f'{OUT}/fig1_mixing_concept.pdf', bbox_inches='tight'); plt.close(fig)
save_xlsx('fig1_mixing_concept', sh)

# ------------------------------------------------------------------ Fig. 2
fig, ax = plt.subplots(1, 4, figsize=(7.2, 2.1), gridspec_kw=dict(wspace=0.38)); rows = []
for (nm, dev, lam), e in EV.items():
    rows.append([nm, LIB[nm][0], dev, lam, e['eta'], e['closed'], e['m0'], e['red_w'], e['red'], e['p'], e['lam2'], e['pinf']])
for k, (j, lab) in enumerate(((5, 'memoryless formula'), (6, 'angle preserving, m = 0'), (7, 'm = 1 − λ₂, q ∝ w'), (8, 'reciprocal, q ∝ w(1 − p)'))):
    for r in rows:
        col, mk, _ = CLS[r[1]]; ax[k].plot(r[4], r[j], mk, color=col, ms=2.4, mew=0, alpha=0.8)
    err = np.array([abs(r[j] - r[4]) for r in rows]) * 100
    ax[k].plot([0.3, 1], [0.3, 1], color=INK2, lw=0.6); ax[k].set_xlim(0.3, 1); ax[k].set_ylim(0.3, 1)
    ax[k].set_title(f"({'abcd'[k]}) {lab}\nmean |err| {err.mean():.2f} pp, max {err.max():.1f} pp", loc='left', fontsize=5.6)
    ax[k].set_xlabel('η_sta, full series'); ax[k].set_ylabel('η_sta, model' if k == 0 else '')
for c, (col, mk, lab) in CLS.items():
    ax[0].plot([], [], mk, color=col, label=lab)
ax[0].legend(fontsize=4.8, frameon=False, loc='lower right')
fig.savefig(f'{OUT}/fig2_reduced_model.png', dpi=300, bbox_inches='tight'); fig.savefig(f'{OUT}/fig2_reduced_model.pdf', bbox_inches='tight'); plt.close(fig)
save_xlsx('fig2_reduced_model', [('cases', ['structure', 'class', 'device', 'lambda_nm', 'eta_full', 'eta_memoryless', 'eta_m0', 'eta_reduced_w', 'eta_reduced_recip', 'p_Lamb', 'lambda2_85', 'p_inf'], rows)])
print(f'Fig 2: {len(rows)} cases; reduced(w) mean {np.mean([abs(r[7]-r[4]) for r in rows])*100:.2f} max {np.max([abs(r[7]-r[4]) for r in rows])*100:.2f}; recip mean {np.mean([abs(r[8]-r[4]) for r in rows])*100:.2f} pp, max {np.max([abs(r[8]-r[4]) for r in rows])*100:.2f}; memoryless mean {np.mean([abs(r[5]-r[4]) for r in rows])*100:.2f}, max {np.max([abs(r[5]-r[4]) for r in rows])*100:.1f}')

# ------------------------------------------------------------------ Fig. 3
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw=dict(wspace=0.42, width_ratios=[1.25, 1, 1])); rows = []
P0 = DV['Al_P'][75].reshape(90, 2).sum(1); R1 = DV['Al_R'][75].reshape(90, 2).mean(1)
pT0 = LIB['MLA hemi_total'][1].sum(1)
pg = np.linspace(0.18, 0.47, 30); Ng = np.logspace(0, 1.4, 30); Z = np.zeros((30, 30))
for i, p in enumerate(pg):
    pT = pT0 * p / (pT0 @ (W1 / W1.sum()))
    for j, Nm in enumerate(Ng):
        Z[j, i] = mixing.reduced(np.clip(pT, 0, 1), W1, R1, P0, 1 / Nm)
cs = ax[0].contour(pg, Ng, Z, levels=[0.6, 0.7, 0.75, 0.78, 0.8, 0.82], colors=INK2, linewidths=0.5); [ax[0].text(sg[np.argmin(np.abs(sg[:, 1] - 15)), 0], 15, f'{lv:.2f}', fontsize=5, color=INK2, ha='center', va='center', rotation=80, bbox=dict(fc='white', ec='none', pad=0.3))
 for lv, segs in zip(cs.levels, cs.allsegs) for sg in segs[:1] if len(sg)]
for nm, (c, BT, BR) in LIB.items():
    e = EV[(nm, 'Al', 550.0)]; Nm = 1 / max(1 - e['lam2'], 1e-9)
    col, mk, _ = CLS[c]; ax[0].plot(e['p'], min(Nm, 25), mk, color=col, ms=4, mec='white', mew=0.3)
    rows.append([nm, c, e['p'], e['lam2'], Nm, e['pinf'], e['pinf'] / e['p85'], e['eta'], EV[(nm, 'Ag', 550.0)]['eta'], EV[(nm, 'LL15', 550.0)]['eta']])
    ax[2].plot(min(Nm, 25), e['pinf'] / e['p85'], mk, color=col, ms=4, mec='white', mew=0.3)
for c, (col, mk, lab) in CLS.items():
    ax[0].plot([], [], mk, color=col, label=lab)
ax[0].set_yscale('log'); ax[0].set_ylim(1, 26); ax[0].set_xlabel('escape probability p (uniform radiance)'); ax[0].set_ylabel('angular mixing number N_mix')
ax[0].text(0.405, 24, '↑ planar: N_mix = ∞', fontsize=5, color=INK2)
ax[0].set_title('(a) structure map (η_sta, Al, 550 nm)', loc='left', fontsize=6.5); ax[0].legend(fontsize=4.8, frameon=True, framealpha=0.9, edgecolor='none', loc='lower left', bbox_to_anchor=(0.0, 0.08))
for pre, lab, col, key in (('A_shape', 'lenslet shape spread s', '#6250d6', 'shape'), ('C_jitter', 'position jitter / pitch', '#2a78d6', 'jitter'),
                           ('W_slope', 'wavy surface, rms slope', '#e87ba4', 'slope'), ('B_aspect', 'height spread', '#eda100', 'aspect')):
    pts = sorted((float(nm.split(key)[1].split('_')[0]), 1 / max(1 - EV[(nm, 'Al', 550.0)]['lam2'], 1e-9)) for nm in LIB if nm.startswith(pre) and not nm.startswith('C2'))
    if pts:
        x, y = zip(*pts); ax[1].plot(x, y, 'o-', color=col, ms=3, label=lab)
ns = 1 / (1 - EV[('scatter g0.8_t1.0', 'Al', 550.0)]['lam2'])
ax[1].axhline(ns, color='#1baf7a', ls=':', lw=0.8); ax[1].text(0.02, ns * 1.06, 'scattering layer (g 0.8, τ* 1)', fontsize=5, color='#1baf7a')
ax[1].axhline(1 / (1 - EV[('MLA hemi_total', 'Al', 550.0)]['lam2']), color='#2a78d6', ls=':', lw=0.8)
ax[1].set_yscale('log'); ax[1].set_yticks([1, 1.5, 2, 3, 5, 8]); ax[1].set_yticklabels(['1', '1.5', '2', '3', '5', '8']); ax[1].minorticks_off(); ax[1].set_xlabel('randomness parameter'); ax[1].set_ylabel('N_mix'); ax[1].set_title('(b) randomness vs angle retention', loc='left')
ax[1].legend(fontsize=4.6, frameon=False, loc='center right')
ax[2].set_xscale('log'); ax[2].set_xlim(1, 30); ax[2].set_ylim(-0.03, 1.08); ax[2].set_xticks([1, 2, 5, 10, 25]); ax[2].set_xticklabels(['1', '2', '5', '10', '∞'])
ax[2].minorticks_off(); ax[2].axhline(1, color=INK2, lw=0.5, ls=':')
ax[2].set_xlabel('N_mix'); ax[2].set_ylabel('late-pass escape ratio ξ = p∞ / p'); ax[2].set_title('(c) late-pass escape', loc='left')
ax[2].text(1.05, 0.03, 'planar (ξ = 0)', fontsize=5, color=INK2); ax[2].text(1.05, 1.01, 'no dead-angle trapping', fontsize=5, color=INK2)
fig.savefig(f'{OUT}/fig3_structure_map.png', dpi=300, bbox_inches='tight'); fig.savefig(f'{OUT}/fig3_structure_map.pdf', bbox_inches='tight'); plt.close(fig)
save_xlsx('fig3_structure_map', [('structures', ['structure', 'class', 'p_Lamb', 'lambda2_85', 'N_mix', 'p_inf', 'xi', 'eta_Al', 'eta_Ag', 'eta_LL15'], rows),
                                 ('contours', ['N_mix \\ p'] + list(pg), [[Ng[j]] + list(Z[j]) for j in range(30)])])

# ------------------------------------------------------------------ Fig. 4
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.6), gridspec_kw=dict(wspace=0.35)); rows = []
P05 = DV['Al_P'][75]; dE = np.diff(sel.S2)
def gbound(Pv, a_, b_):
    need = (np.sin(np.radians(b_)) ** 2 - np.sin(np.radians(a_)) ** 2) / N_SUB ** 2; L = Pv / dE; acc = pw = 0
    for k in np.argsort(-L):
        t = min(dE[k], need - acc)
        if t <= 0: break
        pw += L[k] * t; acc += t
    return pw
GB = [gbound(P05, *b) for b in BANDS]; UB = [(np.sin(np.radians(b)) ** 2 - np.sin(np.radians(a_)) ** 2) / N_SUB ** 2 for a_, b in BANDS]
x = np.arange(4)
ax[0].bar(x - 0.3, UB, 0.15, color='#d9d8d3', label='uniform-radiance bound'); ax[0].bar(x - 0.15, GB, 0.15, color='#7d7c76', label='radiance-ordered bound (cavity)')
for k, (nm, lab, col) in enumerate((('MLA hemi_total', 'hemispherical MLA', '#2a78d6'), ('W_slope1.0', 'wavy surface', '#e87ba4'), ('scatter g0.8_t1.0', 'scattering layer', '#1baf7a'), ('planar', 'planar', '#b5b4ad'))):
    if nm not in LIB: continue
    first = EV[(nm, 'Al', 550.0)]['first']; v = [first[a_:b_].sum() for a_, b_ in BANDS]
    ax[0].bar(x + 0.15 * k, v, 0.15, color=col, label=lab); rows.append([lab] + v + [first.sum()])
ax[0].set_xticks(x); ax[0].set_xticklabels(['0–20°', '20–40°', '40–60°', '60–80°']); ax[0].set_ylabel('first-pass power into band / substrate power')
ax[0].set_title('(a) band-resolved radiance bound (Al, 550 nm)', loc='left'); ax[0].set_ylim(0, 0.32); ax[0].legend(fontsize=4.8, frameon=False, ncol=2, loc='upper left')
rows = [['uniform bound'] + UB + [1 / N_SUB ** 2], ['radiance-ordered bound'] + GB + [gbound(P05, 0, 90)]] + rows
win = []
for nm, (c, BT, BR) in LIB.items():
    e = EV[(nm, 'Al', 550.0)]; col, mk, _ = CLS[c]
    ax[1].plot(e['eta'], e['S'][2], mk, color=col, ms=4, mec='white', mew=0.3); win.append([nm, c, e['eta']] + e['S'])
for c, (col, mk, lab) in CLS.items():
    ax[1].plot([], [], mk, color=col, label=lab)
ax[1].axhspan(0.33, 0.36, color='#eef4fb', zorder=0); ax[1].axhline(0.337, color=INK2, ls=':', lw=0.7)
ax[1].set_ylim(0.28, 0.40); ax[1].set_xlabel('η_sta (Al device, 550 nm)'); ax[1].set_ylabel('selectivity S₄₀₋₆₀ of the output')
ax[1].set_title('(b) composition window of passive structures', loc='left'); ax[1].legend(fontsize=4.8, frameon=False, loc='lower left', ncol=2)
fig.savefig(f'{OUT}/fig4_bound_window.png', dpi=300, bbox_inches='tight'); fig.savefig(f'{OUT}/fig4_bound_window.pdf', bbox_inches='tight'); plt.close(fig)
save_xlsx('fig4_bound_window', [('a_first_pass', ['curve', '0-20', '20-40', '40-60', '60-80', 'total'], rows),
                                ('b_composition', ['structure', 'class', 'eta_sta', 'S0-20', 'S20-40', 'S40-60', 'S60-80'], win)])
S = [w[5] for w in win]
print(f'Fig 4: S40-60 range {min(S):.3f}-{max(S):.3f} over {len(S)} structures; bound total {gbound(P05, 0, 90):.3f}')
json.dump({f'{k[0]}|{k[1]}|{k[2]}': {kk: vv for kk, vv in v.items() if kk != 'first'} for k, v in EV.items()}, open('data/oe_eval.json', 'w'))
