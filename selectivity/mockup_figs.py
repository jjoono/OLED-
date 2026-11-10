"""Mock-ups of Figs. 1-4 for the Optics Express storyline (angular mixing number of passive outcouplers).
Uses whatever data exist now (mixing validation, random sweep so far); final figures will be regenerated."""
import os, json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
import sel, mixing

os.makedirs('figures/mockup', exist_ok=True)
plt.rcParams.update({'font.size': 7, 'axes.titlesize': 7.5, 'axes.spines.top': False, 'axes.spines.right': False})
INK2 = '#52514e'
CLS = {'planar': ('#7d7c76', 's'), 'scatter': ('#1baf7a', 'o'), 'MLA hemi': ('#2a78d6', 'D'), 'MLA free': ('#eb6834', '^'),
       'random lens': ('#6250d6', 'v'), 'wavy': ('#e87ba4', 'P')}
ST = mixing.structures()
DV = np.load('data/devices.npz')
SW = json.load(open('../multipass_reset/random_sweep/summary.json')) if os.path.exists('../multipass_reset/random_sweep/summary.json') else {}
VAL = json.load(open('data/mixing_validation.json'))


def cls_of(name):
    for k in ('planar', 'scatter', 'MLA hemi', 'MLA free'):
        if name.startswith(k):
            return k
    return 'wavy' if name.startswith('W_') else 'random lens'


def kernel(BR):
    return BR / np.maximum(BR.sum(0, keepdims=True), 1e-12)


# ---------------- Fig. 1: concept -----------------------------------------------------------------
fig = plt.figure(figsize=(7.2, 2.6)); gs = fig.add_gridspec(1, 4, width_ratios=[1.35, 1, 1, 1], wspace=0.45)
a = fig.add_subplot(gs[0]); a.axis('off'); a.set_xlim(0, 10); a.set_ylim(0, 10)
a.add_patch(Rectangle((0.5, 0.6), 9, 1.0, color='#d9d8d3')); a.text(5, 1.1, "OLED stack  (R_LED(θ), round-trip loss A')", ha='center', va='center', fontsize=6)
a.add_patch(Rectangle((0.5, 1.6), 9, 4.6, color='#eef4fb')); a.text(9.3, 3.9, 'substrate', rotation=90, ha='center', va='center', fontsize=6, color=INK2)
xs = np.linspace(0.5, 9.5, 10)
for x0 in xs[:-1]:
    t = np.linspace(0, np.pi, 30); a.fill_between(x0 + 0.5 + 0.5 * np.cos(t)[::-1], 6.2, 6.2 + 0.9 * np.sin(t), color='#bcd3ef', lw=0)
a.text(5, 7.6, 'passive outcoupler:  p(θ)  and  B_R(θ_i→θ_r)', ha='center', fontsize=6)
for x, y, x2, y2, c in ((2.2, 1.7, 3.2, 6.1, '#2a78d6'), (3.2, 6.1, 3.9, 9.0, '#2a78d6'), (3.2, 6.1, 4.6, 1.7, '#eb6834'), (4.6, 1.7, 6.0, 6.1, '#eb6834'), (6.0, 6.1, 6.6, 9.0, '#eb6834')):
    a.add_patch(FancyArrowPatch((x, y), (x2, y2), arrowstyle='-|>', mutation_scale=6, color=c, lw=0.9))
a.set_title('(a) multi-pass extraction\nas a Markov chain on angle', loc='left', fontsize=6.5)
for k, (name, lab) in enumerate((('planar', 'planar  (λ₂ = 1)'), ('MLA hemi_total', 'hemispherical MLA'), ('scatter g0.8_t1.0', 'scattering layer'))):
    ax = fig.add_subplot(gs[k + 1]); BR = ST[name][1]; M = kernel(BR); g = ST[name][3]
    ext = [0, 90, 0, 90]
    ax.imshow(np.log10(M + 1e-4), origin='lower', extent=ext, cmap='Blues', vmin=-3, vmax=0, aspect='auto')
    lam2, _ = mixing.mixing_defs(BR, ST[name][2], np.ones(BR.shape[0]))
    ax.set_title(f"({'bcd'[k]}) {lab}\nλ₂ = {lam2:.2f},  N_mix = {1 / max(1 - lam2, 1e-9):.1f}" if lam2 < 0.999 else f"({'bcd'[k]}) {lab}\nN_mix = ∞", loc='left', fontsize=6.5)
    ax.set_xlabel('incident angle θ_i (deg)'); ax.set_ylabel('returned angle θ_r (deg)' if k == 0 else '')
fig.savefig('figures/mockup/fig1_concept.png', dpi=250, bbox_inches='tight'); plt.close(fig)

# ---------------- Fig. 2: validation of the reduced model ----------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4), gridspec_kw=dict(wspace=0.4))
for k, (key, lab) in enumerate((('eta_closed', 'memoryless model, Eq. (2)'), ('eta_m0', 'angle-preserving (m = 0)'), ('eta_red_lam2', 'reduced model, m = 1 − λ₂'))):
    for r in VAL:
        c, mk = CLS[cls_of(r['structure'])]
        ax[k].plot(r['eta_full'], r[key], mk, color=c, ms=2.5, mew=0, alpha=0.8)
    ax[k].plot([0.3, 1], [0.3, 1], color=INK2, lw=0.6); ax[k].set_xlim(0.3, 1); ax[k].set_ylim(0.3, 1)
    e = np.array([abs(r[key] - r['eta_full']) for r in VAL]) * 100
    ax[k].set_title(f"({'abc'[k]}) {lab}\nmean |err| {e.mean():.1f} pp, max {e.max():.0f} pp", loc='left', fontsize=6.5)
    ax[k].set_xlabel('η_sta, full angle-resolved series'); ax[k].set_ylabel('η_sta, model' if k == 0 else '')
for c, (col, mk) in CLS.items():
    if c in ('planar', 'scatter', 'MLA hemi', 'MLA free'):
        ax[0].plot([], [], mk, color=col, label=c)
ax[0].legend(fontsize=5, frameon=False, loc='upper left')
fig.savefig('figures/mockup/fig2_validation.png', dpi=250, bbox_inches='tight'); plt.close(fig)

# ---------------- Fig. 3: (p, N_mix) map ------------------------------------------------------------------
pts = []
for r in VAL:
    if r['device'] == 'Al' and r['lam'] == 550.0:
        pts.append((cls_of(r['structure']), r['p_lamb'], 1 / max(1 - r['lam2'], 1e-9), r['eta_full']))
for name, m in SW.items():
    pts.append((cls_of(name), m['p'], m['N_mix'], m['eta_Al']))
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw=dict(wspace=0.35))
# contours of the reduced model for a flat p(theta) (escape cone transmission scaled) on the Al device
P0 = DV['Al_P'][75].reshape(90, 2).sum(1); R1 = DV['Al_R'][75].reshape(90, 2).mean(1); w1 = sel.WL.reshape(90, 2).sum(1)
Rf, Tf = sel.tmm_RT([], 550.0, np.arange(90) + .5, 1.51)
pg = np.linspace(0.18, 0.46, 30); Ng = np.logspace(0, 1.6, 30); Z = np.zeros((30, 30))
for i, p in enumerate(pg):
    pT = Tf * p / (Tf @ w1)
    for j, Nm in enumerate(Ng):
        Z[j, i] = mixing.reduced(pT, w1, R1, P0, 1 - (1 - 1 / Nm))
cs = ax[0].contour(pg, Ng, Z, levels=[0.5, 0.6, 0.65, 0.7, 0.75, 0.8], colors=INK2, linewidths=0.5)
ax[0].clabel(cs, fontsize=5, fmt='%.2f')
for c, p, Nm, e in pts:
    col, mk = CLS[c]; ax[0].plot(p, min(Nm, 40), mk, color=col, ms=4, mec='white', mew=0.3)
for c, (col, mk) in CLS.items():
    ax[0].plot([], [], mk, color=col, label=c)
ax[0].set_yscale('log'); ax[0].set_xlabel('escape probability p (Lambertian)'); ax[0].set_ylabel('angular mixing number N_mix = 1/(1−λ₂)')
ax[0].set_title('(a) structure map, Al device (contours: η_sta)', loc='left'); ax[0].legend(fontsize=5, frameon=False, loc='upper left')
ax[0].text(0.45, 40, 'planar (∞)', ha='right', va='top', fontsize=5, color=INK2)
sw_s = [(float(k.split('shape')[1]), m['N_mix']) for k, m in SW.items() if k.startswith('A_shape')]
sw_j = [(float(k.split('jitter')[1].split('_')[0]), m['N_mix']) for k, m in SW.items() if k.startswith('C_jitter')]
sw_w = [(float(k.split('slope')[1]), m['N_mix']) for k, m in SW.items() if k.startswith('W_')]
for data, lab, col in ((sw_s, 'lenslet shape spread s', '#6250d6'), (sw_j, 'position jitter (fraction of pitch)', '#2a78d6'), (sw_w, 'wavy surface rms slope', '#e87ba4')):
    if data:
        x, y = zip(*sorted(data)); ax[1].plot(x, y, 'o-', color=col, ms=3, label=lab)
ax[1].axhline(1.2, color='#1baf7a', ls=':', lw=0.8); ax[1].text(0.02, 1.25, 'scattering layer', fontsize=5, color='#1baf7a')
ax[1].set_yscale('log'); ax[1].set_xlabel('randomness parameter'); ax[1].set_ylabel('N_mix')
ax[1].set_title('(b) randomness erases angular memory', loc='left'); ax[1].legend(fontsize=5, frameon=False)
fig.savefig('figures/mockup/fig3_map.png', dpi=250, bbox_inches='tight'); plt.close(fig)

# ---------------- Fig. 4: radiance bound and composition window -----------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.6), gridspec_kw=dict(wspace=0.35))
P05 = DV['Al_P'][75]; dE = np.diff(sel.S2); n = 1.51
def gbound(Pv, a_, b_):
    need = (np.sin(np.radians(b_)) ** 2 - np.sin(np.radians(a_)) ** 2) / n ** 2; L = Pv / dE; acc = pw = 0
    for k in np.argsort(-L):
        t = min(dE[k], need - acc)
        if t <= 0: break
        pw += L[k] * t; acc += t
    return pw
bands = [(0, 20), (20, 40), (40, 60), (60, 80)]
GB = [gbound(P05, *b) for b in bands]; LB = [(np.sin(np.radians(b)) ** 2 - np.sin(np.radians(a_)) ** 2) / n ** 2 for a_, b in bands]
x = np.arange(4)
ax[0].bar(x - 0.3, LB, 0.18, color='#d9d8d3', label='uniform-radiance bound (Yablonovitch)')
ax[0].bar(x - 0.1, GB, 0.18, color='#7d7c76', label='radiance-ordered bound (actual cavity)')
P1 = P05.reshape(90, 2).sum(1)
for k, (name, col, lab) in enumerate((('MLA hemi_total', '#2a78d6', 'hemispherical MLA'), ('scatter g0.8_t1.0', '#1baf7a', 'scattering layer'))):
    BT = ST[name][0]
    if ST[name][3] == '1':
        MB = np.load('../multipass_reset/bsdfs.npz'); first = P1 @ MB['hemi_total_BT']; v = [first[a_:b_].sum() for a_, b_ in bands]
    else:
        BTm, BRm = None, None
        Rf2, Tf2 = sel.tmm_RT([], 550.0, sel.TH, n); SL = np.load('data/slabs.npz'); Ts, Rs = SL['g0.8_t1.0_w1']
        U = np.linalg.solve(np.eye(sel.NB) - Rs * Rf2[None, :], Ts); first = Tf2 * (U @ P05)
        F = [np.clip((np.sin(np.radians(np.minimum(sel.EDG[1:], np.degrees(np.arcsin(np.sin(np.radians(b_)) / n))))) ** 2 - np.maximum(sel.S2[:-1], np.sin(np.radians(np.degrees(np.arcsin(np.sin(np.radians(a_)) / n)))) ** 2)) / dE, 0, 1) for a_, b_ in bands]
        v = [first @ f for f in F]
    ax[0].bar(x + 0.1 + 0.2 * k, v, 0.18, color=col, label=lab)
ax[0].set_xticks(x); ax[0].set_xticklabels(['0–20°', '20–40°', '40–60°', '60–80°']); ax[0].set_ylabel('first-pass power into band / substrate power')
ax[0].set_title('(a) band-resolved radiance bound (Al, 550 nm)', loc='left'); ax[0].legend(fontsize=5, frameon=False)
rows = [(cls_of(r['structure']), r['eta_full'], None) for r in VAL if r['device'] == 'Al' and r['lam'] == 550.0]
for name, m in SW.items():
    col, mk = CLS[cls_of(name)]; ax[1].plot(m['eta_Al'], m['S4060_Al'], mk, color=col, ms=4, mec='white', mew=0.3)
ax[1].axhspan(0.33, 0.36, color='#eef4fb', zorder=0); ax[1].axhline(0.337, color=INK2, ls=':', lw=0.7)
ax[1].text(0.02, 0.03, 'band: 0.33–0.36;  dotted: Lambertian partition 0.337', fontsize=5, color=INK2, transform=ax[1].transAxes)
ax[1].set_ylim(0.30, 0.40)
ax[1].set_xlabel('η_sta (Al device)'); ax[1].set_ylabel('selectivity S_40–60')
ax[1].set_title('(b) composition window of non-selective passive structures', loc='left', fontsize=6.5)
fig.savefig('figures/mockup/fig4_bound_window.png', dpi=250, bbox_inches='tight'); plt.close(fig)
print('mockups saved')
