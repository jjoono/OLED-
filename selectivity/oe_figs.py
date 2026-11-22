"""Final Figs. 1-4 (+ Excel source data) for the angular-mixing manuscript.  Fig. 5 is the benchmark figure of the
freeform study (make_figures.py).  N_mix uses lambda2 on theta < 85 deg (mixing.lam2_cut)."""
import os, sys, json, glob, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from openpyxl import Workbook
import sel, mixing

OUT = '../figures_oe'; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, '..'); import figstyle as FS; FS.apply(); FW = FS.W
ETA = r'$\eta_{\mathrm{sta}}$'; NMIX = r'$N_{\mathrm{mix}}$'; PMIX = r'$p_{\mathrm{mix}}$'
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

# ================================================================== figures (Optics Express layout, see ../figstyle.py)
from matplotlib.lines import Line2D
CAP = lambda t: t[0].upper() + t[1:]
def class_legend(ax, **kw):
    return ax.legend(handles=FS.handles(CLS), **kw)

# ------------------------------------------------------------------ Fig. 1
fig = plt.figure(figsize=(FW, 4.55)); gs = fig.add_gridspec(2, 1, height_ratios=[0.74, 1.3], hspace=0.42)
g1 = gs[0].subgridspec(1, 5, width_ratios=[1.3, 0.3, 1, 1, 1], wspace=0.3)
a = fig.add_subplot(g1[0]); a.axis('off'); a.set_xlim(0, 10); a.set_ylim(0, 10)
a.add_patch(Rectangle((0.3, 0.3), 9.4, 1.5, color='#d9d8d3')); a.text(5, 1.05, r"OLED, $R(\theta)$", ha='center', va='center', fontsize=7.5)
a.add_patch(Rectangle((0.3, 1.8), 9.4, 4.4, color='#eef4fb')); a.text(8.6, 4.0, 'substrate', rotation=90, ha='center', va='center', fontsize=7, color=INK2)
for x0 in np.linspace(0.3, 8.13, 6):
    t = np.linspace(0, np.pi, 30); a.fill_between(x0 + 0.783 + 0.783 * np.cos(t)[::-1], 6.2, 6.2 + 1.3 * np.sin(t), color='#bcd3ef', lw=0)
a.text(5, 9.25, r'structure: $p(\theta)$, $B_R$', ha='center', fontsize=7.5)
for (x, y, x2, y2, c) in ((1.3, 1.9, 2.6, 6.1, '#2a78d6'), (2.6, 6.1, 3.4, 8.9, '#2a78d6'), (2.6, 6.1, 4.3, 1.9, '#eb6834'), (4.3, 1.9, 6.0, 6.1, '#eb6834'), (6.0, 6.1, 6.8, 8.9, '#eb6834')):
    a.add_patch(FancyArrowPatch((x, y), (x2, y2), arrowstyle='-|>', mutation_scale=7, color=c, lw=1.1))
FS.label(a, 'a', y=1.035)
sh = []; mats = []
for k, (nm, lab) in enumerate((('planar', 'Planar'), ('MLA hemi_total', 'Hemisph. MLA'), ('scatter g0.8_t1.0', 'Scattering'))):
    ax = fig.add_subplot(g1[k + 2]); mats.append(ax); BR = LIB[nm][2]; M = kernel = BR / np.maximum(BR.sum(0, keepdims=True), 1e-12)
    im = ax.imshow(np.log10(np.maximum(M, 1e-3)), origin='lower', extent=[0, 90, 0, 90], cmap='Blues', vmin=-3, vmax=0, aspect='equal')
    ax.set_xticks([0, 30, 60, 90]); ax.set_yticks([0, 30, 60, 90]); ax.set_xlabel(r'$\theta_i$ (°)')
    if k == 0: ax.set_ylabel(r'Returned angle $\theta_r$ (°)')
    else: ax.set_yticklabels([])
    FS.label(ax, 'bcd'[k]); l2 = mixing.lam2_cut(BR); print(f'Fig 1{"bcd"[k]} {nm}: lambda2 {l2:.3f} N_mix {1 / max(1 - l2, 1e-9):.1f}')
    sh.append((f'{lab} M', ['row = returned bin, col = incident bin (1 deg)'], M.tolist()))
cax = mats[-1].inset_axes([1.1, 0.0, 0.075, 1.0]); cb = fig.colorbar(im, cax=cax, ticks=[-3, -2, -1, 0]); cb.ax.set_yticklabels([r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', '1']); cb.ax.tick_params(labelsize=7, length=2)
cb.set_label(r'$M(\theta_r \leftarrow \theta_i)$', fontsize=7.5, labelpad=2)
# (e), (f): light incident on the structure at each pass (Al device, 550 nm)
P0 = DV['Al_P'][75].reshape(90, 2).sum(1); P0 = P0 / P0.sum(); Rv = DV['Al_R'][75].reshape(90, 2).mean(1); wn = W1 / W1.sum(); thc = np.arange(90) + 0.5
def passes(nm, K=11):
    c, BT, BR = LIB[nm]; p = BT.sum(1); v = P0.copy(); V = []
    for k in range(K):
        V.append(v.copy()); v = Rv * (BR @ v)
    pb = p @ wn; return p, np.array(V), pb, pb - ((p - pb) ** 2 @ wn) / (1 - pb)
g2 = gs[1].subgridspec(1, 2, wspace=0.34)
ae = fig.add_subplot(g2[0]); p_h, V, pb_h, pm_h = passes('MLA hemi_total'); cm = plt.cm.Blues; rows_e = []
for k in range(11): rows_e.append([k, float(V[k].sum())] + (V[k] / V[k].sum() * 100).tolist())
for k, shade in ((0, 0.3), (1, 0.45), (2, 0.6), (4, 0.78), (10, 1.0)):
    ae.plot(thc, V[k] / V[k].sum() * 100, color=cm(shade), lw=1.3, label=rf'$k$ = {k}')
ae.legend(loc='upper left', fontsize=7, handlelength=1.3)
ae.set_xlim(0, 90); ae.set_ylim(0, 5.8); ae.set_xticks([0, 30, 60, 90]); ae.set_xlabel(r'Incident angle $\theta_i$ (°)'); ae.set_ylabel('Share of incident power (% per degree)')
FS.label(ae, 'e')
af = fig.add_subplot(g2[1]); rows_f = []
for nm, lab, xo in (('scatter g0.8_t1.0', 'Scattering layer', 12.0), ('W_slope1.0', 'Wavy surface', 12.55), ('A_shape1.0', 'Random MLA', 11.45), ('MLA hemi_total', 'Hemisph. MLA', 12.0), ('planar', 'Planar', 12.0)):
    p_, V_, pb_, pm_ = passes(nm); pk = (V_ @ p_) / V_.sum(1); col, mk, _ = CLS[LIB[nm][0]]
    af.plot(range(11), pk / pb_, '-', color=col, lw=1.1, marker=mk, ms=3.6, mec='white', mew=0.3, label=lab)
    af.plot(xo, pm_ / pb_, mk, color=col, mfc='white', ms=4.4, mew=1.0)
    rows_f.append([nm, pb_, pm_] + pk.tolist())
af.plot([-0.6, 10.8], [1, 1], color='k', lw=0.8, ls=':', label='Fully randomized'); af.plot([10.8, 10.8], [0, 1.22], color='k', lw=0.6)
af.text(12.0, 1.08, PMIX + r'$/p$', fontsize=7, ha='center', va='bottom')
af.set_xlim(-0.6, 13.2); af.set_ylim(0, 1.72); af.set_xticks(range(0, 11, 2)); af.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]); af.set_xlabel(r'Pass $k$'); af.set_ylabel(r'Escape probability at pass $k$ / $p$')
af.legend(fontsize=6.2, loc='upper left', ncol=2, labelspacing=0.25, handlelength=1.3, columnspacing=0.6, handletextpad=0.35, borderaxespad=0.3); FS.label(af, 'f')
FS.save(fig, f'{OUT}/fig1_mixing_concept'); plt.close(fig)
save_xlsx('fig1_mixing_concept', sh + [('e_pass_resolved', ['pass k', 'total incident power (first pass = 1)'] + [f'share (%/deg) at {t:.1f} deg' for t in thc], rows_e), ('f_escape_per_pass', ['structure', 'p', 'p_mix'] + [f'p_k (k={k})' for k in range(11)], rows_f)])
print('Fig 1f:', [(r_[0], round(r_[1], 3), round(r_[2], 3), [round(x, 3) for x in r_[3:]]) for r_ in rows_f])
print('Fig 1e total incident power by pass:', [round(r_[1], 3) for r_ in rows_e])

# ------------------------------------------------------------------ Fig. 2 (where the light is mixed)
fig2 = plt.figure(figsize=(FW, 6.7)); g = fig2.add_gridspec(3, 2, wspace=0.36, hspace=0.42)
A, B, C, LG, D, E = [fig2.add_subplot(g[i, j]) for i in range(3) for j in range(2)]
pr = []
for nm, (c, BT, BR) in LIB.items():
    pT = BT.sum(1); w_ = W1 / W1.sum(); p = float(pT @ w_); var = float(((pT - p) ** 2) @ w_); pm = p - var / (1 - p)
    col, mk, _ = CLS[c]; A.plot(p, pm, mk, color=col, ms=4.8, mec='white', mew=0.4); pr.append([nm, c, p, np.sqrt(var), pm])
A.plot([0, 0.5], [0, 0.5], color='k', lw=0.7, ls=':'); A.text(0.262, 0.292, PMIX + r'$ = p$', fontsize=7, rotation=41, ha='center', va='center')
A.set_xlim(0.17, 0.45); A.set_ylim(0, 0.45); A.set_xlabel(r'Escape probability $p$'); A.set_ylabel(r'$p_{\mathrm{mix}} = p - \sigma_p^2/(1-p)$'); FS.label(A, 'a')
PLC = json.load(open('data/placement.json')); rows_c = []
for nm, (c, BT, BR) in LIB.items():
    v = PLC['structures'][f'{nm}|Al|550']; x = v['p'] - v['p_mix']; y = (v['dev_mix'] - v['exit_mix']) * 100
    col, mk, _ = CLS[c]; B.plot(x, y, mk, color=col, ms=4.8, mec='white', mew=0.4)
    rows_c.append([nm, c, v['p'], v['p_mix'], x, v['eta'], v['exit_mix'], v['dev_mix'], v['rand'], y])
rc = np.corrcoef([r[4] for r in rows_c], [r[9] for r in rows_c])[0, 1]; print(f'Fig 2b correlation r = {rc:.3f}')
B.set_xlim(-0.015, 0.37); B.set_ylim(-1.5, 30); B.set_xlabel(r'$p - p_{\mathrm{mix}} = \sigma_p^2/(1-p)$'); B.set_ylabel(r'$\eta_{\mathrm{sta}}$: device $-$ interface mixing (pp)'); FS.label(B, 'b')
B.legend(handles=FS.handles(CLS), loc='upper left', fontsize=6.6, borderaxespad=0.4, labelspacing=0.3, handletextpad=0.3)
rows_d = []
for R, ls in ((0.90, '-'), (0.98, '--')):
    sl = [z for z in PLC['slab'] if z['R'] == R]; nn = [z['n'] for z in sl]
    C.plot(nn, [z['none'] for z in sl], ls, color='#7d7c76'); C.plot(nn, [z['exit_mix'] for z in sl], ls, color='#eb6834'); C.plot(nn, [z['dev_mix'] for z in sl], ls, color='#2a78d6')
    rows_d += [[z['n'], z['R'], z['p'], z['p_mix'], z['none'], z['exit_mix'], z['dev_mix']] for z in sl]
HI = json.load(open('data/highindex.json')); rows_d2 = []; TEX = (('hemiMLA', '^', 'Hemispherical-MLA face', -0.05), ('wavy0.4', 'D', 'Wavy face, slope 0.4', 0.0), ('wavy1.0', 's', 'Wavy face, slope 1.0', 0.05))
for nm, mk, lab, dx in TEX:
    hz = [z for z in HI if z['structure'] == nm and z['R'] == 0.90]
    C.plot([z['n'] + dx for z in hz], [z['none'] for z in hz], mk, color='k', mfc='white', ms=4.2, mew=0.9)
    rows_d2 += [[z['structure'], z['n'], z['R'], z['p'], z['p_mix'], z['none'], z['exit_mix'], z['rand']] for z in HI if z['structure'] == nm]
C.set_xlabel(r'Refractive index $n$ of the slab'); C.set_ylabel('Extraction efficiency'); C.set_ylim(0, 1); C.set_xlim(1.4, 3.6); FS.label(C, 'c')
SV = json.load(open('data/scatter_volume.json')); rows_sv = []
for g_, col in ((0.0, '#9fd8c0'), (0.5, '#5cc29b'), (0.8, '#1baf7a'), (0.9, '#0d7a53')):
    pts = sorted((r_['tau'], r_['Al']['eta']) for r_ in SV['slabs'] if r_['g'] == g_)
    D.plot(*zip(*pts), 'o-', color=col, ms=3.4, lw=1.1, label=rf'$g$ = {g_:g}'); rows_sv += [[g_, t_, e_] for t_, e_ in pts]
RF = SV['refs']
for val, ls, col, lab, below in ((RF['planar']['Al']['rand'], ':', '#2a78d6', 'Ideal transmissive mixing', False), (RF['planar']['Al']['mix'], '-', '#eb6834', 'Reflective mixing at the interface', False), (RF['planar']['Al']['eta'], '-', '#7d7c76', 'Planar, no mixing', False)):
    D.axhline(val, ls=ls, color=col, lw=1.0); D.text(4.3, val + 0.008, lab, fontsize=6.6, color=col, ha='right', va='bottom')
D.set_xscale('log'); D.set_xlim(0.08, 4.8); D.set_ylim(0.33, 0.93); D.set_xticks([0.1, 0.25, 0.5, 1, 2, 4]); D.set_xticklabels(['0.1', '0.25', '0.5', '1', '2', '4']); D.minorticks_off()
D.set_xlabel(r'Reduced optical thickness $\tau^*$'); D.set_ylabel(ETA + ' (Al, 550 nm)'); D.legend(fontsize=6.3, loc='upper center', ncol=4, columnspacing=0.7, handlelength=0.9, handletextpad=0.3, borderaxespad=0.35); FS.label(D, 'd')
hd = []; OFF = []
for nm, (c, BT, BR) in LIB.items():
    e = EV[(nm, 'Al', 550.0)]; Nm = min(1 / max(1 - e['lam2'], 1e-9), 25); col, mk, _ = CLS[c]
    g_mix = (e['m1'] - e['eta']) * 100; g_mem = (e['closed'] - e['eta']) * 100; hd.append([nm, c, Nm, e['eta'], e['m1'], e['closed'], g_mix, g_mem])
    if c == 'planar': continue
    YTOP = 7.2; gm_plot = min(g_mem, YTOP - 0.25)
    E.plot([Nm, Nm], [g_mix, gm_plot], color=col, lw=0.5, alpha=0.6)
    E.plot(Nm, g_mix, mk, color=col, ms=4.6, mec='white', mew=0.4)
    if g_mem > YTOP:
        E.plot(Nm, YTOP - 0.25, marker=r'$\uparrow$', color=col, ms=7, mew=0); OFF.append((Nm, g_mem, col))
    else:
        E.plot(Nm, g_mem, mk, mfc='none', mec=col, ms=4.6, mew=0.8)
pl = [h for h in hd if h[1] == 'planar'][0]; print(f'Fig 2e planar (off scale): reflective mixing +{pl[6]:.1f} pp, fully randomized +{pl[7]:.1f} pp')
E.axhline(0, color='k', lw=0.6); E.set_xscale('log'); E.set_xlim(1, 12); E.set_xticks([1, 2, 5, 10]); E.set_xticklabels(['1', '2', '5', '10']); E.minorticks_off()
for k_, (x_, v_, c_) in enumerate(sorted(OFF)):
    E.text(x_ * (0.93 if k_ % 2 else 1.07), 6.55, f'{v_:.0f}', fontsize=6.2, color=c_, va='center', ha='right' if k_ % 2 else 'left')
E.set_ylim(-1.9, 7.2); E.set_xlabel('Angular mixing number ' + NMIX); E.set_ylabel('Predicted gain in ' + ETA + ' (pp)')
E.legend(handles=[Line2D([], [], ls='none', marker='o', color='k', ms=4.6, label='Reflective mixing'), Line2D([], [], ls='none', marker='o', mfc='none', color='k', ms=4.6, mew=0.8, label='Fully randomized')], fontsize=6.4, loc='lower right', ncol=2, columnspacing=0.8, handletextpad=0.2, borderaxespad=0.25); FS.label(E, 'e')
LG.axis('off')
LG.legend(handles=[Line2D([], [], color='#7d7c76', label='No mixing'), Line2D([], [], color='#eb6834', label='Complete mixing at the exit face'), Line2D([], [], color='#2a78d6', label='Complete mixing at the back\nreflector or in transmission'),
                   Line2D([], [], color='k', ls='-', label=r'$R$ = 0.90'), Line2D([], [], color='k', ls='--', label=r'$R$ = 0.98')] + [Line2D([], [], ls='none', marker=mk, color='k', mfc='white', ms=4.2, mew=0.9, label=lab + r' ($R$ = 0.90)') for _, mk, lab, _ in TEX],
          loc='center left', bbox_to_anchor=(-0.12, 0.5), title='Curves and symbols in (c)', fontsize=6.6, alignment='left', labelspacing=0.4)
FS.save(fig2, f'{OUT}/fig2_complete_mixing'); plt.close(fig2)
save_xlsx('fig2_complete_mixing', [('a_p_mix', ['structure', 'class', 'p', 'sigma_p', 'p_mix'], pr),
                                   ('d_gain', ['structure', 'class', 'N_mix (capped 25)', 'eta_full', 'eta_complete_reciprocal_mixing', 'eta_fully_randomized', 'gain_mixing_pp', 'gain_fully_randomized_pp'], hd),
                                   ('d_HG_volume', ['g', 'tau_star', 'eta_Al_550'], rows_sv), ('b_placement', ['structure', 'class', 'p', 'p_mix', 'p - p_mix', 'eta_actual', 'eta_escape_surface_mixing', 'eta_device_side_mixing', 'eta_fully_randomized', 'device - escape surface (pp)'], rows_c),
                                   ('c_slab', ['n', 'R', 'p', 'p_mix', 'eta_no_mixing', 'eta_exit_face_mixing', 'eta_back_reflector_mixing (= fully randomized)'], rows_d),
                                   ('c_textured', ['exit face', 'n', 'R', 'p', 'p_mix', 'eta_no_mixing', 'eta_exit_face_mixing', 'eta_back_reflector_mixing (= fully randomized)'], rows_d2)])

# ------------------------------------------------------------------ Fig. 3 (reduced model) and Fig. S2 (simpler models)
rows = []
for (nm, dev, lam), e in EV.items():
    rows.append([nm, LIB[nm][0], dev, lam, e['eta'], e['closed'], e['m0'], e['red_w'], e['red'], e['p'], e['lam2'], e['pinf']])
def _scatter(a, j, eq=True):
    for r in rows:
        col, mk, _ = CLS[r[1]]; a.plot(r[4], r[j], mk, color=col, ms=3.2, mew=0, alpha=0.85)
    err = np.array([abs(r[j] - r[4]) for r in rows]) * 100
    a.plot([0.3, 1], [0.3, 1], color='k', lw=0.7); a.set_xlim(0.3, 1); a.set_ylim(0.3, 1); a.set_xticks([0.4, 0.6, 0.8, 1.0]); a.set_yticks([0.4, 0.6, 0.8, 1.0])
    if eq: a.set_aspect('equal')
    a.text(0.96, 0.05, f'Mean |error| {err.mean():.2f} pp\nMax. {err.max():.1f} pp', transform=a.transAxes, fontsize=6.8, ha='right', va='bottom')
    a.set_xlabel(ETA + ', full series'); return err
figS, axS = plt.subplots(1, 3, figsize=(FW, 2.25), gridspec_kw=dict(wspace=0.12), sharey=True)
for k, j in enumerate((5, 6, 7)):
    _scatter(axS[k], j); FS.label(axS[k], 'abc'[k])
axS[0].set_ylabel(ETA + ', model'); axS[2].legend(handles=FS.handles(CLS), fontsize=5.6, loc='upper left', labelspacing=0.2, handletextpad=0.2, borderaxespad=0.2)
FS.save(figS, f'{OUT}/figS_model_comparison'); plt.close(figS)
fig = plt.figure(figsize=(FW, 4.9)); gs2 = fig.add_gridspec(2, 2, width_ratios=[1, 1.3], wspace=0.36, hspace=0.42)
a0 = fig.add_subplot(gs2[0, 0]); _scatter(a0, 8, eq=False); a0.set_ylabel(ETA + ', reduced model'); FS.label(a0, 'a')
ae = fig.add_subplot(gs2[0, 1]); xt = []; DCOL = {'Al': '#000000', 'Ag': '#000000', 'LL15': '#000000'}; DLAB = {'Al': 'Al', 'Ag': 'Ag', 'LL15': 'LL'}
for gi, (dv, (_, lam)) in enumerate([(d, w) for d in DEVS for w in WLS]):
    for r in rows:
        if r[2] != dv or r[3] != lam: continue
        col, mk, _ = CLS[r[1]]; xj = gi + 0.36 * (np.random.default_rng(hash(r[0]) % 2**32).random() - 0.5)
        ae.plot(xj, (r[8] - r[4]) * 100, mk, color=col, ms=2.8, mew=0, alpha=0.85)
    er = np.array([(r[8] - r[4]) * 100 for r in rows if r[2] == dv and r[3] == lam])
    ae.plot([gi - 0.32, gi + 0.32], [np.mean(np.abs(er))] * 2, color=DCOL[dv], lw=1.4)
    xt.append(f'{DLAB[dv]}\n{lam:.0f}')
ae.axhline(0, color='k', lw=0.6); ae.set_xticks(range(len(xt))); ae.set_xticklabels(xt, fontsize=6.6); ae.set_xlim(-0.6, len(xt) - 0.4); ae.tick_params(axis='x', which='both', top=False, length=0)
for b in (2.5, 5.5): ae.axvline(b, color='#b0b0b0', lw=0.6)
ae.set_ylabel('Reduced model $-$ full series (pp)'); ae.set_xlabel('Device and wavelength (nm)'); FS.label(ae, 'b')
ac = fig.add_subplot(gs2[1, 0]); CM = json.load(open('data/cutoff_modes.json'))['modes']; rows_c = []; PCOL = '#b8860b'
ac.axvspan(0.35, 0.62, color='#ececec', lw=0); ac.text(0.485, 24.0, 'Gap > 0.35', fontsize=6.6, ha='center', va='top', color='#4d4d4d')
for nm, m in CM.items():
    if nm == 'planar': continue
    blind = nm.startswith('blind:'); c = LIB[nm][0] if not blind else None
    col, mk = (CLS[c][0], CLS[c][1]) if not blind else ('k', 'o')
    gap = m['lam2'] - m['lam_eff']; red = m['lam2'] > 0.999
    ac.plot(gap, m['err_lam2'][1], 'X' if red else mk, color=col, ms=5 if red else 4, mfc=col if not blind or red else 'white', mew=0.8 if blind else 0, alpha=0.9, clip_on=False, zorder=3)
    rows_c.append([nm, 'out-of-sample' if blind else c, m['lam2'], m['lam_eff'], gap, m['w_slow'], m['err_lam2'][0], m['err_lam2'][1], bool(red)])
PR = json.load(open('data/prospective_results.json'))
for nm, v in PR.items():
    red = v['lam2'] > 0.999
    ac.plot(v['gap'], v['err_max'], 'X' if red else '^', color=PCOL, ms=5 if red else 4.2, mfc=PCOL if red else 'white', mew=0.9, clip_on=False, zorder=3)
    rows_c.append([nm, 'prospective', v['lam2'], v['lam_eff'], v['gap'], v['w_slow'], v['err_mean'], v['err_max'], bool(red)])
ac.set_xlim(0, 0.62); ac.set_ylim(0, 25); ac.set_xlabel(r'$\lambda_2 - \lambda_{\mathrm{eff}}$'); ac.set_ylabel('Maximum |error| (pp)'); FS.label(ac, 'c')
lg = fig.add_subplot(gs2[1, 1]); lg.axis('off')
l1 = lg.legend(handles=FS.handles(CLS), loc='upper left', bbox_to_anchor=(0.02, 1.02), title='Library (38 structures)', fontsize=6.8, alignment='left'); lg.add_artist(l1)
lg.legend(handles=[Line2D([], [], ls='none', marker='o', color='k', mfc='white', ms=4, mew=0.8, label='Out-of-sample (9)'), Line2D([], [], ls='none', marker='^', color=PCOL, mfc='white', ms=4.2, mew=0.9, label='Prospective (12)'),
                   Line2D([], [], ls='none', marker='X', color='k', ms=5, label=r'Reducible chain ($\lambda_2$ = 1)'), Line2D([], [], color='k', lw=1.4, label='Mean |error| in (b)')],
          loc='lower left', bbox_to_anchor=(0.02, -0.12), title='Further symbols', fontsize=6.8, alignment='left')
save_xlsx('fig3c_applicability', [('structures', ['structure', 'class', 'lambda2', 'lambda_eff', 'gap', 'slow_mode_share', 'err_mean_pp', 'err_max_pp', 'closed_set'], rows_c)])
FS.save(fig, f'{OUT}/fig3_reduced_model'); plt.close(fig)
save_xlsx('fig3_reduced_model', [('cases', ['structure', 'class', 'device', 'lambda_nm', 'eta_full', 'eta_fully_randomized', 'eta_m0', 'eta_reduced_w', 'eta_reduced_recip', 'p_Lamb', 'lambda2_85', 'p_inf'], rows)])
print(f'Fig 3: {len(rows)} cases; reduced(w) mean {np.mean([abs(r[7]-r[4]) for r in rows])*100:.2f} max {np.max([abs(r[7]-r[4]) for r in rows])*100:.2f}; recip mean {np.mean([abs(r[8]-r[4]) for r in rows])*100:.2f} pp, max {np.max([abs(r[8]-r[4]) for r in rows])*100:.2f}; memoryless mean {np.mean([abs(r[5]-r[4]) for r in rows])*100:.2f}, max {np.max([abs(r[5]-r[4]) for r in rows])*100:.1f}')

# ------------------------------------------------------------------ Fig. 4 (structure map)
fig = plt.figure(figsize=(FW, 4.9)); g4 = fig.add_gridspec(2, 2, wspace=0.36, hspace=0.42); ax = [fig.add_subplot(g4[0, 0]), fig.add_subplot(g4[0, 1]), fig.add_subplot(g4[1, 0])]; lg = fig.add_subplot(g4[1, 1]); lg.axis('off'); rows = []
P0 = DV['Al_P'][75].reshape(90, 2).sum(1); R1 = DV['Al_R'][75].reshape(90, 2).mean(1)
pT0 = LIB['MLA hemi_total'][1].sum(1)
pg = np.linspace(0.18, 0.47, 30); Ng = np.logspace(0, 1.4, 30); Z = np.zeros((30, 30))
for i, p in enumerate(pg):
    pT = pT0 * p / (pT0 @ (W1 / W1.sum()))
    for j, Nm in enumerate(Ng):
        Z[j, i] = mixing.reduced(np.clip(pT, 0, 1), W1, R1, P0, 1 / Nm)
cs = ax[0].contour(pg, Ng, Z, levels=[0.6, 0.7, 0.75, 0.78, 0.8], colors='#6a6a6a', linewidths=0.6)
for lv, segs in zip(cs.levels, cs.allsegs):
    for sg in segs[:1]:
        if len(sg): ax[0].text(sg[np.argmin(np.abs(sg[:, 1] - 14)), 0], 14, f'{lv:.2f}', fontsize=6.3, color='#4d4d4d', ha='center', va='center', rotation=84, bbox=dict(fc='white', ec='none', pad=0.4))
for nm, (c, BT, BR) in LIB.items():
    e = EV[(nm, 'Al', 550.0)]; Nm = 1 / max(1 - e['lam2'], 1e-9)
    col, mk, _ = CLS[c]
    if c != 'planar': ax[0].plot(e['p'], Nm, mk, color=col, ms=4.6, mec='white', mew=0.4, zorder=3)
    rows.append([nm, c, e['p'], e['lam2'], Nm, e['pinf'], e['pinf'] / e['p85'], e['eta'], EV[(nm, 'Ag', 550.0)]['eta'], EV[(nm, 'LL15', 550.0)]['eta']])
    ax[2].plot(min(Nm, 25), e['pinf'] / e['p85'], mk, color=col, ms=4.6, mec='white', mew=0.4, zorder=3, clip_on=False)
ax[0].set_yscale('log'); ax[0].set_ylim(1, 22); ax[0].set_xlim(0.17, 0.47); ax[0].set_yticks([1, 2, 5, 10, 20]); ax[0].set_yticklabels(['1', '2', '5', '10', '20']); ax[0].minorticks_off()
ax[0].set_xlabel(r'Escape probability $p$'); ax[0].set_ylabel('Angular mixing number ' + NMIX); FS.label(ax[0], 'a')
RLEG = []
for pre, lab, col, key in (('A_shape', 'Lenslet shape spread $s$', '#6250d6', 'shape'), ('C_jitter', 'Position jitter / pitch', '#2a78d6', 'jitter'),
                           ('W_slope', 'Wavy surface, rms slope', '#e87ba4', 'slope'), ('B_aspect', 'Lenslet height spread', '#eda100', 'aspect')):
    pts = sorted((float(nm.split(key)[1].split('_')[0]), 1 / max(1 - EV[(nm, 'Al', 550.0)]['lam2'], 1e-9)) for nm in LIB if nm.startswith(pre) and not nm.startswith('C2'))
    if pre == 'C_jitter':
        pts = [(0.0, 1 / (1 - EV[('A_shape0.0', 'Al', 550.0)]['lam2']))] + [q for q in pts if q[0] > 0]
    if pre == 'B_aspect':
        pts = [(0.0, 1 / (1 - EV[('A_shape0.0', 'Al', 550.0)]['lam2']))] + pts
    if pts:
        x, y = zip(*pts); xm = max(x); ax[1].plot(np.array(x) / xm, y, 'o-', color=col, ms=3.6); RLEG.append(Line2D([], [], color=col, marker='o', ms=3.6, label=f'{lab} (0–{xm:g})'))
ns = 1 / (1 - EV[('scatter g0.8_t1.0', 'Al', 550.0)]['lam2']); nh = 1 / (1 - EV[('MLA hemi_total', 'Al', 550.0)]['lam2'])
ax[1].axhline(ns, color='#1baf7a', ls=':', lw=1.0); ax[1].axhline(nh, color='#2a78d6', ls=':', lw=1.0)
ax[1].set_xlim(-0.04, 1.04); ax[1].set_yscale('log'); ax[1].set_ylim(1, 9); ax[1].set_yticks([1, 2, 3, 5, 8]); ax[1].set_yticklabels(['1', '2', '3', '5', '8']); ax[1].minorticks_off()
ax[1].set_xlabel('Normalized randomness parameter'); ax[1].set_ylabel(NMIX); FS.label(ax[1], 'b')
ax[2].set_xscale('log'); ax[2].set_xlim(1, 32); ax[2].set_ylim(-0.04, 1.08); ax[2].set_xticks([1, 2, 5, 10, 25]); ax[2].set_xticklabels(['1', '2', '5', '10', r'$\infty$']); ax[2].minorticks_off()
ax[2].axhline(1, color='k', lw=0.6, ls=':'); ax[2].set_ylabel(r'Late-pass escape ratio $\xi$'); ax[2].set_xlabel('Angular mixing number ' + NMIX); FS.label(ax[2], 'c')
l1 = lg.legend(handles=FS.handles(CLS), loc='upper left', bbox_to_anchor=(0.02, 1.04), title='Structures in (a), (c)', fontsize=6.6, alignment='left', labelspacing=0.3); lg.add_artist(l1)
lg.legend(handles=RLEG + [Line2D([], [], color='#2a78d6', ls=':', lw=1.0, label='Hemispherical MLA'), Line2D([], [], color='#1baf7a', ls=':', lw=1.0, label=r'Scattering layer ($g$ = 0.8, $\tau^*$ = 1)')], loc='lower left', bbox_to_anchor=(0.02, -0.16), title='Randomization in (b)', fontsize=6.6, alignment='left', labelspacing=0.3)
FS.save(fig, f'{OUT}/fig4_structure_map'); plt.close(fig)
save_xlsx('fig4_structure_map', [('structures', ['structure', 'class', 'p_Lamb', 'lambda2_85', 'N_mix', 'p_inf', 'xi', 'eta_Al', 'eta_Ag', 'eta_LL15'], rows),
                                 ('contours', ['N_mix \\ p'] + list(pg), [[Ng[j]] + list(Z[j]) for j in range(30)])])

# ------------------------------------------------------------------ Fig. S4 (radiance bound and angular composition)
fig, ax = plt.subplots(1, 2, figsize=(FW, 2.5), gridspec_kw=dict(wspace=0.36)); rows = []
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
ax[0].bar(x - 0.3, UB, 0.15, color='#d9d8d3', label='Uniform-radiance bound'); ax[0].bar(x - 0.15, GB, 0.15, color='#7d7c76', label='Radiance-ordered bound')
for k, (nm, lab, col) in enumerate((('MLA hemi_total', 'Hemispherical MLA', '#2a78d6'), ('W_slope1.0', 'Wavy surface', '#e87ba4'), ('scatter g0.8_t1.0', 'Scattering layer', '#1baf7a'), ('planar', 'Planar', '#b5b4ad'))):
    if nm not in LIB: continue
    first = EV[(nm, 'Al', 550.0)]['first']; v = [first[a_:b_].sum() for a_, b_ in BANDS]
    ax[0].bar(x + 0.15 * k, v, 0.15, color=col, label=lab); rows.append([lab] + v + [first.sum()])
ax[0].set_xticks(x + 0.075); ax[0].set_xticklabels(['0–20°', '20–40°', '40–60°', '60–80°']); ax[0].tick_params(axis='x', which='both', top=False, length=0); ax[0].set_ylabel('First-pass power into band')
ax[0].set_ylim(0, 0.38); ax[0].legend(fontsize=5.6, ncol=2, loc='upper left', columnspacing=0.6, handlelength=0.8, handletextpad=0.4, borderaxespad=0.3); FS.label(ax[0], 'a')
rows = [['uniform bound'] + UB + [1 / N_SUB ** 2], ['radiance-ordered bound'] + GB + [gbound(P05, 0, 90)]] + rows
win = []
for nm, (c, BT, BR) in LIB.items():
    e = EV[(nm, 'Al', 550.0)]; col, mk, _ = CLS[c]
    ax[1].plot(e['eta'], e['S'][2], mk, color=col, ms=4.6, mec='white', mew=0.4); win.append([nm, c, e['eta']] + e['S'])
ax[1].axhspan(0.33, 0.36, color='#eef4fb', zorder=0); ax[1].axhline(0.337, color='k', ls=':', lw=0.7)
ax[1].set_ylim(0.28, 0.40); ax[1].set_xlabel(ETA + ' (Al, 550 nm)'); ax[1].set_ylabel('Fraction of the output in 40–60°')
ax[1].legend(handles=FS.handles(CLS), fontsize=5.8, loc='lower left', ncol=2, columnspacing=0.6, handletextpad=0.2); FS.label(ax[1], 'b')
FS.save(fig, f'{OUT}/fig5_bound_window'); plt.close(fig)
save_xlsx('fig5_bound_window', [('a_first_pass', ['curve', '0-20', '20-40', '40-60', '60-80', 'total'], rows),
                                ('b_composition', ['structure', 'class', 'eta_sta', 'S0-20', 'S20-40', 'S40-60', 'S60-80'], win)])
S = [w[5] for w in win]
print(f'Fig S4: S40-60 range {min(S):.3f}-{max(S):.3f} over {len(S)} structures; bound total {gbound(P05, 0, 90):.3f}')
json.dump({f'{k[0]}|{k[1]}|{k[2]}': {kk: vv for kk, vv in v.items() if kk != 'first'} for k, v in EV.items()}, open('data/oe_eval.json', 'w'))
