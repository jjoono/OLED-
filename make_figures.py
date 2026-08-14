"""
make_figures.py  --  manuscript figure set (Fig. 1-5)

Regenerates every figure in the manuscript from the archived result files, so
that the numbers in the captions and the numbers in the panels cannot drift
apart. Run from the repository root:

    python3 make_figures.py

Inputs
  pareto_front_result.mat        weighted-sum sweep, 691 evaluations   -> Fig 2b, Fig 3
  opt_4band_result_25by25.mat    convex per-band campaign              -> Fig 2a, Fig 2c, Fig 5
  opt_hemisphere_result.mat      hemispherical reference               -> Fig 2c
  opt_4band_inverted_result.mat  inverted (concave) family             -> Fig 5
  stress_random_result.mat       randomly assembled family             -> Fig 5
  warmstart_hemisphere_result.mat  hemisphere-seeded control run       -> Fig 2e,f
  angular_recycling_result.npz   Markov recycling / DBR angular data   -> Fig 4
  angular_recycling_bandwidth.npz  DBR bandwidth sweep                 -> Fig 4

EVAL_LOG column layout (shared by every campaign):
  1..13  design vector | 14 EQE_total | 15..18 band EQE (0-20/20-40/40-60/60-80)
  19 phase | 20 weight

Outputs: fig1_platform.png/.pdf ... fig5_families.png/.pdf, plus a console
summary of every number quoted in the captions.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Wedge, Rectangle, Polygon
from scipy.io import loadmat
from scipy.interpolate import PchipInterpolator

plt.rcParams.update({
    'font.size': 8.5,
    'axes.titlesize': 9,
    'axes.labelsize': 8.5,
    'legend.fontsize': 7.2,
    'xtick.labelsize': 7.5,
    'ytick.labelsize': 7.5,
    'axes.linewidth': 0.8,
    'figure.dpi': 110,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

BANDS = ['0-20', '20-40', '40-60', '60-80']
BAND_TEX = [r'0-20$^\circ$', r'20-40$^\circ$', r'40-60$^\circ$', r'60-80$^\circ$']
S_LAMB = np.array([0.117, 0.296, 0.337, 0.220])
CB = ['#3B6EA5', '#C6512F', '#4F8F5B', '#8A5FA8']   # per-band colour
N_TOP = 20

SUMMARY = []


def note(s):
    SUMMARY.append(s)
    print(s)


def save(fig, stem):
    fig.savefig(stem + '.png')
    fig.savefig(stem + '.pdf')
    plt.close(fig)
    print(f'  saved -> {stem}.png / .pdf')


def clean_log(L):
    """Rows with a usable total EQE and four finite band values."""
    T, B = L[:, 13], L[:, 14:18]
    ok = np.isfinite(T) & (T > 0) & np.isfinite(B).all(axis=1)
    return T[ok], B[ok, :]


def natural_composition(T, B, n_top=N_TOP, pct=(10, 90)):
    """Median / percentile band selectivity of the n_top most efficient designs.

    This is the internal no-steering reference: the composition a design ends up
    with when it is optimized for efficiency rather than for any one band.
    """
    top = np.argsort(T)[::-1][:min(n_top, T.size)]
    S = B[top, :] / T[top, None]
    return np.median(S, axis=0), np.percentile(S, pct[0], axis=0), np.percentile(S, pct[1], axis=0)


def band_corr(T, B):
    """Pearson R between total EQE and each band selectivity."""
    return np.array([np.corrcoef(T, B[:, j] / T)[0, 1] for j in range(4)])


def profile(xv, n=400):
    """Lens profile from a 13-variable design vector.

    Control points are (0,1), (x2..x6, y2..y6), (1,0) in units of the lens
    radius; the height axis is then scaled by stretch_Z. Returned as
    (radial coordinate, height), both normalized to the lens radius.
    """
    px = np.concatenate(([0.0], np.sort(xv[0:5]), [1.0]))
    py = np.concatenate(([1.0], xv[5:10], [0.0]))
    px, uniq = np.unique(px, return_index=True)
    py = py[uniq]
    t = np.linspace(0, 1, n)
    return t, PchipInterpolator(px, py)(t) * xv[12]


# ============================================================
#  load
# ============================================================
P = loadmat('pareto_front_result.mat')
C = loadmat('opt_4band_result_25by25.mat')
H = loadmat('opt_hemisphere_result.mat')
I = loadmat('opt_4band_inverted_result.mat')
Rn = loadmat('stress_random_result.mat')
W  = loadmat('warmstart_hemisphere_result.mat')
d1 = np.load('angular_recycling_result.npz')
d3 = np.load('angular_recycling_bandwidth.npz')

Tp, Bp = clean_log(P['EVAL_LOG'])
phase_p = P['EVAL_LOG'][:, 18]
phase_p = phase_p[np.isfinite(P['EVAL_LOG'][:, 13]) & (P['EVAL_LOG'][:, 13] > 0)
                  & np.isfinite(P['EVAL_LOG'][:, 14:18]).all(axis=1)]
wgt_p = P['EVAL_LOG'][:, 19]
wgt_p = wgt_p[np.isfinite(P['EVAL_LOG'][:, 13]) & (P['EVAL_LOG'][:, 13] > 0)
              & np.isfinite(P['EVAL_LOG'][:, 14:18]).all(axis=1)]

Tc, Bc = clean_log(C['EVAL_LOG'])
band_eqe = C['band_eqe_hi'].ravel()
band_tot = C['band_tot_hi'].ravel()
band_x = C['band_x']

hemi_val = H['hemi_val'].ravel()
hemi_tot = H['hemi_tot'].ravel()
hemi_x = H['hemi_x']
HEMI_X, HEMI_Y = H['HEMI_X'].ravel(), H['HEMI_Y'].ravel()

# E_star is the best *high-precision* total EQE (50,000 rays, 151 lambda, x3),
# i.e. the best of the five weighted-sum optima -- not the maximum of the
# coarse-fidelity search log, which reaches 0.5591 at 10,000 rays.
E_star = float(np.max(P['pareto_tot'].ravel()))
R_full = band_corr(Tp, Bp)
S_nat, S_lo, S_hi = natural_composition(Tc, Bc)

note('=' * 62)
note(f'best high-precision total EQE                E* = {E_star:.4f}')
note(f'coarse search log maximum                         {np.max(Tp):.4f}')
note(f'total EQE range                              {Tp.min():.4f} - {Tp.max():.4f}')
note('full-population selectivity-efficiency R     ' +
     ' '.join(f'{r:+.2f}' for r in R_full))
note('natural composition (convex, top-20 median)  ' +
     ' '.join(f'{s:.3f}' for s in S_nat))
note('=' * 62)


# ============================================================
#  Fig. 1 | Platform and benchmark design
# ============================================================
def fig1():
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.26)

    # ---- (a) architecture ----
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(-0.05, 1.20); ax.set_ylim(-0.16, 1.10); ax.axis('off')
    layers = [(0.030, 0.052, '#8E8E93', 'Al cathode'),
              (0.082, 0.052, '#B8CBE0', 'ETL  (10-150 nm)'),
              (0.134, 0.052, '#E8A33D', 'EML  (dipole)'),
              (0.186, 0.052, '#B8CBE0', 'HTL  (10-150 nm)'),
              (0.238, 0.052, '#9FD1C7', 'ITO')]
    for y0, h, col, lab in layers:
        ax.add_patch(Rectangle((0.10, y0), 0.58, h, fc=col, ec='k', lw=0.5))
        ax.annotate(lab, xy=(0.68, y0 + h / 2), xytext=(0.73, y0 + h / 2),
                    va='center', fontsize=6.2,
                    arrowprops=dict(arrowstyle='-', lw=0.4, color='0.5'))
    ax.add_patch(Rectangle((0.02, 0.300), 0.72, 0.36, fc='#DCEAF5', ec='k', lw=0.6))
    ax.text(0.06, 0.47, 'glass substrate\n$n$ = 1.51,  1.295 mm', fontsize=6.6, va='center')

    # lenslets
    for i in range(14):
        xc = 0.033 + i * 0.0517
        ax.add_patch(Wedge((xc, 0.660), 0.0258, 0, 180, fc='#F2D9A8', ec='k', lw=0.4))
    ax.annotate('', xy=(0.02, 0.745), xytext=(0.74, 0.745),
                arrowprops=dict(arrowstyle='<->', lw=0.8))
    ax.text(0.38, 0.765, 'textured patch  25 $\\times$ 25 mm', ha='center', fontsize=6.6)
    ax.text(0.76, 0.665, 'MLA\n$r\\approx$10 $\\mu$m', fontsize=6.2, va='center')

    # emission
    for ang in (-38, -14, 14, 38):
        r = np.radians(ang)
        ax.annotate('', xy=(0.39 + 0.36 * np.sin(r), 0.165 + 0.36 * np.cos(r)),
                    xytext=(0.39, 0.165),
                    arrowprops=dict(arrowstyle='->', lw=0.8, color='#C6512F'))
    ax.annotate('', xy=(0.10, -0.012), xytext=(0.68, -0.012),
                arrowprops=dict(arrowstyle='<->', lw=0.8))
    ax.text(0.39, -0.085, 'OLED disc,  $r_{\\rm OLED}$ = 1 mm', ha='center', fontsize=6.6)
    ax.set_title('(a) architecture', loc='left')

    # ---- (b) lens classes ----
    ax = fig.add_subplot(gs[0, 1])
    xs = np.linspace(0, 1, 200)
    ax.plot(np.concatenate((-xs[::-1], xs)),
            np.concatenate((np.sqrt(1 - xs[::-1] ** 2), np.sqrt(1 - xs ** 2))),
            color='#3B6EA5', lw=2, label='hemisphere  (3 free vars)')
    t, h = profile(band_x[2, :])
    ax.plot(np.concatenate((-t[::-1], t)), np.concatenate((h[::-1], h)),
            color='#C6512F', lw=2, label='freeform  (13 free vars)')
    cx = np.sort(band_x[2, 0:5]); cy = band_x[2, 5:10] * band_x[2, 12]
    ax.plot(cx, cy, 'o', ms=3.6, mfc='w', mec='#C6512F', mew=1.0)
    ax.plot(-cx, cy, 'o', ms=3.6, mfc='w', mec='#C6512F', mew=1.0)
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xlim(-1.25, 1.25); ax.set_ylim(-0.12, 2.55)
    ax.set_xlabel('radial coordinate  $r / r_{\\rm lens}$')
    ax.set_ylabel('height  $z / r_{\\rm lens}$')
    ax.legend(loc='upper center', framealpha=0.9)
    ax.text(0.0, 0.10, 'open circles: spline control points;\nendpoints fixed at (0,1) and (1,0)',
            fontsize=6.2, va='bottom', ha='center')
    ax.set_title('(b) lens classes', loc='left')
    ax.grid(alpha=0.25)

    # ---- (c) angular bands ----
    ax = fig.add_subplot(gs[1, 0])
    ax.set_xlim(-1.15, 1.15); ax.set_ylim(-0.22, 1.18); ax.axis('off')
    for j, (lo, hi) in enumerate([(0, 20), (20, 40), (40, 60), (60, 80)]):
        ax.add_patch(Wedge((0, 0), 1.0, 90 - hi, 90 - lo, fc=CB[j], alpha=0.55, ec='w', lw=0.8))
        ax.add_patch(Wedge((0, 0), 1.0, 90 + lo, 90 + hi, fc=CB[j], alpha=0.55, ec='w', lw=0.8))
        am = np.radians(90 - (lo + hi) / 2)
        ax.text(1.06 * np.cos(am), 1.06 * np.sin(am), BANDS[j], fontsize=6.6,
                ha='left', va='center', color=CB[j])
    ax.plot([-1.15, 1.15], [0, 0], 'k-', lw=0.9)
    ax.plot([0, 0], [0, 1.12], 'k--', lw=0.7)
    ax.text(0.02, 1.13, r'$\theta = 0$', fontsize=6.8)
    ax.add_patch(Rectangle((-0.30, -0.10), 0.60, 0.10, fc='#F2D9A8', ec='k', lw=0.5))
    ax.text(0, -0.175, 'emitter + MLA (far field over full azimuth)',
            ha='center', fontsize=6.8)
    ax.text(-1.12, 0.92, r'$S_j = {\rm EQE}_{{\rm band},j}\,/\,{\rm EQE}_{\rm total}$',
            fontsize=7.4, va='top')
    ax.set_title('(c) polar bands', loc='left')

    # ---- (d) workflow ----
    ax = fig.add_subplot(gs[1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
    steps = [('CPS microcavity dipole model\n$I_{\\rm sub}(\\theta,\\lambda)$, 453-753 nm', '#DCEAF5'),
             ('LightTools 3D ray trace\n(MATLAB COM)', '#E7F0E2'),
             ('surrogateopt global search\n10,000 rays / 31 $\\lambda$', '#FBEBD7'),
             ('patternsearch polish', '#FBEBD7'),
             ('high-ray re-evaluation\n50,000 rays / 151 $\\lambda$ / $\\times$3', '#F5DCDC')]
    yb, hgt, gap = 0.93, 0.155, 0.035
    for i, (txt, col) in enumerate(steps):
        y = yb - i * (hgt + gap) - hgt
        ax.add_patch(FancyBboxPatch((0.06, y), 0.88, hgt,
                                    boxstyle='round,pad=0.012,rounding_size=0.02',
                                    fc=col, ec='k', lw=0.6))
        ax.text(0.50, y + hgt / 2, txt, ha='center', va='center', fontsize=6.6)
        if i:
            ax.add_patch(FancyArrowPatch((0.50, y + hgt + gap), (0.50, y + hgt),
                                         arrowstyle='-|>', mutation_scale=8, lw=0.8, color='k'))
    ax.set_title('(d) workflow', loc='left')

    save(fig, 'fig1_platform')


# ============================================================
#  Fig. 2 | Achievable region and Pareto collapse
# ============================================================
def fig2():
    fig = plt.figure(figsize=(7.4, 5.6))
    gs = fig.add_gridspec(2, 6, hspace=0.48, wspace=1.25)

    S_win = band_eqe / band_tot
    E_max_c = float(np.max(Tc))
    gain_sel = S_win / S_nat
    gain_tot = band_tot / E_max_c
    gain_net = band_eqe / (S_nat * E_max_c)

    note('band-dedicated optima (convex, patch 25):')
    note(f'{"band":>8} {"S_win":>7} {"S_nat":>7} {"sel gain":>9} {"tot ratio":>10} {"net":>6}')
    for j in range(4):
        note(f'{BANDS[j]:>8} {S_win[j]:7.3f} {S_nat[j]:7.3f} '
             f'{gain_sel[j]:9.3f} {gain_tot[j]:10.3f} {gain_net[j]:6.3f}')

    # ---- (a) selectivity vs natural spread ----
    ax = fig.add_subplot(gs[0, 0:2])
    xx = np.arange(1, 5)
    for j in range(4):
        ax.add_patch(Rectangle((xx[j] - 0.30, S_lo[j]), 0.60, S_hi[j] - S_lo[j],
                               fc='0.82', ec='0.45', lw=0.7, zorder=1))
        ax.plot([xx[j] - 0.36, xx[j] + 0.36], [S_LAMB[j]] * 2, 'k--', lw=1.0, zorder=2)
    ax.plot(xx, S_win, 'o', ms=6.5, mfc='#C6512F', mec='k', mew=0.7, ls='none', zorder=3,
            label='band-dedicated optimum')
    ax.set_xticks(xx); ax.set_xticklabels(BANDS, fontsize=7.0)
    ax.set_xlim(0.5, 4.5); ax.set_ylim(0, 0.50)
    ax.set_ylabel(r'selectivity  $S_j$')
    ax.set_title('(a) dedicated vs natural', loc='left')
    ax.legend(loc='upper left', framealpha=0.92)
    ax.text(0.60, 0.415, 'grey: 10-90% of top-20\ndashed: Lambertian', fontsize=6.2, va='top')
    ax.grid(alpha=0.25)

    # ---- (b) gain decomposition ----
    ax = fig.add_subplot(gs[0, 2:4])
    w = 0.26
    ax.bar(xx - w, gain_sel, w, color='#3B6EA5', ec='k', lw=0.4, label='selectivity gain')
    ax.bar(xx, gain_tot, w, color='#9AA7B5', ec='k', lw=0.4, label='total-EQE ratio')
    ax.bar(xx + w, gain_net, w, color='#C6512F', ec='k', lw=0.4, label='net band gain')
    ax.axhline(1, color='k', lw=0.9)
    for j in range(4):
        ax.text(xx[j] + w, gain_net[j] + 0.012, f'{gain_net[j]:.2f}', ha='center', fontsize=6.0)
    ax.set_xticks(xx); ax.set_xticklabels(BANDS, fontsize=7.0)
    ax.set_xlim(0.5, 4.5); ax.set_ylim(0.82, 1.50)
    ax.set_ylabel('ratio')
    ax.set_title(r'(b) gain decomposition', loc='left')
    ax.legend(loc='upper right', framealpha=0.95, fontsize=6.6)
    ax.grid(alpha=0.25, axis='y')

    # ---- (c) Pareto collapse ----
    ax = fig.add_subplot(gs[0, 4:6])
    rnd = phase_p == 1
    ax.scatter(Tp[rnd], Bp[rnd, 2], s=9, c='0.65', alpha=0.55, lw=0,
               label=f'random feasible ($n$={rnd.sum()})')
    ax.scatter(Tp[~rnd], Bp[~rnd, 2], s=9, c='#4F8F5B', alpha=0.55, lw=0,
               label=f'weighted-sum search ($n$={(~rnd).sum()})')
    pf = np.polyfit(Tp, Bp[:, 2], 1)
    xs = np.linspace(Tp.min(), Tp.max(), 50)
    r2 = np.corrcoef(Tp, Bp[:, 2])[0, 1] ** 2
    ax.plot(xs, np.polyval(pf, xs), 'k-', lw=1.3,
            label=f'linear fit  $R^2$ = {r2:.3f}')
    ax.plot(P['pareto_tot'].ravel(), P['pareto_band'].ravel(), 'D', ms=6,
            mfc='#C6512F', mec='k', mew=0.7, ls='none',
            label='weighted optima (5 weights)')
    ax.set_xlabel(r'EQE$_{\rm total}$')
    ax.set_ylabel(r'EQE$_{40-60^\circ}$')
    ax.set_title('(c) a collapse, not a front', loc='left')
    ax.set_ylim(0.02, 0.245)
    ax.legend(loc='lower right', framealpha=1.0, fontsize=5.9, handlelength=1.2,
              borderpad=0.35, labelspacing=0.3)
    ax.grid(alpha=0.25)
    note(f'Fig 2b linear collapse R^2 = {r2:.3f}, '
         f'weighted optima total EQE spread = '
         f'{P["pareto_tot"].ravel().min():.4f}-{P["pareto_tot"].ravel().max():.4f}')

    # ---- (d) best profiles ----
    ax = fig.add_subplot(gs[1, 0:2])
    hx = np.concatenate(([0.0], HEMI_X, [1.0]))
    hy = np.concatenate(([1.0], HEMI_Y, [0.0]))
    th = np.linspace(0, 1, 300)
    k_ht = int(np.argmax(hemi_tot))
    ax.plot(th, PchipInterpolator(hx, hy)(th) * hemi_x[k_ht, 2], 'k-', lw=2.2,
            label=f'hemisphere ref. ($\\times${hemi_x[k_ht, 2]:.2f})')
    for j in range(4):
        t, h = profile(band_x[j, :])
        ax.plot(t, h, lw=1.5, color=CB[j], label=f'{BAND_TEX[j]} optimum')
    ax.axhline(0, color='k', lw=0.7)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.1, 4.1)
    ax.set_xlabel('radial coordinate  $r / r_{\\rm lens}$')
    ax.set_ylabel('height  $z / r_{\\rm lens}$')
    ax.set_title('(d) best freeform profiles', loc='left')
    ax.legend(loc='upper right', ncol=1, framealpha=0.92, fontsize=5.8)
    ax.grid(alpha=0.25)

    # ---- (e) freeform vs hemisphere ----
    #  The freeform entry is the best over BOTH freeform campaigns: the original
    #  per-band search and the hemisphere-seeded control run. Taking only the
    #  first would report a search artifact (two bands came out below the
    #  hemisphere there purely because a 13-variable search from random seeds
    #  had not converged) rather than what the design class can reach.
    ax = fig.add_subplot(gs[1, 2:4])
    ws_val = W['ws_val'].ravel()
    ff_orig = np.concatenate((band_eqe, [E_star]))
    ff = np.fmax(ff_orig, np.nan_to_num(ws_val, nan=-np.inf))
    hm = hemi_val[:5]
    G = ff / hm
    note('freeform best per objective (orig / warm-start / adopted):')
    for j, nm in enumerate(BANDS + ['total']):
        note(f'  {nm:>6} {ff_orig[j]:.5f} / '
             f'{ws_val[j] if np.isfinite(ws_val[j]) else float("nan"):.5f} -> {ff[j]:.5f}')
    xx5 = np.arange(1, 6)
    ax.bar(xx5 - 0.19, ff, 0.38, color='#C6512F', ec='k', lw=0.5, label='freeform (13 vars)')
    ax.bar(xx5 + 0.19, hm, 0.38, color='#3B6EA5', ec='k', lw=0.5, label='hemisphere (3 vars)')
    for j in range(5):
        ax.text(xx5[j], max(ff[j], hm[j]) + 0.013, f'{G[j]:.3f}',
                ha='center', fontsize=6.4,
                color=('#2E6B3C' if G[j] >= 1 else '#9A3324'))
    ax.set_xticks(xx5); ax.set_xticklabels(BANDS + ['total'], fontsize=6.2)
    ax.set_ylim(0, 0.72)
    ax.set_ylabel('EQE at the matching objective')
    ax.set_title(r'(e) gain $G_j$ over hemisphere', loc='left')
    ax.legend(loc='upper left', framealpha=0.92, fontsize=6.2)
    ax.grid(alpha=0.25, axis='y')
    note('G_j (best freeform / hemisphere): ' + ' '.join(f'{g:.3f}' for g in G))

    # ---- (f) hemisphere-seeded control ----
    #  Answers the objection that G_j < 1 in the original campaign meant an
    #  under-converged search rather than an absent gain. Each arm's search is
    #  restarted AT that arm's hemisphere optimum. Where a gain exists the same
    #  procedure finds it at enormous significance; where the original campaign
    #  fell short, restarting recovers the deficit and then stops.
    ax = fig.add_subplot(gs[1, 4:6])
    base = W['base_val'].ravel()
    sd_b, sd_w = W['base_sd'].ravel(), W['ws_sd'].ravel()
    se = np.sqrt(sd_b ** 2 + sd_w ** 2) / np.sqrt(3)
    tval = (ws_val - base) / se
    rel = 100 * (ws_val - base) / base
    T_CRIT = 2.132
    ok = np.isfinite(rel)
    idx = np.arange(5)[ok]
    cols = ['#4F8F5B' if tval[j] > T_CRIT else '#9AA7B5' for j in idx]
    ax.bar(np.arange(len(idx)), rel[idx], 0.6, color=cols, ec='k', lw=0.5)
    for i, j in enumerate(idx):
        ax.text(i, rel[j] + 0.16, f'$t$={tval[j]:.0f}' if tval[j] > 10 else f'$t$={tval[j]:.1f}',
                ha='center', fontsize=6.2)
    ax.axhline(0, color='k', lw=0.9)
    ax.set_xticks(np.arange(len(idx)))
    ax.set_xticklabels([(BANDS + ['total'])[j] for j in idx], fontsize=6.2)
    ax.set_ylim(-0.4, 6.2)
    ax.set_ylabel('gain over hemisphere (%)')
    ax.set_title('(f) restarted at the hemisphere', loc='left')
    ax.text(0.04, 0.97, 'green: exceeds noise ($t>2.13$)\ngrey: does not',
            transform=ax.transAxes, fontsize=5.8, va='top')
    ax.grid(alpha=0.25, axis='y')
    note('warm start vs hemisphere: ' +
         ' '.join(f'{(BANDS+["total"])[j]} {rel[j]:+.2f}% (t={tval[j]:.1f})' for j in idx))

    save(fig, 'fig2_achievable_region')


# ============================================================
#  Fig. 3 | Selectivity map across the design space
# ============================================================
def fig3():
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.0), sharex=True)
    for j, ax in enumerate(axs.ravel()):
        S = Bp[:, j] / Tp
        rnd = phase_p == 1
        ax.scatter(Tp[rnd], S[rnd], s=8, c='0.68', alpha=0.55, lw=0, label='random feasible')
        ax.scatter(Tp[~rnd], S[~rnd], s=8, c=CB[j], alpha=0.55, lw=0, label='optimized')
        pf = np.polyfit(Tp, S, 1)
        xs = np.linspace(Tp.min(), Tp.max(), 50)
        ax.plot(xs, np.polyval(pf, xs), 'k-', lw=1.2)
        ax.axhline(S_LAMB[j], color='r', ls='--', lw=1.0,
                   label=f'Lambertian {S_LAMB[j]:.3f}')
        ax.set_title(f'{BAND_TEX[j]}   $R$ = {R_full[j]:+.2f}', loc='left')
        ax.set_ylabel(r'$S_j$')
        if j >= 2:
            ax.set_xlabel(r'EQE$_{\rm total}$')
        ax.grid(alpha=0.25)
        if j == 0:
            ax.legend(loc='upper left', framealpha=0.92)
    fig.suptitle('band selectivity drifts with efficiency, but the band ordering never inverts',
                 fontsize=8.6, y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, 'fig3_selectivity_map')


# ============================================================
#  Fig. 4 | The recycling route and the design-route map
# ============================================================
def fig4():
    WALL = float(d1['wall']) * 100
    fig = plt.figure(figsize=(7.2, 5.2))
    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.30)

    # ---- (a) loss dependence of the Markov recycling model ----
    ax = fig.add_subplot(gs[0, 0])
    a_list = np.array([0.0, 0.01, 0.02, 0.05, 0.10, 0.20]) * 100
    ideal_b = np.array([100.0, 94.2, 89.1, 76.6, 62.1, 45.0])
    flat_b = np.array([33.8, 33.3, 32.8, 31.3, 29.1, 25.5])
    ax.plot(a_list, ideal_b, 'o-', lw=1.8, color='#3B6EA5', ms=4.5,
            label='ideal angular filter')
    ax.plot(a_list, flat_b, 's-', lw=1.8, color='#9AA7B5', ms=4.5,
            label='planar / MLA / scatterer')
    ax.axhline(WALL, color='#C6512F', ls='--', lw=1.2,
               label=f'single-pass wall  {WALL:.1f}%')
    ax.plot([10], [62.1], 'o', ms=9, mfc='none', mec='#3B6EA5', mew=1.4)
    ax.annotate('62.1% at $a$ = 10%', xy=(10, 62.1), xytext=(12.2, 68),
                fontsize=6.4, arrowprops=dict(arrowstyle='->', lw=0.7))
    ax.annotate('29.1%', xy=(10, 29.1), xytext=(13.0, 17), fontsize=6.4,
                arrowprops=dict(arrowstyle='->', lw=0.7))
    ax.set_xlabel('round-trip loss per bounce  $a$ (%)')
    ax.set_ylabel(r'into 40-60$^\circ$ band (% of generated)')
    ax.set_title('(a) recycling gain is capped by loss', loc='left')
    ax.set_ylim(5, 112)
    ax.legend(loc='lower left', framealpha=0.95, fontsize=6.4)
    ax.grid(alpha=0.25)

    # ---- (b) bandwidth dependence of a real DBR ----
    ax = fig.add_subplot(gs[0, 1])
    dl = d3['dlam']; bd = d3['band'] * 100; sel = d3['sel'] * 100
    ax.plot(dl, bd, 'o-', lw=1.8, color='#8A5FA8', ms=4.5, label='8-pair DBR, band delivery')
    ax.plot(dl, sel, '^--', lw=1.4, color='#8A5FA8', ms=4.5, alpha=0.6,
            label='DBR selectivity (of escaped)')
    ax.axhline(WALL, color='#C6512F', ls='--', lw=1.2, label=f'single-pass wall {WALL:.1f}%')
    ax.annotate(f'{bd[0]:.1f}% monochromatic', xy=(dl[0], bd[0]), xytext=(11, 54),
                fontsize=6.4, arrowprops=dict(arrowstyle='->', lw=0.7))
    ax.annotate(f'{bd[-1]:.1f}% at 100 nm', xy=(dl[-1], bd[-1]), xytext=(52, 24),
                fontsize=6.6, arrowprops=dict(arrowstyle='->', lw=0.7))
    ax.set_xlabel(r'source bandwidth  $\Delta\lambda$ (nm)')
    ax.set_ylabel('into band (% of generated)')
    ax.set_title('(b) a real filter needs a narrow source', loc='left')
    ax.set_ylim(0, 80)
    ax.legend(loc='upper right', framealpha=0.92)
    ax.grid(alpha=0.25)
    note(f'Fig 4b DBR band delivery: {bd[0]:.1f}% (mono) -> {bd[-1]:.1f}% (100 nm)')

    # ---- (c) single-pass angular response and the wall ----
    ax = fig.add_subplot(gs[1, 0])
    th = d1['th']
    lo, hi = float(d1['th_sub_lo']), float(d1['th_sub_hi'])
    ax.axvspan(lo, hi, color='0.87', label=r'40-60$^\circ$ in air')
    ax.plot(th, d1['T_flat'], lw=1.8, color='#9AA7B5', label='planar interface')
    ax.plot(th, d1['T_dbr'], lw=1.8, color='#8A5FA8', label='optimized 8-pair DBR')
    Ti = np.zeros_like(th); Ti[(th >= lo) & (th <= hi)] = 1
    ax.plot(th, Ti, 'k--', lw=1.1, label='ideal angular filter')
    ax.set_xlim(0, 50); ax.set_ylim(0, 1.06)
    ax.set_xlabel(r'$\theta$ inside the substrate (deg)')
    ax.set_ylabel('transmittance')
    ax.set_title('(c) single-pass angular response', loc='left')
    ax.text(1.5, 0.13, f'single-pass wall\n{WALL:.1f}% analytic / 33.8% numeric\n'
            r'(independent of $g$)', fontsize=6.4, va='bottom')
    ax.legend(loc='upper right', framealpha=0.95, fontsize=6.4)
    ax.grid(alpha=0.25)

    # ---- (d) design-route map ----
    ax = fig.add_subplot(gs[1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
    ax.add_patch(FancyBboxPatch((0.06, 0.845), 0.88, 0.125,
                                boxstyle='round,pad=0.012,rounding_size=0.02',
                                fc='#EDEDF2', ec='k', lw=0.7))
    ax.text(0.50, 0.907, 'hemispherical benchmark already reproduced\n'
                         '(shape freedom exhausted)',
            ha='center', va='center', fontsize=6.6)
    ax.add_patch(FancyArrowPatch((0.50, 0.845), (0.50, 0.795),
                                 arrowstyle='-|>', mutation_scale=8, lw=0.8, color='k'))
    routes = [('raise total extraction',
               'source / cavity engineering\n'
               r'($\eta_{\rm rad}$, dipole orientation, cavity order)', '#DCEAF5'),
              ('raise absolute band power',
               'aperture expansion\n(extraction area, not lens shape)', '#E7F0E2'),
              ('raise angular selectivity',
               'angular-selective recycling\n'
               r'(needs $a\lesssim$ 2%, $\Delta\lambda\lesssim$ 30 nm)', '#F5DCDC')]
    for i, (fom, lever, col) in enumerate(routes):
        y = 0.575 - i * 0.205
        ax.add_patch(FancyBboxPatch((0.06, y), 0.88, 0.165,
                                    boxstyle='round,pad=0.010,rounding_size=0.02',
                                    fc=col, ec='k', lw=0.6))
        ax.text(0.50, y + 0.126, 'if the target is to ' + fom + ':',
                ha='center', va='center', fontsize=6.3, style='italic', color='0.25')
        ax.text(0.50, y + 0.055, lever, ha='center', va='center', fontsize=6.5)
    ax.text(0.50, -0.005, 'lens-shape optimization is on none of these routes',
            ha='center', fontsize=6.5, style='italic', color='0.35')
    ax.set_title('(d) design-route map', loc='left')

    save(fig, 'fig4_recycling_routes')


# ============================================================
#  Fig. 5 | Generality across MLA families
# ============================================================
def fig5():
    # Every family is reduced by the SAME statistic from its OWN log: the
    # top-20 median composition and the own-population correlation. Mixing
    # statistics between families (e.g. a top-20 median for one and a
    # population mean for another) produces apparent family differences as
    # large as the real ones, so the population mean is shown alongside as an
    # open marker rather than substituted for the median.
    fams = [(C, 'convex freeform (reference)', '#C6512F'),
            (I, 'inverted (concave)', '#3B6EA5'),
            (Rn, 'randomly assembled', '#4F8F5B')]
    fig, axs = plt.subplots(3, 3, figsize=(7.2, 6.6))

    S_conv_ref, R_conv_ref = None, None
    for r, (D, name, col) in enumerate(fams):
        T, B = clean_log(D['EVAL_LOG'])
        Snat, _, _ = natural_composition(T, B)
        Smean = np.mean(B / T[:, None], axis=0)
        Rsel = band_corr(T, B)
        if r == 0:
            S_conv_ref, R_conv_ref = Snat.copy(), Rsel.copy()
        note(f'{name:<28} E*={T.max():.4f}  top20med=' + ' '.join(f'{s:.3f}' for s in Snat)
             + '  popmean=' + ' '.join(f'{s:.3f}' for s in Smean)
             + '  R=' + ' '.join(f'{v:+.2f}' for v in Rsel))

        # (1) collapse
        ax = axs[r, 0]
        ax.scatter(T, B[:, 2], s=10, c=col, alpha=0.5, lw=0)
        pf = np.polyfit(T, B[:, 2], 1)
        xs = np.linspace(T.min(), T.max(), 50)
        ax.plot(xs, np.polyval(pf, xs), 'k-', lw=1.2)
        r2 = np.corrcoef(T, B[:, 2])[0, 1] ** 2
        ax.set_xlabel(r'EQE$_{\rm total}$'); ax.set_ylabel(r'EQE$_{40-60^\circ}$')
        ax.set_title(f'({"abc"[r]}1) {name}\n$R^2$ = {r2:.3f},  $n$ = {T.size}', loc='left')
        ax.grid(alpha=0.25)

        # (2) natural composition
        ax = axs[r, 1]
        xx = np.arange(1, 5)
        ax.bar(xx - 0.26, Snat, 0.26, color=col, ec='k', lw=0.4, label='this family')
        ax.bar(xx, S_conv_ref, 0.26, color='0.62', ec='k', lw=0.4, label='convex ref.')
        ax.bar(xx + 0.26, S_LAMB, 0.26, color='none', ec='k', lw=0.8, ls='--',
               label='Lambertian')
        ax.plot(xx - 0.26, Smean, 'o', ms=3.8, mfc='w', mec='k', mew=0.8,
                ls='none', label='same family, population mean')
        ax.set_xticks(xx); ax.set_xticklabels(BANDS, fontsize=6.6)
        ax.set_ylim(0, 0.58)
        ax.set_ylabel(r'natural $S_j$')
        ax.set_title(f'({"abc"[r]}2) composition', loc='left')
        if r == 0:
            ax.legend(loc='upper center', fontsize=5.6, framealpha=0.95, ncol=2,
                      handlelength=1.1, columnspacing=0.8, borderpad=0.3)
        ax.grid(alpha=0.25, axis='y')

        # (3) drift
        ax = axs[r, 2]
        ax.bar(xx - 0.17, Rsel, 0.34, color=col, ec='k', lw=0.4, label='this family')
        ax.bar(xx + 0.17, R_conv_ref, 0.34, color='0.62', ec='k', lw=0.4, label='convex ref.')
        ax.axhline(0, color='k', lw=0.8)
        ax.set_xticks(xx); ax.set_xticklabels(BANDS, fontsize=6.6)
        ax.set_ylim(-1.05, 1.05)
        ax.set_ylabel(r'$R\,($EQE$_{\rm total},\,S_j)$')
        ax.set_title(f'({"abc"[r]}3) drift', loc='left')
        if r == 0:
            ax.legend(loc='lower left', fontsize=6.0, framealpha=0.92)
        ax.grid(alpha=0.25, axis='y')

    fig.tight_layout()
    save(fig, 'fig5_families')


if __name__ == '__main__':
    print('\n[Fig. 1]'); fig1()
    print('\n[Fig. 2]'); fig2()
    print('\n[Fig. 3]'); fig3()
    print('\n[Fig. 4]'); fig4()
    print('\n[Fig. 5]'); fig5()
    print('\n' + '=' * 62)
    print('all figures written.')
