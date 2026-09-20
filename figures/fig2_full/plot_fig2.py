# -*- coding: utf-8 -*-
"""Figure 2, composed: parasitic absorption is the binding loss.

(a) round-trip loss against ITO thickness, split into the mirror and TCO parts, and the
    resulting extraction efficiency
(b) the cost of a high-index outcoupling structure: p falls, the two eta_ext curves peel apart
(c) the result: substrate-delivered power keeps rising, the shaded gap never escapes
(d) round-trip loss of the candidate electrode structures, resolved in angle and wavelength

Sized for a full journal width (180 mm).  Colour is the reflector throughout:
vermillion = Al, blue = Ag, grey = the TCO part of the loss."""
import numpy as np, matplotlib, os, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'sim', 'design_rule4'))
import fig2d as F

C_AL, C_AG, C_TCO = '#D55E00', '#0072B2', '#9C9C9C'
C_AG2, C_THIN, C_DBR = '#56B4E9', '#009E73', '#9467BD'
DASH = (0, (3.2, 1.8))
FS_LAB, FS_TICK, FS_NOTE, FS_LET = 6.6, 6.0, 5.7, 8.5
plt.rcParams.update({'font.size': FS_LAB, 'axes.linewidth': 0.7,
                     'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
                     'xtick.major.size': 2.6, 'ytick.major.size': 2.6})

A = os.path.join(ROOT, 'figures', 'fig2a_mock')
B = os.path.join(ROOT, 'figures', 'fig2b_mock')
def a_load(m):
    T = np.loadtxt(f'{A}/kn2_{m}_15.csv', delimiter=','); Z = np.loadtxt(f'{A}/kn3_{m}_15_k0.csv', delimiter=',')
    return T[:, 0], Z[:, 2], T[:, 2] - Z[:, 2], T[:, 8]
al, ag = a_load('al'), a_load('ag')
D = np.genfromtxt(f'{B}/fig2b_curves_konig.csv', delimiter=',', names=True)
nb, pb = D['n_sub'], D['p']
eAl, eAg, sAl, sAg, qAl, qAg = (D['eta_ext_Al'], D['eta_ext_Ag'], D['eta_sub_Al'],
                               D['eta_sub_Ag'], D['EQE_Al'], D['EQE_Ag'])

fig = plt.figure(figsize=(7.2, 5.0))
R1, H1 = 0.585, 0.340
R2, H2 = 0.100, 0.340
w1, gap = 0.186, 0.052
xs = [0.075 + i*(w1 + gap) for i in range(4)]
ax_a1 = fig.add_axes([xs[0], R1, w1, H1])
ax_a2 = fig.add_axes([xs[1], R1, w1, H1])
ax_b  = fig.add_axes([xs[2], R1, w1, H1])
ax_c  = fig.add_axes([xs[3], R1, w1, H1])
ax_d1 = fig.add_axes([0.075, R2, 0.130, H2])
ax_d2 = fig.add_axes([0.215, R2, 0.130, H2])
ax_cb = fig.add_axes([0.355, R2, 0.008, H2])
ax_d3 = fig.add_axes([0.450, R2, 0.245, H2])
ax_leg = fig.add_axes([0.722, R2, 0.270, H2]); ax_leg.axis('off')

def clean(ax):
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=FS_TICK, pad=1.5)
    ax.grid(axis='y', color='0.93', lw=0.6)
def letter(ax, t, dx=-0.30):
    ax.text(dx, 1.16, t, transform=ax.transAxes, fontsize=FS_LET, fontweight='bold', va='top')

# ---------------- (a) left: the round-trip loss, split --------------------
DS = [30, 50, 100, 150, 200]
idx = [int(np.argmin(abs(al[0] - t))) for t in DS]
x = np.arange(len(DS)); bw = 0.36
for off, (dd, Am, At, _), c, tag in ((-0.5*bw, al, C_AL, 'Al'), (0.5*bw, ag, C_AG, 'Ag')):
    m = np.array([Am[i] for i in idx]); t = np.array([At[i] for i in idx])
    ax_a1.bar(x + off, m, bw, color=c, lw=0)
    ax_a1.bar(x + off, t, bw, bottom=m, color=C_TCO, lw=0)
ax_a1.set_xticks(x); ax_a1.set_xticklabels([str(t) for t in DS])
ax_a1.set_xlim(-0.55, len(DS) - 0.45); ax_a1.set_ylim(0, 0.19)
ax_a1.set_yticks([0, 0.05, 0.10, 0.15])
ax_a1.set_xlabel('ITO thickness (nm)'); ax_a1.set_ylabel("round-trip loss  A′")
ax_a1.legend(handles=[Patch(facecolor=C_AL, label='mirror, Al'),
                      Patch(facecolor=C_AG, label='mirror, Ag'),
                      Patch(facecolor=C_TCO, label='TCO')],
             loc='upper left', fontsize=FS_NOTE, frameon=False, handlelength=1.1,
             handleheight=0.9, labelspacing=0.28, borderpad=0.1)
clean(ax_a1); letter(ax_a1, 'a')

# ---------------- (a) right: extraction efficiency ------------------------
for (dd, _, _, ext), c in ((al, C_AL), (ag, C_AG)):
    ax_a2.plot(dd, ext, color=c, lw=1.8, solid_capstyle='round')
ax_a2.plot(150, 0.916, '*', ms=8, mfc=C_AL, mec='w', mew=0.8, zorder=6)
ax_a2.annotate('measured\n(green device)', (150, 0.916), xytext=(-5, -24),
               textcoords='offset points', ha='right', fontsize=FS_NOTE, color='0.4',
               arrowprops=dict(arrowstyle='->', color='0.6', lw=0.7))
ax_a2.text(196, ag[3][-1] + 0.012, 'Ag', color=C_AG, fontsize=FS_LAB, fontweight='bold', ha='right')
ax_a2.text(196, al[3][-1] + 0.012, 'Al', color=C_AL, fontsize=FS_LAB, fontweight='bold', ha='right')
ax_a2.set_xlim(30, 200); ax_a2.set_ylim(0.70, 1.0)
ax_a2.set_xlabel('ITO thickness (nm)'); ax_a2.set_ylabel('$\\eta_{ext}$')
clean(ax_a2)

# ---------------- (b) ------------------------------------------------------
ax_b.fill_between(nb, eAl, eAg, color=C_AG, alpha=0.10, lw=0)
ax_b.plot(nb, pb, color='0.55', lw=1.1, ls=DASH)
ax_b.plot(nb, eAl, color=C_AL, lw=1.8); ax_b.plot(nb, eAg, color=C_AG, lw=1.8)
ax_b.annotate('', (1.965, eAl[-2]), (1.965, eAg[-2]),
              arrowprops=dict(arrowstyle='<->', color='0.45', lw=0.7, shrinkA=0, shrinkB=0))
ax_b.text(1.945, 0.5*(eAl[-2] + eAg[-2]), '21 %p', ha='right', va='center',
          fontsize=FS_NOTE, color='0.35')
ax_b.text(1.70, eAg[8] + 0.02, 'Ag', color=C_AG, fontsize=FS_LAB, fontweight='bold', ha='center')
ax_b.text(1.70, eAl[8] - 0.055, 'Al', color=C_AL, fontsize=FS_LAB, fontweight='bold', ha='center')
ax_b.text(1.88, pb[12] + 0.03, '$p$', color='0.45', fontsize=FS_LAB, ha='center')
ax_b.set_xlim(1.3, 2.0); ax_b.set_ylim(0, 1.0)
ax_b.set_xlabel('$n_{sub}$'); ax_b.set_ylabel('$\\eta_{ext}$,   $p$')
clean(ax_b); letter(ax_b, 'b')

# ---------------- (c) ------------------------------------------------------
ax_c.fill_between(nb, qAl, sAl, color=C_AL, alpha=0.16, lw=0)
ax_c.fill_between(nb, qAg, sAg, color=C_AG, alpha=0.16, lw=0)
ax_c.plot(nb, sAl, color=C_AL, lw=0.8, alpha=0.7); ax_c.plot(nb, sAg, color=C_AG, lw=0.8, alpha=0.7)
ax_c.plot(nb, qAl, color=C_AL, lw=1.9); ax_c.plot(nb, qAg, color=C_AG, lw=1.9)
for q, c in ((qAg, C_AG), (qAl, C_AL)):
    i = int(np.argmax(q))
    ax_c.plot(nb[i], q[i], 'o', ms=3.6, mfc=c, mec='w', mew=0.8, zorder=5)
    ax_c.annotate(f'{q[i]:.2f}', (nb[i], q[i]), xytext=(5, -9), textcoords='offset points',
                  fontsize=FS_NOTE, fontweight='bold', color=c)
ax_c.axvline(1.8, color='0.85', lw=0.7, ls=':', zorder=1)
ax_c.annotate('never escapes', (1.62, 0.5*(sAl[6] + qAl[6])), xytext=(1.315, 0.05),
              textcoords='data', ha='left', va='bottom', fontsize=FS_NOTE, color=C_AL,
              arrowprops=dict(arrowstyle='->', color=C_AL, lw=0.7, connectionstyle='arc3,rad=-0.2'))
ax_c.text(1.83, 0.20, '$n_{sub}$=$n_{EML}$', fontsize=FS_NOTE, color='0.45', ha='left')
ax_c.set_xlim(1.3, 2.0); ax_c.set_ylim(0, 1.0)
ax_c.set_xlabel('$n_{sub}$'); ax_c.set_ylabel('$\\eta_{sub}^{(0)}$,   EQE')
clean(ax_c); letter(ax_c, 'c')

# ---------------- (d) maps and angle plot ---------------------------------
VIS = np.arange(430.0, 701.0)
NSUB = 1.5
MAPPED = [('Al / ITO 150 nm', ax_d1), ('Ag / IZO 50 nm', ax_d2)]
for name, ax in MAPPED:
    th, lam, R, T = F.maps(name, n_sub=NSUB, lam=VIS)
    im = ax.pcolormesh(np.degrees(th), lam, 100*(1 - R), cmap='inferno_r', vmin=0, vmax=30,
                       shading='auto', rasterized=True)
    ax.set_title(name, fontsize=FS_NOTE + 0.4, pad=2.5)
    ax.set_xlabel('$\\theta$ in substrate (°)')
    ax.set_xticks([0, 30, 60, 90]); ax.tick_params(labelsize=FS_TICK, pad=1.5)
ax_d1.set_ylabel('wavelength (nm)'); ax_d2.set_yticklabels([])
ax_d1.set_xticklabels(['0', '30', '60', ''])
letter(ax_d1, 'd')
cb = fig.colorbar(im, cax=ax_cb); cb.set_label('1 − R  (%)', fontsize=FS_LAB, labelpad=2)
cb.ax.tick_params(labelsize=FS_TICK, pad=1.5); cb.outline.set_linewidth(0.6)

STRUCT = [('Al / ITO 150 nm', C_AL), ('DBR / IZO 50 nm', C_DBR), ('Ag / Ag 10 nm', C_THIN),
          ('Ag / ITO 150 nm', C_AG2), ('Ag / IZO 50 nm', C_AG)]
th = np.linspace(0, np.pi/2, 452); th = 0.5*(th[1:] + th[:-1])
w = np.cos(th)*np.sin(th); w /= w.max()
ax_d3.fill_between(np.degrees(th), 0, 7*w, color='0.91', lw=0, zorder=0)
ax_d3.text(46, 0.4, '$\\cos\\theta\\sin\\theta$ weight', fontsize=FS_NOTE, color='0.5', ha='center')
W = {}
for name, c in STRUCT:
    t2, lam, R, T = F.maps(name, n_sub=NSUB, lam=VIS)
    sw = np.clip(np.interp(lam, F.LAM, F.GREEN), 0, None); sw /= sw.sum()
    ax_d3.plot(np.degrees(t2), 100*(1 - (R*sw[:, None]).sum(0)), color=c, lw=1.5)
    W[name] = F.weighted(name, n_sub=NSUB, spectrum=F.GREEN, lam=VIS)
ax_d3.axvline(np.degrees(np.arcsin(1/NSUB)), color=C_DBR, lw=0.7, ls=':', zorder=1)
ax_d3.text(88, 31.2, 'the dielectric mirror leaks only\ninside the escape cone to air',
           fontsize=FS_NOTE, color=C_DBR, ha='right', va='top')
ax_d3.set_xlim(0, 90); ax_d3.set_ylim(0, 32); ax_d3.set_xticks([0, 30, 60, 90])
ax_d3.set_xlabel('$\\theta$ in substrate (°)'); ax_d3.set_ylabel('round-trip loss 1 − R  (%)')
clean(ax_d3)

ax_leg.set_xlim(0, 1); ax_leg.set_ylim(0, 1)
ax_leg.text(0.0, 1.02, 'flux- and spectrum-weighted\nround-trip loss', fontsize=FS_NOTE + 0.4,
            va='top', color='0.25')
for i, (name, c) in enumerate(STRUCT):
    y = 0.70 - i*0.108
    ax_leg.plot([0.0, 0.10], [y, y], color=c, lw=1.8, clip_on=False)
    lab = f"{name}   {100*W[name]['loss']:.1f} %"
    ax_leg.text(0.14, y, lab, fontsize=FS_NOTE, va='center', color='0.2')
ax_leg.text(0.14, 0.70 - 1*0.108 - 0.068, '(3.0 absorbed + 8.5 leaked)',
            fontsize=FS_NOTE - 0.6, va='center', color=C_DBR)
ax_leg.text(0.0, 0.02, 'generic stack:  substrate / electrode /\n420 nm organics (n = 1.8) / reflector,\nglass substrate, green emitter',
            fontsize=FS_NOTE - 0.6, va='bottom', color='0.45')

fig.savefig('fig2_full.png', dpi=320)
fig.savefig('fig2_full.pdf')
print('saved', {k: round(100*v['loss'], 1) for k, v in W.items()})
