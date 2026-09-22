# -*- coding: utf-8 -*-
"""Fig. 1: (a) full-width ray schematics, lossy vs low-loss, with the 3D renders as small insets;
(b) round-trip decay and (c) master curve side by side below, vertical axes aligned."""
import numpy as np, matplotlib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
C_LOSSY, C_LOW = '#6E6E6E', '#C00000'
FS, FSN = 7.4, 6.4
plt.rcParams.update({'font.size': FS, 'axes.linewidth': 0.8})
fig = plt.figure(figsize=(7.2, 5.0))
gs = fig.add_gridspec(2, 2, height_ratios=[0.62, 1.0], left=0.10, right=0.985, top=0.97, bottom=0.095, hspace=0.35, wspace=0.28)

# ---------- (a) ----------
axa = fig.add_subplot(gs[0, :]); axa.axis('off')
lossy = np.asarray(Image.open(os.path.join(HERE, 'rays_lossy.png')).convert('RGBA'))
low = np.asarray(Image.open(os.path.join(HERE, 'rays_lowloss.png')).convert('RGBA'))
lossy = lossy[:, :int(lossy.shape[1] * 0.745)]           # drop the layer labels on the right
render = np.asarray(Image.open(os.path.join(HERE, 'render_both.png')).convert('RGBA'))
h, w = render.shape[:2]
ren_l, ren_r = render[:, :w // 2], render[:, w // 2:]
axa.set_xlim(0, 1); axa.set_ylim(0, 1)
for k, (img, title, col) in enumerate(((lossy, 'conventional OLED: large round-trip loss', C_LOSSY),
                                         (low, 'this work: small round-trip loss', C_LOW))):
    x0 = 0.0 if k == 0 else 0.52
    ia = axa.inset_axes([x0 + 0.125, 0.02, 0.37, 0.86]); ia.imshow(img); ia.axis('off')
    ir = axa.inset_axes([x0 - 0.005, 0.12, 0.135, 0.66]); ir.imshow(ren_l if k == 0 else ren_r); ir.axis('off')
    axa.text(x0 + 0.245, 1.02, title, ha='center', va='top', fontsize=FS, color=col, fontweight='bold')
# legend
leg = [Line2D([], [], color='#1A9E7E', lw=1.6, label='light path'),
       Line2D([], [], color='#E0472B', lw=1.2, label='Joule loss at the reflector  /  absorption in the TCO')]
axa.legend(handles=leg, loc='lower center', bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False, fontsize=FSN, handlelength=1.6)
axa.text(-0.075, 1.04, '(a)', transform=axa.transAxes, fontsize=9, fontweight='bold', va='top')

# ---------- (b) ----------
p = 0.4
n = np.arange(0, 11)
axb = fig.add_subplot(gs[1, 0])
for A, col, ls in ((0.02, C_LOW, '-'), (0.10, '#B0B0B0', '-'), (0.25, C_LOSSY, '-')):
    q = (1 - p) * (1 - A)
    cum = p * (1 - q ** (n + 1)) / (1 - q)
    axb.plot(n, cum, color=col, lw=1.7, ls=ls, marker='o', ms=3, mew=0, label="A′ = %.2f  →  %.2f" % (A, p / (p + (1 - p) * A)))
axb.set_xlim(-0.3, 10.3); axb.set_ylim(0.2, 1.0)
axb.set_xlabel('number of round trips'); axb.set_ylabel('cumulative escaped fraction')
axb.text(0.03, 0.97, '$p$ = 0.4', transform=axb.transAxes, va='top', color='0.3')
axb.legend(loc='lower right', frameon=False, fontsize=FSN, title='round-trip loss  →  η$_{\\rm ext}$', title_fontsize=FSN)
axb.text(-0.22, 1.0, '(b)', transform=axb.transAxes, fontsize=9, fontweight='bold', va='top')

# ---------- (c) ----------
axc = fig.add_subplot(gs[1, 1], sharey=axb)
R = np.linspace(0, 1, 400)
for pp, col, ls, lab in ((0.6, '0.55', ':', '$p$ = 0.6'), (0.4, C_LOW, '-', '$p$ = 0.4'), (0.25, '0.55', '--', '$p$ = 0.25')):
    A = 1 - R
    axc.plot(R, pp / (pp + (1 - pp) * A), color=col, lw=1.6 if pp == 0.4 else 1.1, ls=ls)
    axc.text(0.04, pp + 0.03, lab, color=col, fontsize=FSN)
for A, col, lab, dy in ((0.25, C_LOSSY, 'large loss (A′ = 0.25)', -0.13), (0.02, C_LOW, 'low-loss reflector (A′ = 0.02)', 0.045)):
    y = p / (p + (1 - p) * A)
    axc.plot([1 - A], [y], marker='o', ms=6, color=col, mew=0, zorder=5)
    axc.annotate(lab, (1 - A, y), xytext=(1 - A - 0.02, y + dy), ha='right', va='center', fontsize=FSN, color=col,
                 arrowprops=dict(arrowstyle='-', lw=0.6, color=col))
axc.set_xlim(0, 1.0); axc.set_xlabel('round-trip reflectance  $R_{\\rm LED}$ = 1 − $A$′')
axc.set_ylabel('extraction efficiency  η$_{\\rm ext}$')
axc.text(-0.22, 1.0, '(c)', transform=axc.transAxes, fontsize=9, fontweight='bold', va='top')
for ax in (axb, axc):
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
axc.tick_params(labelleft=True)
for e in ('png', 'pdf'):
    fig.savefig(os.path.join(HERE, 'fig1.' + e), dpi=400)
print('fig1 written')
