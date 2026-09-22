# -*- coding: utf-8 -*-
"""Fig. 1 as one row: (a) the two ray schematics stacked, (b) decay, (c) master curve."""
import numpy as np, matplotlib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__))
C_LOSSY, C_LOW = '#6E6E6E', '#C00000'
FS, FSN = 7.0, 6.0
plt.rcParams.update({'font.size': FS, 'axes.linewidth': 0.8})
fig = plt.figure(figsize=(7.2, 2.45))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1], left=0.005, right=0.99, top=0.9, bottom=0.19, wspace=0.42)

axa = fig.add_subplot(gs[0, 0]); axa.axis('off'); axa.set_xlim(0, 1); axa.set_ylim(0, 1)
lossy = np.asarray(Image.open(os.path.join(HERE, 'rays_lossy.png')).convert('RGBA'))
lossy = lossy[:, :int(lossy.shape[1] * 0.745)]
low = np.asarray(Image.open(os.path.join(HERE, 'rays_lowloss.png')).convert('RGBA'))
for k, (img, title, col) in enumerate(((lossy, 'conventional: large round-trip loss', C_LOSSY),
                                         (low, 'this work: small round-trip loss', C_LOW))):
    y0 = 0.56 if k == 0 else 0.08
    ia = axa.inset_axes([0.06, y0, 0.92, 0.36]); ia.imshow(img); ia.axis('off')
    axa.text(0.52, y0 + 0.40, title, ha='center', va='bottom', fontsize=FSN + 0.4, color=col, fontweight='bold')
leg = [Line2D([], [], color='#1A9E7E', lw=1.4, label='light path'),
       Line2D([], [], color='#E0472B', lw=1.0, label='Joule loss / TCO absorption')]
axa.legend(handles=leg, loc='lower center', bbox_to_anchor=(0.52, -0.2), ncol=2, frameon=False, fontsize=FSN - 0.4, handlelength=1.4, columnspacing=1.2)
axa.text(0.0, 1.06, '(a)', transform=axa.transAxes, fontsize=8.5, fontweight='bold', va='bottom')

p = 0.4; n = np.arange(0, 11)
axb = fig.add_subplot(gs[0, 1])
for A, col in ((0.02, C_LOW), (0.10, '#B0B0B0'), (0.25, C_LOSSY)):
    q = (1 - p) * (1 - A)
    axb.plot(n, p * (1 - q ** (n + 1)) / (1 - q), color=col, lw=1.5, marker='o', ms=2.6, mew=0, label="A′ = %.2f → %.2f" % (A, p / (p + (1 - p) * A)))
axb.set_xlim(-0.3, 10.3); axb.set_ylim(0.2, 1.0)
axb.set_xlabel('number of round trips'); axb.set_ylabel('cumulative escaped fraction')
axb.text(0.04, 0.97, '$p$ = 0.4', transform=axb.transAxes, va='top', color='0.3', fontsize=FSN)
axb.legend(loc='lower right', frameon=False, fontsize=FSN - 0.4, title='round-trip loss → η$_{\\rm ext}$', title_fontsize=FSN - 0.4, handlelength=1.6)
axb.text(-0.3, 1.06, '(b)', transform=axb.transAxes, fontsize=8.5, fontweight='bold', va='bottom')

axc = fig.add_subplot(gs[0, 2], sharey=axb)
R = np.linspace(0, 1, 400)
for pp, col, ls, lab in ((0.6, '0.55', ':', '$p$ = 0.6'), (0.4, C_LOW, '-', '$p$ = 0.4'), (0.25, '0.55', '--', '$p$ = 0.25')):
    axc.plot(R, pp / (pp + (1 - pp) * (1 - R)), color=col, lw=1.5 if pp == 0.4 else 1.0, ls=ls)
    axc.text(0.04, pp + 0.03, lab, color=col, fontsize=FSN)
for A, col, lab, xy in ((0.25, C_LOSSY, 'large loss\n(A′ = 0.25)', (0.80, 0.60)), (0.02, C_LOW, 'low-loss reflector\n(A′ = 0.02)', (0.50, 1.0))):
    y = p / (p + (1 - p) * A)
    axc.plot([1 - A], [y], marker='o', ms=5, color=col, mew=0, zorder=5)
    axc.annotate(lab, (1 - A, y), xytext=xy, ha='center', va='top', fontsize=FSN - 0.4, color=col,
                 arrowprops=dict(arrowstyle='-', lw=0.6, color=col, shrinkB=3))
axc.set_xlim(0, 1.0); axc.set_xlabel('round-trip reflectance  $R_{\\rm LED}$ = 1 − $A$′')
axc.set_ylabel('extraction efficiency  η$_{\\rm ext}$')
axc.tick_params(labelleft=True)
axc.text(-0.3, 1.06, '(c)', transform=axc.transAxes, fontsize=8.5, fontweight='bold', va='bottom')
for ax in (axb, axc):
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=FSN)
for e in ('png', 'pdf'):
    fig.savefig(os.path.join(HERE, 'fig1_row.' + e), dpi=400)
print('fig1_row written')
