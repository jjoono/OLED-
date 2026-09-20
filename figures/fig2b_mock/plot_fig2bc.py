"""Fig. 2(b) and 2(c) as two small square panels, sharing the substrate-index axis.
(b) the cost:   p falls with n_sub, so the two eta_ext curves peel apart.
(c) the result: eta_sub keeps rising, but the product turns over for Al and not for Ag.
Data: fig2b_curves.csv."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import textwrap, os

D = np.genfromtxt(os.environ.get('DATA', 'fig2b_curves_konig.csv'), delimiter=',', names=True)
n, p = D['n_sub'], D['p']
eAl, eAg = D['eta_ext_Al'], D['eta_ext_Ag']
sAl, sAg = D['eta_sub_Al'], D['eta_sub_Ag']
qAl, qAg = D['EQE_Al'], D['EQE_Ag']
C_AL, C_AG, C_P = '#D55E00', '#0072B2', '#8C8C8C'
DASH = (0, (4, 2.2))

fig = plt.figure(figsize=(7.6, 4.9))
axb = fig.add_axes([0.095, 0.275, 0.375, 0.600])
axc = fig.add_axes([0.600, 0.275, 0.375, 0.600])
for ax, tag in ((axb, '(b)'), (axc, '(c)')):
    ax.set_xlim(1.3, 2.0); ax.set_ylim(0, 1.0)
    ax.set_xticks([1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0])
    ax.set_xticklabels(['1.3', '', '1.5', '', '1.7', '', '1.9', ''])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel('$n_{sub}$  (substrate / outcoupling structure)', fontsize=9.5)
    ax.tick_params(labelsize=9)
    ax.grid(axis='y', color='0.93', lw=0.8)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    ax.set_box_aspect(1.0)
    ax.text(-0.26, 1.10, tag, transform=ax.transAxes, fontsize=13, fontweight='bold', va='top')

# ---- (b) the cost: p down, the eta_ext curves peel apart --------------------
axb.fill_between(n, eAl, eAg, color=C_AG, alpha=0.09, lw=0)
axb.plot(n, p,   color=C_P,  lw=1.6, ls=DASH, zorder=3)
axb.plot(n, eAl, color=C_AL, lw=2.4, zorder=4, solid_capstyle='round')
axb.plot(n, eAg, color=C_AG, lw=2.4, zorder=4, solid_capstyle='round')
axb.set_ylabel('$\\eta_{ext}$,   $p$', fontsize=10.5)
axb.annotate('', (1.965, eAl[-2]), (1.965, eAg[-2]),
             arrowprops=dict(arrowstyle='<->', color='0.45', lw=0.9, shrinkA=0, shrinkB=0))
axb.annotate(f'{100*(eAg[-2]-eAl[-2]):.0f} %p', (1.965, 0.5*(eAl[-2]+eAg[-2])), xytext=(-5, 0),
             textcoords='offset points', ha='right', va='center', fontsize=8.5, color='0.35')
axb.annotate('Ag', (1.70, eAg[8]), xytext=(0, 5), textcoords='offset points',
             color=C_AG, fontsize=10, fontweight='bold', ha='center')
axb.annotate('Al', (1.70, eAl[8]), xytext=(0, -14), textcoords='offset points',
             color=C_AL, fontsize=10, fontweight='bold', ha='center')
axb.annotate('$p$', (1.90, p[12]), xytext=(0, 7), textcoords='offset points',
             color='0.35', fontsize=10.5, ha='center')
axb.text(1.325, 0.075, 'fewer escapes per pass\n→ more round trips\n→ the mirror is tested harder',
         fontsize=8.2, color='0.35', va='bottom')

# ---- (c) the result: the payoff and the turnover ----------------------------
axc.plot(n, sAg, color=C_AG, lw=1.5, ls=DASH, zorder=3)
axc.plot(n, sAl, color=C_AL, lw=1.5, ls=DASH, zorder=3)
axc.plot(n, qAg, color=C_AG, lw=2.6, zorder=4, solid_capstyle='round')
axc.plot(n, qAl, color=C_AL, lw=2.6, zorder=4, solid_capstyle='round')
axc.set_ylabel('$\\eta_{sub}^{(0)}$,   EQE', fontsize=10.5)
for q, c, dy in ((qAg, C_AG, -16), (qAl, C_AL, -16)):
    i = int(np.argmax(q))
    axc.plot(n[i], q[i], 'o', ms=6, mfc=c, mec='w', mew=1.3, zorder=5)
    axc.annotate(f'{q[i]:.2f}', (n[i], q[i]), xytext=(9, dy), textcoords='offset points',
                 fontsize=9, fontweight='bold', color=c, ha='left')
axc.legend(handles=[Line2D([], [], color='0.45', lw=1.5, ls=DASH, label='$\\eta_{sub}^{(0)}$'),
                    Line2D([], [], color='0.45', lw=2.4, label='EQE')],
           loc='lower right', bbox_to_anchor=(1.0, 0.015), fontsize=9, frameon=False,
           handlelength=2.0, labelspacing=0.35)
axc.text(1.325, 0.03, 'Al: the gain in $\\eta_{sub}^{(0)}$\nis spent on the loss\nin $\\eta_{ext}$ — EQE turns\nover near $n_{sub}$ = 1.8',
         fontsize=8.2, color='0.35', va='bottom')

fig.legend(handles=[Line2D([], [], color=C_AG, lw=2.4, label='Ag reflector  (low-loss)'),
                    Line2D([], [], color=C_AL, lw=2.4, label='Al reflector  (conventional)')],
           loc='lower center', bbox_to_anchor=(0.535, 0.125), ncol=2, fontsize=9, frameon=False,
           handlelength=2.4, columnspacing=2.6)
foot = ('550 nm, isotropic dipole, PLQY = 1; reflector 100 nm / ETL 200 nm / EML 20 nm / HTL 200 nm / ITO 50 nm (n = 1.864 + 0.0032i at 550 nm, Koenig 2014) / substrate, '
        'n_sub index-matched to the outcoupling structure.  η_ext = p/[p + (1−p)A′], EQE = η_sub^(0) η_ext.')
fig.text(0.085, 0.02, '\n'.join(textwrap.wrap(foot, 118)), fontsize=7.8, va='bottom', ha='left', color='0.3')
fig.savefig(os.environ.get('OUT', 'fig2bc_squares.png'), dpi=170)
print('saved')
