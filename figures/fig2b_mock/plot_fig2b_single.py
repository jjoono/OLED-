"""Fig. 2(b), single-panel variants.
A: the three curves as drafted by the author -- p, eta_ext(Al), eta_ext(Ag).
B: the same panel with the two EQE curves added, so the payoff of a high index also shows.
Data: fig2b_curves.csv."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import textwrap

D = np.genfromtxt('fig2b_curves.csv', delimiter=',', names=True)
n, p = D['n_sub'], D['p']
eAl, eAg = D['eta_ext_Al'], D['eta_ext_Ag']
qAl, qAg = D['EQE_Al'], D['EQE_Ag']
C_AL, C_AG, C_P = '#D55E00', '#0072B2', '#8C8C8C'

fig, axs = plt.subplots(1, 2, figsize=(11.8, 5.6), sharey=True)
fig.subplots_adjust(left=0.072, right=0.988, top=0.855, bottom=0.315, wspace=0.09)

for ax in axs:
    ax.fill_between(n, eAl, eAg, color=C_AG, alpha=0.09, lw=0)
    ax.plot(n, p,   color=C_P,  lw=1.7, ls=(0, (4, 2.2)), zorder=3)
    ax.plot(n, eAl, color=C_AL, lw=2.6, zorder=4, solid_capstyle='round')
    ax.plot(n, eAg, color=C_AG, lw=2.6, zorder=4, solid_capstyle='round')
    ax.set_xlabel('Refractive index of the substrate / MLA,  $n_{sub}$')
    ax.set_xlim(1.3, 2.0); ax.set_ylim(0, 1.05)
    ax.grid(axis='y', color='0.93', lw=0.8)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
axs[0].set_ylabel('Efficiency  /  escape probability')

# --- A: the author's three curves ------------------------------------------
a = axs[0]
a.set_title('A   three curves, as drafted', fontsize=11.5, loc='left', color='0.25')
for x, txt in ((1.33, None), (1.97, None)):
    pass
i0, i1 = 1, len(n) - 2
for i, dy in ((i0, 8), (i1, 8)):
    a.annotate('', (n[i], eAl[i]), (n[i], eAg[i]),
               arrowprops=dict(arrowstyle='<->', color='0.45', lw=1.0, shrinkA=0, shrinkB=0))
a.annotate(f'{100*(eAg[i0]-eAl[i0]):.0f} %p', (n[i0], 0.5*(eAl[i0]+eAg[i0])),
           xytext=(7, 0), textcoords='offset points', va='center', fontsize=9, color='0.35')
a.annotate(f'{100*(eAg[i1]-eAl[i1]):.0f} %p', (n[i1], 0.5*(eAl[i1]+eAg[i1])),
           xytext=(-7, 0), textcoords='offset points', va='center', ha='right', fontsize=9, color='0.35')
a.annotate('fewer photons escape per pass\n→ more round trips\n→ the mirror is tested harder',
           (1.45, 0.425), xytext=(-2, -66), textcoords='offset points', fontsize=9, color='0.3', ha='left',
           arrowprops=dict(arrowstyle='->', color='0.55', lw=1.0))
a.annotate('this widening gap is the whole point:\nthe cost of a lossy mirror grows with $n_{sub}$',
           (1.86, 0.725), xytext=(0, -74), textcoords='offset points', fontsize=9, color='0.3',
           ha='center', arrowprops=dict(arrowstyle='->', color='0.55', lw=1.0))

# --- B: + the two EQE curves ------------------------------------------------
b = axs[1]
b.set_title('B   the same panel, with EQE = $\\eta_{sub}^{(0)}\\,\\eta_{ext}$ added', fontsize=11.5, loc='left', color='0.25')
b.plot(n, qAl, color=C_AL, lw=2.8, ls=(0, (1.2, 1.6)), zorder=5, dash_capstyle='round')
b.plot(n, qAg, color=C_AG, lw=2.8, ls=(0, (1.2, 1.6)), zorder=5, dash_capstyle='round')
for q, c in ((qAg, C_AG), (qAl, C_AL)):
    i = int(np.argmax(q))
    b.plot(n[i], q[i], 'o', ms=7.5, mfc=c, mec='w', mew=1.4, zorder=6)
    b.annotate(f'{q[i]:.2f}', (n[i], q[i]), xytext=(4, -13), textcoords='offset points',
               fontsize=9.5, fontweight='bold', color=c)
b.annotate('Al: the gain in $\\eta_{sub}^{(0)}$ is spent on\nthe loss in $\\eta_{ext}$ — EQE turns over',
           (1.86, 0.535), xytext=(-6, -96), textcoords='offset points', fontsize=9, color='0.3',
           ha='right', arrowprops=dict(arrowstyle='->', color='0.55', lw=1.0))
b.annotate('Ag: no turnover —\nthe ceiling is gone', (1.93, 0.739), xytext=(-46, 34),
           textcoords='offset points', fontsize=9, color='0.3', ha='center',
           arrowprops=dict(arrowstyle='->', color='0.55', lw=1.0))
fig.legend(handles=[Line2D([], [], color=C_AG, lw=2.6, label='$\\eta_{ext}$,  Ag reflector'),
                    Line2D([], [], color=C_AL, lw=2.6, label='$\\eta_{ext}$,  Al reflector'),
                    Line2D([], [], color=C_P, lw=1.7, ls=(0, (4, 2.2)), label='$p$,  single-pass escape'),
                    Line2D([], [], color=C_AG, lw=2.8, ls=(0, (1.2, 1.6)), label='EQE, Ag  (panel B only)'),
                    Line2D([], [], color=C_AL, lw=2.8, ls=(0, (1.2, 1.6)), label='EQE, Al  (panel B only)')],
           loc='lower center', bbox_to_anchor=(0.53, 0.175), ncol=5, fontsize=9.5, frameon=True,
           framealpha=0.95, edgecolor='0.85', borderpad=0.6, handlelength=2.3, columnspacing=1.6)

fig.suptitle('Fig. 2(b) — one panel or two?  what each version can and cannot say',
             x=0.072, ha='left', fontsize=12.5, y=0.955)
foot = ('Same model as before (550 nm, isotropic dipole, PLQY = 1, reflector 100 nm / ETL 200 nm / EML 20 nm / HTL 200 nm / ITO 50 nm (n = 1.9 + 0.02i) / substrate, '
        'n_sub swept 1.30–2.00 and index-matched to the outcoupling structure).  '
        'A says: as n_sub rises, p falls, so each photon makes more round trips and the penalty for a lossy mirror grows — the two η_ext curves peel apart, from a few %p at n_sub = 1.3 to 15 %p at 2.0.  '
        'That is a complete and honest statement, but every curve in it falls, so the reader is left with "a high index only costs you".  '
        'B adds the missing half: η_sub^(0) rises steeply over the same range, so the product turns over at 0.54 for Al near n_sub = 1.8 while Ag keeps climbing to 0.75.')
fig.text(0.072, 0.018, '\n'.join(textwrap.wrap(foot, 190)), fontsize=8.3, va='bottom', ha='left', color='0.28')
fig.savefig('fig2b_single_vs_eqe.png', dpi=150)
print('saved')
