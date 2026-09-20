"""Fig. 2(b), third variant: the author's three curves as the hero panel,
with the two EQE curves tucked into an inset so the payoff of a high index still shows."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import textwrap

D = np.genfromtxt('fig2b_curves.csv', delimiter=',', names=True)
n, p = D['n_sub'], D['p']
eAl, eAg = D['eta_ext_Al'], D['eta_ext_Ag']
qAl, qAg = D['EQE_Al'], D['EQE_Ag']
C_AL, C_AG, C_P = '#D55E00', '#0072B2', '#8C8C8C'

fig = plt.figure(figsize=(7.6, 6.0))
ax = fig.add_axes([0.105, 0.265, 0.875, 0.615])
ax.fill_between(n, eAl, eAg, color=C_AG, alpha=0.09, lw=0)
ax.plot(n, p,   color=C_P,  lw=1.7, ls=(0, (4, 2.2)), zorder=3)
ax.plot(n, eAl, color=C_AL, lw=2.6, zorder=4, solid_capstyle='round')
ax.plot(n, eAg, color=C_AG, lw=2.6, zorder=4, solid_capstyle='round')
ax.set_xlim(1.3, 2.0); ax.set_ylim(0, 1.05)
ax.set_xlabel('Refractive index of the substrate / MLA,  $n_{sub}$')
ax.set_ylabel('Efficiency  /  escape probability')
ax.grid(axis='y', color='0.93', lw=0.8)
for s in ('top', 'right'): ax.spines[s].set_visible(False)

ax.annotate('', (1.965, eAl[-2]), (1.965, eAg[-2]),
            arrowprops=dict(arrowstyle='<->', color='0.45', lw=1.0, shrinkA=0, shrinkB=0))
ax.annotate(f'{100*(eAg[-2]-eAl[-2]):.0f} %p', (1.965, 0.5*(eAl[-2]+eAg[-2])), xytext=(-6, 0),
            textcoords='offset points', ha='right', va='center', fontsize=9, color='0.35')
k = 8   # n_sub = 1.70
ax.annotate('$\\eta_{ext}$,  Ag', (n[k], eAg[k]), xytext=(0, 7), textcoords='offset points',
            color=C_AG, fontsize=10.5, fontweight='bold', ha='center')
ax.annotate('$\\eta_{ext}$,  Al', (n[k], eAl[k]), xytext=(0, -17), textcoords='offset points',
            color=C_AL, fontsize=10.5, fontweight='bold', ha='center')
ax.annotate('$p$,  single-pass escape\nfewer escapes per pass → more round trips',
            (1.84, 0.293), xytext=(0, 40), textcoords='offset points',
            color='0.3', fontsize=9.2, ha='center',
            arrowprops=dict(arrowstyle='->', color='0.55', lw=1.0))

ins = ax.inset_axes([0.045, 0.045, 0.395, 0.265])
ins.plot(n, qAg, color=C_AG, lw=2.0)
ins.plot(n, qAl, color=C_AL, lw=2.0)
for q, c in ((qAg, C_AG), (qAl, C_AL)):
    i = int(np.argmax(q))
    ins.plot(n[i], q[i], 'o', ms=5, mfc=c, mec='w', mew=1.2, zorder=5)
    ins.annotate(f'{q[i]:.2f}', (n[i], q[i]), xytext=(4, -12), textcoords='offset points',
                 fontsize=8.5, fontweight='bold', color=c)
ins.set_xlim(1.3, 2.0); ins.set_ylim(0, 1.05)
ins.set_yticks([0, 0.4, 0.8]); ins.set_xticks([1.4, 1.6, 1.8, 2.0])
ins.tick_params(labelsize=8, length=3, pad=2)
ins.text(0.5, 0.97, 'EQE = $\\eta_{sub}^{(0)}\\,\\eta_{ext}$', transform=ins.transAxes,
         fontsize=9, ha='center', va='top')
ins.grid(axis='y', color='0.93', lw=0.7)
for s in ('top', 'right'): ins.spines[s].set_visible(False)
ins.set_facecolor('white')

fig.suptitle('Fig. 2(b) — variant C: three curves, EQE in an inset', x=0.105, ha='left',
             fontsize=12, y=0.965)
foot = ('The main panel carries the message as drafted — p falls with n_sub, so the two η_ext curves peel apart and the cost of a lossy mirror grows. '
        'The inset keeps the other half on the page without competing for attention: the Al device turns over at 0.54, the Ag device climbs to 0.75 and stays there.')
fig.text(0.105, 0.02, '\n'.join(textwrap.wrap(foot, 105)), fontsize=8.3, va='bottom', ha='left', color='0.28')
fig.savefig('fig2b_single_inset.png', dpi=150)
print('saved')
