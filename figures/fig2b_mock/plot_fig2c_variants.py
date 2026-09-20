"""Three readings of Fig. 2(c), same data (fig2b_curves_konig.csv).
C1 the current four-curve panel, C2 the gap drawn as an area, C3 one reference curve for
the substrate-delivered power and the two EQE curves in colour."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import textwrap

D = np.genfromtxt('fig2b_curves_konig.csv', delimiter=',', names=True)
n = D['n_sub']
sAl, sAg, qAl, qAg = D['eta_sub_Al'], D['eta_sub_Ag'], D['EQE_Al'], D['EQE_Ag']
C_AL, C_AG, GREY = '#D55E00', '#0072B2', '#8C8C8C'
DASH = (0, (4, 2.2))

fig, axs = plt.subplots(1, 3, figsize=(13.2, 5.0))
fig.subplots_adjust(left=0.055, right=0.99, top=0.84, bottom=0.30, wspace=0.22)
for ax, tag in zip(axs, ('C1   as drawn now', 'C2   the gap as an area', 'C3   one reference curve')):
    ax.set_xlim(1.3, 2.0); ax.set_ylim(0, 1.0)
    ax.set_xticks([1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0])
    ax.set_xticklabels(['1.3', '', '1.5', '', '1.7', '', '1.9', ''])
    ax.set_xlabel('$n_{sub}$', fontsize=10)
    ax.grid(axis='y', color='0.93', lw=0.8)
    ax.tick_params(labelsize=9)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    ax.set_box_aspect(1.0)
    ax.set_title(tag, fontsize=11, loc='left', color='0.25')

def peak(ax, q, c, dy=-16):
    i = int(np.argmax(q))
    ax.plot(n[i], q[i], 'o', ms=6, mfc=c, mec='w', mew=1.3, zorder=6)
    ax.annotate(f'{q[i]:.2f}', (n[i], q[i]), xytext=(9, dy), textcoords='offset points',
                fontsize=9, fontweight='bold', color=c)

# --- C1 ---------------------------------------------------------------------
a = axs[0]
a.plot(n, sAg, color=C_AG, lw=1.5, ls=DASH); a.plot(n, sAl, color=C_AL, lw=1.5, ls=DASH)
a.plot(n, qAg, color=C_AG, lw=2.6); a.plot(n, qAl, color=C_AL, lw=2.6)
peak(a, qAg, C_AG); peak(a, qAl, C_AL)
a.set_ylabel('$\\eta_{sub}^{(0)}$,   EQE', fontsize=10.5)
a.legend(handles=[Line2D([], [], color='0.45', lw=1.5, ls=DASH, label='$\\eta_{sub}^{(0)}$'),
                  Line2D([], [], color='0.45', lw=2.4, label='EQE')],
         loc='lower right', bbox_to_anchor=(1.0, 0.015), fontsize=9, frameon=False, handlelength=2.0)

# --- C2 ---------------------------------------------------------------------
b = axs[1]
b.fill_between(n, qAl, sAl, color=C_AL, alpha=0.16, lw=0)
b.fill_between(n, qAg, sAg, color=C_AG, alpha=0.16, lw=0)
b.plot(n, sAl, color=C_AL, lw=1.0, alpha=0.65); b.plot(n, sAg, color=C_AG, lw=1.0, alpha=0.65)
b.plot(n, qAl, color=C_AL, lw=2.8); b.plot(n, qAg, color=C_AG, lw=2.8)
peak(b, qAg, C_AG); peak(b, qAl, C_AL)
b.set_ylabel('power fraction', fontsize=10.5)
b.annotate('delivered to the substrate,\nnever got out', (1.66, 0.5*(sAl[7]+qAl[7])),
           xytext=(-4, -72), textcoords='offset points', ha='center', fontsize=8.6, color=C_AL,
           arrowprops=dict(arrowstyle='->', color=C_AL, lw=0.9, alpha=0.8))
b.annotate('', (1.95, 0.5*(sAg[-2]+qAg[-2])), xytext=(-34, 26), textcoords='offset points',
           arrowprops=dict(arrowstyle='->', color=C_AG, lw=0.9, alpha=0.8))
b.text(1.905, 0.988, 'the same gap,\nwith a low-loss stack', ha='right', va='top',
       fontsize=8.6, color=C_AG)

# --- C3 ---------------------------------------------------------------------
c = axs[2]
c.plot(n, sAg, color=GREY, lw=1.8, ls=DASH, zorder=3)
c.plot(n, qAg, color=C_AG, lw=2.8); c.plot(n, qAl, color=C_AL, lw=2.8)
peak(c, qAg, C_AG); peak(c, qAl, C_AL)
c.set_ylabel('power fraction', fontsize=10.5)
c.annotate('$\\eta_{sub}^{(0)}$  delivered\nto the substrate', (1.63, sAg[6]), xytext=(4, 14),
           textcoords='offset points', fontsize=8.8, color='0.35',
           arrowprops=dict(arrowstyle='->', color='0.6', lw=0.9))

for ax in axs:
    ax.axvline(1.8, color='0.8', lw=1.0, ls=':', zorder=1)
axs[1].annotate('$n_{sub}$ = $n_{EML}$:\nthe waveguide\nmode vanishes', (1.8, 0.035),
                xytext=(-8, 0), textcoords='offset points', ha='right', va='bottom',
                fontsize=8.2, color='0.45')

fig.legend(handles=[Line2D([], [], color=C_AG, lw=2.8, label='Ag reflector  (low-loss)'),
                    Line2D([], [], color=C_AL, lw=2.8, label='Al reflector  (conventional)'),
                    Patch(facecolor='0.6', alpha=0.25, label='extraction loss,  $\\eta_{sub}^{(0)}$ − EQE  (C2)')],
           loc='lower center', bbox_to_anchor=(0.52, 0.145), ncol=3, fontsize=9.5, frameon=False,
           handlelength=2.4, columnspacing=2.4)
foot = ('Same data in all three. C1 asks the reader to hold two line styles and two colours at once, and above n_sub = 1.8 the Ag pair nearly touches (0.96 against 0.89).  '
        'C2 keeps both quantities but draws the difference — what reached the substrate and never escaped — as an area; the two bands never overlap, so the thin blue sliver against '
        'the thick orange band is the message.  C3 drops the Al η_sub curve and keeps one grey reference, which is the cleanest read but mixes two effects in the orange gap, '
        'since the Al stack also absorbs before the light reaches the substrate (η_sub 0.84 against 0.96 at n_sub = 1.8).')
fig.text(0.055, 0.02, '\n'.join(textwrap.wrap(foot, 168)), fontsize=8.3, va='bottom', ha='left', color='0.3')
fig.suptitle('Fig. 2(c) — three arrangements of the same numbers', x=0.055, ha='left', fontsize=12.5, y=0.955)
fig.savefig('fig2c_variants.png', dpi=160)
print('saved')
