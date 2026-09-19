"""Fig. 1(b) and 1(c): the recycling process and the law it sums to.
Everything here is eq. (2); the only inputs are p and A'."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import textwrap

P = 0.40
CASES = [(0.02, '#0072B2', 'lossless reflector'),
         (0.10, '#7FB1DC', None),
         (0.25, '#C8553D', 'Al electrode')]
BAR = {0.02: '#4A3B2A', 0.10: '#9C7B52', 0.25: '#E3CBA5'}
N = 10
eta = lambda a: P/(P + (1-P)*a)

fig, (axb, axc) = plt.subplots(1, 2, figsize=(12.4, 4.9))
fig.subplots_adjust(left=0.065, right=0.945, top=0.90, bottom=0.30, wspace=0.42)

# ---------------- (b) process ----------------
n = np.arange(N+1); w = 0.26
axr = axb.twinx()
for k, (a, c, _) in enumerate(CASES):
    q = (1-P)*(1-a)
    axb.bar(n + (k-1)*w, P*q**n, width=w, color=BAR[a], edgecolor='0.35', lw=0.4, zorder=3)
    axr.plot(n, P*(1-q**(n+1))/(1-q), color=c, lw=2.4, zorder=4, solid_capstyle='round')
    axr.annotate(f'{eta(a):.2f}', (N-0.4, eta(a)), xytext=(0, 7), textcoords='offset points',
                 color=c, fontsize=9.5, fontweight='bold', va='bottom', ha='right')
axb.set_xlabel('Number of round trips'); axb.set_xlim(-0.6, N+0.6)
axb.set_ylabel('Escaped fraction per pass'); axb.set_ylim(0, 0.44)
axr.set_ylabel('Cumulative escaped fraction  ( → $\\eta_{ext}$ )', color=CASES[0][1])
axr.set_ylim(0, 1.0); axr.tick_params(axis='y', colors=CASES[0][1])
axr.spines['right'].set_color(CASES[0][1]); axb.spines['top'].set_visible(False); axr.spines['top'].set_visible(False)
axb.set_xticks(n)
axb.annotate('n = 0: first pass,\nbefore any round trip', (0, P/2), xytext=(40, 2),
             textcoords='offset points', fontsize=8.6, color='0.35',
             arrowprops=dict(arrowstyle='-', color='0.6', lw=0.8))
axb.legend(handles=[Line2D([], [], marker='s', ls='none', mfc=BAR[a], mec='0.35', ms=9,
                           label=f"A′ = {a}") for a, _, _ in CASES],
           loc='center right', bbox_to_anchor=(1.0, 0.52), fontsize=9, frameon=False,
           handletextpad=0.5, title='per pass', title_fontsize=8.5)
axb.set_title('(b)  the recycling process', fontsize=11.5, loc='left')

# ---------------- (c) the law ----------------
R = np.linspace(0, 1, 500)
axc.plot(R, P/(P + (1-P)*(1-R)), color='#B00', lw=2.6, zorder=4)
for pp, sty in ((0.25, (0, (4, 2))), (0.60, (0, (1.2, 2)))):
    axc.plot(R, pp/(pp + (1-pp)*(1-R)), color='0.65', lw=1.2, ls=sty, zorder=2)
    xl = 0.32 if pp > 0.4 else 0.55
    axc.annotate(f'p = {pp}', (xl, pp/(pp+(1-pp)*(1-xl))), xytext=(0, -9 if pp > 0.4 else 10),
                 textcoords='offset points', ha='center', va='top' if pp > 0.4 else 'bottom',
                 fontsize=8.6, color='0.5')
for a, c, lab in CASES:
    axc.plot(1-a, eta(a), 'o', ms=9, mfc=c, mec='w', mew=1.6, zorder=6)
    if lab:
        # park both labels in the empty upper-left wedge, well clear of the curves
        tx, ty = (0.30, 0.86) if a > 0.1 else (0.30, 0.99)
        axc.annotate(f"A′ = {a}   {lab}", (1-a, eta(a)), xytext=(tx, ty), textcoords='data',
                     ha='left', va='center', fontsize=9, color=c, fontweight='bold',
                     arrowprops=dict(arrowstyle='-', color=c, lw=1.0,
                                     connectionstyle='angle3,angleA=0,angleB=60'))
axc.set_xlabel('Round-trip reflectance,  $R_{LED} = 1 - A′$')
axc.set_ylabel('Extraction efficiency  ($\\eta_{ext}$)')
axc.set_xlim(0, 1.02); axc.set_ylim(0, 1.12)
for s in ('top', 'right'): axc.spines[s].set_visible(False)
axc.set_title('(c)  what it sums to:  $\\eta_{ext} = p\\,/\\,[\\,p + (1-p)A′\\,]$', fontsize=11.5, loc='left')

axb.text(0.03, 0.97, f'p = {P}', transform=axb.transAxes, ha='left', va='top',
         fontsize=10, color='0.35', style='italic')
axc.text(0.03, 0.70, f'p = {P}', transform=axc.transAxes, ha='left', va='top',
         fontsize=10, color='0.35', style='italic')
foot = ('Both panels are eq. (2) with no other input than the single-pass escape probability p and the round-trip loss A′.  '
        '(b) bars are the fraction escaping at each pass, p[(1−p)(1−A′)]ⁿ, read on the left axis; lines are the running total, read on the right.  '
        'n = 0 is the first pass, before any round trip.  The value each line saturates at is η_ext, and those three values are the three points in (c).  '
        '(c) the grey curves show that p only shifts the family; the marked points are the same three A′ as in (b).  p = 0.4 is the escape probability of ordinary glass.')
fig.text(0.065, 0.02, '\n'.join(textwrap.wrap(foot, 190)), fontsize=8.3, va='bottom', ha='left', color='0.28')
fig.savefig('fig1bc_mock.png', dpi=150)
print('saved')
