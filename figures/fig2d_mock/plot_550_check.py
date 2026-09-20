"""Why the 550 nm cut and the spectrum-averaged curve disagree about the dielectric mirror."""
import numpy as np, matplotlib, os, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import textwrap
sys.path.insert(0, os.path.join('..', '..', 'sim', 'design_rule4'))
import fig2d as F

VIS = np.arange(430.0, 701.0)
C = {'Al': '#D55E00', 'Ag': '#0072B2', 'DBR, re-optimised': '#9467BD'}
LAB = {'Al': 'Al', 'Ag': 'Ag', 'DBR, re-optimised': 'DBR (optimised)'}
fig, axs = plt.subplots(1, 3, figsize=(13.2, 4.3))
fig.subplots_adjust(left=0.055, right=0.99, top=0.85, bottom=0.30, wspace=0.26)
sw = np.clip(np.interp(VIS, F.LAM, F.GREEN), 0, None); sw /= sw.sum()
for name, c in C.items():
    th, lam, R, T = F.maps(name, n_sub=1.5, lam=VIS)
    d = np.degrees(th); i550 = int(np.argmin(abs(lam - 550)))
    axs[0].plot(d, 100*(1 - R[i550]), color=c, lw=1.8)
    axs[1].plot(d, 100*(1 - (R*sw[:, None]).sum(0)), color=c, lw=1.8)
    w = np.cos(th)*np.sin(th); w /= w.sum()
    axs[2].plot(lam, 100*((1 - R)*w[None, :]).sum(1), color=c, lw=1.8)
axs[2].fill_between(VIS, 0, 32*sw/sw.max(), color='0.9', lw=0, zorder=0)
axs[2].text(640, 1.0, 'emission spectrum', fontsize=8.5, color='0.5', ha='center')
for ax, t in ((axs[0], 'at 550 nm only'), (axs[1], 'averaged over the emission spectrum')):
    ax.set_xlim(0, 90); ax.set_ylim(0, 25); ax.set_xticks([0, 30, 60, 90])
    ax.set_xlabel('$\\theta$ in substrate (°)'); ax.set_ylabel('round-trip loss 1 − R  (%)')
    ax.set_title(t, fontsize=11, loc='left')
    ax.axvline(41.81, color='0.8', lw=0.8, ls=':')
axs[2].set_xlim(440, 690); axs[2].set_ylim(0, 32)
axs[2].set_xlabel('wavelength (nm)'); axs[2].set_ylabel('flux-weighted 1 − R  (%)')
axs[2].set_title('flux-weighted, per wavelength', fontsize=11, loc='left')
for ax in axs:
    ax.grid(axis='y', color='0.93', lw=0.8)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
axs[0].legend(handles=[Line2D([], [], color=c, lw=1.8, label=LAB[n]) for n, c in C.items()],
              loc='upper left', fontsize=9, frameon=False)
axs[0].text(0.97, 0.95, 'flux-weighted at 550 nm\nAl 15.3 %   Ag 4.29 %   DBR 3.96 %',
            transform=axs[0].transAxes, ha='right', va='top', fontsize=8.6, color='0.3')
axs[1].text(0.97, 0.95, 'over the spectrum\nAl 15.6 %   Ag 4.50 %   DBR 5.47 %',
            transform=axs[1].transAxes, ha='right', va='top', fontsize=8.6, color='0.3')
fig.suptitle('The dielectric mirror does beat Ag at 550 nm — and loses it over the band',
             x=0.055, ha='left', fontsize=12.5, y=0.955)
foot = ('Left: the 550 nm row of the map, which is what a single-wavelength cut of the raw data gives — the optimised DBR sits at 2–3 % from normal incidence to the escape cone, below Ag throughout, '
        'and the narrow spikes at 41.8°, 72° and 82° are the escape-cone edge, the angle where LiF stops propagating (n_LiF = 1.409 against an in-plane index of 1.5 sin θ) and a further stack resonance.  '
        'Centre: the same curves averaged over the green emission spectrum; the ranking reverses.  Right: why — the metal mirrors are flat in wavelength while the dielectric stack is only good inside its '
        'stopband, and the emitter runs from 511 to 597 nm at 10–90 % of its integral.  The optimiser placed the band where the emission and the cos·sin weight are largest, which costs normal incidence, '
        'where the flux weight is zero.')
fig.text(0.055, 0.02, '\n'.join(textwrap.wrap(foot, 190)), fontsize=8.2, va='bottom', ha='left', color='0.3')
fig.savefig('dbr_550_vs_band.png', dpi=160)
print('saved')
