"""Where the beyond-cone loss of the dielectric mirror comes from.
Beyond the escape cone the light is totally reflected only at the last interface, so it still
propagates through every layer on the way there; the comb is the stack's own resonances."""
import numpy as np, matplotlib, os, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import textwrap
sys.path.insert(0, os.path.join('..', '..', 'sim', 'design_rule4'))
import fig2d as F

VIS = np.arange(430.0, 701.0); SEL = np.isin(F.LAM, VIS)
TH = np.radians(np.linspace(0.05, 89.75, 900))
SW = np.clip(np.interp(VIS, F.LAM, F.GREEN), 0, None); SW /= SW.sum()
THC = np.degrees(np.arcsin(1/1.5))
ITO0 = F.ITO.real.astype(complex); ZNS0 = F.ZNS.real.astype(complex); LIF0 = F.LIF.real.astype(complex)
DBR = lambda z, l: [(z if i % 2 == 0 else l, d) for i, (m, d) in enumerate(F.dbr_chirp())]

def curve(layers, lam=None):
    sel = np.isin(F.LAM, VIS) if lam is None else np.isin(F.LAM, lam)
    g = VIS if lam is None else lam
    R, T = F.RT([(m[sel], d) for m, d in layers], np.full(g.shape, 1.5, dtype=complex),
                np.ones(g.shape, dtype=complex), TH, g)
    return R, T

base = [(F.ITO, 150.0), (F.ORG, 420.0)]
sets = [('every layer as measured', base + F.dbr_chirp(), '#9467BD', '-'),
        ('ITO made lossless', [(ITO0, 150.0), (F.ORG, 420.0)] + F.dbr_chirp(), '#D55E00', (0, (4, 2))),
        ('LiF made lossless', base + DBR(F.ZNS, LIF0), '#009E73', (0, (1.5, 1.5))),
        ('everything lossless', [(ITO0, 150.0), (F.ORG, 420.0)] + DBR(ZNS0, LIF0), '0.35', (0, (6, 2, 1, 2)))]

fig, axs = plt.subplots(1, 2, figsize=(11.6, 4.5))
fig.subplots_adjust(left=0.065, right=0.985, top=0.85, bottom=0.30, wspace=0.22)
L550 = np.array([550.0])
for lab, lay, c, ls in sets:
    R, T = curve(lay, L550)
    axs[0].plot(np.degrees(TH), 100*(1 - R[0]), color=c, lw=1.7, ls=ls, label=lab)
    R, T = curve(lay)
    axs[1].plot(np.degrees(TH), 100*((1 - R)*SW[:, None]).sum(0), color=c, lw=1.7, ls=ls, label=lab)
for ax, t in ((axs[0], 'at 550 nm'), (axs[1], 'averaged over the emission spectrum')):
    ax.axvspan(0, THC, color='0.94', lw=0)
    ax.text(THC/2, 9.4, 'escape cone\nto air', ha='center', va='top', fontsize=8.5, color='0.45')
    ax.set_xlim(0, 90); ax.set_ylim(0, 10); ax.set_xticks([0, 30, 60, 90])
    ax.set_xlabel('$\\theta$ in substrate (°)'); ax.set_ylabel('round-trip loss 1 − R  (%)')
    ax.set_title(t, fontsize=11, loc='left')
    ax.grid(axis='y', color='0.93', lw=0.8)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
axs[0].legend(loc='upper left', fontsize=8.8, frameon=False)
axs[1].text(0.97, 0.95, 'beyond the cone, mean of the spectrum average:\nas measured 4.03 %   ·   ITO lossless 1.24 %\nLiF lossless 2.87 %   ·   all lossless 0.000 %',
            transform=axs[1].transAxes, ha='right', va='top', fontsize=8.5, color='0.3')
fig.suptitle('Beyond the escape cone the light still crosses every layer — that is what the comb is',
             x=0.065, ha='left', fontsize=12.5, y=0.955)
foot = ('Chirped ZnS/LiF stack on the generic Fig. 2 stack, n_sub = 1.5.  Beyond 41.8° the wave is evanescent in air, so nothing escapes — but it is still propagating in the ITO (n = 1.86), the organics (1.8), '
        'ZnS (2.36) and LiF (1.41, until 70.1°), and the total reflection happens at the last interface rather than at the entrance.  The light therefore samples every absorbing layer on the way there and back, '
        'and at the angle-wavelength pairs where the stack is resonant the field builds up inside it and the absorption is enhanced: that is the comb.  Remove the absorption and the comb goes with it — with every '
        'layer lossless the loss beyond the cone is identically zero, as pure total internal reflection requires.  What is left is 2.9 %p from the 150 nm ITO and 1.2 %p from the ten LiF layers, whose k is only '
        '2e-4 but which add up to about a micrometre of material.')
fig.text(0.065, 0.02, '\n'.join(textwrap.wrap(foot, 190)), fontsize=8.2, va='bottom', ha='left', color='0.3')
fig.savefig('dbr_fringes.png', dpi=160)
print('saved')
