"""Fig. 2(d) mock: angle- and wavelength-resolved round-trip loss for the candidate
electrode structures, with the flux-weighted averages the text quotes."""
import numpy as np, matplotlib, sys, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import textwrap
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'sim', 'design_rule4'))
import fig2d as F

NSUB = float(os.environ.get('NSUB', 1.5))
SPEC = F.GREEN if NSUB < 1.65 else F.ORANGE
VIS = np.arange(430.0, 701.0)
NAMES = list(F.stacks())
MAPPED = ['Al / ITO 150 nm', 'Ag / IZO 50 nm', 'DBR / IZO 50 nm']
COL = {'Al / ITO 150 nm': '#D55E00', 'Ag / ITO 150 nm': '#56B4E9',
       'Ag / IZO 50 nm': '#0072B2', 'Ag / Ag 10 nm': '#009E73', 'DBR / IZO 50 nm': '#9467BD'}

W = {k: F.weighted(k, n_sub=NSUB, spectrum=SPEC, lam=VIS) for k in NAMES}

fig = plt.figure(figsize=(14.0, 5.3))
axs = [fig.add_axes([0.045 + i*0.185, 0.375, 0.150, 0.47]) for i in range(3)]
cax = fig.add_axes([0.605, 0.375, 0.009, 0.47])
axl = fig.add_axes([0.715, 0.375, 0.265, 0.47])

for ax, name in zip(axs, MAPPED):
    th, lam, R, T = F.maps(name, n_sub=NSUB, lam=VIS)
    im = ax.pcolormesh(np.degrees(th), lam, 100*(1 - R), cmap='inferno_r', vmin=0, vmax=40,
                       shading='auto', rasterized=True)
    ax.set_xlabel('angle in the substrate (°)', fontsize=8.8)
    ax.set_title(name, fontsize=9.5, loc='left')
    ax.set_xticks([0, 30, 60, 90]); ax.tick_params(labelsize=8.5)
    if ax is axs[0]: ax.set_ylabel('wavelength (nm)', fontsize=9.5)
    else: ax.set_yticklabels([])
cb = fig.colorbar(im, cax=cax); cb.set_label('round-trip loss 1 − R  (%)', fontsize=8.8)
cb.ax.tick_params(labelsize=8)

th = np.linspace(0, np.pi/2, 452); th = 0.5*(th[1:] + th[:-1])
w = np.cos(th)*np.sin(th); w /= w.max()
axl.fill_between(np.degrees(th), 0, 9*w, color='0.90', lw=0, zorder=0)
axl.text(45, 0.4, 'weight of each angle,  $\\cos\\theta\\,\\sin\\theta$',
         fontsize=7.6, color='0.5', ha='center', va='bottom', zorder=1)
for name in NAMES:
    t2, lam, R, T = F.maps(name, n_sub=NSUB, lam=VIS)
    sw = np.clip(np.interp(lam, F.LAM, SPEC), 0, None); sw /= sw.sum()
    axl.plot(np.degrees(t2), 100*(1 - (R*sw[:, None]).sum(0)), color=COL[name], lw=2.2)
axl.set_xlim(0, 90); axl.set_ylim(0, 40); axl.set_xticks([0, 30, 60, 90])
axl.set_xlabel('angle in the substrate (°)', fontsize=8.8)
axl.set_ylabel('round-trip loss 1 − R  (%)', fontsize=9.5)
axl.tick_params(labelsize=8.5)
axl.grid(axis='y', color='0.93', lw=0.8)
for s in ('top', 'right'): axl.spines[s].set_visible(False)
axl.set_title('spectrum-averaged, all five structures', fontsize=9.5, loc='left')
order = sorted(NAMES, key=lambda k: -W[k]['loss'])
handles = []
for k in order:
    lab = f"{k}:  {100*W[k]['loss']:.1f} %"
    if W[k]['transmitted'] > 0.005:
        lab += f"  ({100*W[k]['absorbed']:.1f} absorbed + {100*W[k]['transmitted']:.1f} leaked)"
    handles.append(Line2D([], [], color=COL[k], lw=2.4, label=lab))
fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.52, 0.155), ncol=3,
           fontsize=8.6, frameon=False, handlelength=2.0, columnspacing=2.0,
           title='flux- and spectrum-weighted round-trip loss', title_fontsize=8.6)
axl.annotate('beyond the escape cone to air\nthe dielectric mirror cannot leak',
             (np.degrees(np.arcsin(1/NSUB)), 21), xytext=(6, 30), textcoords='offset points',
             fontsize=7.6, color='#9467BD', ha='left',
             arrowprops=dict(arrowstyle='->', color='#9467BD', lw=0.9))
axl.axvline(np.degrees(np.arcsin(1/NSUB)), color='#9467BD', lw=0.8, ls=':', zorder=1)

fig.suptitle(f'Fig. 2(d)  Round-trip loss of the candidate electrode structures  '
             f'(n$_{{sub}}$ = {NSUB}, {"green" if NSUB < 1.65 else "orange"} emitter)',
             x=0.045, ha='left', fontsize=12.5, y=0.955)
foot = ('Light arriving from the substrate on the generic stack of Fig. 2(a)–(c): substrate / transparent electrode / 420 nm of non-absorbing organics (n = 1.8) / reflector, '
        'dispersive n,k for every electrode and mirror (McPeak Ag, Johnson–Christy-type Al, Koenig ITO, measured IZO, ZnS and LiF).  '
        '1 − R is what one round trip loses: for the metal mirrors all of it is absorption, for the ZnS/LiF stack (4.5 pairs, 70/115 nm as deposited on the orange device) most of it is '
        'light leaking through outside the stopband, which is why it is listed separately.  The percentages are averaged over the emission spectrum and over angle with the cos·sin weight '
        'that an angularly randomising outcoupling structure enforces, so they are the A′ that eq. (2) takes.')
fig.text(0.045, 0.015, '\n'.join(textwrap.wrap(foot, 205)), fontsize=8.0, va='bottom', ha='left', color='0.28')
fig.savefig(f'fig2d_mock_{str(NSUB).replace(".","")}.png', dpi=150)
print('saved', NSUB, {k: round(100*v['loss'], 1) for k, v in W.items()})
