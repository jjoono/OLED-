# -*- coding: utf-8 -*-
"""Figure 2 with the derived panels folded away.

eta_ext is a deterministic function of A' through eq. (2), so the two panels that plotted it
on its own carried no data the loss panels did not already have.  Here the ITO-thickness
panels get a second axis in eta_ext instead, and the substrate-index sweep is one panel with
p, A' and eta_ext in three separate bands.  Two rows of four."""
import numpy as np, matplotlib, os, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'sim', 'design_rule4'))
import fig2d as F

C_AL, C_AG, C_TCO, C_DBR = '#D55E00', '#0072B2', '#9C9C9C', '#9467BD'
DASH = (0, (3.2, 1.8))
FS_LAB, FS_TICK, FS_NOTE, FS_LET = 6.8, 6.2, 5.9, 8.5
plt.rcParams.update({'font.size': FS_LAB, 'axes.linewidth': 0.7,
                     'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
                     'xtick.major.size': 2.4, 'ytick.major.size': 2.4})
P_AB = 0.38                     # escape probability on glass, where (a) and (b) are computed

A = os.path.join(ROOT, 'figures', 'fig2a_mock')
Bd = os.path.join(ROOT, 'figures', 'fig2b_mock')
def a_load(m):
    T = np.loadtxt(f'{A}/kn2_{m}_15.csv', delimiter=','); Z = np.loadtxt(f'{A}/kn3_{m}_15_k0.csv', delimiter=',')
    return T[:, 0], Z[:, 2], T[:, 2] - Z[:, 2]
al, ag = a_load('al'), a_load('ag')
D = np.genfromtxt(f'{Bd}/fig2b_curves_konig.csv', delimiter=',', names=True)
nb, pb = D['n_sub'], D['p']
eAl, eAg, sAl, sAg, qAl, qAg = (D['eta_ext_Al'], D['eta_ext_Ag'], D['eta_sub_Al'],
                               D['eta_sub_Ag'], D['EQE_Al'], D['EQE_Ag'])
aAl, aAg = D['Aprime_Al'], D['Aprime_Ag']

fig = plt.figure(figsize=(7.2, 5.1))
H = 0.325
Y = [0.605, 0.175]
H2 = 0.275
x4 = [0.075, 0.330, 0.585, 0.820]; w4 = 0.150
ax = {}
for k, i in zip('abcd', range(4)): ax[k] = fig.add_axes([x4[i], Y[0], w4, H])
xm = [0.075, 0.228, 0.381]; wm = 0.135
for k, i in zip('efg', range(3)): ax[k] = fig.add_axes([xm[i], Y[1], wm, H2])
ax_cb = fig.add_axes([0.524, Y[1], 0.008, H2])
ax['h'] = fig.add_axes([0.655, Y[1], 0.315, H2])

def clean(a_):
    for s in ('top', 'right'): a_.spines[s].set_visible(False)
    a_.tick_params(labelsize=FS_TICK, pad=1.5)
    a_.grid(axis='y', color='0.93', lw=0.6)
def letter(k, dx=-0.34, dy=1.16):
    ax[k].text(dx, dy, k, transform=ax[k].transAxes, fontsize=FS_LET, fontweight='bold', va='top')

# ---- (a), (b): the loss split, with eta_ext as a second axis --------------
for k, (d, Am, At), c, name, ymax, eticks in (
        ('a', al, C_AL, 'Al reflector', 0.19, [1.0, 0.95, 0.90, 0.85, 0.80]),
        ('b', ag, C_AG, 'Ag reflector', 0.062, [1.0, 0.98, 0.96, 0.94, 0.92])):
    a_ = ax[k]
    a_.fill_between(d, 0, Am, color=c, alpha=0.9, lw=0)
    a_.fill_between(d, Am, Am + At, color=C_TCO, alpha=0.9, lw=0)
    a_.plot(d, Am + At, color='k', lw=0.8)
    a_.text(196, Am[-1]/2, 'mirror', ha='right', va='center', fontsize=FS_NOTE, color='w', fontweight='bold')
    a_.text(196, Am[-1] + At[-1]/2, 'TCO', ha='right', va='center', fontsize=FS_NOTE, color='0.2', fontweight='bold')
    a_.text(0.04, 0.95, name, transform=a_.transAxes, va='top', fontsize=FS_NOTE + 0.4, color=c, fontweight='bold')
    a_.text(0.96, 0.95, '$p$ = 0.38', transform=a_.transAxes, va='top', ha='right',
            fontsize=FS_NOTE, color='0.35')
    a_.set_xlim(30, 200); a_.set_ylim(0, ymax)
    a_.set_xlabel('ITO thickness (nm)'); a_.set_ylabel('round-trip loss  A′')
    clean(a_); letter(k)
    r = a_.twinx()
    r.set_ylim(0, ymax)
    r.set_yticks([(P_AB/e - P_AB)/(1 - P_AB) for e in eticks])
    r.set_yticklabels([f'{e:.2f}' for e in eticks])
    r.set_ylabel('$\\eta_{ext}$', labelpad=1)
    r.tick_params(labelsize=FS_TICK, pad=1.5)
    r.spines['top'].set_visible(False)

# ---- (c): the substrate-index sweep, three bands in one panel -------------
c_ = ax['c']
c_.plot(nb, eAg, color=C_AG, lw=1.9); c_.plot(nb, eAl, color=C_AL, lw=1.9)
c_.plot(nb, pb, color='0.45', lw=1.5, ls=DASH)
c_.plot(nb, aAl, color=C_AL, lw=1.0); c_.plot(nb, aAg, color=C_AG, lw=1.0)
c_.text(1.315, 0.865, '$\\eta_{ext}$', color='0.25', fontsize=FS_NOTE + 1.0, va='center')
c_.text(1.315, 0.470, '$p$', color='0.45', fontsize=FS_NOTE + 1.0, va='center')
c_.text(1.315, 0.195, "A′", color='0.25', fontsize=FS_NOTE + 1.0, va='center')
c_.text(1.72, eAl[9] - 0.05, 'Al', color=C_AL, fontsize=FS_NOTE + 0.4, fontweight='bold', ha='center')
c_.text(1.72, eAg[9] + 0.022, 'Ag', color=C_AG, fontsize=FS_NOTE + 0.4, fontweight='bold', ha='center')
c_.set_xlim(1.3, 2.0); c_.set_ylim(0, 1.0)
c_.set_xlabel('$n_{sub}$'); c_.set_ylabel('$\\eta_{ext}$,   $p$,   A′')
clean(c_); letter('c')

# ---- (d): substrate-delivered power and the product -----------------------
d_ = ax['d']
d_.fill_between(nb, qAl, sAl, color=C_AL, alpha=0.16, lw=0)
d_.fill_between(nb, qAg, sAg, color=C_AG, alpha=0.16, lw=0)
d_.plot(nb, sAl, color=C_AL, lw=0.8, alpha=0.7); d_.plot(nb, sAg, color=C_AG, lw=0.8, alpha=0.7)
d_.plot(nb, qAl, color=C_AL, lw=1.9); d_.plot(nb, qAg, color=C_AG, lw=1.9)
for q, c in ((qAg, C_AG), (qAl, C_AL)):
    i = int(np.argmax(q))
    d_.plot(nb[i], q[i], 'o', ms=3.6, mfc=c, mec='w', mew=0.8, zorder=5)
    d_.annotate(f'{q[i]:.2f}', (nb[i], q[i]), xytext=(5, -9), textcoords='offset points',
                fontsize=FS_NOTE, fontweight='bold', color=c)
d_.axvline(1.8, color='0.85', lw=0.7, ls=':', zorder=1)
d_.annotate('never escapes', (1.62, 0.5*(sAl[6] + qAl[6])), xytext=(1.315, 0.04),
            textcoords='data', ha='left', va='bottom', fontsize=FS_NOTE, color=C_AL,
            arrowprops=dict(arrowstyle='->', color=C_AL, lw=0.7, connectionstyle='arc3,rad=-0.2'))
d_.text(1.83, 0.17, '$n_{sub}$=$n_{EML}$', fontsize=FS_NOTE, color='0.45', ha='left')
d_.legend(handles=[Line2D([], [], color='0.45', lw=0.9, label='$\\eta_{sub}^{(0)}$'),
                   Line2D([], [], color='0.45', lw=1.9, label='EQE')],
          loc='upper left', fontsize=FS_NOTE, frameon=False, handlelength=1.4, labelspacing=0.25)
d_.set_xlim(1.3, 2.0); d_.set_ylim(0, 1.0)
d_.set_xlabel('$n_{sub}$'); d_.set_ylabel('$\\eta_{sub}^{(0)}$,   EQE')
clean(d_); letter('d')

# ---- (e)-(h): the three reflectors ---------------------------------------
VIS = np.arange(430.0, 701.0); NSUB = 1.5
THC = np.degrees(np.arcsin(1/NSUB))
NAMES = [('Al', C_AL), ('Ag', C_AG), ('TCO + DBR', C_DBR)]
TITLE = {'TCO + DBR': 'ITO 50 nm + DBR'}
W = {n: F.weighted(n, n_sub=NSUB, spectrum=F.GREEN, lam=VIS) for n, _ in NAMES}
for k, (name, c) in zip('efg', NAMES):
    th, lam, R, T = F.maps(name, n_sub=NSUB, lam=VIS)
    im = ax[k].pcolormesh(np.degrees(th), lam, 100*(1 - R), cmap='inferno_r', vmin=0, vmax=30,
                          shading='auto', rasterized=True)
    ax[k].set_title(TITLE.get(name, name), fontsize=FS_NOTE + 0.5, pad=2.5, color=c)
    ax[k].set_xlabel('$\\theta$ in substrate (°)')
    ax[k].set_xticks([0, 30, 60, 90])
    ax[k].set_xticklabels(['0', '30', '60', ''] if k != 'g' else ['0', '30', '60', '90'])
    ax[k].tick_params(labelsize=FS_TICK, pad=1.5)
    if k != 'e': ax[k].set_yticklabels([])
    w = W[name]
    txt = (f"absorbed {100*w['absorbed']:.1f} %" if w['transmitted'] < 0.002 else
           f"absorbed {100*w['absorbed']:.1f} %\nleaked {100*w['transmitted']:.1f} %")
    ax[k].text(0.5, -0.33, txt, transform=ax[k].transAxes, ha='center', va='top',
               fontsize=FS_NOTE, color=c)
ax['g'].axvline(THC, color='w', lw=0.8, ls=(0, (2.5, 1.5)))
ax['e'].set_ylabel('wavelength (nm)')
for k in 'efg': letter(k, dx=-0.40 if k == 'e' else -0.12)
cb = fig.colorbar(im, cax=ax_cb); cb.set_label('1 − R  (%)', fontsize=FS_LAB, labelpad=2)
cb.ax.tick_params(labelsize=FS_TICK, pad=1.5); cb.outline.set_linewidth(0.6)

h_ = ax['h']
th = np.linspace(0, np.pi/2, 452); th = 0.5*(th[1:] + th[:-1])
wt = np.cos(th)*np.sin(th); wt /= wt.max()
h_.fill_between(np.degrees(th), 0, 6*wt, color='0.91', lw=0, zorder=0)
h_.text(46, 0.3, '$\\cos\\theta\\sin\\theta$ weight', fontsize=FS_NOTE, color='0.5', ha='center')
for name, c in NAMES:
    t2, lam, R, T = F.maps(name, n_sub=NSUB, lam=VIS)
    sw = np.clip(np.interp(lam, F.LAM, F.GREEN), 0, None); sw /= sw.sum()
    h_.plot(np.degrees(t2), 100*(1 - (R*sw[:, None]).sum(0)), color=c, lw=1.5)
h_.axvline(THC, color=C_DBR, lw=0.7, ls=':', zorder=1)
h_.text(89, 12.5, 'the dielectric mirror leaks only\ninside the escape cone to air',
        fontsize=FS_NOTE, color=C_DBR, ha='right', va='top')
h_.set_xlim(0, 90); h_.set_ylim(0, 31); h_.set_xticks([0, 30, 60, 90])
h_.set_xlabel('$\\theta$ in substrate (°)'); h_.set_ylabel('round-trip loss 1 − R  (%)')
clean(h_); letter('h', dx=-0.18)
h_.legend(handles=[Line2D([], [], color=c, lw=1.5,
                          label=f"{TITLE.get(n, n)}   {100*W[n]['loss']:.1f} %") for n, c in NAMES],
          loc='upper left', bbox_to_anchor=(0.01, 0.99), fontsize=FS_NOTE, frameon=False,
          handlelength=1.6, labelspacing=0.3)

fig.savefig('fig2_compact.png', dpi=320)
fig.savefig('fig2_compact.pdf')
print('saved')
