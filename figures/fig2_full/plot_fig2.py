# -*- coding: utf-8 -*-
"""Figure 2, composed as three rows of panels.

a-c  ITO thickness: the round-trip loss split into mirror and TCO, for Al and for Ag, and
     the resulting extraction efficiency
d-f  substrate index: the two inputs of eq. (2), the extraction efficiency they give, and
     the product with the substrate-delivered power
g-j  the three reflectors at a fixed 150 nm ITO, resolved in angle and wavelength

Colour is the reflector throughout: vermillion Al, blue Ag, purple DBR; grey is reserved
for the TCO part of the loss."""
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

A = os.path.join(ROOT, 'figures', 'fig2a_mock')
Bd = os.path.join(ROOT, 'figures', 'fig2b_mock')
def a_load(m):
    T = np.loadtxt(f'{A}/kn2_{m}_15.csv', delimiter=','); Z = np.loadtxt(f'{A}/kn3_{m}_15_k0.csv', delimiter=',')
    return T[:, 0], Z[:, 2], T[:, 2] - Z[:, 2], T[:, 8]
al, ag = a_load('al'), a_load('ag')
D = np.genfromtxt(f'{Bd}/fig2b_curves_konig.csv', delimiter=',', names=True)
nb, pb = D['n_sub'], D['p']
eAl, eAg, sAl, sAg, qAl, qAg = (D['eta_ext_Al'], D['eta_ext_Ag'], D['eta_sub_Al'],
                               D['eta_sub_Ag'], D['EQE_Al'], D['EQE_Ag'])
aAl, aAg = D['Aprime_Al'], D['Aprime_Ag']

fig = plt.figure(figsize=(7.2, 6.7))
H = 0.225
Y = [0.740, 0.435, 0.120]
x3 = [0.078, 0.398, 0.718]; w3 = 0.252
ax = {}
for k, i in zip('abc', range(3)): ax[k] = fig.add_axes([x3[i], Y[0], w3, H])
for k, i in zip('def', range(3)): ax[k] = fig.add_axes([x3[i], Y[1], w3, H])
xm = [0.078, 0.246, 0.414]; wm = 0.150
for k, i in zip('ghi', range(3)): ax[k] = fig.add_axes([xm[i], Y[2], wm, H])
ax_cb = fig.add_axes([0.572, Y[2], 0.008, H])
ax['j'] = fig.add_axes([0.700, Y[2], 0.272, H])

def clean(a_, grid=True):
    for s in ('top', 'right'): a_.spines[s].set_visible(False)
    a_.tick_params(labelsize=FS_TICK, pad=1.5)
    if grid: a_.grid(axis='y', color='0.93', lw=0.6)
def letter(k, dx=-0.235, dy=1.17):
    ax[k].text(dx, dy, k, transform=ax[k].transAxes, fontsize=FS_LET, fontweight='bold', va='top')

# ================= row 1: ITO thickness ===================================
for k, (d, Am, At, _), c, name, ymax in (('a', al, C_AL, 'Al reflector', 0.19),
                                         ('b', ag, C_AG, 'Ag reflector', 0.062)):
    a_ = ax[k]
    a_.fill_between(d, 0, Am, color=c, alpha=0.9, lw=0)
    a_.fill_between(d, Am, Am + At, color=C_TCO, alpha=0.9, lw=0)
    a_.plot(d, Am + At, color='k', lw=0.8)
    a_.text(196, Am[-1]/2, 'mirror', ha='right', va='center', fontsize=FS_NOTE, color='w', fontweight='bold')
    a_.text(196, Am[-1] + At[-1]/2, 'TCO', ha='right', va='center', fontsize=FS_NOTE,
            color='0.2', fontweight='bold')
    j, kk = np.argmin(abs(d - 50)), np.argmin(abs(d - 150))
    a_.text(0.04, 0.95, f'{name}\nTCO {At[j]/(Am[j]+At[j]):.0%} → {At[kk]/(Am[kk]+At[kk]):.0%} of A′',
            transform=a_.transAxes, va='top', fontsize=FS_NOTE, color='0.25')
    a_.set_xlim(30, 200); a_.set_ylim(0, ymax)
    a_.set_xlabel('ITO thickness (nm)'); a_.set_ylabel('round-trip loss  A′')
    clean(a_); letter(k)
ax['b'].text(0.97, 0.97, 'note the axis:\nthree times finer', transform=ax['b'].transAxes,
             ha='right', va='top', fontsize=FS_NOTE, color='0.45')

for (dd, _, _, ext), c in ((al, C_AL), (ag, C_AG)):
    ax['c'].plot(dd, ext, color=c, lw=1.8, solid_capstyle='round')
ax['c'].plot(150, 0.916, '*', ms=8, mfc=C_AL, mec='w', mew=0.8, zorder=6)
ax['c'].annotate('measured\n(green device)', (150, 0.916), xytext=(-5, -22),
                 textcoords='offset points', ha='right', fontsize=FS_NOTE, color='0.4',
                 arrowprops=dict(arrowstyle='->', color='0.6', lw=0.7))
ax['c'].text(196, ag[3][-1] + 0.012, 'Ag', color=C_AG, fontsize=FS_LAB, fontweight='bold', ha='right')
ax['c'].text(196, al[3][-1] + 0.012, 'Al', color=C_AL, fontsize=FS_LAB, fontweight='bold', ha='right')
ax['c'].set_xlim(30, 200); ax['c'].set_ylim(0.70, 1.0)
ax['c'].set_xlabel('ITO thickness (nm)'); ax['c'].set_ylabel('$\\eta_{ext}$')
clean(ax['c']); letter('c')

# ================= row 2: substrate index =================================
d_ = ax['d']
d_.plot(nb, pb, color='0.45', lw=1.7, ls=DASH)
d_.plot(nb, aAl, color=C_AL, lw=1.7); d_.plot(nb, aAg, color=C_AG, lw=1.7)
d_.set_xlim(1.3, 2.0); d_.set_ylim(0, 0.72)
d_.set_xlabel('$n_{sub}$'); d_.set_ylabel('$p$   and   A′')
d_.text(1.40, 0.545, '$p$', color='0.4', fontsize=FS_LAB + 1.2)
d_.text(1.97, aAl[-1] + 0.030, "Al,  A′", color=C_AL, fontsize=FS_NOTE + 0.4,
        fontweight='bold', ha='right')
d_.text(1.97, aAg[-1] + 0.020, "Ag,  A′", color=C_AG, fontsize=FS_NOTE + 0.4,
        fontweight='bold', ha='right')
d_.text(0.97, 0.97, 'the escape probability falls\nby 2.7×; the round-trip loss\nbarely moves',
        transform=d_.transAxes, ha='right', va='top', fontsize=FS_NOTE, color='0.3')
clean(d_); letter('d')

e_ = ax['e']
e_.fill_between(nb, eAl, eAg, color=C_AG, alpha=0.10, lw=0)
e_.plot(nb, eAl, color=C_AL, lw=1.9); e_.plot(nb, eAg, color=C_AG, lw=1.9)
e_.annotate('', (1.965, eAl[-2]), (1.965, eAg[-2]),
            arrowprops=dict(arrowstyle='<->', color='0.45', lw=0.7, shrinkA=0, shrinkB=0))
e_.text(1.945, 0.5*(eAl[-2] + eAg[-2]), '21 %p', ha='right', va='center', fontsize=FS_NOTE, color='0.35')
e_.text(1.66, eAg[7] + 0.02, 'Ag', color=C_AG, fontsize=FS_LAB, fontweight='bold', ha='center')
e_.text(1.66, eAl[7] - 0.05, 'Al', color=C_AL, fontsize=FS_LAB, fontweight='bold', ha='center')
e_.text(0.5, 0.06, 'the same loss costs more\nwhen $p$ is small', transform=e_.transAxes,
        ha='center', va='bottom', fontsize=FS_NOTE, color='0.3')
e_.set_xlim(1.3, 2.0); e_.set_ylim(0.6, 1.0)
e_.set_xlabel('$n_{sub}$'); e_.set_ylabel('$\\eta_{ext}$')
clean(e_); letter('e')

f_ = ax['f']
f_.fill_between(nb, qAl, sAl, color=C_AL, alpha=0.16, lw=0)
f_.fill_between(nb, qAg, sAg, color=C_AG, alpha=0.16, lw=0)
f_.plot(nb, sAl, color=C_AL, lw=0.8, alpha=0.7); f_.plot(nb, sAg, color=C_AG, lw=0.8, alpha=0.7)
f_.plot(nb, qAl, color=C_AL, lw=1.9); f_.plot(nb, qAg, color=C_AG, lw=1.9)
for q, c in ((qAg, C_AG), (qAl, C_AL)):
    i = int(np.argmax(q))
    f_.plot(nb[i], q[i], 'o', ms=3.6, mfc=c, mec='w', mew=0.8, zorder=5)
    f_.annotate(f'{q[i]:.2f}', (nb[i], q[i]), xytext=(5, -9), textcoords='offset points',
                fontsize=FS_NOTE, fontweight='bold', color=c)
f_.axvline(1.8, color='0.85', lw=0.7, ls=':', zorder=1)
f_.annotate('never escapes', (1.62, 0.5*(sAl[6] + qAl[6])), xytext=(1.315, 0.04),
            textcoords='data', ha='left', va='bottom', fontsize=FS_NOTE, color=C_AL,
            arrowprops=dict(arrowstyle='->', color=C_AL, lw=0.7, connectionstyle='arc3,rad=-0.2'))
f_.text(1.83, 0.17, '$n_{sub}$=$n_{EML}$', fontsize=FS_NOTE, color='0.45', ha='left')
f_.legend(handles=[Line2D([], [], color='0.45', lw=0.9, label='$\\eta_{sub}^{(0)}$'),
                   Line2D([], [], color='0.45', lw=1.9, label='EQE')],
          loc='upper left', fontsize=FS_NOTE, frameon=False, handlelength=1.4, labelspacing=0.25)
f_.set_xlim(1.3, 2.0); f_.set_ylim(0, 1.0)
f_.set_xlabel('$n_{sub}$'); f_.set_ylabel('$\\eta_{sub}^{(0)}$,   EQE')
clean(f_); letter('f')

# ================= row 3: the three reflectors at 150 nm ITO ==============
VIS = np.arange(430.0, 701.0); NSUB = 1.5
THC = np.degrees(np.arcsin(1/NSUB))
NAMES = [('Al', C_AL), ('Ag', C_AG), ('TCO + DBR', C_DBR)]
TITLE = {'TCO + DBR': 'ITO 50 nm + DBR\n(10 chirped pairs)'}
WMAP = {n: F.weighted(n, n_sub=NSUB, spectrum=F.GREEN, lam=VIS) for n, _ in NAMES}
for k, (name, c) in zip('ghi', NAMES):
    th, lam, R, T = F.maps(name, n_sub=NSUB, lam=VIS)
    im = ax[k].pcolormesh(np.degrees(th), lam, 100*(1 - R), cmap='inferno_r', vmin=0, vmax=30,
                          shading='auto', rasterized=True)
    ax[k].set_title(TITLE.get(name, name), fontsize=FS_NOTE + 0.5, pad=2.5, color=c)
    ax[k].set_xlabel('$\\theta$ in substrate (°)')
    ax[k].set_xticks([0, 30, 60, 90])
    ax[k].set_xticklabels(['0', '30', '60', ''] if k != 'i' else ['0', '30', '60', '90'])
    ax[k].tick_params(labelsize=FS_TICK, pad=1.5)
    if k != 'g': ax[k].set_yticklabels([])
    W = WMAP[name]
    txt = (f"absorbed {100*W['absorbed']:.1f} %" if W['transmitted'] < 0.002 else
           f"absorbed {100*W['absorbed']:.1f} %\nleaked {100*W['transmitted']:.1f} %")
    ax[k].text(0.5, -0.30, txt, transform=ax[k].transAxes, ha='center', va='top',
               fontsize=FS_NOTE, color=c)
ax['i'].axvline(THC, color='w', lw=0.8, ls=(0, (2.5, 1.5)))
ax['i'].annotate('leaks to air', (THC/2, 690), fontsize=FS_NOTE - 0.3, color='w',
                 ha='center', va='top')
ax['i'].annotate('absorption only', (THC + (90 - THC)/2, 690), fontsize=FS_NOTE - 0.3,
                 color='0.15', ha='center', va='top')
ax['g'].set_ylabel('wavelength (nm)')
for k in 'ghi': letter(k, dx=-0.36 if k == 'g' else -0.12)
cb = fig.colorbar(im, cax=ax_cb); cb.set_label('1 − R  (%)', fontsize=FS_LAB, labelpad=2)
cb.ax.tick_params(labelsize=FS_TICK, pad=1.5); cb.outline.set_linewidth(0.6)

j_ = ax['j']
th = np.linspace(0, np.pi/2, 452); th = 0.5*(th[1:] + th[:-1])
w = np.cos(th)*np.sin(th); w /= w.max()
j_.fill_between(np.degrees(th), 0, 6*w, color='0.91', lw=0, zorder=0)
j_.text(46, 0.3, '$\\cos\\theta\\sin\\theta$ weight', fontsize=FS_NOTE, color='0.5', ha='center')
W = {}
CURVES = NAMES
for name, c in CURVES:
    t2, lam, R, T = F.maps(name, n_sub=NSUB, lam=VIS)
    sw = np.clip(np.interp(lam, F.LAM, F.GREEN), 0, None); sw /= sw.sum()
    ls = '-'
    j_.plot(np.degrees(t2), 100*(1 - (R*sw[:, None]).sum(0)), color=c, lw=1.5, ls=ls)
    W[name] = F.weighted(name, n_sub=NSUB, spectrum=F.GREEN, lam=VIS)
j_.axvline(np.degrees(np.arcsin(1/NSUB)), color=C_DBR, lw=0.7, ls=':', zorder=1)
j_.text(89, 11.5, 'the dielectric mirror leaks only\ninside the escape cone to air',
        fontsize=FS_NOTE, color=C_DBR, ha='right', va='top')
j_.set_xlim(0, 90); j_.set_ylim(0, 31); j_.set_xticks([0, 30, 60, 90])
j_.set_xlabel('$\\theta$ in substrate (°)'); j_.set_ylabel('round-trip loss 1 − R  (%)')
clean(j_); letter('j', dx=-0.20)
lab = [('Al', f"Al   {100*W['Al']['loss']:.1f} %", C_AL, '-'),
       ('Ag', f"Ag   {100*W['Ag']['loss']:.1f} %", C_AG, '-'),
       ('TCO + DBR', f"ITO + DBR   {100*W['TCO + DBR']['loss']:.1f} %", C_DBR, '-')]
j_.legend(handles=[Line2D([], [], color=c, lw=1.5, ls=s, label=t) for _, t, c, s in lab],
          loc='upper left', bbox_to_anchor=(0.01, 0.99), fontsize=FS_NOTE, frameon=False,
          handlelength=1.6, labelspacing=0.3)

fig.savefig('fig2_full.png', dpi=320)
fig.savefig('fig2_full.pdf')
print('saved', {k: round(100*v['loss'], 1) for k, v in W.items()})
