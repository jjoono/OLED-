# -*- coding: utf-8 -*-
"""Fig 3(c): the surface plasmon has to be pushed past the substrate light line.

ci   where the dissipated power sits in the in-plane effective index k_x/k0.
     Left of n_sub it can reach the substrate, right of it it is bound.  A low
     out-of-plane index in the ETL walks the plasmon leftwards.
cii  the substrate-delivered power against the plasmon index, at three ETL
     thicknesses.  A step, not a slope -- and the measured ETLs sit on the
     wrong side of it or barely on the right side."""
import numpy as np, matplotlib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import nspp, materials as M

HERE = os.path.dirname(os.path.abspath(__file__))
D_SHOW, N_SUB = 150.0, 1.80
CF = [('n_e=1.80', 1.80, '#BDD7E7'), ('n_e=1.70', 1.70, '#6BAED6'),
      ('n_e=1.60', 1.60, '#2171B5'), ('n_e=1.50', 1.50, '#08306B')]
C_B3, C_B4 = '#D55E00', '#009E73'
CD = {60.0: '#C6C6C6', 100.0: '#7F7F7F', 150.0: '#1A1A1A'}
FS_LAB, FS_TICK, FS_NOTE, FS_LET = 7.2, 6.6, 6.1, 9.0
plt.rcParams.update({'font.size': FS_LAB, 'axes.linewidth': 0.7,
                     'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
                     'xtick.major.size': 2.4, 'ytick.major.size': 2.4})

fig, ax = plt.subplots(1, 2, figsize=(6.9, 2.7))

# ---------------- ci : where the power sits in k-space ----------------------
S = np.genfromtxt(os.path.join(HERE, 'fig3c_spectrum.csv'), delimiter=',',
                  names=True, dtype=None, encoding='utf-8')
S = S[S['d_ETL_nm'] == D_SHOW]
# The EML light line at k_x/k0 = n_EML = 1.80 is a branch point of the kernel
# (an integrable 1/sqrt divergence).  It is blanked out so that it does not
# compete visually with the plasmon peaks; the dashed line marks it instead.
S = S[np.abs(S['n_eff'] - 1.80) > 0.009]
a = ax[0]
a.axvspan(N_SUB, 2.6, color='#EFEFEF', zorder=0, lw=0)
a.axvline(N_SUB, color='0.35', lw=0.9, dashes=(2.6, 1.6), zorder=3)
for key, ne, c in CF:
    m = S[S['series'] == key]
    a.semilogy(m['n_eff'], m['TM'], color=c, lw=1.3, zorder=4)
for key, c in (('B3PyMPM', C_B3), ('B4PyMPM', C_B4)):
    m = S[S['series'] == key]
    a.semilogy(m['n_eff'], m['TM'], color=c, lw=1.2, dashes=(3.0, 1.5), zorder=5)
a.set_xlim(1.45, 2.35)
a.set_ylim(3e-3, 22)
a.set_xlabel(r'in-plane effective index   $k_x/k_0$')
a.set_ylabel('TM power dissipation density')
a.text(N_SUB - 0.04, 13, 'reaches the\nsubstrate', fontsize=FS_NOTE,
       color='0.4', ha='right', va='top', linespacing=1.15)
a.text(N_SUB + 0.04, 13, 'bound\n(SPP)', fontsize=FS_NOTE, color='0.4',
       ha='left', va='top', linespacing=1.15)
a.text(N_SUB - 0.02, 4.5e-3, r'$n_{\rm sub}$ = $n_{\rm EML}$ = 1.80',
       fontsize=FS_NOTE, color='0.35', ha='right', va='bottom', rotation=90)
a.legend(handles=[Line2D([], [], color=c, lw=1.3, label=r'$n_e$ = %.2f' % ne)
                  for _, ne, c in CF]
         + [Line2D([], [], color=C_B3, lw=1.2, dashes=(3, 1.5), label='B3PyMPM'),
            Line2D([], [], color=C_B4, lw=1.2, dashes=(3, 1.5), label='B4PyMPM')],
         loc='lower left', frameon=False, fontsize=FS_NOTE - 0.3,
         handlelength=1.7, labelspacing=0.26, borderaxespad=0.25)
a.set_title(r'$d_{\rm ETL}$ = %.0f nm,  $n_o$ = 1.80 (families)' % D_SHOW,
            fontsize=FS_NOTE, pad=3.5)

# ---------------- cii : the step ------------------------------------------
C = np.genfromtxt(os.path.join(HERE, 'fig3c_collapse.csv'), delimiter=',',
                  names=True, dtype=None, encoding='utf-8')
b = ax[1]
b.axvspan(N_SUB, 2.3, color='#EFEFEF', zorder=0, lw=0)
b.axvline(N_SUB, color='0.35', lw=0.9, dashes=(2.6, 1.6), zorder=3)
for d in (60.0, 100.0, 150.0):
    m = C[(C['series'] == 'n_o=1.80') & (C['d_ETL_nm'] == d)]
    o = np.argsort(m['n_SPP_Ag'])
    b.plot(m['n_SPP_Ag'][o], 100 * m['eta_sub'][o], color=CD[d], lw=1.4, zorder=4)
for key, c, dy in (('B3PyMPM', C_B3, 10), ('B4PyMPM', C_B4, 10)):
    m = C[(C['series'] == key) & (C['d_ETL_nm'] == D_SHOW)]
    b.plot(m['n_SPP_Ag'], 100 * m['eta_sub'], marker='o', ms=4.4, color=c,
           mew=0, ls='none', zorder=6)
    b.annotate('%s\n$n_e$ = %.2f' % (key, m['n_e'][0]),
               (m['n_SPP_Ag'][0], 100 * m['eta_sub'][0]),
               textcoords='offset points',
               xytext=(-5 if key == 'B4PyMPM' else 6, -3 if key == 'B4PyMPM' else 7),
               va='top' if key == 'B4PyMPM' else 'bottom',
               ha='right' if key == 'B4PyMPM' else 'left',
               fontsize=FS_NOTE - 0.4, color=c, linespacing=1.15)
b.set_xlim(1.56, 2.14)
b.set_ylim(56, 101)
b.set_xlabel(r'surface-plasmon index   $n_{\rm SPP}(n_o, n_e)$')
b.set_ylabel(r'into the substrate mode  $\eta_{\rm sub}$  (%)')
b.legend(handles=[Line2D([], [], color=CD[d], lw=1.4,
                         label=r'$d_{\rm ETL}$ = %.0f nm' % d)
                  for d in (60.0, 100.0, 150.0)],
         loc='lower right', frameon=False, fontsize=FS_NOTE,
         handlelength=1.7, labelspacing=0.26, borderaxespad=0.3,
         title='markers: measured ETLs, 150 nm')
b.get_legend().get_title().set_fontsize(FS_NOTE - 0.6)
b.get_legend().get_title().set_color('0.35')
sec = b.secondary_xaxis('top', functions=(
    lambda x: np.interp(x, *(lambda t: (nspp.n_spp(M.AG, 1.8, t), t))(
        np.linspace(1.30, 1.80, 400))),
    lambda e: nspp.n_spp(M.AG, 1.8, np.clip(e, 1.30, 1.80))))
sec.set_xlabel(r'$n_e$  at  $n_o$ = 1.80', fontsize=FS_NOTE, labelpad=2)
sec.tick_params(labelsize=FS_TICK - 0.3, width=0.7, size=2.2)

for k, a_ in enumerate(ax):
    a_.tick_params(labelsize=FS_TICK)
    a_.text(-0.20, 1.14, 'c' + ('i' if k == 0 else 'ii'), transform=a_.transAxes,
            fontsize=FS_LET, fontweight='bold', va='top')

fig.tight_layout(pad=0.4, w_pad=2.0)
for e in ('png', 'pdf'):
    fig.savefig(os.path.join(HERE, 'fig3c_mock.' + e), dpi=400)
print('fig3c_mock.png written')
