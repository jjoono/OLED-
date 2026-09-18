"""Fig 3(c) mock: 2-D map of the substrate-delivered power eta_sub^(0)(d_ETL, n_sub)
for an Ag reflector, from fig3c_Ag_ne18.csv / fig3c_Ag_ne15.csv (columns:
n_sub, d_ETL, air, sub_confined, wg, spp, abs; produced by fig3c_Ag.m)."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
import textwrap

def load(fn):
    d = np.loadtxt(fn, delimiter=',')
    ns, dE = np.unique(d[:, 0]), np.unique(d[:, 1])
    eta = np.full((len(ns), len(dE)), np.nan); spp = eta.copy(); wg = eta.copy(); ab = eta.copy()
    for r in d:
        i, j = np.searchsorted(ns, r[0]), np.searchsorted(dE, r[1])
        eta[i, j] = r[2] + r[3]; wg[i, j] = r[4]; spp[i, j] = r[5]; ab[i, j] = r[6]
    return ns, dE, eta, wg, spp, ab

cases = [('fig3c_Ag_ne18.csv', '(c1) isotropic ETL:  n_e,ETL = n_EML = 1.8'),
         ('fig3c_Ag_ne15.csv', '(c2) low-index ETL:  n_e,ETL = 1.5  (n_o = 1.8)')]
N_EML = 1.8
D_ELEC = 200          # nm, start of the driving-voltage band (undoped ETL)
devices = [  # (d_ETL, n_sub, marker, label)  -- positions indicative, see footer
    (55, 1.50, '*', 'Ag device (glass + n≈1.5 MLA film,\nd_ETL ≈ 55 nm): η_ext = 91.6 % measured'),
    (50, 1.77, 'o', 'DBR device (n = 1.77 MLA,\nd_ETL = [[x]], Al cathode)')]

fig, axes = plt.subplots(1, 2, figsize=(13.6, 6.2), sharey=True)
fig.subplots_adjust(left=0.06, right=0.885, top=0.90, bottom=0.25, wspace=0.08)
norm = Normalize(0.1, 1.0)
levels = np.linspace(0.1, 1.0, 37)
for ax, (fn, title) in zip(axes, cases):
    ns, dE, eta, wg, spp, ab = load(fn)
    cf = ax.contourf(dE, ns, np.clip(eta, 0.1, 1.0), levels=levels, cmap='viridis', norm=norm, extend='min')
    cs = ax.contour(dE, ns, eta, levels=[0.5, 0.7, 0.8, 0.9], colors='w', linewidths=[0.9, 0.9, 0.9, 1.8])
    # label the 0.5/0.7/0.8 contours on their horizontal branches (right-hand side); the thick line is 0.9
    pos = []
    for v, x in zip((0.5, 0.7, 0.8), (255, 325, 385)):          # spread the labels along x
        j = np.argmin(abs(dE - x)); col = eta[:, j]
        i = np.argmax(col >= v) - 1                                  # first crossing from below (eta rises with n_sub here)
        y = ns[i] + (v - col[i]) / (col[i + 1] - col[i]) * (ns[i + 1] - ns[i]) if 0 <= i < len(ns) - 1 else ns[np.argmin(abs(col - v))]
        pos.append((x, y))
    ax.clabel(cs, levels=[0.5, 0.7, 0.8], fmt=lambda v: f'{v:.1f}', fontsize=8.5, inline=True, inline_spacing=3,
              manual=pos)
    # --- horizontal transition: n_sub = n_EML (waveguide cut-off) ---
    ax.axhline(N_EML, color='w', ls=':', lw=1.5)
    i18 = np.searchsorted(ns, 1.8)
    thr = dE[np.argmax(spp[i18] < 0.05)]
    ax.text(8 if thr > 150 else thr + 8, N_EML + 0.012, 'n_sub = n_EML : WG modes vanish above this line', color='w',
            ha='left', va='bottom', fontsize=8.6, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.15', fc='k', ec='none', alpha=0.35))
    # --- vertical transition: d_ETL where the u > 1 (SPP bin) power < 5 % at n_sub = 1.8 ---
    ax.axvline(thr, color='w', ls=':', lw=1.5)
    ax.text(thr + 6, 1.395, f'SPP < 5 %  at\nd_ETL ≈ {thr:.0f} nm', color='w', fontsize=8.6,
            va='bottom', ha='left', fontweight='bold')
    # --- electrical band ---
    ax.add_patch(Rectangle((D_ELEC, 1.375), 400 - D_ELEC, 0.65, facecolor='none', hatch='///',
                           edgecolor='k', lw=0, alpha=0.35))
    ax.text(398, 1.395, 'driving-voltage\npenalty (undoped)', color='k', fontsize=8.6, ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.25', fc='w', ec='none', alpha=0.75))
    # --- design region (top-right corner) ---
    ax.add_patch(Rectangle((thr, N_EML), 400 - thr, 2.025 - N_EML, facecolor='none', edgecolor='#ff7f0e',
                           lw=2.0, ls='--'))
    eta_corner = eta[ns >= 1.8][:, dE >= max(thr, D_ELEC)]
    ax.text(398, 2.005, f'design region\nη_sub ≈ {eta_corner.min():.2f}–{eta_corner.max():.2f}',
            color='#ff7f0e', fontsize=9, ha='right', va='top', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.25', fc='k', ec='none', alpha=0.45))
    # --- measured devices ---
    for x, y, m, lab in devices:
        ax.plot(x, y, marker=m, ms=13 if m == '*' else 9, mfc='w' if m == 'o' else '#ff7f0e', mec='k', mew=1.2,
                ls='none', zorder=5)
    ax.set_title(title, fontsize=11, loc='left')
    ax.set_xlabel('d_ETL (nm)   (dipole–metal distance = d_ETL + 10 nm)')
    ax.set_xlim(5, 400); ax.set_ylim(1.375, 2.025)
    ax.set_xticks([50, 100, 150, 200, 250, 300, 350, 400])
axes[0].set_ylabel('n_sub  (substrate / MLA refractive index)')
# device labels on the left panel only
axes[0].annotate(devices[0][3], xy=(55, 1.50), xytext=(120, 1.455), fontsize=8.3, color='k',
                 arrowprops=dict(arrowstyle='-', color='k', lw=0.8),
                 bbox=dict(boxstyle='round,pad=0.25', fc='w', ec='none', alpha=0.85))
axes[0].annotate(devices[1][3], xy=(50, 1.77), xytext=(110, 1.66), fontsize=8.3, color='k',
                 arrowprops=dict(arrowstyle='-', color='k', lw=0.8),
                 bbox=dict(boxstyle='round,pad=0.25', fc='w', ec='none', alpha=0.85))
cax = fig.add_axes([0.90, 0.25, 0.017, 0.65])
cb = fig.colorbar(cf, cax=cax)
cb.set_label('η_sub^(0):  power delivered to the substrate\n(air + substrate-confined; PLQY = 1)', fontsize=9.5)
cb.set_ticks([0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0])
fig.suptitle('(c)  Combining the two conditions: substrate-delivered power versus (d_ETL, n_sub), Ag reflector',
             x=0.06, ha='left', fontsize=12.5, y=0.97)
foot = ('Generic model, 550 nm, isotropic dipole, PLQY = 1: Ag (McPeak n,k; 100 nm) / ETL (n_o = 1.8, n_e as stated, '
        'thickness d_ETL) / EML 1.8 (20 nm, dipole at centre) / HTL 1.8 (50 nm) / TCO 1.8+0.02i (50 nm, index-matched) / '
        'substrate n_sub.  Colour: η_sub^(0) = 1 − WG − SPP − absorption; white contours 0.5/0.7/0.8 and (thick) 0.9.  '
        'Dotted lines: horizontal transition n_sub = n_EML (WG cut-off) and vertical transition d_ETL(SPP < 5 % at n_sub = 1.8).  '
        'Hatched: d_ETL > 200 nm (driving-voltage penalty for an undoped ETL).  '
        'Orange dashed box: design region (index ladder n_e,ETL < n_EML ≲ n_sub completed AND d_ETL above threshold).  '
        'Markers: measured devices, positions indicative only ([[confirm d_ETL; the DBR device has an Al cathode]]).')
fig.text(0.06, 0.02, '\n'.join(textwrap.wrap(foot, 190)), fontsize=8.3, va='bottom', ha='left', color='0.25')
fig.savefig('fig3c_Ag_mock.png', dpi=150)
print('saved c')
