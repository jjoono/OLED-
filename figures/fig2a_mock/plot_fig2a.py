"""Fig. 2(a) mock: ITO thickness vs round-trip loss (mirror vs TCO) and extraction efficiency, Al vs Ag.
Data: tco2_{al,ag}_15.csv (k_TCO = 0.02) and tco3_{al,ag}_15_k0.csv (k_TCO = 0) — n_sub = 1.5, p = 0.38."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import textwrap
C_MIR, C_TCO = '#6E6E6E', '#E69F00'
C_AL, C_AG = '#555555', '#0072B2'
def load(m):
    F = np.loadtxt(f'tco2_{m}_15.csv', delimiter=','); Z = np.loadtxt(f'tco3_{m}_15_k0.csv', delimiter=',')
    d = F[:, 0]; Am = Z[:, 2]; At = F[:, 2] - Z[:, 2]; return d, Am, At, F[:, 8]
al = load('al'); ag = load('ag')

fig, axs = plt.subplots(1, 3, figsize=(14.4, 4.9), gridspec_kw={'width_ratios': [1, 1, 1.15]})
fig.subplots_adjust(left=0.055, right=0.985, top=0.83, bottom=0.30, wspace=0.30)
for ax, (d, Am, At, _), name in zip(axs[:2], (al, ag), ('Al', 'Ag')):
    ax.fill_between(d, 0, Am, color=C_MIR, alpha=0.85, lw=0)
    ax.fill_between(d, Am, Am + At, color=C_TCO, alpha=0.85, lw=0)
    ax.plot(d, Am + At, color='k', lw=1.2)
    i = -1
    if Am[i] > 0.05:
        ax.text(d[i] - 4, Am[i] / 2, f'mirror ohmic\n{Am[i]:.3f}', ha='right', va='center', fontsize=9, color='w', fontweight='bold')
    else:
        ax.text(34, Am[i] + 0.012, f'mirror ohmic  {Am[i]:.3f}', ha='left', va='bottom', fontsize=9, color=C_MIR, fontweight='bold')
    ax.text(d[i] - 4, Am[i] + At[i] / 2, f'TCO absorption\n{At[i]:.3f}', ha='right', va='center', fontsize=9, color='k', fontweight='bold')
    j = np.argmin(abs(d - 50)); k = np.argmin(abs(d - 150))
    ax.annotate(f'TCO share {At[j]/(Am[j]+At[j]):.0%} → {At[k]/(Am[k]+At[k]):.0%}\n(50 → 150 nm)',
                xy=(0.04, 0.96), xycoords='axes fraction', ha='left', va='top', fontsize=9,
                bbox=dict(boxstyle='round,pad=0.3', fc='w', ec='0.8'))
    ax.set_title(f'{name} reflector:  round-trip loss A′ = 1 − ⟨R_LED⟩', fontsize=10.5, loc='left')
    ax.set_xlabel('ITO thickness (nm)'); ax.set_ylabel("round-trip loss  A′")
    ax.set_xlim(30, 200); ax.set_ylim(0, 0.36)
    ax.grid(axis='y', color='0.9', lw=0.8)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax = axs[2]
for (d, Am, At, ext), name, c in ((al, 'Al', C_AL), (ag, 'Ag', C_AG)):
    ax.plot(d, ext, color=c, lw=2.8, solid_capstyle='round')
    j = np.argmin(abs(d - 50)); k = np.argmin(abs(d - 150))
    ax.text(0.97, 0.94 if name == 'Ag' else 0.06, f'{name}:  {ext[j]:.2f} → {ext[k]:.2f}  ({100*(ext[k]-ext[j]):+.1f} pp)',
            transform=ax.transAxes, ha='right', va='top' if name == 'Ag' else 'bottom', fontsize=9.5, color=c, fontweight='bold')
for x in (50, 150): ax.axvline(x, color='0.85', lw=1, ls=':')
ax.set_title('Substrate-to-air extraction efficiency', fontsize=10.5, loc='left')
ax.set_xlabel('ITO thickness (nm)'); ax.set_ylabel('$\\eta_{ext}$')
ax.set_xlim(30, 200); ax.set_ylim(0.6, 0.95); ax.grid(axis='y', color='0.9', lw=0.8)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
fig.suptitle('Fig. 2(a)  Once the mirror is low-loss, the transparent electrode is the only loss left',
             x=0.055, ha='left', fontsize=12.5, y=0.955)
foot = ('Generic stack, 550 nm, isotropic dipole, PLQY = 1: reflector (Al, JO n,k / Ag, McPeak n,k; 100 nm) / ETL 200 nm / EML 20 nm / HTL 200 nm / '
        'ITO (n = 1.9 + 0.02i, thickness swept) / glass n = 1.5 with a microlens film, p = 0.38.  '
        'A′ is the flux-weighted reflectance deficit of the OLED stack seen from the substrate; its mirror part is A′ recomputed with k_TCO = 0 and the TCO part is the remainder.  '
        'η_ext = p/[p + (1−p)A′].  Companion numbers at n_sub = 1.8 (p = 0.30): η_ext falls 0.675 → 0.600 for Al and 0.819 → 0.694 for Ag over the same 50 → 150 nm.')
fig.text(0.055, 0.02, '\n'.join(textwrap.wrap(foot, 200)), fontsize=8.2, va='bottom', ha='left', color='0.28')
fig.savefig('fig2a_mock.png', dpi=150); print('saved fig2a')
