"""Fig. 2(a): ITO thickness vs round-trip loss (mirror vs TCO) and extraction efficiency, Al vs Ag.
Data: kn2_{al,ag}_15.csv (Konig ITO) and kn3_{al,ag}_15_k0.csv (same stack at k_TCO = 0),
n_sub = 1.5, p = 0.38.  Set PRE = 'tco' to redraw the earlier k_TCO = 0.02 version."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import textwrap, os
PRE2, PRE3 = os.environ.get('PRE2', 'kn2'), os.environ.get('PRE3', 'kn3')
C_MIR, C_TCO = '#6E6E6E', '#E69F00'
C_AL, C_AG = '#555555', '#0072B2'

def load(m):
    F = np.loadtxt(f'{PRE2}_{m}_15.csv', delimiter=','); Z = np.loadtxt(f'{PRE3}_{m}_15_k0.csv', delimiter=',')
    return F[:, 0], Z[:, 2], F[:, 2] - Z[:, 2], F[:, 8]
al, ag = load('al'), load('ag')
YMAX = 0.19

fig, axs = plt.subplots(1, 3, figsize=(14.4, 5.1), gridspec_kw={'width_ratios': [1, 1, 1.15]})
fig.subplots_adjust(left=0.055, right=0.985, top=0.83, bottom=0.31, wspace=0.30)
for ax, (d, Am, At, _), name in zip(axs[:2], (al, ag), ('Al', 'Ag')):
    ax.fill_between(d, 0, Am, color=C_MIR, alpha=0.85, lw=0)
    ax.fill_between(d, Am, Am + At, color=C_TCO, alpha=0.85, lw=0)
    ax.plot(d, Am + At, color='k', lw=1.2)
    if Am[-1] > 0.05:
        ax.text(d[-1] - 4, Am[-1]/2, f'mirror ohmic\n{Am[-1]:.3f}', ha='right', va='center',
                fontsize=9, color='w', fontweight='bold')
        ax.text(d[-1] - 4, Am[-1] + At[-1]/2, f'TCO absorption\n{At[-1]:.3f}', ha='right', va='center',
                fontsize=9, color='k', fontweight='bold')
    j, k = np.argmin(abs(d - 50)), np.argmin(abs(d - 150))
    ax.annotate(f'TCO share of A′  {At[j]/(Am[j]+At[j]):.0%} → {At[k]/(Am[k]+At[k]):.0%}\n(50 → 150 nm)',
                xy=(0.04, 0.96), xycoords='axes fraction', ha='left', va='top', fontsize=9,
                bbox=dict(boxstyle='round,pad=0.3', fc='w', ec='0.8'))
    ax.set_title(f'{name} reflector:  round-trip loss A′ = 1 − ⟨R_LED⟩', fontsize=10.5, loc='left')
    ax.set_xlabel('ITO thickness (nm)'); ax.set_ylabel('round-trip loss  A′')
    ax.set_xlim(30, 200); ax.set_ylim(0, YMAX)
    ax.grid(axis='y', color='0.9', lw=0.8)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)

# the Ag stack is a fifth of the Al stack, so the split needs a zoom to stay legible
d, Am, At, _ = ag
ins = axs[1].inset_axes([0.46, 0.30, 0.50, 0.42])
ins.fill_between(d, 0, Am, color=C_MIR, alpha=0.85, lw=0)
ins.fill_between(d, Am, Am + At, color=C_TCO, alpha=0.85, lw=0)
ins.plot(d, Am + At, color='k', lw=1.0)
ins.set_xlim(30, 200); ins.set_ylim(0, 0.07); ins.set_yticks([0, 0.03, 0.06])
ins.tick_params(labelsize=7.5, length=3, pad=2)
ins.text(0.5, 0.94, '×2.7 zoom', transform=ins.transAxes, ha='center', va='top', fontsize=7.5, color='0.35')
ins.text(196, Am[-1]/2, f'mirror  {Am[-1]:.3f}', ha='right', va='center', fontsize=7.5, color='w', fontweight='bold')
ins.text(196, Am[-1] + At[-1]/2, f'TCO  {At[-1]:.3f}', ha='right', va='center', fontsize=7.5, color='k', fontweight='bold')
for s in ('top', 'right'): ins.spines[s].set_visible(False)
ins.set_facecolor('white')

ax = axs[2]
for (d, Am, At, ext), name, c in ((al, 'Al', C_AL), (ag, 'Ag', C_AG)):
    ax.plot(d, ext, color=c, lw=2.8, solid_capstyle='round')
    j, k = np.argmin(abs(d - 50)), np.argmin(abs(d - 150))
    ax.text(0.97, 0.95 if name == 'Ag' else 0.06,
            f'{name}:  {ext[j]:.3f} → {ext[k]:.3f}  ({100*(ext[k]-ext[j]):+.1f} pp)',
            transform=ax.transAxes, ha='right', va='top' if name == 'Ag' else 'bottom',
            fontsize=9.5, color=c, fontweight='bold')
for x in (50, 150): ax.axvline(x, color='0.85', lw=1, ls=':')
ax.set_title('Substrate-to-air extraction efficiency', fontsize=10.5, loc='left')
ax.set_xlabel('ITO thickness (nm)'); ax.set_ylabel('$\\eta_{ext}$')
ax.set_xlim(30, 200); ax.set_ylim(0.70, 1.0); ax.grid(axis='y', color='0.9', lw=0.8)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.plot(150, 0.916, '*', ms=15, mfc='#D55E00', mec='w', mew=1.2, zorder=6)
ax.annotate('measured, green device\n(150 nm ITO, Ag):  0.916', (150, 0.916),
            xytext=(-12, -40), textcoords='offset points', ha='right', fontsize=8.5, color='0.35',
            arrowprops=dict(arrowstyle='->', color='0.6', lw=1.0))

fig.suptitle('Fig. 2(a)  With Al the mirror is the whole loss; once it is low-loss, the transparent electrode takes over',
             x=0.055, ha='left', fontsize=12.5, y=0.955)
foot = ('Generic stack, 550 nm, isotropic dipole, PLQY = 1: reflector (Al, JO n,k / Ag, McPeak n,k; 100 nm) / ETL 200 nm / EML 20 nm / HTL 200 nm / '
        'ITO (n = 1.864 + 0.0032i at 550 nm, Koenig et al. 2014; thickness swept) / glass n = 1.5 with a microlens film, p = 0.38.  '
        'A′ is the flux-weighted reflectance deficit of the OLED stack seen from the substrate; its mirror part is A′ recomputed with k_TCO = 0 and the TCO part is the remainder.  '
        'η_ext = p/[p + (1−p)A′].  The Ag mirror loss is 0.018 against 0.134 for Al, so the TCO — the same film in both stacks — goes from a rounding error to the main loss path.  '
        'Companion numbers at n_sub = 1.8 (p = 0.30): η_ext falls 0.741 → 0.721 for Al and 0.933 → 0.894 for Ag over the same 50 → 150 nm.')
fig.text(0.055, 0.02, '\n'.join(textwrap.wrap(foot, 200)), fontsize=8.2, va='bottom', ha='left', color='0.28')
fig.savefig('fig2a_mock.png', dpi=150); print('saved fig2a')
