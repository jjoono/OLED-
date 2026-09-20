"""Fig. 2(b): why the substrate index alone does not buy EQE, and what removes the ceiling.
Data: fig2b_curves.csv (n_sub, p, then eta_sub / A' / eta_ext / EQE for Al and for Ag)."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import textwrap

D = np.genfromtxt('fig2b_curves.csv', delimiter=',', names=True)
n, p = D['n_sub'], D['p']
C_SUB, C_EXT, C_EQE, C_P = '#E69F00', '#0072B2', '#9467BD', '#999999'

fig, axs = plt.subplots(1, 2, figsize=(11.6, 5.0), sharey=True)
fig.subplots_adjust(left=0.075, right=0.985, top=0.86, bottom=0.30, wspace=0.10)
for ax, tag, name in ((axs[0], 'Al', 'Al electrode  (lossy mirror)'),
                      (axs[1], 'Ag', 'Ag electrode  (low-loss mirror)')):
    es, ee, eq = D[f'eta_sub_{tag}'], D[f'eta_ext_{tag}'], D[f'EQE_{tag}']
    ax.plot(n, es, color=C_SUB, lw=2.4, solid_capstyle='round')
    ax.plot(n, ee, color=C_EXT, lw=2.4, solid_capstyle='round')
    ax.plot(n, eq, color=C_EQE, lw=3.2, solid_capstyle='round')
    ax.plot(n, p, color=C_P, lw=1.4, ls=(0, (4, 2)))
    i = int(np.argmax(eq))
    ax.plot(n[i], eq[i], 'o', ms=8, mfc=C_EQE, mec='w', mew=1.5, zorder=5)
    ax.annotate(f'{eq[i]:.2f}', (n[i], eq[i]), xytext=(-6, -16), textcoords='offset points',
                ha='right', va='top', color=C_EQE, fontsize=10, fontweight='bold')
    ax.set_title(name, fontsize=11.5, loc='left')
    ax.set_xlabel('Refractive index of the substrate / MLA')
    ax.set_xlim(1.3, 2.0); ax.set_ylim(0, 1.05)
    ax.grid(axis='y', color='0.92', lw=0.8)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
axs[0].set_ylabel('Efficiency  /  power fraction')

from matplotlib.lines import Line2D
handles = [Line2D([], [], color=C_SUB, lw=2.4, label='$\\eta_{sub}^{(0)}$   delivered to the substrate'),
           Line2D([], [], color=C_EXT, lw=2.4, label='$\\eta_{ext}$   substrate → air'),
           Line2D([], [], color=C_EQE, lw=3.2, label='EQE = $\\eta_{sub}^{(0)}\\eta_{ext}$'),
           Line2D([], [], color=C_P, lw=1.4, ls=(0, (4, 2)), label='$p$   single-pass escape')]
axs[0].legend(handles=handles, loc='lower right', fontsize=9, frameon=True, framealpha=0.95,
              edgecolor='0.85', borderpad=0.6, handlelength=2.2)
axs[0].annotate('η$_{sub}$ rises but η$_{ext}$ falls:\nEQE turns over', (1.80, 0.543),
                xytext=(-104, 30), textcoords='offset points', ha='left', va='bottom',
                fontsize=9, color='0.3',
                arrowprops=dict(arrowstyle='->', color='0.55', lw=1.0))
axs[1].annotate('with the round-trip loss suppressed,\nthe ceiling is gone', (1.86, 0.742),
                xytext=(-24, -70), textcoords='offset points', ha='center', fontsize=9, color='0.3',
                arrowprops=dict(arrowstyle='->', color='0.55', lw=1.0))
fig.suptitle('(b)  raising the substrate index only pays off once the round-trip loss is small',
             x=0.075, ha='left', fontsize=12.5, y=0.955)
foot = ('Full model at 550 nm, isotropic dipole, PLQY = 1: reflector 100 nm / ETL 200 nm / EML 20 nm / HTL 200 nm / ITO 50 nm (n = 1.9 + 0.02i) / substrate, '
        'substrate index swept 1.30–2.00 and index-matched to the outcoupling structure.  '
        'η_ext = p/[p + (1−p)A′] with A′ = 1 − ⟨R_LED⟩ from the model and p the single-pass escape probability of the outcoupling structure (grey dashed), '
        'which falls roughly as 1/n_sub².  EQE = η_sub^(0) η_ext.  '
        'The Al device peaks at 0.54 near n_sub = 1.8 and then flattens, because the gain in η_sub is spent on the loss in η_ext; the Ag device reaches 0.75 and holds it.')
fig.text(0.075, 0.02, '\n'.join(textwrap.wrap(foot, 185)), fontsize=8.3, va='bottom', ha='left', color='0.28')
fig.savefig('fig2b_mock.png', dpi=150); print('saved')
