"""Companion checks: (a,b) the TCO real index, (c,d) the ETL extraordinary index.
Reads p_ito_n{18,19,20}.csv and p_netl_*.csv (param, eta_sub, A', wg, spp, abs, eta_ext x2, EQE x2)."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import textwrap

C = ['#0072B2', '#E69F00', '#9467BD']
KMIN = 0.0025      # the k = 0 limit is not resolved for n_TCO = 2.0 (lossless TCO-guided pole)
def nspp(ne, no=1.8, nm=0.044, km=3.819):
    em = (nm+1j*km)**2
    return np.sqrt(em*ne**2*(em-no**2)/(em**2-no**2*ne**2)).real
GRID = np.linspace(1.30, 1.95, 800)
def ne_threshold(ns): return np.interp(ns, [nspp(x) for x in GRID], GRID)

fig, axs = plt.subplots(1, 4, figsize=(16.6, 5.4))
fig.subplots_adjust(left=0.05, right=0.99, top=0.855, bottom=0.315, wspace=0.30)

D = {n: np.loadtxt(f'p_ito_n{n}.csv', delimiter=',') for n in ('18', '19', '20')}
for c, n in zip(C, ('18', '19', '20')):
    d = D[n][D[n][:, 0] >= KMIN]; lab = f'$n_{{TCO}}$ = {n[0]}.{n[1]}'
    axs[0].plot(d[:, 0], d[:, 1], color=c, lw=2.4, label=lab)
    axs[0].plot(d[:, 0], d[:, 6], color=c, lw=1.5, ls=(0, (5, 2)))
    axs[1].plot(d[:, 0], d[:, 4], color=c, lw=2.4, label=lab)
axs[0].set_ylabel('Efficiency'); axs[0].set_ylim(0.55, 0.98)
axs[0].set_title('(a)  a higher-index TCO costs substrate power\n      while flattering $\\eta_{ext}$', fontsize=10.3, loc='left')
axs[0].legend(loc='upper right', fontsize=8.8, frameon=False)
axs[0].text(0.077, 0.585, 'solid  $\\eta_{sub}^{(0)}$\ndashed  $\\eta_{ext}$', ha='right', va='bottom', fontsize=9)
axs[1].set_ylabel('Power in the  $u>1$  bin'); axs[1].set_ylim(0, 0.095)
axs[1].set_title('(b)  above $n_{TCO}=n_{sub}$ that bin also holds\n      TCO-guided light, not only SPP', fontsize=10.3, loc='left')
axs[1].legend(loc='center left', fontsize=8.8, frameon=False)
for ax in axs[:2]:
    ax.set_xlim(0, 0.08); ax.set_xlabel('Extinction coefficient of the TCO,  $k_{TCO}$')

N = {1.8: np.loadtxt('p_netl_ito.csv', delimiter=','), 1.9: np.loadtxt('p_netl_ito_sub19.csv', delimiter=',')}
for c, ns in zip(C[:2], (1.8, 1.9)):
    d = N[ns]; thr = ne_threshold(ns)
    axs[2].plot(d[:, 0], d[:, 1], color=c, lw=2.4, label=f'$n_{{sub}}$ = {ns}')
    axs[3].plot(d[:, 0], d[:, 4], color=c, lw=2.4, label=f'$n_{{sub}}$ = {ns}')
    for k, ax in ((2, axs[2]), (3, axs[3])):
        ax.axvline(thr, color=c, lw=1.2, ls=':')
    axs[3].annotate(f'$n_{{SPP}}=n_{{sub}}$\nat $n_e$ = {thr:.2f}', (thr, 0.031 if ns == 1.8 else 0.016),
                    color=c, fontsize=8.6, ha='right' if ns == 1.8 else 'left', va='center',
                    xytext=(-5 if ns == 1.8 else 5, 0), textcoords='offset points')
for ax in axs[2:]:
    ax.axvline(1.6, color='#B00', lw=1.6)
    ax.set_xlim(1.40, 1.80); ax.set_xlabel('Extraordinary index of the ETL,  $n_e$   ($n_o$ = 1.8)')
axs[2].annotate('requested\n$n_e$ = 1.6', (1.6, 0.9505), color='#B00', fontsize=9, fontweight='bold',
                ha='left', va='top', xytext=(6, 0), textcoords='offset points')
axs[2].set_ylabel('$\\eta_{sub}^{(0)}$'); axs[2].set_ylim(0.885, 0.952)
axs[2].set_title('(c)  $n_e$ = 1.6 lands in the dip:  the plasmon is\n      bound again, costing 1.7 %p of substrate power', fontsize=10.3, loc='left')
axs[2].legend(loc='lower left', fontsize=8.8, frameon=False)
axs[3].set_ylabel('Power in the  $u>1$  bin'); axs[3].set_ylim(0, 0.055)
axs[3].set_title('(d)  the plasmon switches on exactly where\n      $n_{SPP}$ crosses $n_{sub}$', fontsize=10.3, loc='left')
axs[3].legend(loc='upper left', fontsize=8.8, frameon=False)
for ax in axs:
    ax.grid(axis='y', color='0.9', lw=0.8)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)

fig.suptitle('Two modelling choices, settled by calculation: the TCO real index and the ETL birefringence',
             x=0.05, ha='left', fontsize=12.5, y=0.955)
foot = ('Device A: Ag 100 nm (McPeak) / ETL 200 nm / EML 20 nm / HTL 200 nm / TCO 50 nm / substrate; 550 nm, isotropic dipole, PLQY = 1, u grid 3000 points, five-channel closure 1e-6 or better.  '
        '(a,b) ETL at n_o = 1.8, n_e = 1.6, k_TCO swept; curves start at k_TCO = 0.0025 because at exactly k_TCO = 0 the TCO-guided mode of the n = 2.0 film is a lossless pole that the u grid under-samples.  '
        '(c,d) TCO fixed at 1.9 + 0.02i, n_e swept.  The SPP index follows k_SPP² = k0² ε_m ε_e (ε_m − ε_o)/(ε_m² − ε_o ε_e) for a uniaxial dielectric on Ag.  '
        'Below the threshold the plasmon is leaky and its power reappears as substrate light; above it the plasmon is bound and the power is lost. The threshold tracks n_sub, moving from n_e = 1.60 to 1.68 as the substrate goes 1.8 → 1.9.')
fig.text(0.05, 0.015, '\n'.join(textwrap.wrap(foot, 215)), fontsize=8.2, va='bottom', ha='left', color='0.28')
fig.savefig('dr4f_mock.png', dpi=150)
print('saved fig2')
