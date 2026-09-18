"""Design rule #4 (parasitic absorption), simplified stack, 550 nm.
Reads dr4b_*.csv (columns: param, eta_sub, A', eta_ext(p=0.30), eta_ext(p=0.40), EQE(0.30), EQE(0.40))."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import textwrap

PRIME = '\u2032'
C_EXT, C_SUB, C_EQE = '#0072B2', '#E69F00', '#9467BD'
P_USE, P_ALT = 0.30, 0.40          # p read off the blue curve of Fig. 1c at n_sub = 1.8 / value used on the slide

kito = np.loadtxt('dr4b_kito.csv', delimiter=',')
nagb = np.loadtxt('dr4b_nagb.csv', delimiter=',')
nag2 = np.loadtxt('dr4b_nagb_k02.csv', delimiter=',')

fig = plt.figure(figsize=(15.2, 6.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.85], wspace=0.28,
                      left=0.055, right=0.985, top=0.885, bottom=0.275)
axes = [fig.add_subplot(gs[0, i]) for i in range(3)]

bands = {
    0: [(0.000, 0.020, '#E8E4B8', 'ITO\nw/ RTP'), (0.020, 0.040, '#DDD0E8', 'IZO\nw/o RTP'),
        (0.040, 0.080, '#CBE3EE', 'ITO\nw/o RTP')],
    1: [(0.00, 0.10, '#E8E4B8', 'Thick\nAg'), (0.10, 0.35, '#DDD0E8', 'Thin Ag w/\nideal seed'),
        (0.35, 0.50, '#CBE3EE', 'Mg:Ag, Yb:Ag,\nAu, Al…')],
}
for k, ax in enumerate(axes[:2]):
    d = [kito, nagb][k]
    for x0, x1, col, lab in bands[k]:
        ax.axvspan(x0, x1, color=col, alpha=0.55, lw=0, zorder=0)
        ax.text((x0+x1)/2, 0.245, lab, ha='center', va='bottom', fontsize=8.6, color='#333', zorder=1)
    ax.plot(d[:, 0], d[:, 3], color=C_EXT, lw=3.0, zorder=4, solid_capstyle='round')
    ax.plot(d[:, 0], d[:, 1], color=C_SUB, lw=2.0, zorder=3, solid_capstyle='round')
    ax.plot(d[:, 0], d[:, 5], color=C_EQE, lw=2.0, ls=(0, (5, 2)), zorder=3)
    ax.set_ylim(0.22, 1.0); ax.set_yticks(np.arange(0.3, 1.01, 0.1))
    ax.grid(axis='y', color='0.9', lw=0.8); ax.set_axisbelow(False)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    ax.set_ylabel('Efficiency')

axes[0].set_xlim(0, 0.08); axes[0].set_xlabel('Extinction coefficient of the TCO,  $k_{TCO}$')
axes[0].set_title('(a)  transparent electrode:  $\\bar n_{TCO}=1.8+k_{TCO}i$', fontsize=11, loc='left')
axes[1].set_xlim(0, 0.50); axes[1].set_xlabel('Refractive index of the thin bottom Ag,  $n_{Ag}$')
axes[1].set_title('(b)  10-nm bottom Ag:  $\\bar n_{Ag}=n_{Ag}+3.5i$   ($k_{TCO}=0$)', fontsize=11, loc='left')

# direct labels instead of a legend box
for ax, d in ((axes[0], kito), (axes[1], nagb)):
  for col, xf, dy, lab, cc in ((3, 0.55, -18, '$\\eta_{ext}$  substrate\u2192air', C_EXT),
                               (1, 0.85, 10, '$\\eta_{sub}^{(0)}$', C_SUB),
                               (5, 0.55, -19, 'EQE $=\\eta_{sub}^{(0)}\\eta_{ext}$', C_EQE)):
    i = int(xf*(len(d)-1))
    ax.annotate(lab, (d[i, 0], d[i, col]), textcoords='offset points', xytext=(0, dy),
                color=cc, fontsize=10, fontweight='bold', ha='center',
                va='top' if dy < 0 else 'bottom')

# --- panel (c): sensitivity to p ---
ax = axes[2]
pp = np.linspace(0.15, 0.9, 200)
for A, lab, sty in ((kito[0, 2], f'$k_{{TCO}}=0$   (A' + PRIME + f' = {kito[0,2]:.3f})', '-'),
                    (kito[8, 2], f'$k_{{TCO}}=0.02$   (A' + PRIME + f' = {kito[8,2]:.3f})', '--'),
                    (kito[-1, 2], f'$k_{{TCO}}=0.08$   (A' + PRIME + f' = {kito[-1,2]:.3f})', ':')):
    ax.plot(pp, pp/(pp+(1-pp)*A), color=C_EXT, ls=sty, lw=2.2, label=lab)
ax.legend(loc='lower right', fontsize=8.8, frameon=True, framealpha=0.9, edgecolor='0.85',
          borderpad=0.6, handlelength=2.6)
ax.axvline(P_USE, color='#B00', lw=1.4)
ax.text(P_USE-0.012, 0.243, f'p = {P_USE:.2f}', color='#B00', fontsize=9.2, fontweight='bold',
        va='bottom', ha='right')
ax.axvline(P_ALT, color='0.45', lw=1.0, ls=':')
ax.text(P_ALT+0.014, 0.988, f'p = {P_ALT:.2f} (slide)', color='0.35', fontsize=8.6, va='top', ha='left')
ax.set_xlim(0.15, 0.9); ax.set_ylim(0.22, 1.0); ax.set_yticks(np.arange(0.3, 1.01, 0.1))
ax.grid(axis='y', color='0.9', lw=0.8)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.set_xlabel('MLA single-pass escape probability,  p'); ax.set_ylabel('$\\eta_{ext}$')
ax.set_title('(c)  how much p matters:  $\\eta_{ext}=p/[p+(1-p)A′]$', fontsize=11, loc='left')

fig.suptitle('Design rule: parasitic absorption sets the substrate-to-air extraction efficiency',
             x=0.055, ha='left', fontsize=13, y=0.965)
foot = ('Simplified stack, 550 nm, isotropic dipole, PLQY = 1:  Ag 100 nm (McPeak n,k = 0.044 + 3.819i, fixed) / ETL 200 nm / EML 20 nm (dipole at the centre) / '
        'HTL 200 nm / TCO 50 nm / Ag 10 nm / substrate.  All organics n = 1.8 with k = 0; substrate n = 1.8, index-matched, so the waveguide channel is identically zero.  '
        'A′ = 1 − ⟨R_LED⟩, the reflectance of the OLED stack seen from the substrate, flux-weighted (cosθ sinθ) over substrate angles and averaged over p and s.  '
        'η_ext = p/[p + (1−p)A′]; η_sub^(0) = air + substrate-confined from the five-channel budget (closure to 3e-15).  '
        'In (b) the 10-nm Ag carries the swept index while the 100-nm reflector stays at McPeak; in (a) both Ag layers are McPeak.  '
        'Shaded bands reproduce the material ranges of the original slide.  '
        'p = 0.30 is read off the blue single-pass-escape-probability curve of Fig. 1c at n_sub = 1.8 (±0.01 from the digitisation); '
        'the grey line in (c) marks p = 0.40, the value that reproduces the original slide. Because η_ext depends on p only through the formula above, '
        'a different p needs no re-run — only a re-plot from the stored A′.')
fig.text(0.055, 0.015, '\n'.join(textwrap.wrap(foot, 205)), fontsize=8.2, va='bottom', ha='left', color='0.28')
fig.savefig('dr4b_mock.png', dpi=150)
print('saved')
