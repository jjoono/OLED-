"""Design rule #4 (parasitic absorption): two ALTERNATIVE bottom electrodes, 550 nm.
(a) reads f3a_ito_konig.csv, (b) reads f3a_ag_konig.csv.
CSV columns: param, eta_sub, A', wg, spp, abs, abs_top, abs_bottom, eta_ext(0.30), eta_ext(0.40), EQE(0.30), EQE(0.40)."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import textwrap

C_EXT, C_SUB, C_EQE = '#0072B2', '#E69F00', '#9467BD'
P_USE, P_ALT = 0.30, 0.40
ito = np.loadtxt('f3a_ito_konig.csv', delimiter=',')
ag  = np.loadtxt('f3a_ag_konig.csv',  delimiter=',')

fig = plt.figure(figsize=(15.2, 6.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.9], wspace=0.28,
                      left=0.055, right=0.985, top=0.875, bottom=0.275)
axes = [fig.add_subplot(gs[0, i]) for i in range(3)]

bands = {0: [(0.000, 0.020, '#E8E4B8', 'ITO\nw/ RTP'), (0.020, 0.040, '#DDD0E8', 'IZO\nw/o RTP'),
             (0.040, 0.080, '#CBE3EE', 'ITO\nw/o RTP')],
         1: [(0.00, 0.10, '#E8E4B8', 'Thick\nAg'), (0.10, 0.35, '#DDD0E8', 'Thin Ag w/\nideal seed'),
             (0.35, 0.50, '#CBE3EE', 'Mg:Ag, Yb:Ag,\nAu, Al…')]}
for k, (ax, d) in enumerate(zip(axes[:2], (ito, ag))):
    for x0, x1, col, lab in bands[k]:
        ax.axvspan(x0, x1, color=col, alpha=0.55, lw=0, zorder=0)
        ax.text((x0+x1)/2, 0.245, lab, ha='center', va='bottom', fontsize=8.6, color='#333')
    ax.plot(d[:, 0], d[:, 8], color=C_EXT, lw=3.0, zorder=4, solid_capstyle='round')
    ax.plot(d[:, 0], d[:, 1], color=C_SUB, lw=2.0, zorder=3, solid_capstyle='round')
    ax.plot(d[:, 0], d[:, 10], color=C_EQE, lw=2.0, ls=(0, (5, 2)), zorder=3)
    for col, xf, dy, lab, cc in ((8, 0.55, -18, '$\\eta_{ext}$  substrate→air', C_EXT),
                                 (1, 0.85, 10, '$\\eta_{sub}^{(0)}$', C_SUB),
                                 (10, 0.55, -19, 'EQE $=\\eta_{sub}^{(0)}\\eta_{ext}$', C_EQE)):
        i = int(xf*(len(d)-1))
        ax.annotate(lab, (d[i, 0], d[i, col]), textcoords='offset points', xytext=(0, dy),
                    color=cc, fontsize=10, fontweight='bold', ha='center',
                    va='top' if dy < 0 else 'bottom')
    ax.set_ylim(0.22, 1.0); ax.set_yticks(np.arange(0.3, 1.01, 0.1))
    ax.grid(axis='y', color='0.9', lw=0.8)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    ax.set_ylabel('Efficiency')

axes[0].set_xlim(0, 0.08); axes[0].set_xlabel('Extinction coefficient of the TCO,  $k_{TCO}$')
axes[0].set_title('(a)  device A — TCO electrode, 50 nm:  $\\bar n=1.864+k_{TCO}i$', fontsize=11, loc='left')
axes[1].set_xlim(0, 0.50); axes[1].set_xlabel('Refractive index of the thin Ag,  $n_{Ag}$')
axes[1].set_title('(b)  device B — thin-Ag electrode, 10 nm:  $\\bar n=n_{Ag}+3.819i$', fontsize=11, loc='left')
axes[1].axvline(0.044, color='#B00', lw=1.2, ls='-')
axes[1].text(0.052, 0.965, 'bulk Ag\n(McPeak, 0.044)', color='#B00', fontsize=8.5, va='top')

# --- (c) every measured electrode on the same stack ---
import json as _json
FILMS = _json.load(open('/tmp/films.json')) if __import__('os').path.exists('/tmp/films.json') \
        else _json.load(open('films.json'))
C_T, C_A = '#0072B2', '#009E73'
ax = axes[2]
ax.plot(ito[1:, 7], ito[1:, 10], color=C_T, lw=1.3, alpha=0.4, zorder=1)
ax.plot(ag[1:, 7], ag[1:, 10], color=C_A, lw=1.3, alpha=0.4, zorder=1)
SHOW = {'l_ITO': ('ITO used here\n(n 1.86, k 0.0032)', (10, 7), 'left'),
        'etri_ITO': ('ITO, as-deposited\n(k 0.048)', (-9, 3), 'right'),
        'Ag_bulk': ('Ag 10 nm, bulk n', (10, 4), 'left'),
        'Ag_SNU': ('Ag 10 nm, poor seed', (-9, 3), 'right')}
for f in FILMS:
    c = C_T if f['family'] == 'TCO' else C_A
    ax.plot(f['abs_electrode'], f['EQE'], 'o', ms=6.5, mfc=c, mec='w', mew=1.0, zorder=4)
    if f['film'] in SHOW:
        lab, (dx, dy), ha = SHOW[f['film']]
        ax.annotate(lab, (f['abs_electrode'], f['EQE']), xytext=(dx, dy), textcoords='offset points',
                    fontsize=8.0, color=c, ha=ha, va='center')
iz = [f for f in FILMS if f['film'] == 'IZO'][0]
ax.annotate('IZO,  n 2.06, k 0.0012\nthe cleanest film, and still below\nthe ITO: its index exceeds the substrate',
            (iz['abs_electrode'], iz['EQE']), xytext=(0.0032, 0.575), textcoords='data',
            fontsize=8.0, color=C_T, ha='left', va='center',
            arrowprops=dict(arrowstyle='->', color=C_T, lw=0.8, alpha=0.8))
ax.set_xscale('log'); ax.set_xlim(1.5e-3, 0.4); ax.set_ylim(0.35, 0.95)
ax.set_xlabel('Power absorbed in the transparent electrode')
ax.set_ylabel('EQE $=\\eta_{sub}^{(0)}\\eta_{ext}$')
ax.grid(color='0.9', lw=0.8, which='both')
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
ax.legend(handles=[Line2D([], [], marker='o', ls='', mfc=C_T, mec='w', ms=6.5, label='TCO, 50 nm'),
                   Line2D([], [], marker='o', ls='', mfc=C_A, mec='w', ms=6.5, label='thin Ag, 10 nm')],
          loc='lower left', fontsize=8.8, frameon=False, handletextpad=0.4)
ax.set_title('(c)  measured electrodes, same stack', fontsize=11, loc='left')

fig.suptitle('Design rule: parasitic absorption sets the substrate-to-air extraction efficiency',
             x=0.055, ha='left', fontsize=13, y=0.96)
foot = ('Two ALTERNATIVE devices, not one stack. Common part, 550 nm, isotropic dipole, PLQY = 1:  '
        'Ag 100 nm (McPeak n,k = 0.044 + 3.819i, fixed) / ETL 200 nm / EML 20 nm (dipole at the centre) / HTL 200 nm / [bottom electrode] / substrate.  '
        'Device A takes a 50-nm TCO whose real index is 1.8636, the Koenig et al. (2014) ITO at 550 nm, so the sweep passes through the operating point of Fig. 2; device B a 10-nm thin Ag.  The 100-nm reflector is never swept.  '
        'EML and HTL isotropic at n = 1.8 with k = 0; the ETL is isotropic at n = 1.8, as in Fig. 2; substrate n = 1.8; u grid 3000 points; five-channel closure 1e-15.  '
        'A′ = 1 − ⟨R_LED⟩, the reflectance of the OLED stack seen from the substrate, flux-weighted (cosθ sinθ) over substrate angles and averaged over p and s.  '
        'η_ext = p/[p + (1−p)A′]; η_sub^(0) = air + substrate-confined.  Shaded bands reproduce the material ranges of the original slide.  '
        'p = 0.30 throughout, read off the single-pass escape probability at n_sub = 1.8; because η_ext depends on p only through the formula above, a different p needs no re-run.  '
        'Panel (c) puts every transparent electrode in the project library on this same stack — each point is a measured film at 550 nm, the faint lines are the sweeps of (a) and (b) on that axis.  '
        'The two families fall on one trend: what costs EQE is how much of the light the electrode absorbs, whatever it is made of.  The oxides span a factor of thirty in absorption, from a research-grade '
        'film to an as-deposited one, so the electrode is a process variable rather than a material constant.')
fig.text(0.055, 0.012, '\n'.join(textwrap.wrap(foot, 205)), fontsize=8.2, va='bottom', ha='left', color='0.28')
fig.savefig('dr4e_mock.png', dpi=150)
print('saved')
