"""Fig. S9: finite lateral size (data from finite_size.py, Al device, 550 nm)."""
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
src = open('oe_figs.py').read(); exec(src[:src.index('DV = np.load')])          # rcParams, OUT, CLS, save_xlsx
exec(src[src.index('def save_xlsx'):src.index('\n\n', src.index('def save_xlsx'))])
D = json.load(open('data/finite_size.json'))
NM = {'planar': ('planar', 'planar interface'), 'MLA hemi_total': ('MLA hemi', 'hemispherical MLA'),
      'W_slope0.2': ('wavy', 'wavy, rms slope 0.2'), 'scatter g0.8_t1.0': ('scatter', 'scattering layer')}
fig, ax = plt.subplots(1, 3, figsize=(7.1, 2.5), gridspec_kw=dict(wspace=0.38)); sheets = []
# (a) substrate size
rows = []
for nm, (c, lab) in NM.items():
    col = CLS[c][0]
    for var, ls in (('actual', '-'), ('exit_mix', '--')):
        z = sorted([r for r in D['size'] if r['structure'] == nm and r['variant'] == var], key=lambda r: r['W_over_t'])
        ax[0].plot([r['W_over_t'] for r in z], [r['eta'] for r in z], ls, color=col, lw=1, label=lab if var == 'actual' else None)
        ax[0].axhline(z[0]['eta_inf'], color=col, lw=0.4, ls=':') if var == 'actual' else None
        rows += [[nm, var, r['az'], r['W_over_t'], r['eta'], r['edge'], r['eta_inf']] for r in z]
ax[0].axvline(12.5 / 1.295, color=INK2, lw=0.6); ax[0].text(12.5 / 1.295 * 1.08, 0.33, 'benchmark', fontsize=5.5, color=INK2)
ax[0].set_xscale('log'); ax[0].set_xlabel('substrate half-width / thickness'); ax[0].set_ylabel('η_sta')
ax[0].set_ylim(0.3, 0.85); ax[0].set_title('(a) substrate size (glass/air sides)', loc='left')
ax[0].legend(fontsize=5, frameon=False, loc='lower right')
sheets.append(('a_size', ['structure', 'variant', 'azimuth', 'W/t', 'eta', 'side loss', 'eta_infinite'], rows))
# (b) stack reflectance only under the emitting disc
rows = []
for nm in ('planar', 'MLA hemi_total', 'W_slope0.2'):
    col = CLS[NM[nm][0]][0]
    for outside, ls in (('glass/air', '-'), ('black', '--')):
        z = sorted([r for r in D['pixel'] if r['structure'] == nm and r['variant'] == 'actual' and r['outside'] == outside], key=lambda r: r['a_over_t'])
        ax[1].plot([r['a_over_t'] for r in z], [r['eta'] for r in z], ls, color=col, lw=1)
    rows += [[r['structure'], r['variant'], r['outside'], r['a_over_t'], r['eta'], r['edge'], r['eta_inf']] for r in D['pixel'] if r['structure'] == nm]
ax[1].plot([], [], 'k-', lw=1, label='glass/air outside the pixel'); ax[1].plot([], [], 'k--', lw=1, label='absorbing outside the pixel')
ax[1].axvline(1 / 1.295, color=INK2, lw=0.6)
ax[1].set_xscale('log'); ax[1].set_xlabel('pixel radius / thickness'); ax[1].set_ylabel('η_sta'); ax[1].set_ylim(0.3, 0.85)
ax[1].set_title('(b) stack reflectance only under the pixel', loc='left'); ax[1].legend(fontsize=5, frameon=False, loc='upper left')
sheets.append(('b_pixel', ['structure', 'variant', 'outside', 'a/t', 'eta', 'side loss', 'eta_infinite'], rows))
# (c) gain from complete reciprocal mixing at the escape interface: infinite vs finite
rows = []; labs = []; xs = np.arange(4)
for k, nm in enumerate(NM):
    g_inf = [r for r in D['size'] if r['structure'] == nm and r['W_over_t'] == 256]
    ginf = next(r['eta_inf'] for r in g_inf if r['variant'] == 'exit_mix') - next(r['eta_inf'] for r in g_inf if r['variant'] == 'actual')
    b = {v: next(r['eta'] for r in D['size'] if r['structure'] == nm and r['variant'] == v and abs(r['W_over_t'] - 9.65) < 1e-6) for v in ('actual', 'exit_mix')}
    p = {v: next(r['eta'] for r in D['pixel'] if r['structure'] == nm and r['variant'] == v and r['outside'] == 'black' and abs(r['a_over_t'] - 1 / 1.295) < 1e-6) for v in ('actual', 'exit_mix')}
    vals = [ginf * 100, (b['exit_mix'] - b['actual']) * 100, (p['exit_mix'] - p['actual']) * 100]
    for j, (v, hatch) in enumerate(zip(vals, ('', '////', '....'))):
        ax[2].bar(k + (j - 1) * 0.26, v, 0.25, color=CLS[NM[nm][0]][0], alpha=1 - 0.3 * j, hatch=hatch, edgecolor='white', lw=0.3)
    rows.append([nm] + vals)
for j, (lab, hatch) in enumerate((('infinite', ''), ('benchmark size', '////'), ('benchmark size, absorbing outside pixel', '....'))):
    ax[2].bar([], [], color='#9c9b96', alpha=1 - 0.3 * j, hatch=hatch, edgecolor='white', label=lab)
ax[2].axhline(0, color=INK2, lw=0.5); ax[2].set_xticks(xs); ax[2].set_xticklabels(['planar', 'hemi. MLA', 'wavy 0.2', 'scattering'], fontsize=6)
ax[2].set_ylabel('gain from complete mixing at the escape interface (pp)'); ax[2].set_title('(c) mixing gain (Al, 550 nm)', loc='left')
ax[2].legend(fontsize=5, frameon=False, loc='upper right')
sheets.append(('c_mixing_gain', ['structure', 'infinite (pp)', 'benchmark size (pp)', 'benchmark size, absorbing outside pixel (pp)'], rows))
fig.savefig(f'{OUT}/figS_finite_size.png', dpi=300, bbox_inches='tight'); fig.savefig(f'{OUT}/figS_finite_size.pdf', bbox_inches='tight')
save_xlsx('figS_finite_size', sheets)
