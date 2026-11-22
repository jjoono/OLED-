"""Fig. S9: finite lateral size (data from finite_size.py, Al device, 550 nm)."""
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
src = open('oe_figs.py').read(); exec(src[:src.index('DV = np.load')])          # style, OUT, CLS, save_xlsx
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
exec(src[src.index('def save_xlsx'):src.index('\n\n', src.index('def save_xlsx'))])
D = json.load(open('data/finite_size.json'))
NM = {'planar': ('planar', 'Planar'), 'MLA hemi_total': ('MLA hemi', 'Hemispherical MLA'),
      'W_slope0.2': ('wavy', 'Wavy surface, slope 0.2'), 'scatter g0.8_t1.0': ('scatter', 'Scattering layer')}
fig = plt.figure(figsize=(FW, 4.7)); g_ = fig.add_gridspec(2, 2, wspace=0.36, hspace=0.42); ax = [fig.add_subplot(g_[0, 0]), fig.add_subplot(g_[0, 1]), fig.add_subplot(g_[1, 0])]; lg = fig.add_subplot(g_[1, 1]); lg.axis('off'); sheets = []
# (a) substrate size
rows = []
for nm, (c, lab) in NM.items():
    col = CLS[c][0]
    for var, ls in (('actual', '-'), ('exit_mix', '--')):
        z = sorted([r for r in D['size'] if r['structure'] == nm and r['variant'] == var], key=lambda r: r['W_over_t'])
        ax[0].plot([r['W_over_t'] for r in z], [r['eta'] for r in z], ls, color=col, lw=1.2)
        ax[0].axhline(z[0]['eta_inf'], color=col, lw=0.6, ls=':') if var == 'actual' else None
        rows += [[nm, var, r['az'], r['W_over_t'], r['eta'], r['edge'], r['eta_inf']] for r in z]
ax[0].axvline(12.5 / 1.295, color='k', lw=0.6); ax[0].text(12.5 / 1.295 * 1.1, 0.315, 'Benchmark', fontsize=6.6)
ax[0].set_xscale('log'); ax[0].set_xlabel('Substrate half-width / thickness'); ax[0].set_ylabel(ETA)
ax[0].set_ylim(0.3, 0.85); FS.label(ax[0], 'a')
sheets.append(('a_size', ['structure', 'variant', 'azimuth', 'W/t', 'eta', 'side loss', 'eta_infinite'], rows))
# (b) stack reflectance only under the emitting disc
rows = []
for nm in ('planar', 'MLA hemi_total', 'W_slope0.2'):
    col = CLS[NM[nm][0]][0]
    for outside, ls in (('glass/air', '-'), ('black', '--')):
        z = sorted([r for r in D['pixel'] if r['structure'] == nm and r['variant'] == 'actual' and r['outside'] == outside], key=lambda r: r['a_over_t'])
        ax[1].plot([r['a_over_t'] for r in z], [r['eta'] for r in z], ls if outside == 'glass/air' else '-.', color=col, lw=1.2)
    rows += [[r['structure'], r['variant'], r['outside'], r['a_over_t'], r['eta'], r['edge'], r['eta_inf']] for r in D['pixel'] if r['structure'] == nm]
ax[1].axvline(1 / 1.295, color='k', lw=0.6)
ax[1].set_xscale('log'); ax[1].set_xlabel('Pixel radius / thickness'); ax[1].set_ylabel(ETA); ax[1].set_ylim(0.3, 0.85); FS.label(ax[1], 'b')
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
ax[2].axhline(0, color='k', lw=0.6); ax[2].set_xticks(xs); ax[2].set_xticklabels(['Planar', 'Hemisph.\nMLA', 'Wavy\nsurface', 'Scatt.\nlayer'], fontsize=6.6); ax[2].tick_params(axis='x', which='both', top=False, length=0)
ax[2].set_ylabel('Gain from interface mixing (pp)'); FS.label(ax[2], 'c')
l1 = lg.legend(handles=[Line2D([], [], color=CLS[c][0], lw=1.6, label=lab) for c, lab in NM.values()], loc='upper left', bbox_to_anchor=(0.0, 1.04), title='Structure (color)', fontsize=6.6, alignment='left', labelspacing=0.3); lg.add_artist(l1)
lg.legend(handles=[Line2D([], [], color='k', ls='-', lw=1.2, label='(a), (b): no extra mixing; glass/air outside'), Line2D([], [], color='k', ls='--', lw=1.2, label='(a): complete mixing at the interface'), Line2D([], [], color='k', ls='-.', lw=1.2, label='(b): absorbing outside the pixel'),
                   Line2D([], [], color='k', ls=':', lw=0.6, label='(a): infinite lateral extent')] + [Patch(fc='#7d7c76', alpha=1 - 0.3 * j, hatch=h, ec='white', label='(c): ' + lab) for j, (lab, h) in enumerate((('infinite', ''), ('benchmark size', '////'), ('benchmark size, absorbing outside', '....')))],
          loc='lower left', bbox_to_anchor=(0.0, -0.16), fontsize=6.6, labelspacing=0.3)
sheets.append(('c_mixing_gain', ['structure', 'infinite (pp)', 'benchmark size (pp)', 'benchmark size, absorbing outside pixel (pp)'], rows))
FS.save(fig, f'{OUT}/figS_finite_size')
save_xlsx('figS_finite_size', sheets)
