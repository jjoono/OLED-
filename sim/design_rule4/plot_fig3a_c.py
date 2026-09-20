# -*- coding: utf-8 -*-
"""Replacement for the p-sensitivity panel of Fig. 3(a): every transparent electrode in the
project's own ellipsometry library, simulated on the same stack.  Every point is a measured
film, so nothing here rests on an assumed k."""
import numpy as np, json, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import textwrap

F = json.load(open('/tmp/films.json'))
C_TCO, C_AG = '#0072B2', '#009E73'
fig, ax = plt.subplots(figsize=(5.6, 5.0))
fig.subplots_adjust(left=0.135, right=0.975, top=0.90, bottom=0.315)

A = np.loadtxt('f3a_ito_konig.csv', delimiter=',')
B = np.loadtxt('f3a_ag_konig.csv', delimiter=',')
ax.plot(A[1:, 7], A[1:, 10], color=C_TCO, lw=1.3, alpha=0.4, zorder=1)
ax.plot(B[1:, 7], B[1:, 10], color=C_AG, lw=1.3, alpha=0.4, zorder=1)

SHOW = {'l_ITO': ('ITO used here\n(n 1.86, k 0.0032)', (10, 7), 'left'),
        
        'etri_ITO': ('ITO, as-deposited\n(k 0.048)', (-9, 3), 'right'),
        'Ag_bulk': ('Ag 10 nm, bulk n', (10, 4), 'left'),
        'Ag_SNU': ('Ag 10 nm, poor seed\n(n 0.53)', (-9, 3), 'right')}
for f in F:
    c = C_TCO if f['family'] == 'TCO' else C_AG
    ax.plot(f['abs_electrode'], f['EQE'], 'o', ms=6.5, mfc=c, mec='w', mew=1.0, zorder=4)
    if f['film'] in SHOW:
        lab, (dx, dy), ha = SHOW[f['film']]
        ax.annotate(lab, (f['abs_electrode'], f['EQE']), xytext=(dx, dy),
                    textcoords='offset points', fontsize=7.8, color=c, ha=ha, va='center')
iz = [f for f in F if f['film'] == 'IZO'][0]
ax.annotate('IZO,  n 2.06, k 0.0012\nthe cleanest film, and still\nbelow the ITO: its index\nexceeds the substrate',
            (iz['abs_electrode'], iz['EQE']), xytext=(0.0032, 0.61), textcoords='data',
            fontsize=7.8, color=C_TCO, ha='left', va='center',
            arrowprops=dict(arrowstyle='->', color=C_TCO, lw=0.8, alpha=0.8))
ax.set_xscale('log'); ax.set_xlim(1.5e-3, 0.4); ax.set_ylim(0.35, 0.95)
ax.set_xlabel('power absorbed in the transparent electrode')
ax.set_ylabel('EQE  $=\\eta_{sub}^{(0)}\\,\\eta_{ext}$')
ax.grid(color='0.93', lw=0.8, which='both')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.set_title('(c)  every electrode in the library, on the same stack', fontsize=11, loc='left')
ax.legend(handles=[Line2D([], [], marker='o', ls='', mfc=C_TCO, mec='w', ms=6.5, label='TCO, 50 nm'),
                   Line2D([], [], marker='o', ls='', mfc=C_AG, mec='w', ms=6.5, label='thin Ag, 10 nm')],
          loc='lower left', fontsize=8.6, frameon=False, handletextpad=0.4)
foot = ('Each point is a measured film from the project\'s own library at 550 nm, run on the stack of (a) and (b): Ag 100 nm / ETL 200 / EML 20 / HTL 200 / electrode / substrate 1.80, p = 0.30.  '
        'The faint lines are the sweeps of (a) and (b) on the same axis.  Two things follow.  The families collapse: whatever the electrode is made of, what costs EQE is how much of the light it absorbs, '
        'and a 10 nm silver film with a bulk-like index is worth the same as a clean 50 nm oxide.  And the spread within the oxides is a factor of thirty in absorption, from 0.2 % to 10 %, which is the whole '
        'range from a research-grade film to an as-deposited one — the electrode is a process variable, not a material constant.')
fig.text(0.02, 0.02, '\n'.join(textwrap.wrap(foot, 98)), fontsize=7.4, va='bottom', ha='left', color='0.3')
fig.savefig('fig3a_c_films.png', dpi=200)
print('saved')
