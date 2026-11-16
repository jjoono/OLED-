import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
P = (2, 3, 4, 5, 7, 10, 15, 25)
ld = lambda p: json.load(open('results.json' if p == 10 else f'results_patch{p:g}.json'))
D = {p: ld(p) for p in P}; C = {'bare': '#7d7c76', 'hemisphere': '#2a78d6', 'freeform': '#eb6834'}
M = {'bare': (32.2, 1.01, [30.6, 34.1], [0.95, 1.08]), 'hemisphere': (25.3, 0.66, [23.8, 27.7], [0.63, 0.71]), 'freeform': (30.3, 1.065, [25.9, 33.7], [0.99, 1.13])}
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4), gridspec_kw=dict(wspace=0.3))
for nm, c in C.items():
    for o, ls in (('Al', '-'), ('noAl', '--')):
        k = f'{nm}|{o}'
        ax[0].plot(P, [D[p][k]['EQE'] * 100 for p in P], ls, color=c, marker='o', ms=3, lw=1.1, label=f'{nm}, ' + ('Al outside' if o == 'Al' else 'no Al outside'))
        ax[1].plot(P, [D[p][k]['EQE'] / D[p][k]['EQE_lambertian'] for p in P], ls, color=c, marker='o', ms=3, lw=1.1)
    m, pe, rr, pr = M[nm]
    ax[0].errorbar(3, m, yerr=[[m - rr[0]], [rr[1] - m]], fmt='s', color=c, mec='k', ms=6, capsize=3, zorder=5)
    ax[1].errorbar(3, pe, yerr=[[pe - pr[0]], [pr[1] - pe]], fmt='s', color=c, mec='k', ms=6, capsize=3, zorder=5)
for a in ax:
    a.set_xscale('log'); a.set_xticks(P); a.set_xticklabels([str(p) for p in P]); a.set_xlabel('lens patch width (mm, square, centred)'); a.spines[['top', 'right']].set_visible(False)
    a.axvline(3, color='k', lw=0.4, ls=':')
ax[0].set_ylabel('EQE (%)'); ax[0].set_title('(a) EQE vs lens patch size (squares: measured, 3 mm)', loc='left', fontsize=8); ax[0].legend(fontsize=6, frameon=False, loc='lower right')
ax[1].set_ylabel('PE factor = EQE / (π I(0))'); ax[1].set_ylim(0.5, 1.5); ax[1].set_title('(b) angular factor (squares: measured)', loc='left', fontsize=8)
fig.savefig('patch_sweep.png', dpi=200, bbox_inches='tight')
