import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d = json.load(open('spot.json')); r = np.array(d['r_mm'])
C = {'bare': '#7d7c76', 'hemisphere': '#2a78d6', 'freeform': '#eb6834'}; MEAS = {'bare': (30.1, [29.7, 27.7, 30.7, 31.4, 30.8]), 'hemisphere': (24.4, [22.9, 27.1, 22.5, 25.3, 24.1]), 'freeform': (25.5, [22.9, 26.7, 23.7, 28.7])}
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4), gridspec_kw=dict(wspace=0.3))
for nm, c in C.items():
    for o, ls in (('Al', '-'), ('noAl', '--')):
        ax[0].plot(r, np.array(d[f'{nm}|{o}']['EQE_within']) * 100, ls, color=c, lw=1.2, label=f'{nm}, ' + ('Al outside' if o == 'Al' else 'no Al outside'))
    m, v = MEAS[nm]; ax[0].axhspan(min(v), max(v), color=c, alpha=0.12, lw=0); ax[0].axhline(m, color=c, lw=0.6, ls=':')
ax[0].set_xlim(0, 6); ax[0].set_ylim(0, 50); ax[0].axvline(1.0, color='k', lw=0.5); ax[0].text(1.03, 2, 'pixel edge', fontsize=6)
ax[0].set_xlabel('radius of the collected area on the substrate (mm)'); ax[0].set_ylabel('EQE collected within radius (%)')
ax[0].set_title('(a) simulated EQE vs collection radius; bands = measured range', loc='left', fontsize=8); ax[0].legend(fontsize=6, frameon=False, loc='lower right')
b = np.array(d['bare|Al']['EQE_within'])
for nm in ('hemisphere', 'freeform'):
    for o, ls in (('Al', '-'), ('noAl', '--')):
        ax[1].plot(r, np.array(d[f'{nm}|{o}']['EQE_within']) / np.array(d[f'bare|{o}']['EQE_within']), ls, color=C[nm], lw=1.2)
    ax[1].axhline(MEAS[nm][0] / 30.1, color=C[nm], lw=0.8, ls=':')
ax[1].axhline(1, color='k', lw=0.4); ax[1].set_xlim(0.8, 6); ax[1].set_ylim(0.4, 1.5)
ax[1].set_xlabel('radius of the collected area (mm)'); ax[1].set_ylabel('MLA / bare'); ax[1].set_title('(b) enhancement vs collection radius (dotted: measured)', loc='left', fontsize=8)
for a in ax: a.spines[['top', 'right']].set_visible(False)
fig.savefig('spot_comparison.png', dpi=200, bbox_inches='tight')
