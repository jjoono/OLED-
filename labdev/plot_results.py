import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R = json.load(open('results.json')); A = np.load('results.npz'); th = np.arange(90) + 0.5
C = {'bare': '#7d7c76', 'hemisphere': '#2a78d6', 'freeform': '#eb6834'}
fig, ax = plt.subplots(1, 3, figsize=(10, 3.1), gridspec_kw=dict(wspace=0.35))
for nm, c in C.items():
    for o, ls in (('Al', '-'), ('noAl', '--')):
        I = A[f'{nm}_{o}']; ax[0].plot(th, I * 100, ls, color=c, lw=1.1, label=f'{nm}, ' + ('Al outside' if o == 'Al' else 'no Al outside'))
        ax[1].plot(th, I / (I[:5] * np.diff(1 - np.cos(np.radians(np.arange(6))))).sum() * np.diff(1 - np.cos(np.radians(np.arange(6)))).sum(), ls, color=c, lw=1.1)
ax[1].plot(th, np.cos(np.radians(th)), ':', color='k', lw=0.8, label='Lambertian')
ax[0].set_xlabel('air angle (°)'); ax[0].set_ylabel('photons / sr / electron (%)'); ax[0].set_title('(a) angular intensity', loc='left'); ax[0].legend(fontsize=6, frameon=False)
ax[1].set_xlabel('air angle (°)'); ax[1].set_ylabel('normalized to 0°'); ax[1].set_title('(b) normalized', loc='left'); ax[1].legend(fontsize=6, frameon=False)
bands = ['0-20', '20-40', '40-60', '60-80', '80-90']; x = np.arange(len(bands)); k = 0
for nm, c in C.items():
    for o, hatch in (('Al', ''), ('noAl', '////')):
        ax[2].bar(x + (k - 2.5) * 0.13, [R[f'{nm}|{o}']['bands'][b] * 100 for b in bands], 0.13, color=c, hatch=hatch, edgecolor='white', lw=0.3); k += 1
ax[2].set_xticks(x); ax[2].set_xticklabels([b + '°' for b in bands]); ax[2].set_ylabel('EQE in band (%)'); ax[2].set_title('(c) band EQE (hatched: no Al outside)', loc='left')
for a in ax: a.spines[['top', 'right']].set_visible(False)
fig.savefig('lab_comparison.png', dpi=200, bbox_inches='tight')
