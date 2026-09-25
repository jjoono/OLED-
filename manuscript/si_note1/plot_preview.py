import csv, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
rows = list(csv.DictReader(open('/home/user/OLED-/sim/pass_resolved/pass_resolved.csv')))
C = {'Ag n_sub 1.5': ('#0571b0', '-'), 'Ag n_sub 1.8': ('#0571b0', '--'), 'Al n_sub 1.5': ('#d95f02', '-'), 'Al n_sub 1.8': ('#d95f02', '--')}
fig, ax = plt.subplots(1, 3, figsize=(11, 3.4), dpi=200); K = 15
for c, (col, ls) in C.items():
    r = [x for x in rows if x['case'] == c]; k = np.array([int(x['k']) for x in r])
    p = np.array([float(x['p_k']) for x in r]); A = np.array([float(x["A'_k"]) for x in r]); cum = np.array([float(x['cumulative eta_ext']) for x in r])
    pc, Ac, e3 = float(r[0]['p (cos.sin)']), float(r[0]["A' (cos.sin)"]), float(r[0]['eq3'])
    m = k <= K
    ax[0].plot(k[m], p[m], 'o', color=col, ls=ls, ms=3, lw=1.2, label=c); ax[0].axhline(pc, color=col, lw=0.6, ls=':')
    ax[1].plot(k[m], 100*A[m], 'o', color=col, ls=ls, ms=3, lw=1.2); ax[1].axhline(100*Ac, color=col, lw=0.6, ls=':')
    ax[2].plot(k[m], cum[m], 'o', color=col, ls=ls, ms=3, lw=1.2); ax[2].axhline(e3, color=col, lw=0.6, ls=':')
ax[0].set(xlabel='pass k', ylabel='escape probability $p_k$', title='(a)'); ax[0].legend(fontsize=7, frameon=False)
ax[1].set(xlabel='pass k', ylabel="round-trip loss $A'_k$ (%)", title='(b)')
ax[2].set(xlabel='pass k', ylabel=r'cumulative $\eta_{ext}$', title='(c)', ylim=(0.2, 1.0))
fig.tight_layout(); fig.savefig('pass_resolved_preview.png')
