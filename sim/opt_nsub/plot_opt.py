import csv, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
r = list(csv.DictReader(open('/home/user/OLED-/sim/opt_nsub/opt_nsub.csv')))
fig, ax = plt.subplots(1, 3, figsize=(12, 3.6), dpi=200)
for m, c in (('Ag', '#0571b0'), ('Al', '#d95f02')):
    for meth, ls in (('eq3', '--'), ('series', '-')):
        x = [q for q in r if q['reflector'] == m and q['method'] == meth]; n = [float(q['n_sub']) for q in x]; g = lambda k: [float(q[k]) for q in x]
        lab = f"{m}, {'eq. (3)' if meth == 'eq3' else 'series'}"
        ax[0].plot(n, g('EQE'), color=c, ls=ls, marker='o', ms=2.5, label=lab)
        ax[1].plot(n, g('eta_ext'), color=c, ls=ls, marker='o', ms=2.5)
        ax[2].plot(n, g('max_eta_ext_any_d'), color=c, ls=ls, marker='o', ms=2.5)
ax[0].set(xlabel='$n_{sub}$', ylabel='max EQE (thickness-optimised)', title='(a) EQE at each optimum'); ax[0].legend(fontsize=7, frameon=False)
ax[1].set(xlabel='$n_{sub}$', ylabel=r'$\eta_{ext}$ at the EQE optimum', title=r'(b) $\eta_{ext}$ of that optimum')
ax[2].set(xlabel='$n_{sub}$', ylabel=r'max $\eta_{ext}$ over all thicknesses', title=r'(c) max $\eta_{ext}$ alone')
fig.tight_layout(); fig.savefig('/home/user/OLED-/sim/opt_nsub/opt_nsub.png')
