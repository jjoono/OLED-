"""Preview of Supplementary Fig. 1 from figure_update_rawdata.xlsx (sheet SI_Fig1_nsub)."""
import openpyxl, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
s = openpyxl.load_workbook('/home/user/OLED-/figures/update/figure_update_rawdata.xlsx', data_only=True)['SI_Fig1_nsub']
h = [c.value for c in s[1]]; d = np.array([[c.value for c in r] for r in s.iter_rows(min_row=2) if r[0].value is not None], float)
col = lambda k: d[:, h.index(k)]; n = col('n_sub')
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6), dpi=200)
for m, c, cut in (('Ag', '#0571b0', (1.46, 1.59, 1.76, 1.78)), ('Al', '#d95f02', (1.42, 1.58, 1.75, 1.78))):
    ax[0].plot(n, col('eta_ext_series_' + m), color=c, lw=1.5, label=f'{m}: matrix series')
    ax[0].plot(n, col('eta_ext_eq3_' + m), color=c, lw=1.0, ls='--', label=f'{m}: eq. (3)')
    ax[1].plot(n, col('Psub_beyond70deg_' + m), color=c, lw=1.5, label=m)
    for x in cut:
        for a in ax: a.axvline(x, color=c, lw=0.5, ls=':')
ax[0].set(xlabel='$n_{sub}$ (= $n_{MLA}$)', ylabel=r'$\eta_{ext}$', title='(a)'); ax[0].legend(fontsize=7, frameon=False)
ax[1].set(xlabel='$n_{sub}$ (= $n_{MLA}$)', ylabel='fraction of $P_{sub}$ beyond 70°', title='(b)'); ax[1].legend(fontsize=7, frameon=False)
fig.text(0.5, 0.005, 'dotted: cut-off indices at which guided light starts to enter the substrate (Ag blue, Al orange); 550 nm, Fig. 2(d)-(f) stack', ha='center', fontsize=6.5)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig('/home/user/OLED-/figures/si_fig1/si_fig1_preview.png')
