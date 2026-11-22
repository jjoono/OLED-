"""Fig. S5: best-so-far objective versus evaluation count for the hemispherical campaign, the original freeform
per-band campaign, and the restart control (search fidelity: 10,000 rays, 31 wavelengths), with the final
re-evaluated values.  Output: figures_oe/figS5_optimizer_convergence.{png,pdf,xlsx}"""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.io import loadmat
from openpyxl import Workbook
from matplotlib.lines import Line2D
import figstyle as FS
FS.apply()
H = loadmat('opt_hemisphere_result.mat'); O = loadmat('opt_4band_result_25by25.mat'); W = loadmat('warmstart_hemisphere_result.mat')
F = loadmat('freeform_EQEtotal_result.mat', squeeze_me=True)
arms = [('0–20°', 0, 14), ('20–40°', 20, 15), ('40–60°', 40, 16), ('60–80°', 60, 17), ('Total EQE', -1, 13)]
fig, axs = plt.subplots(2, 3, figsize=(FS.W, 3.5), gridspec_kw=dict(wspace=0.42, hspace=0.5)); axs = axs.ravel(); axs[5].axis('off'); wb = Workbook(); wb.remove(wb.active)
for a, (lab, wtag, col) in enumerate(arms):
    ax = axs[a]; ws = wb.create_sheet(lab.replace('°', 'deg')); ws.append(['campaign', 'evaluation', 'best_so_far'])
    for name, D, c in (('Hemisphere (3 variables)', H, '#2a78d6'), ('Freeform, original', O, '#eb6834'), ('Freeform, restart', W, '#1baf7a')):
        L = D['EVAL_LOG']; L = L[L[:, 19] == wtag]
        if len(L) == 0:
            continue
        v = np.where(np.isfinite(L[:, col]), L[:, col], -np.inf); b = np.maximum.accumulate(v); ok = np.isfinite(b)
        ax.plot(np.arange(1, len(b) + 1)[ok], b[ok], color=c, lw=1.1, label=name)
        for i, x in enumerate(b): ws.append([name, i + 1, float(x) if np.isfinite(x) else None])
    k = [0, 20, 40, 60, -1].index(wtag)
    ax.axhline(H['hemi_val'].ravel()[k], color='#2a78d6', ls=':', lw=0.8)
    if wtag == -1:
        for s in F['start_log']:
            ax.plot(140, float(s['bestEQE']), 'o', color='#eb6834', ms=3)
        ax.axhline(float(W['ws_val'].ravel()[k]), color='#1baf7a', ls=':', lw=0.8)
    else:
        fin = O['band_eqe_hi'].ravel()[k]; ax.axhline(fin, color='#eb6834', ls=':', lw=0.8); ax.axhline(W['ws_val'].ravel()[k], color='#1baf7a', ls=':', lw=0.8)
    ax.text(0.95, 0.07, f"({'abcde'[a]}) {lab}", transform=ax.transAxes, fontsize=7, ha='right', va='bottom'); ax.tick_params(labelsize=7)
    if a >= 3 or a == 2: ax.set_xlabel('Evaluation')
    if a in (0, 3): ax.set_ylabel('Best EQE so far')
h, l = axs[0].get_legend_handles_labels(); axs[5].legend(h + [Line2D([], [], color='k', ls=':', lw=0.8), Line2D([], [], ls='none', marker='o', color='#eb6834', ms=3)], l + ['Final re-evaluated values', 'Independent starts (total)'], fontsize=6.6, loc='center left', bbox_to_anchor=(-0.1, 0.5))
FS.save(fig, 'figures_oe/figS5_optimizer_convergence')
wb.save('figures_oe/figS5_optimizer_convergence.xlsx')
for k, (lab, wtag, col) in enumerate(arms):
    print(lab, 'n evals hemi %d, original %d, restart %d' % ((H['EVAL_LOG'][:, 19] == wtag).sum(), (O['EVAL_LOG'][:, 19] == wtag).sum(), (W['EVAL_LOG'][:, 19] == wtag).sum()))
