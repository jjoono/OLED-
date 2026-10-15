"""Figures D1-D2: scatterer thickness rule and pixel-level confinement (Excel alongside)."""
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from openpyxl import Workbook
plt.rcParams.update({'font.size': 7, 'axes.titlesize': 7.5})
D = json.load(open('pixel_results.json')); ref = D['ref']; rows = [r for r in D['rows'] if r['kind'] == 'scatter' and r['d'] > 0]
B = json.load(open('pitch_best.json'))


def save(fig, name, sheets):
    fig.savefig(name + '.png', dpi=300, bbox_inches='tight'); fig.savefig(name + '.pdf', bbox_inches='tight')
    wb = Workbook(); wb.remove(wb.active)
    for t, h, data in sheets:
        ws = wb.create_sheet(t[:31]); ws.append(h)
        for d in data: ws.append(d)
    wb.save(name + '.xlsx'); print('saved', name)


fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw=dict(wspace=0.3))
mk = {0.0: 'o', 0.6: 's', 0.8: '^'}; col = {0.5: '#2E7D5B', 1.0: '#3B6EA5', 2.0: '#C6512F'}
sh = []
for g in (0.0, 0.6, 0.8):
    for ls in (0.5, 1.0, 2.0):
        rr = sorted([r for r in rows if r['g'] == g and r['ls'] == ls and r['alb'] == 1.0], key=lambda r: r['tau_tr'])
        x = [r['tau_tr'] for r in rr]
        ax[0].plot(x, [r['open_beta'] / ref for r in rr], marker=mk[g], color=col[ls], ms=3, lw=0.8, label=f'g={g}, $\\ell_s$={ls} µm')
        sh.append((f'g{g}_ls{ls}', ['d_um', 'tau_tr', 'gain_open', 'gain_TIRwall_P6', 'crosstalk_open'],
                   [[r['d'], r['tau_tr'], r['open_beta'] / ref, r['tir_air_beta'] / ref, r['crosstalk']] for r in rr]))
ax[0].axhline(1, color='k', lw=0.6); ax[0].set_xscale('log')
ax[0].set_xlabel(r'transport optical thickness $\tau^* = d(1-g)/\ell_s$'); ax[0].set_ylabel('gain vs strong-cavity micro-OLED')
ax[0].set_title('(a) scatterer thickness rule (no walls, absorption-free)', loc='left'); ax[0].legend(fontsize=4.8, ncol=3); ax[0].grid(alpha=.25)
# (b) pixel pitch with optimized scatterer, spacer 0.1 um
pr = sorted([r for r in B['rows'] if r['s'] == 0.1], key=lambda r: r['P'])
P = [r['P'] for r in pr]
for key, lab, c, ls_ in (('open_beta', 'no walls (blurred, crosstalk)', '0.6', ':'), ('bank0.995_beta', 'mirror wall R=0.995', '#2E7D5B', '-'),
                         ('bank0.98_beta', 'mirror wall R=0.98', '#3B6EA5', '-'), ('tir_air_beta', 'air-gap TIR wall + Ag', '#8E5BA6', '--'),
                         ('bank0.95_beta', 'Ag wall R=0.95', '#C6512F', '-'), ('bank0.9_beta', 'Al wall R=0.90', '#7a4a2a', '-')):
    ax[1].plot(P, [r[key] / ref for r in pr], marker='o', ms=3, color=c, ls=ls_, label=lab)
ax[1].axhline(1, color='k', lw=0.6)
ax2 = ax[1].twinx(); ax2.plot(P, [r['crosstalk'] for r in pr], 'k:', lw=0.8); ax2.set_ylabel('crosstalk without walls', fontsize=6); ax2.set_ylim(0, 1)
ax[1].set_xlabel('pixel pitch (µm)'); ax[1].set_ylabel('gain vs strong-cavity micro-OLED')
ax[1].set_title('(b) pixel confinement (scatterer 0.5 µm, spacer 0.1 µm)', loc='left'); ax[1].legend(fontsize=5, loc='lower right'); ax[1].grid(alpha=.25)
sh.append(('pitch_scan', ['P_um', 'spacer_um', 'gain_open', 'crosstalk_open', 'gain_Rw0.90', 'gain_Rw0.95', 'gain_Rw0.98', 'gain_Rw0.995', 'gain_TIRair', 'gain_TIRn1.25', 'mean_wall_hits'],
           [[r['P'], r['s'], r['open_beta'] / ref, r['crosstalk'], r['bank0.9_beta'] / ref, r['bank0.95_beta'] / ref, r['bank0.98_beta'] / ref,
             r['bank0.995_beta'] / ref, r['tir_air_beta'] / ref, r['tir_n125_beta'] / ref, r['mean_hits']] for r in B['rows']]))
save(fig, 'figD_pixel_scatterer', sh)
