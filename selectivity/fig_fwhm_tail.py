"""figSel6: emitter requirement map (FWHM x red-tail fraction) for on-axis cd/A, low-loss n 1.8 + scatterer + filter."""
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from openpyxl import Workbook
R = json.load(open('data/fwhm_tail.json'))
FW = sorted({r['fwhm'] for r in R}); TL = sorted({r['tail'] for r in R})
get = lambda k, w, t: next(r for r in R if r['fwhm'] == w and r['tail'] == t)[k]['cdA']
plt.rcParams.update({'font.size': 7})
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw=dict(wspace=0.35))
wb = Workbook(); wb.remove(wb.active)
for a, key, title in ((ax[0], 'dbr', '(a) quarter-wave edge DBR'), (ax[1], 'cavity', '(b) thin-Ag strong cavity (reference)')):
    Z = np.array([[get(key, w, t) for w in FW] for t in TL])
    im = a.imshow(Z, origin='lower', cmap='Blues', vmin=0, vmax=1700, aspect='auto')
    a.set_xticks(range(len(FW))); a.set_xticklabels(FW); a.set_yticks(range(len(TL))); a.set_yticklabels([f'{int(t * 100)}' for t in TL])
    for i in range(len(TL)):
        for j in range(len(FW)):
            a.text(j, i, f'{Z[i, j]:.0f}', ha='center', va='center', fontsize=5.5, color='white' if Z[i, j] > 1000 else '#0b0b0b')
    if key == 'dbr':
        a.contour(Z, levels=[1000], colors='#eb6834', linewidths=1.2)
    a.set_xlabel('main-band FWHM (nm, 530 nm)'); a.set_ylabel('red-tail photon fraction (%)'); a.set_title(title, loc='left')
    ws = wb.create_sheet(key); ws.append(['tail_fraction \\ FWHM_nm'] + FW)
    for t, row in zip(TL, Z): ws.append([t] + [float(v) for v in row])
fig.colorbar(im, ax=ax, label='on-axis current efficiency (cd/A), PLQY 1', shrink=0.9)
ws = wb.create_sheet('all'); ws.append(['FWHM', 'tail', 'config', 'cdA', 'EQE', 'duv30', 'CIE_x', 'CIE_y', 'dbr_design'])
for r in R:
    for k in ('nofilter', 'dbr', 'ideal20', 'ideal10'):
        ws.append([r['fwhm'], r['tail'], k, r[k]['cdA'], r[k]['EQE'], r[k]['duv30'], r[k]['xy'][0], r[k]['xy'][1], r['dbr_design']])
    for k in ('cavity', 'MLA'):
        ws.append([r['fwhm'], r['tail'], k, r[k]['cdA'], r[k]['EQE'], None, None, None, None])
fig.savefig('figures/figSel6_fwhm_tail.png', dpi=300, bbox_inches='tight'); fig.savefig('figures/figSel6_fwhm_tail.pdf', bbox_inches='tight')
wb.save('figures/figSel6_fwhm_tail.xlsx'); print('saved')
