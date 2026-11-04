import json, numpy as np, matplotlib, collim as C
matplotlib.use('Agg'); import matplotlib.pyplot as plt
ref = {'MLA, max-EQE design': C.evaluate(C.load_dev('dev/Pt_IZO_n1.8_dbr4_mla.npz'), 'mla', 740, 80, 0.8),
       'Planar ITO/glass': C.evaluate(C.load_dev('dev/Pt_ITO_n1.5_dbr0_flat.npz'), 'flat', 740, 80, 0.8)}
res = {k: np.array(json.load(open('results/%s.json' % t))['D']) for k, t in
       [('Planar strong cavity (4-pair DBR)', 'Pt_planar_IZO_n1.5_N4_Leq10'), ('Collimator, 8-pair DBR (FWHM 80 nm)', 'Pt_N8_Leq10'),
        ('Collimator, FWHM 40 nm', 'FW40_N8_Leq10'), ('Collimator, FWHM 20 nm', 'FW20_N8_Leq10')]}
allD = {**ref, **res}; cols = ['#2a9d8f', '#888888', '#bbbbbb', '#d1495b', '#e07a5f', '#7b2cbf']
sr = 2 * np.pi * (np.cos(np.radians(np.arange(90))) - np.cos(np.radians(np.arange(1, 91))))
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5})
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6), dpi=200)
for (k, D), c in zip(allD.items(), cols):
    ax[0].plot(C.TH, np.pi * D / sr, color=c, lw=1.6, label='%s (EQE %.0f%%)' % (k, 100 * D.sum()))
ax[0].set(xlabel='Angle in air (°)', ylabel='Radiance, Lambertian-equivalent EQE', xlim=(0, 90), ylim=(0, None)); ax[0].legend(fontsize=6.8, frameon=False)
ax[0].set_title('Angular emission (Pt(fprpz)$_2$-like, PLQY 0.8)', fontsize=9)
names = list(allD); x = np.arange(len(names)); w = 0.2
ms = [C.metrics(allD[k]) for k in names]
for i, (key, lab) in enumerate([('EQE', 'EQE (total)'), ('Leq10', 'On-axis radiance (±10°)'), ('E_z1_a1', 'Irradiance @1 mm ×1'), ('E_z3_a1', 'Irradiance @3 mm ×3')]):
    v = np.array([m[key] for m in ms]) * (3 if key == 'E_z3_a1' else 1)
    ax[1].bar(x + (i - 1.5) * w, v, w, label=lab)
ax[1].set_xticks(x); ax[1].set_xticklabels(['MLA\nmax-EQE', 'Planar\nITO', 'Planar\ncavity', 'Collim.\n80 nm', 'Collim.\n40 nm', 'Collim.\n20 nm'], fontsize=7.5)
ax[1].set_ylabel('photons per electron (or Lambertian-eq.)'); ax[1].legend(fontsize=7, frameon=False); ax[1].set_title('Metrics (source radius 1 mm)', fontsize=9)
fig.tight_layout(); fig.savefig('collimator_summary.png')
