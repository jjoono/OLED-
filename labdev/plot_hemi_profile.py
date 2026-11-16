"""Measured (260328, hemisphere devices, values relayed from the local analysis session: median of 4 points, negative-angle
side, normalized to 0 deg) versus simulated normalized angular intensity (actual substrate layout, layout_results.json)."""
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
A = np.arange(0, 90, 10)
MEAS = {'S3D2': [1.00, 0.97, 0.86, 0.71, 0.49, 0.32, 0.28, 0.13, 0.08], 'S2D4 (#2)': [1.00, 0.96, 0.84, 0.69, 0.52, 0.37, 0.27, 0.20, 0.13],
        'S1D4': [1.00, 0.97, 0.85, 0.65, 0.44, 0.29, 0.21, 0.17, 0.15], 'S4D2 (#3)': [1.00, 0.95, 0.84, 0.69, 0.53, 0.42, 0.33, 0.22, 0.10],
        'S5D4': [1.00, 0.96, 0.83, 0.64, 0.44, 0.32, 0.26, 0.23, 0.15]}
R = json.load(open('layout_results.json')); E = np.arange(91); dOm = 2 * np.pi * (np.cos(np.radians(E[:-1])) - np.cos(np.radians(E[1:])))
def norm(prof):
    p = np.array(prof); i0 = (p[:5] * dOm[:5]).sum() / dOm[:5].sum()
    out = [1.0] + [(p[a - 1:a + 1] * dOm[a - 1:a + 1]).sum() / dOm[a - 1:a + 1].sum() / i0 for a in A[1:]]
    return np.array(out), p / i0
fig, ax = plt.subplots(figsize=(4.6, 3.3))
for k, v in MEAS.items(): ax.plot(A, v, 'o', ms=3.5, mfc='none', color='#2a78d6', alpha=0.8)
ax.plot([], [], 'o', ms=3.5, mfc='none', color='#2a78d6', label='measured, 5 hemisphere devices')
th = np.arange(90) + 0.5
for key, ls, lab in (('hemisphere|patch3mm|rectAl', '-', 'simulated, 3 × 3 mm patch'), ('hemisphere|patch10mm|rectAl', '--', 'simulated, 10 × 10 mm patch'),
                     ('bare|patch0mm|rectAl', ':', 'simulated, bare')):
    s, full = norm(R[key]['profile']); ax.plot(th, np.convolve(full, np.ones(3) / 3, 'same'), ls, color='k' if 'bare' not in key else '#7d7c76', lw=1.1, label=lab + f" (PE {R[key]['PE']:.2f})")
ax.plot(th, np.cos(np.radians(th)), color='#bbb', lw=0.6, label='Lambertian')
ax.set_xlim(0, 85); ax.set_ylim(0, 1.15); ax.set_xlabel('angle (°)'); ax.set_ylabel('I(θ) / I(0)'); ax.legend(fontsize=6, frameon=False, loc='lower left')
ax.spines[['top', 'right']].set_visible(False); ax.set_title('Hemispherical MLA: measured vs simulated (no fitting)', fontsize=8, loc='left')
fig.savefig('hemi_profile_compare.png', dpi=220, bbox_inches='tight')
s3, _ = norm(R['hemisphere|patch3mm|rectAl']['profile']); m = np.mean(list(MEAS.values()), 0)
print('sim 3mm :', np.round(s3, 2)); print('meas mean:', np.round(m, 2)); print('RMSE 10-80:', round(float(np.sqrt(np.mean((s3[1:] - m[1:]) ** 2))), 3))
