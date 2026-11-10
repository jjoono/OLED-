"""Markdown table and summary for direct ray tracing and blind structures (reads blind/results.json)."""
import json, numpy as np
R = json.load(open('blind/results.json'))
LAB = {'hemi_total': 'hemispherical MLA (library)', 'free_40-60': 'freeform MLA, 40–60° (library)', 'W_slope0.2': 'wavy, rms slope 0.2 (library)',
       'W_slope1.0': 'wavy, rms slope 1.0 (library)', 'cap60': 'spherical caps, 60° contact', 'cap35': 'spherical caps, 35° contact',
       'paraboloid': 'paraboloids', 'cone45': 'cones, 45° side', 'cone60': 'cones, 60° side', 'gauss_bump': 'Gaussian bumps',
       'mixed_cone_hemi': 'random mix of cones and hemispheres', 'ring_wave0.6': 'ring-spectrum wavy, slope 0.6', 'ring_wave1.2': 'ring-spectrum wavy, slope 1.2'}
print('| Structure | $p$ | $N_{\\mathrm{mix}}$ | $\\xi$ | $\\eta_{\\mathrm{sta}}$ full series | direct, $(\\theta,\\phi)$ | direct, $\\phi$ re-randomized | reduced − full (pp) | memoryless − full (pp) |')
print('|---|---|---|---|---|---|---|---|---|')
for k, r in R.items():
    a = r['Al_550']; d = r.get('direct')
    ds = (f"{d['eta']:.4f} ± {d['se']:.4f}", f"{d['eta_randphi']:.4f}") if d else ('—', '—')
    print(f"| {LAB[k]} | {r['p']:.3f} | {r['N_mix']:.1f} | {r['xi']:.2f} | {a['full']:.4f} | {ds[0]} | {ds[1]} | {(a['red'] - a['full']) * 100:+.2f} | {(a['memoryless'] - a['full']) * 100:+.2f} |")
bl = [r for r in R.values() if r['tag'] == 'blind']
e = np.array([abs(r[f'{dv}_{l}']['red'] - r[f'{dv}_{l}']['full']) for r in bl for dv in ('Al', 'Ag', 'LL15') for l in (500, 550, 600)]) * 100
ew = np.array([abs(r[f'{dv}_{l}']['red_w'] - r[f'{dv}_{l}']['full']) for r in bl for dv in ('Al', 'Ag', 'LL15') for l in (500, 550, 600)]) * 100
em = np.array([abs(r[f'{dv}_{l}']['memoryless'] - r[f'{dv}_{l}']['full']) for r in bl for dv in ('Al', 'Ag', 'LL15') for l in (500, 550, 600)]) * 100
print(f'\nblind: {len(bl)} structures, {e.size} cases; reciprocal mean {e.mean():.2f} max {e.max():.2f}; q~w mean {ew.mean():.2f} max {ew.max():.2f}; memoryless mean {em.mean():.2f} max {em.max():.2f}')
dd = [(k, r['direct']['eta'] - r['Al_550']['full'], r['direct']['eta_randphi'] - r['Al_550']['full'], r['direct']['se']) for k, r in R.items() if 'direct' in r]
for k, a, b, se in dd:
    print(f'direct {k}: azimuth-tracked - series {a * 100:+.2f} pp, phi-randomized - series {b * 100:+.2f} pp (se {se * 100:.2f})')
