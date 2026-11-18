"""Realistic transmissive mixing: non-absorbing Henyey-Greenstein volumes (g = 0-0.9, tau* = 0.1-4) under a planar
glass/air exit face (adding method, as the library scattering layers), compared with the ideal lossless transmissive
randomizer (= fully randomized model), complete reflective mixing at the planar interface, the planar interface without
mixing, and the hemispherical MLA. Al, Ag and low-loss devices at 550 nm. Writes data/scatter_volume.json."""
import json, numpy as np
src = open('oe_figs.py').read(); exec(src[:src.index('# ------------------------------------------------------------------ Fig. 1')])
I = np.eye(sel.NB); Rf, Tf = sel.tmm_RT([], 550.0, sel.TH, N_SUB)
def slab(Ts, Rs):
    U = np.linalg.solve(I - Rs * Rf[None, :], Ts)
    return to_1deg_air(Tf[:, None] * U), Rs + (Ts * Rf[None, :]) @ U
def coarse(BTair, BR):
    w = sel.WL; W2 = w.reshape(90, 2); f = (W2 / W2.sum(1, keepdims=True)).ravel()
    return (BTair * f[:, None]).reshape(90, 2, 90).sum(1), (BR * f[None, :]).reshape(90, 2, 90, 2).sum(axis=(1, 3))
out = {'slabs': [], 'refs': {}}
wl = sel.WL / sel.WL.sum()
for k in SL.files:
    if not k.endswith('_w1'): continue
    g = float(k.split('_')[0][1:]); t = float(k.split('_')[1][1:]); Ts, Rs = SL[k]
    back = float((Rs.sum(0) * wl).sum())                       # diffuse reflectance of the bare slab, Lambertian incidence
    BT, BR = coarse(*slab(Ts, Rs)); PINF[id(BR)] = mixing.pinf_cut(BR); row = dict(g=g, tau=t, slab_backscatter=back, p=float(BT.sum(1) @ (W1 / W1.sum())))
    for dev in DEVS:
        e = evaluate(BT, BR, dev, 75); row[dev] = dict(eta=e['eta'], rand=e['closed'], mix=e['m1'])
    out['slabs'].append(row)
for nm in ('planar', 'MLA hemi_total'):
    c, BT, BR = LIB[nm]; out['refs'][nm] = {dev: dict(eta=EV[(nm, dev, 550.0)]['eta'], rand=EV[(nm, dev, 550.0)]['closed'], mix=EV[(nm, dev, 550.0)]['m1']) for dev in DEVS}
json.dump(out, open('data/scatter_volume.json', 'w'), indent=1)
pl = out['refs']['planar']
for dev in DEVS:
    print(dev, 'planar', round(pl[dev]['eta'], 3), 'reflective mix', round(pl[dev]['mix'], 3), 'ideal transmissive', round(pl[dev]['rand'], 3), 'hemi MLA', round(out['refs']['MLA hemi_total'][dev]['eta'], 3))
    best = sorted(out['slabs'], key=lambda r: -r[dev]['eta'])[:4]
    print('   best HG volumes:', [(r['g'], r['tau'], round(r[dev]['eta'], 3), round(r['slab_backscatter'], 3)) for r in best])
    for g in (0.0, 0.5, 0.8, 0.9):
        print('   g', g, [(r['tau'], round(r[dev]['eta'], 3)) for r in out['slabs'] if r['g'] == g])
