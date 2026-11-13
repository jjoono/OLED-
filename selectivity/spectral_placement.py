"""Spectrally weighted comparison of the mixing modes, and a partially diffuse electrode.
Weighting as in the main text: sum_lam s(lam) eta_sub(lam) eta(lam) / sum_lam s(lam) eta_sub(lam), 400-700 nm every 10 nm,
Ir(dmppy-ph)2tmd spectrum, device recomputed at each wavelength, structure matrices fixed.
Modes: actual (specular device), reflective mixing by the structure (eta_mix), reflective mixing by the device
(rank-one kernel D), transmissive randomizer (= fully randomized model).
Partially diffuse electrode: D_m = (1 - m) diag R + m D (reciprocal for every m).
Writes data/spectral_placement.json."""
import json, numpy as np, sel, mixing
src = open('oe_figs.py').read(); exec(src[:src.index('# ------------------------------------------------------------------ Fig. 1')])
w = W1 / W1.sum(); I = np.eye(90)
J = np.arange(0, 301, 10); LAMs = 400.0 + J; K2 = J // 2           # devices.npz is on a 2 nm grid
s = sel.D['s_Ir(dmppy-ph)2tmd'][J] * LAMs; s /= s.sum()

def modes(BT, BR, P, R, ms=()):
    pT = BT.sum(1)
    eta = float(pT @ np.linalg.solve(I - R[:, None] * BR, P))
    exit_ = mixing.reduced(pT, W1, R, P, 1.0)
    D = np.outer(R * w, R) / (w * R).sum()
    dev = float(pT @ np.linalg.solve(I - D @ BR, P))
    rand = mixing.closed(pT, W1, R)
    part = [float(pT @ np.linalg.solve(I - ((1 - m) * np.diag(R) + m * D) @ BR, P)) for m in ms]
    return eta, exit_, dev, rand, part

MS = (0.0, 0.1, 0.25, 0.5, 1.0); SEL = ('planar', 'MLA hemi_total', 'A_shape0.7', 'W_slope0.2', 'W_slope1.0', 'scatter g0.8_t1.0')
out = {'structures': {}, 'partial': {}, 'm': list(MS)}
for dev in DEVS:
    wt = s * DV[f'{dev}_eta_sub'][K2]; wt = wt / wt.sum()
    for nm, (c, BT, BR) in LIB.items():
        acc = np.zeros(4); accp = np.zeros(len(MS))
        for k, kk in enumerate(K2):
            P = DV[f'{dev}_P'][kk].reshape(90, 2).sum(1); R = DV[f'{dev}_R'][kk].reshape(90, 2).mean(1)
            e, x, d, r, part = modes(BT, BR, P, R, MS if nm in SEL else ())
            acc += wt[k] * np.array([e, x, d, r])
            if nm in SEL: accp += wt[k] * np.array(part)
        out['structures'][f'{nm}|{dev}'] = dict(cls=c, eta=acc[0], exit_mix=acc[1], dev_mix=acc[2], rand=acc[3])
        if nm in SEL: out['partial'][f'{nm}|{dev}'] = accp.tolist()
json.dump(out, open('data/spectral_placement.json', 'w'), indent=1)
S = out['structures']
for dev in DEVS:
    r = [v for k, v in S.items() if k.endswith('|' + dev) and v['cls'] != 'planar']
    ge = np.array([v['exit_mix'] - v['eta'] for v in r]) * 100; gd = np.array([v['dev_mix'] - v['eta'] for v in r]) * 100
    gr = np.array([v['rand'] - v['eta'] for v in r]) * 100; dr = np.array([v['dev_mix'] - v['rand'] for v in r]) * 100
    pl = S['planar|' + dev]
    print(f"{dev}: escape-interface gain {ge.min():+.1f}..{ge.max():+.1f} | device gain {gd.min():+.1f}..{gd.max():+.1f} (median {np.median(gd):.1f}) | transmissive gain {gr.min():+.1f}..{gr.max():+.1f} | device-rand {dr.min():+.1f}..{dr.max():+.1f}")
    print(f"   planar: {pl['eta']:.3f} / {pl['exit_mix']:.3f} / {pl['dev_mix']:.3f} / {pl['rand']:.3f}")
for k, v in out['partial'].items():
    print(k, ' '.join(f'{x:.3f}' for x in v))
