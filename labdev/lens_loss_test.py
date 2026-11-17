"""Which loss in the attached lens glass is consistent with the measured lens devices? Three hypotheses, each scaled to
lower the orange hemispherical EQE by 14% (3 x 3 mm patch, azimuthally averaged far-field model, device.npz):
 (1) exit loss: a fraction of the light transmitted by the lens surface is lost (angle independent);
 (2) interface loss: each crossing of the oil/lens-glass interface transmits tau (angle independent);
 (3) bulk loss in the 0.7 mm lens glass, exp(-alpha * 0.7 mm / cos theta) per crossing.
Prints the EQE ratio to the lossless model and the angular factor for the hemispherical and freeform arrays."""
import json, numpy as np, simulate as SM, layout as LY
from scipy.optimize import brentq
TAN = SM.TAN; COS = np.cos(np.radians(np.arange(90) + 0.5))
dv = np.load('device.npz'); st = np.load('layout_stacks.npz'); wq = SM.Q * dv['s'] * dv['F'] * dv['f_sub'] / (SM.Q * np.sum(dv['s'] * dv['F']) + 1 - SM.Q)
PL = SM.planar_kernel(); K = {'hemisphere': SM.lens_kernel('hemisphere'), 'freeform': SM.lens_kernel('freeform')}
BTn = lambda k: k['BT'] / np.maximum(k['p'][:, None], 1e-15)
E = np.arange(91); dOm = 2 * np.pi * (np.cos(np.radians(E[:-1])) - np.cos(np.radians(E[1:])))
KS = range(0, 31, 3)


def run(nm, first, rt, N=40000):
    """first[theta]: attenuation before the first encounter; rt[theta]: extra attenuation per round trip."""
    air = 0
    for k in KS:
        R = {'pix': dv['R_pix'][k] * rt, **{n: st[n][k] * rt for n in st.files}}
        P = dv['P'][k] * first
        eS, eP, sd = LY.mc(K[nm], PL, R, P / P.sum(), 1.4, 1.5, True, N=N, seed=k)
        air = air + wq[k] * P.sum() / dv['P'][k].sum() * (eS @ BTn(K[nm]) + eP @ BTn(PL))
    eqe = air.sum(); i0 = (air[:5] / dOm[:5]).mean() if False else air[:5].sum() / dOm[:5].sum()
    return eqe, eqe / (np.pi * i0)


ref = {nm: run(nm, np.ones(90), np.ones(90)) for nm in K}
out = {'lossless': {nm: dict(EQE=ref[nm][0], PE=ref[nm][1]) for nm in K}}
print('lossless', {nm: (round(v[0] * 100, 2), round(v[1], 3)) for nm, v in ref.items()})
target = 0.86
hyp = {'exit': None,
       'interface': lambda x: (np.full(90, x), np.full(90, x ** 2)),
       'bulk': lambda x: (np.exp(-x * 0.7 / COS), np.exp(-2 * x * 0.7 / COS))}
for h, f in hyp.items():
    if h == 'exit':
        out[h] = {nm: dict(EQE=ref[nm][0] * target, PE=ref[nm][1], param=target) for nm in K}; print(h, 'ratio', target, 'PE unchanged'); continue
    g = lambda x: run('hemisphere', *f(x), N=15000)[0] / ref['hemisphere'][0] - target
    lo, hi = (0.5, 1.0) if h == 'interface' else (0.0, 2.0)
    x = brentq(g, lo, hi, xtol=1e-3) if h == 'interface' else brentq(g, 1e-4, 1.0, xtol=1e-4)
    out[h] = {}
    for nm in K:
        e, pe = run(nm, *f(x)); out[h][nm] = dict(EQE=e, PE=pe, param=x)
        print(h, f'param {x:.4f}', nm, f'EQE ratio {e / ref[nm][0]:.3f}  PE {pe:.3f} (lossless {ref[nm][1]:.3f})', flush=True)
json.dump(out, open('lens_loss_test.json', 'w'), indent=1)
