"""Lens-patch width sweep in the actual substrate layout (square patch centred on the pixel, upper rectangle Al)."""
import json, numpy as np
import simulate as SM, layout as LY
dv = np.load('device.npz'); st = np.load('layout_stacks.npz'); L = dv['lam']; s = dv['s']; F = dv['F']; fsub = dv['f_sub']
Q = SM.Q; wq = Q * s * F * fsub / (Q * np.sum(s * F) + 1 - Q)
PL = SM.planar_kernel(); E = np.arange(91); dOm = 2 * np.pi * (np.cos(np.radians(E[:-1])) - np.cos(np.radians(E[1:])))
BTn = lambda K: K['BT'] / np.maximum(K['p'][:, None], 1e-15); res = {}
for nm in ('hemisphere', 'freeform'):
    K = SM.lens_kernel(nm)
    for w in (2, 3, 4, 5, 7, 10, 14, 20):
        air = np.zeros((len(L), 90))
        for k in range(len(L)):
            R = {'pix': dv['R_pix'][k], **{n: st[n][k] for n in st.files}}
            eS, eP, sd = LY.mc(K, PL, R, dv['P'][k], 1.4, w / 2, True, N=60000, seed=k)
            air[k] = eS @ BTn(K) + eP @ BTn(PL)
        A = wq[:, None] * air; prof = A.sum(0) / dOm; i0 = (prof[:5] * dOm[:5]).sum() / dOm[:5].sum()
        res[f'{nm}|{w}'] = dict(EQE=float(A.sum()), front_EQE=float(np.pi * i0), PE=float(A.sum() / (np.pi * i0)))
        print(nm, w, f"EQE {A.sum()*100:.1f} PE {res[f'{nm}|{w}']['PE']:.2f}", flush=True)
json.dump(res, open('layout_sweep.json', 'w'))
