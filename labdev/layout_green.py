"""Green lab devices (251120) in the substrate layout: bare (0.7 mm), hemispherical and freeform (cup) MLAs on a
3 x 3 mm patch (1.4 mm glass), plus 10 x 10 mm for reference. Same model as layout.py with the green device optics."""
import json, numpy as np
import simulate as SM, layout as LY
import os
SFX = os.environ.get('GREEN_SFX', '')
dv = np.load(f'device_green{SFX}.npz'); st = np.load(f'layout_stacks_green{SFX}.npz'); L = dv['lam']; s = dv['s']; F = dv['F']; fsub = dv['f_sub']
Q = SM.Q; wq = Q * s * F * fsub / (Q * np.sum(s * F) + 1 - Q)
PL = SM.planar_kernel(); E = np.arange(91); dOm = 2 * np.pi * (np.cos(np.radians(E[:-1])) - np.cos(np.radians(E[1:])))
BTn = lambda K: K['BT'] / np.maximum(K['p'][:, None], 1e-15)
KER = {'bare': PL, 'hemisphere': SM.lens_kernel('hemisphere'), 'freeform': SM.lens_kernel('green_freeform')}
res = {}
for nm, t, hp in [('bare', 0.7, 0.0)] + [(n, 1.4, h) for n in ('hemisphere', 'freeform') for h in (1.5, 5.0)]:
    air = np.zeros((len(L), 90)); side = np.zeros(len(L))
    for k in range(len(L)):
        R = {'pix': dv['R_pix'][k], **{n: st[n][k] for n in st.files}}
        eS, eP, sd = LY.mc(KER[nm], PL, R, dv['P'][k], t, hp, True, seed=k)
        air[k] = eS @ BTn(KER[nm]) + eP @ BTn(PL); side[k] = sd
    A = wq[:, None] * air; prof = A.sum(0) / dOm; i0 = (prof[:5] * dOm[:5]).sum() / dOm[:5].sum()
    eqe = float(A.sum()); fr = float(np.pi * i0)
    rel = [1.0] + [float((prof[a - 3:a + 3] * dOm[a - 3:a + 3]).sum() / dOm[a - 3:a + 3].sum() / i0) for a in range(10, 90, 10)]
    key = f'{nm}|patch{2 * hp:g}mm'
    res[key] = dict(EQE=eqe, front_EQE=fr, PE=eqe / fr, side=float((wq * side).sum() / wq.sum()), rel=rel, profile=prof.tolist())
    print(f"{key:22s} EQE {eqe*100:.1f}%  front {fr*100:.1f}%  PE {eqe/fr:.2f}  side {res[key]['side']:.3f}  I/I0: " + ' '.join(f'{v:.2f}' for v in rel), flush=True)
json.dump(res, open(f'layout_green_results{SFX}.json', 'w'))
