"""Layout simulation (3 x 3 mm patch, actual substrate layout) for the MLA fabrication variants."""
import sys, json, numpy as np
import simulate as SM, layout as LY
dv = np.load('device.npz'); st = np.load('layout_stacks.npz'); L = dv['lam']; s = dv['s']; F = dv['F']; fsub = dv['f_sub']
Q = SM.Q; wq = Q * s * F * fsub / (Q * np.sum(s * F) + 1 - Q)
PL = SM.planar_kernel(); E = np.arange(91); dOm = 2 * np.pi * (np.cos(np.radians(E[:-1])) - np.cos(np.radians(E[1:])))
BTn = lambda K: K['BT'] / np.maximum(K['p'][:, None], 1e-15)
A10 = np.arange(10, 90, 10); res = {}
for nm in sys.argv[1:]:
    K = SM.lens_kernel(nm); air = np.zeros((len(L), 90))
    for k in range(len(L)):
        R = {'pix': dv['R_pix'][k], **{n: st[n][k] for n in st.files}}
        eS, eP, sd = LY.mc(K, PL, R, dv['P'][k], 1.4, 1.5, True, seed=k)
        air[k] = eS @ BTn(K) + eP @ BTn(PL)
    A = wq[:, None] * air; prof = A.sum(0) / dOm; i0 = (prof[:5] * dOm[:5]).sum() / dOm[:5].sum()
    eqe = float(A.sum()); fr = float(np.pi * i0)
    rel = [float((prof[a - 1:a + 1] * dOm[a - 1:a + 1]).sum() / dOm[a - 1:a + 1].sum() / i0) for a in A10]
    res[nm] = dict(EQE=eqe, front_EQE=fr, PE=eqe / fr, rel=rel, profile=prof.tolist())
    print(f"{nm:12s} EQE {eqe*100:.1f}% (x0.955 {eqe*95.5:.1f})  front {fr*100:.1f}%  PE {eqe/fr:.2f}  I/I0 10-80: " + ' '.join(f'{v:.2f}' for v in rel), flush=True)
json.dump(res, open('layout_variants_' + '_'.join(sys.argv[1:]) + '.json', 'w'))
