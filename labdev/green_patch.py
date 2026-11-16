"""Green devices at the bare-calibrated cavity (B3PyMPM 65 nm, TAPC 80 nm): azimuth-averaged EQE, Lambertian-equivalent
on-axis EQE and angular factor for 3 x 3 mm and 10 x 10 mm lens patches. Writes green_patch.json."""
import json, numpy as np, simulate as SM, layout as LY, green_joint_fit as GJ
PL = SM.planar_kernel(); KER = {'bare': PL, 'hemisphere': SM.lens_kernel('hemisphere'), 'freeform': SM.lens_kernel('green_freeform')}
BTn = lambda K: K['BT'] / np.maximum(K['p'][:, None], 1e-15)
ST, wq = GJ.stacks(65.0, 80.0); res = {}
for nm, t, hp in (('bare', 0.7, 0.0), ('hemisphere', 1.4, 1.5), ('hemisphere', 1.4, 5.0), ('freeform', 1.4, 1.5), ('freeform', 1.4, 5.0)):
    air = np.zeros((len(ST), 90))
    for j, (P, R) in enumerate(ST):
        eS, eP, sd = LY.mc(KER[nm], PL, R, P, t, hp, True, N=150000, seed=j); air[j] = eS @ BTn(KER[nm]) + eP @ BTn(PL)
    A = wq[:, None] * air; rl, fr = GJ.rel(A.sum(0) / GJ.dOm)
    res[f'{nm}|{2*hp:g}mm'] = dict(EQE=float(A.sum()), front=float(fr), PE=float(A.sum() / fr), rel=rl.tolist())
    print(nm, 2 * hp, f'EQE {A.sum()*100:.1f} front {fr*100:.1f} PE {A.sum()/fr:.3f}', flush=True)
json.dump(res, open('green_patch.json', 'w'), indent=1)
