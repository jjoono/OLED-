"""Monte Carlo uncertainty of the spectral metrics (Table S14): ray noise (same geometry, other ray seeds) and
realization spread (new geometry, same statistics).  Inputs: ../multipass_reset/random_sweep_unc/, data/slabs_unc.npz."""
import os, glob, json, numpy as np
src = open('oe_figs.py').read(); exec(src[:src.index('LIB = library()')])
LIB = library()


def metrics(BT, BR):
    pT = BT.sum(1); l2 = mixing.lam2_cut(BR); pinf, p85 = mixing.pinf_cut(BR); out = dict(p=float(pT @ (W1 / W1.sum())), lam2=l2,
        Nmix=1 / (1 - l2), pinf=pinf, xi=pinf / p85)
    for dev in DEVS:
        P = DV[f'{dev}_P'][75].reshape(90, 2).sum(1); R = DV[f'{dev}_R'][75].reshape(90, 2).mean(1)
        e = mixing.full(pT, BR, R, P); out[f'eta_{dev}'] = e
        if dev == 'Al': out['err_Al'] = (mixing.reduced(pT, W1, R, P, 1 - l2) - e) * 100
    return out


rows = {}
U = '../multipass_reset/random_sweep_unc'
for nm in sorted({os.path.basename(f).split('__')[0] for f in glob.glob(f'{U}/*.npz')}):
    base = metrics(LIB[nm][1], LIB[nm][2])
    ray = [base] + [metrics(z['BT'], z['BR'].T) for z in (np.load(f) for f in sorted(glob.glob(f'{U}/{nm}__ray*.npz')))]
    real = [base] + [metrics(z['BT'], z['BR'].T) for z in (np.load(f) for f in sorted(glob.glob(f'{U}/{nm}__real*.npz')))]
    rows[nm] = (base, ray, real)
SLU = np.load('data/slabs_unc.npz')
for key in ('g0.8_t1.0', 'g0.0_t0.25'):
    nm = 'scatter ' + key; base = metrics(LIB[nm][1], LIB[nm][2]); ray = [base]
    for sd in (8, 9):
        SL0 = dict(SL); SL0[f'{key}_w1'] = SLU[f'{key}_w1__seed{sd}']
        globals()['SL'] = SL0; L2 = library(); ray.append(metrics(L2[nm][1], L2[nm][2]))
    globals()['SL'] = np.load('data/slabs.npz'); rows[nm] = (base, ray, [base])
KEYS = ('p', 'lam2', 'Nmix', 'pinf', 'xi', 'eta_Al', 'eta_Ag', 'eta_LL15', 'err_Al')
res = {}
for nm, (base, ray, real) in rows.items():
    res[nm] = {k: dict(value=base[k], sd_ray=float(np.std([r[k] for r in ray], ddof=1)) if len(ray) > 1 else None,
                       sd_real=float(np.std([r[k] for r in real], ddof=1)) if len(real) > 1 else None, n_ray=len(ray), n_real=len(real)) for k in KEYS}
json.dump(res, open('data/uncertainty.json', 'w'), indent=1)
fmt = lambda d, k, nd: f"{d[k]['value']:.{nd}f} ± {d[k]['sd_ray']:.{nd}f}" + (f" ({d[k]['sd_real']:.{nd}f})" if d[k]['sd_real'] is not None else '')
print('| Structure | n (ray / real.) | p | λ₂ | N_mix | p∞ | ξ | η Al | η Ag | η LL | err. Al (pp) |'); print('|---|---|---|---|---|---|---|---|---|---|---|')
for nm, d in res.items():
    print(f"| {nm} | {d['p']['n_ray']} / {d['p']['n_real'] if d['p']['n_real'] > 1 else '–'} | " + ' | '.join(fmt(d, k, n) for k, n in (('p', 3), ('lam2', 3), ('Nmix', 2), ('pinf', 3), ('xi', 3), ('eta_Al', 4), ('eta_Ag', 4), ('eta_LL15', 4), ('err_Al', 2))) + ' |')
