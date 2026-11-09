"""Randomization sweep at the best edge DBR of each stack (20-deg cone, Gaussian 530/30 nm and
v-DABNA spectra): idealized partial mixing m, and physical scattering slabs (g, transport
optical thickness tau* = S(1-g), albedo)."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, json, time
from multiprocessing import Pool
import sel, common

SW = common.Sweep()
DV = dict(np.load('data/devices.npz')); SL = dict(np.load('data/slabs.npz'))
SLABS = {k: tuple(v) for k, v in SL.items()}
MS = (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 1.0)
CASES = []
for st in sel.STACKS:
    for sp in ('G30', 'v-DABNA (HF)'):
        d, _ = SW.best_edge(st, 'full', common.SPECTRA[sp], 20.0, mat='TiO2/SiO2', N=16)
        CASES.append((st, sp, d))
        CASES.append((st, sp, ('ideal', 20.0)))
        CASES.append((st, sp, ('fresnel',)))


def one(case):
    st, sp, d = case
    dev = dict(eta_sub=DV[f'{st}_eta_sub'], P=DV[f'{st}_P'], R=DV[f'{st}_R'])
    ns = float(DV[f'{st}_n_sub']); spec = common.SPECTRA[sp]
    Rf, Tf = sel.filter_RT('fresnel' if d[0] == 'fresnel' else d, ns)
    out = {}
    for m in MS:
        Ts, Rs = sel.mixer(m)
        out[f'm{m}'] = sel.metrics(dev['eta_sub'], sel.cones(sel.solve(dev, Rf, Tf, Ts, Rs), ns), spec)
    for k, (Ts, Rs) in SLABS.items():
        out[k] = sel.metrics(dev['eta_sub'], sel.cones(sel.solve(dev, Rf, Tf, Ts, Rs), ns), spec)
    return dict(stack=st, spec=sp, design=str(d), res=out)


if __name__ == '__main__':
    t0 = time.time()
    with Pool(4) as p:
        res = p.map(one, CASES)
    json.dump(res, open('data/random.json', 'w'))
    print('done', len(res), time.time() - t0)
