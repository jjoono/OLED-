"""Main sweep: every stack x randomizer x angular filter (ideal cones, bare Fresnel exit, and
quarter-wave edge DBRs over material, pair count and band-edge wavelength).  Stores, per
configuration, the power inside each air cone per wavelength (per unit substrate-delivered power)."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, time, itertools
from multiprocessing import Pool
import sel

DV = np.load('data/devices.npz'); SL = np.load('data/slabs.npz')
DEV = {k: dict(eta_sub=DV[f'{k}_eta_sub'], P=DV[f'{k}_P'], R=DV[f'{k}_R'], n_sub=float(DV[f'{k}_n_sub'])) for k in sel.STACKS}
RAND = {'full': sel.mixer(1.0), 'none': sel.mixer(0.0), 'scat': tuple(SL['g0.8_t1.0_w1'])}
NSUBS = sorted({d['n_sub'] for d in DEV.values()})
MATS = list(sel.DBR_MATS); NS = (4, 6, 8, 12, 16); LB = np.arange(420.0, 821.0, 5.0)
DESIGNS = [('fresnel',)] + [('ideal', b) for b in sel.BETAS] + [('edge', m, n, lb) for m in MATS for n in NS for lb in LB]


def one(idx):
    d = DESIGNS[idx]; out = {}
    for ns in NSUBS:
        spec = 'fresnel' if d[0] == 'fresnel' else d
        Rf, Tf = sel.filter_RT(spec, ns)
        for k, dev in DEV.items():
            if dev['n_sub'] != ns:
                continue
            for r, (Ts, Rs) in RAND.items():
                out[(k, r)] = sel.cones(sel.solve(dev, Rf, Tf, Ts, Rs), ns).astype(np.float32)
    return idx, out


if __name__ == '__main__':
    t0 = time.time(); res = {}
    with Pool(4) as p:
        for c, (i, o) in enumerate(p.imap_unordered(one, range(len(DESIGNS)), chunksize=4)):
            for key, v in o.items():
                res.setdefault(key, np.zeros((len(DESIGNS), sel.LAM.size, len(sel.BETAS)), np.float32))[i] = v
            if c % 50 == 0:
                print(f'{c}/{len(DESIGNS)} [{time.time() - t0:.0f}s]', flush=True)
    np.savez_compressed('data/sweep.npz', designs=np.array([str(d) for d in DESIGNS]),
                        **{f'{k}|{r}': v for (k, r), v in res.items()})
    print('done', f'[{time.time() - t0:.0f}s]')
