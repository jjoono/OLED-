"""Ray-noise replicas of two scattering layers (default seed 7 is the library; seeds 8, 9 here)."""
import numpy as np, slab
out = {}
for g, tau in ((0.8, 1.0), (0.0, 0.25)):
    for sd in (8, 9):
        Ts, Rs = slab.slab(tau / (1 - g), g, 1.0, N=20000, seed=sd); out[f'g{g}_t{tau}_w1__seed{sd}'] = np.stack([Ts, Rs])
np.savez('data/slabs_unc.npz', **out); print('done')
