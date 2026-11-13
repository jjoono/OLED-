"""Textured exit faces on high-index slabs: hemispherical MLA and Gaussian wavy surfaces (rms slope 0.4, 1.0)
traced at n = 1.51, 2.4, 3.5 with the height-field tracers (index-matched texture on the slab).
Saves highindex/<name>_n<n>.npz with BT[in, air], BR[in, sub], lost."""
import os, time, numpy as np
import hftrace as T, hftrace_gen as G
OUT = 'highindex'; os.makedirs(OUT, exist_ok=True)
rho = np.linspace(0, 1, 2001); hemi = np.sqrt(np.clip(1 - rho ** 2, 0, None))
for n in (1.51, 2.4, 3.5):
    mev = 80 if n < 2 else 400
    for name, fn in (('hemiMLA', lambda: T.bsdf(rho, hemi, n=n, N=20000, max_ev=mev)),
                     ('wavy0.4', lambda: G.bsdf_wave(0.4, n=n, N=10000, max_ev=mev)),
                     ('wavy1.0', lambda: G.bsdf_wave(1.0, n=n, N=10000, max_ev=mev))):
        f = f'{OUT}/{name}_n{n}.npz'
        if os.path.exists(f): continue
        t = time.time(); BT, BR, lost = fn(); np.savez(f, BT=BT, BR=BR, lost=lost)
        print(f'{name} n={n}: {time.time()-t:.0f}s, lost max {lost.max():.4f}, escape(normal) {BT[0].sum():.3f}', flush=True)
