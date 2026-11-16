"""Lattice-resolved (absolute-azimuth) kernels for the hemispherical and green freeform MLAs (revtrace.kernel_abs)."""
import sys, time, numpy as np
sys.path.insert(0, '../multipass_reset'); import revtrace as RT
a = np.linspace(0, np.pi / 2, 400)
P = np.load('green_profile_rz.npy'); P = np.r_[P, [[10.2, 0.0]]] / 10.2
from profile import profile
rho, hh = profile(); ii = np.unique(np.r_[np.linspace(0, 3990, 360).astype(int), np.arange(3990, 4001)])   # orange freeform, ~370 nodes
cases = {'hemisphere': (np.sin(a), np.cos(a)), 'green_freeform': (P[:, 0], P[:, 1]), 'freeform': (rho[ii], hh[ii])}
w = np.diff(np.sin(np.radians(np.arange(91))) ** 2)
for nm in sys.argv[1:] or cases:
    t = time.time(); KTa, Ka = RT.kernel_abs(*cases[nm], N=int(60000))
    np.savez_compressed(f'kernelabs_{nm}.npz', KTa=KTa, Ka=Ka)
    p = KTa.sum((2, 3)).mean(1); ref = np.load(f'kernel_{nm}.npz')['BT'].sum(1)
    print(f'{nm}: {time.time()-t:.0f}s p {np.sum(p*w):.4f} (azimuth-averaged kernel {np.sum(ref*w):.4f}); '
          f'max spread of p over a-bins {np.abs(KTa.sum((2,3)) - p[:, None]).max():.3f}', flush=True)
