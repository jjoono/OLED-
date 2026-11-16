"""Azimuth-resolved kernels of the hexagonal close-packed hemispherical and freeform lenslet arrays (n = 1.51,
index-matched). Writes labdev/kernel_<name>.npz with BT[in, air], K[in, sub, dphi], lost, KT[in, air, dphi]."""
import sys, time, numpy as np
sys.path.insert(0, '../multipass_reset'); sys.path.insert(0, '.')
import hftrace_az as HA
from profile import profile
rho = np.linspace(0, 1, 4001)
cases = {'hemisphere': (rho, np.sqrt(np.clip(1 - rho ** 2, 0, None))), 'freeform': profile()}
for nm in sys.argv[1:] or cases:
    t = time.time(); BT, K, lost, KT = HA.kernel(*cases[nm], N=20000, max_ev=300, transmitted=True)
    np.savez(f'kernel_{nm}.npz', BT=BT, K=K, lost=lost, KT=KT)
    w = np.diff(np.sin(np.radians(np.arange(91))) ** 2)
    print(f'{nm}: {time.time()-t:.0f}s, p {np.sum(BT.sum(1) * w):.4f}, lost max {lost.max():.4f}', flush=True)
