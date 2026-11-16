"""Fabrication variants of the 'hemispherical' MLA: spherical caps of height ratio a = h / r_base (0.5, 0.7) and
hemispheres with a reduced fill factor (base radius 0.85 of the cell radius, flat between lenses)."""
import sys, time, numpy as np
sys.path.insert(0, '../multipass_reset'); import hftrace_az as HA
rho = np.linspace(0, 1, 4001)
def cap(a, rb=1.0):
    R = (rb ** 2 + a ** 2) / (2 * a); h = np.sqrt(np.clip(R ** 2 - rho ** 2, 0, None)) + (a - R)
    return np.where(rho <= rb, np.clip(h, 0, None), 0.0)
V = {'cap0.5': cap(0.5), 'cap0.7': cap(0.7), 'hemiFF0.85': np.where(rho <= 0.85, np.sqrt(np.clip(0.85 ** 2 - rho ** 2, 0, None)), 0.0)}
w = np.diff(np.sin(np.radians(np.arange(91))) ** 2)
for nm in sys.argv[1:] or V:
    t = time.time(); BT, K, lost, KT = HA.kernel(rho, V[nm], N=20000, max_ev=300, transmitted=True)
    np.savez(f'kernel_{nm}.npz', BT=BT, K=K, lost=lost, KT=KT)
    print(f'{nm}: {time.time()-t:.0f}s, p {np.sum(BT.sum(1) * w):.4f}', flush=True)
