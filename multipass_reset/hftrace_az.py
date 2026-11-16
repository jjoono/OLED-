"""Azimuth-resolved return kernel of the hexagonal lenslet array (same physics as hftrace.py).
For each incident substrate polar bin i (azimuth uniform over the cell), records the transmitted fraction and the
joint distribution of the returned polar bin r and the change of the lateral propagation azimuth
dphi = phi_out - phi_in (0: lateral motion continues, pi: reversed).
Returns BT[i] (escape fraction), K[i, r, a] (returned fraction into polar bin r and azimuth-change bin a, NA bins)."""
import inspect, numpy as np
from numba import njit, prange
import hftrace as T

src = inspect.getsource(T._trace.py_func)
src = src[src.index('def _trace('):].replace('def _trace(', 'def _trace_az(')
src = src.replace('return 1, np.arccos(min(1.0, -dz))', 'return 1, np.arccos(min(1.0, -dz)), np.arctan2(dy, dx)')
src = src.replace('return 0, np.arccos(min(1.0, dz))', 'return 0, np.arccos(min(1.0, dz)), 0.0')
src = src.replace('return 2, 0.0', 'return 2, 0.0, 0.0')
ns = dict(np=np, _surf=T._surf)
exec(src, ns)
_trace_az = njit(cache=False)(ns['_trace_az'])


@njit(parallel=True)
def _batch_az(th, N, n, htab, dtab, hmax, lip, max_ev, seed):
    oc = np.empty(N, np.int64); to = np.empty(N); dp = np.empty(N)
    for i in prange(N):
        np.random.seed(seed + i)
        u = np.random.random(); w = np.random.random()
        sx = u * 2.0 + w * 1.0; sy = w * T.SQ3
        phi = 2 * np.pi * np.random.random()
        a, b, po = _trace_az(th, phi, sx, sy, n, htab, dtab, hmax, lip, max_ev)
        oc[i] = a; to[i] = b; dp[i] = (po - phi) % (2 * np.pi)
    return oc, to, dp


def kernel(rho, h, n=1.51, N=20000, max_ev=80, NA=36, seed=1):
    htab, dtab = T.tables(rho, h)
    hmax = float(htab.max()) + 1e-6; lip = float(np.abs(dtab).max())
    BT = np.zeros(90); K = np.zeros((90, 90, NA)); lost = np.zeros(90)
    for i in range(90):
        oc, to, dp = _batch_az(np.radians(i + 0.5), N, n, htab, dtab, hmax, lip, max_ev, seed * 7919 + i * 100003)
        BT[i] = np.mean(oc == 0); lost[i] = np.mean(oc == 2)
        m = oc == 1
        r = np.minimum((np.degrees(to[m])).astype(int), 89); a = np.minimum((dp[m] / (2 * np.pi) * NA).astype(int), NA - 1)
        np.add.at(K[i], (r, a), 1.0 / N)
    return BT, K, lost


if __name__ == '__main__':
    import time
    rho = np.linspace(0, 1, 2001); hemi = np.sqrt(np.clip(1 - rho ** 2, 0, None))
    t = time.time(); BT, K, lost = kernel(rho, hemi)
    np.savez('azkernel_hemiMLA.npz', BT=BT, K=K, lost=lost)
    print(f'{time.time()-t:.0f}s, escape-weighted p {np.sum(BT*np.diff(np.sin(np.radians(np.arange(91)))**2)):.4f}, lost max {lost.max():.4f}')
    ret = K.sum(1); w = ret.sum(1, keepdims=True)
    for i in (10, 30, 45, 60, 75):
        d = ret[i] / max(w[i, 0], 1e-12); print(i, 'dphi share  0±30°: %.2f  ±(30-150°): %.2f  180±30°: %.2f' % (d[[0, 1, 2, 33, 34, 35]].sum(), d[3:15].sum() + d[21:33].sum(), d[15:21].sum()))
