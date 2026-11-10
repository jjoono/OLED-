"""Height-field BSDF tracer for two random surface families, sharing the validated ray routine of hftrace.py:

  'lens'  : union of axisymmetric lenslets on a jittered hexagonal K x K supercell.  Lenslet k has centre
            (lattice site + jitter), radius r_k and profile H[p_k] (height in units of its own radius); where
            lenslets overlap the higher surface is taken, gaps are flat (z = 0).
  'wave'  : periodic Gaussian-correlated random surface  h = h0 + sum a_m cos(g_m . x + phi_m)  on an L x L box.

Ray optics is scale free, so only shapes and slopes matter.  Output as hftrace.bsdf: BT[i_in, j_air],
BR[i_in, j_sub], lost[i_in] (1 deg bins), index-matched (n) medium below the surface."""
import inspect, numpy as np
from numba import njit, prange
import hftrace
from hftrace import NT, SQ3

A = 2.0


@njit(cache=True, inline='always')
def _surf_lens(x, y, P0, P1, P2, P3, P4, P5):
    # P0: H[n_prof, NT]; P1: D[n_prof, NT]; P2: site table [K, K, 4] = (jx, jy, radius, profile id); P3: [K]; P4: [neighbour ring]
    H = P0; D = P1; S = P2; K = P3.shape[0]
    j = 2.0 * y / (A * SQ3); i = x / A - j / 2.0
    i0 = int(np.floor(i + 0.5)); j0 = int(np.floor(j + 0.5))
    bh = 0.0; bd = 0.0; bux = 0.0; buy = 0.0
    ring = int(P4[0])
    for di in range(-ring, ring + 1):
        for dj in range(-ring, ring + 1):
            ii = i0 + di; jj = j0 + dj
            si = ii % K; sj = jj % K
            cx = ii * A + jj * A / 2.0 + S[si, sj, 0]; cy = jj * A * SQ3 / 2.0 + S[si, sj, 1]
            r = S[si, sj, 2]
            ux = x - cx; uy = y - cy
            rho = np.sqrt(ux * ux + uy * uy) / r
            if rho >= 1.0:
                continue
            p = int(S[si, sj, 3])
            f = rho * (NT - 1); k = int(f)
            if k >= NT - 1:
                k = NT - 2
            w = f - k
            h = r * (H[p, k] * (1 - w) + H[p, k + 1] * w)
            if h > bh:
                d = D[p, k] * (1 - w) + D[p, k + 1] * w
                nr = rho * r
                bh = h; bd = d
                if nr > 1e-12:
                    bux = ux / nr; buy = uy / nr
                else:
                    bux = 0.0; buy = 0.0
    return bh, bd, bux, buy


@njit(cache=True, inline='always')
def _surf_wave(x, y, P0, P1, P2, P3, P4, P5):
    # P0: amplitudes, P1: gx, P2: gy, P3: phases, P4: [h0], P5 unused
    h = P4[0]; hx = 0.0; hy = 0.0
    for m in range(P0.shape[0]):
        a = P1[m] * x + P2[m] * y + P3[m]
        h += P0[m] * np.cos(a)
        s = -P0[m] * np.sin(a)
        hx += s * P1[m]; hy += s * P2[m]
    g = np.sqrt(hx * hx + hy * hy)
    if g > 1e-14:
        return h, g, hx / g, hy / g
    return h, 0.0, 0.0, 0.0


def _make_trace(surf):
    src = inspect.getsource(hftrace._trace.py_func)
    src = src[src.index('def _trace'):]
    src = src.replace('def _trace(th, phi, sx, sy, n, htab, dtab, hmax, lip, max_ev):',
                      'def _trace(th, phi, sx, sy, n, P0, P1, P2, P3, P4, P5, hmax, lip, max_ev):')
    for a in ('px, py', 'qx, qy', 'px + mid * dx, py + mid * dy'):
        src = src.replace(f'_surf({a}, htab, dtab)', f'_SURF({a}, P0, P1, P2, P3, P4, P5)')
    assert 'htab' not in src
    ns = {'np': np, '_SURF': surf}
    exec(src, ns)
    return njit(ns['_trace'])


_trace_lens = _make_trace(_surf_lens)
_trace_wave = _make_trace(_surf_wave)


def _make_batch(trace):
    @njit(parallel=True)
    def batch(th, N, n, P0, P1, P2, P3, P4, P5, hmax, lip, max_ev, seed, bx1, by1, bx2, by2):
        oc = np.empty(N, np.int64); to = np.empty(N)
        for i in prange(N):
            np.random.seed(seed + i)
            u = np.random.random(); w = np.random.random()
            sx = u * bx1 + w * bx2; sy = u * by1 + w * by2
            phi = 2 * np.pi * np.random.random()
            a, b = trace(th, phi, sx, sy, n, P0, P1, P2, P3, P4, P5, hmax, lip, max_ev)
            oc[i] = a; to[i] = b
        return oc, to
    return batch


_batch_lens = _make_batch(_trace_lens)
_batch_wave = _make_batch(_trace_wave)


def _run(batch, P, box, hmax, lip, n, N, max_ev, seed):
    edges = np.linspace(0, 90, 91)
    BT = np.zeros((90, 90)); BR = np.zeros((90, 90)); lost = np.zeros(90)
    for i in range(90):
        oc, to = batch(np.radians(i + 0.5), N, n, *P, hmax, lip, max_ev, seed * 7919 + i * 100003, *box)
        d = np.degrees(to)
        BT[i] = np.histogram(d[oc == 0], edges)[0] / N
        BR[i] = np.histogram(d[oc == 1], edges)[0] / N
        lost[i] = np.mean(oc == 2)
    return BT, BR, lost


def bsdf_lens(H, sites, n=1.51, N=20000, max_ev=80, seed=1):
    """H[n_prof, NT] profiles (height / radius); sites[K, K, 4] = (jitter x, jitter y, radius, profile id) in units of
    the nominal radius (pitch 2)."""
    H = np.ascontiguousarray(np.atleast_2d(H), float); D = np.ascontiguousarray(np.gradient(H, np.linspace(0, 1, NT), axis=1))
    S = np.ascontiguousarray(sites, float); K = S.shape[0]
    hmax = float((H.max(1)[S[..., 3].astype(int)] * S[..., 2]).max()) + 1e-6
    lip = float(np.abs(D).max())
    reach = float((S[..., 2] + np.hypot(S[..., 0], S[..., 1])).max())
    ring = 1 if reach <= 1.6 else 2                                   # lattice neighbours that can cover a point
    P = (H, D, S, np.zeros(K), np.array([float(ring)]), np.zeros(1))
    box = (K * A, 0.0, K * A / 2.0, K * A * SQ3 / 2.0)
    return _run(_batch_lens, P, box, hmax, lip, n, N, max_ev, seed)


def wave_surface(rms_slope, L=12.0, lc=1.0, kmax=3.0, seed=0):
    """Gaussian-correlated random surface (correlation length lc) scaled to the requested rms slope."""
    rng = np.random.default_rng(seed); g = []
    M = int(np.ceil(kmax * lc * L / (2 * np.pi)))
    for m in range(0, M + 1):
        for k in range(-M, M + 1):
            if m == 0 and k <= 0:
                continue
            gx, gy = 2 * np.pi * m / L, 2 * np.pi * k / L
            if np.hypot(gx, gy) * lc <= kmax:
                g.append((gx, gy))
    g = np.array(g); kk = np.hypot(g[:, 0], g[:, 1])
    amp = np.exp(-(kk * lc) ** 2 / 8.0) * rng.rayleigh(1.0, len(g))
    ph = rng.uniform(0, 2 * np.pi, len(g))
    slope_rms = np.sqrt(0.5 * np.sum((amp * kk) ** 2))
    amp *= rms_slope / slope_rms
    return amp, g[:, 0], g[:, 1], ph, L


def bsdf_wave(rms_slope, n=1.51, N=20000, max_ev=80, seed=1, surf_seed=0, **kw):
    amp, gx, gy, ph, L = wave_surface(rms_slope, seed=surf_seed, **kw)
    h0 = float(np.abs(amp).sum()) + 1e-6
    xs = np.linspace(0, L, 400); X, Y = np.meshgrid(xs, xs)
    a = gx[:, None, None] * X + gy[:, None, None] * Y + ph[:, None, None]
    hx = -(amp[:, None, None] * gx[:, None, None] * np.sin(a)).sum(0); hy = -(amp[:, None, None] * gy[:, None, None] * np.sin(a)).sum(0)
    lip = 1.3 * float(np.hypot(hx, hy).max())                      # sampled maximum slope with margin
    P = (amp, gx, gy, ph, np.array([h0]), np.zeros(1))
    box = (L, 0.0, 0.0, L)
    return _run(_batch_wave, P, box, 2 * h0, lip, n, N, max_ev, seed)
