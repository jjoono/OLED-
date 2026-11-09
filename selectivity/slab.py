"""Monte-Carlo transmission / reflection matrices of an index-matched scattering slab
(Henyey-Greenstein phase function g, scattering optical thickness S = mu_s d, single-scattering
albedo w).  No refractive-index steps at its faces: the slab only redistributes angles, and the
filter and the device supply all reflections.  Column i = input bin (0.5 deg), row = exit bin."""
import numpy as np
from numba import njit, prange
import sel

NB = sel.NB


@njit(cache=True)
def _hg(g, u):
    if abs(g) < 1e-6:
        return 2 * u - 1
    s = (1 - g * g) / (1 - g + 2 * g * u)
    return (1 + g * g - s * s) / (2 * g)


@njit(cache=True, parallel=True)
def _run(S, g, w, N, seed):
    Ts = np.zeros((NB, NB)); Rs = np.zeros((NB, NB))
    for i in prange(NB):
        np.random.seed(seed + i)
        for k in range(N):
            th = np.radians((i + np.random.random()) * 0.5)
            ph = 2 * np.pi * np.random.random()
            ux = np.sin(th) * np.cos(ph); uy = np.sin(th) * np.sin(ph); uz = np.cos(th)
            z = 0.0; wt = 1.0
            while True:
                d = -np.log(np.random.random() + 1e-300) / S if S > 0 else 1e9
                zn = z + d * uz
                if zn >= 1.0:
                    b = int(np.degrees(np.arccos(min(1.0, uz))) / 0.5)
                    Ts[min(b, NB - 1), i] += wt; break
                if zn <= 0.0:
                    b = int(np.degrees(np.arccos(min(1.0, -uz))) / 0.5)
                    Rs[min(b, NB - 1), i] += wt; break
                z = zn; wt *= w
                if wt < 1e-3:
                    break
                ct = _hg(g, np.random.random()); st = np.sqrt(max(0.0, 1 - ct * ct)); p2 = 2 * np.pi * np.random.random()
                if abs(uz) > 0.99999:
                    ux, uy, uz = st * np.cos(p2), st * np.sin(p2), ct * np.sign(uz)
                else:
                    dd = np.sqrt(1 - uz * uz)
                    nx = st * (ux * uz * np.cos(p2) - uy * np.sin(p2)) / dd + ux * ct
                    ny = st * (uy * uz * np.cos(p2) + ux * np.sin(p2)) / dd + uy * ct
                    nz = -st * np.cos(p2) * dd + uz * ct
                    ux, uy, uz = nx, ny, nz
    return Ts / N, Rs / N


def slab(S, g, w=1.0, N=20000, seed=7):
    return _run(float(S), float(g), float(w), int(N), int(seed))
