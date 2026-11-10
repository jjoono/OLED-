"""BSDF of a hexagonal close-packed microlens array in which every lenslet may have its own axisymmetric
profile: a periodic K x K supercell of lenslets, profile index assign[i mod K, j mod K].  Same conventions,
tracer and output format as hftrace.py (BT[i_in, j_air], BR[i_in, j_sub], 1 deg bins, index-matched lens).

Profiles are given as an array H[n_prof, NT] on rho = linspace(0, 1, NT) (height in lens radii)."""
import numpy as np
from numba import njit, prange
from hftrace import NT, SQ3

A = 2.0                                   # pitch in lens radii (close-packed)


@njit(cache=True, inline='always')
def _cell(x, y):
    j = 2.0 * y / (A * SQ3)
    i = x / A - j / 2.0
    i0 = round(i); j0 = round(j)
    best = 1e18; bx = 0.0; by = 0.0; bi = 0; bj = 0
    for di in range(-1, 2):
        for dj in range(-1, 2):
            ii = i0 + di; jj = j0 + dj
            cx = ii * A + jj * A / 2.0; cy = jj * A * SQ3 / 2.0
            d2 = (x - cx) ** 2 + (y - cy) ** 2
            if d2 < best:
                best = d2; bx = cx; by = cy; bi = int(ii); bj = int(jj)
    return bx, by, bi, bj


@njit(cache=True, inline='always')
def _surf(x, y, H, D, assign, K):
    cx, cy, ci, cj = _cell(x, y)
    ux = x - cx; uy = y - cy
    rho = np.sqrt(ux * ux + uy * uy)
    if rho >= 1.0:
        return 0.0, 0.0, 0.0, 0.0
    p = assign[ci % K, cj % K]
    f = rho * (NT - 1); k = int(f)
    if k >= NT - 1:
        k = NT - 2
    w = f - k
    h = H[p, k] * (1 - w) + H[p, k + 1] * w
    d = D[p, k] * (1 - w) + D[p, k + 1] * w
    if rho > 1e-12:
        ux /= rho; uy /= rho
    else:
        ux = 0.0; uy = 0.0
    return h, d, ux, uy


@njit(cache=True)
def _trace(th, phi, sx, sy, n, H, D, assign, K, hmax, lip, max_ev):
    dx = np.sin(th) * np.cos(phi); dy = np.sin(th) * np.sin(phi); dz = np.cos(th)
    px = sx; py = sy; pz = -1e-6
    glass = True
    kstep = 1.0 / np.sqrt(1.0 + lip * lip)
    for ev in range(max_ev):
        # ---- march to the next interface ----
        hit = False
        for it in range(400000):
            if glass and dz < 0.0 and pz < 0.0:
                return 1, np.arccos(min(1.0, -dz))          # back into substrate
            if (not glass) and dz > 0.0 and pz > hmax:
                return 0, np.arccos(min(1.0, dz))           # out into air
            if glass and pz < -1e-9 and dz > 0.0:           # jump up to z=0
                t = (-1e-9 - pz) / dz
                px += t * dx; py += t * dy; pz += t * dz
            if (not glass) and pz > hmax + 1e-9 and dz < 0.0:
                t = (hmax + 1e-9 - pz) / dz
                px += t * dx; py += t * dy; pz += t * dz
            h, d, ux, uy = _surf(px, py, H, D, assign, K)
            f = pz - h
            step = abs(f) * kstep
            if step < 2e-5:
                step = 2e-5
            qx = px + step * dx; qy = py + step * dy; qz = pz + step * dz
            h2, d2, ux2, uy2 = _surf(qx, qy, H, D, assign, K)
            f2 = qz - h2
            inside_now = f < 0.0
            if (f2 < 0.0) != inside_now:
                lo = 0.0; hi = step
                for b in range(40):
                    mid = 0.5 * (lo + hi)
                    hm, dm, uxm, uym = _surf(px + mid * dx, py + mid * dy, H, D, assign, K)
                    fm = pz + mid * dz - hm
                    if (fm < 0.0) == inside_now:
                        lo = mid
                    else:
                        hi = mid
                px += hi * dx; py += hi * dy; pz += hi * dz
                hit = True
                break
            px = qx; py = qy; pz = qz
        if not hit:
            return 2, 0.0
        h, d, ux, uy = _surf(px, py, H, D, assign, K)
        nx = -d * ux; ny = -d * uy; nz = 1.0
        nn = np.sqrt(nx * nx + ny * ny + nz * nz)
        nx /= nn; ny /= nn; nz /= nn
        n1 = n if glass else 1.0
        n2 = 1.0 if glass else n
        ci = -(dx * nx + dy * ny + dz * nz)
        if ci < 0.0:
            nx = -nx; ny = -ny; nz = -nz; ci = -ci
        eta = n1 / n2
        k = 1.0 - eta * eta * (1.0 - ci * ci)
        if k < 0.0:
            Rf = 1.0; ct = 0.0
        else:
            ct = np.sqrt(k)
            rs = (n1 * ci - n2 * ct) / (n1 * ci + n2 * ct)
            rp = (n1 * ct - n2 * ci) / (n1 * ct + n2 * ci)
            Rf = 0.5 * (rs * rs + rp * rp)
        if np.random.random() < Rf:
            dx += 2 * ci * nx; dy += 2 * ci * ny; dz += 2 * ci * nz
        else:
            dx = eta * dx + (eta * ci - ct) * nx
            dy = eta * dy + (eta * ci - ct) * ny
            dz = eta * dz + (eta * ci - ct) * nz
            glass = not glass
        m = np.sqrt(dx * dx + dy * dy + dz * dz)
        dx /= m; dy /= m; dz /= m
        px += 1e-7 * dx; py += 1e-7 * dy; pz += 1e-7 * dz
    return 2, 0.0



@njit(cache=True, parallel=True)
def _batch(th, N, n, H, D, assign, K, hmax, lip, max_ev, seed):
    oc = np.empty(N, np.int64); to = np.empty(N)
    for i in prange(N):
        np.random.seed(seed + i)
        u = np.random.random() * K; w = np.random.random() * K       # uniform over the K x K supercell
        sx = u * A + w * A / 2.0; sy = w * A * SQ3 / 2.0
        phi = 2 * np.pi * np.random.random()
        a, b = _trace(th, phi, sx, sy, n, H, D, assign, K, hmax, lip, max_ev)
        oc[i] = a; to[i] = b
    return oc, to


def profiles_to_tables(H):
    r = np.linspace(0, 1, NT)
    D = np.gradient(H, r, axis=1)
    return np.ascontiguousarray(H), np.ascontiguousarray(D)


def bsdf(H, assign, n=1.51, N=40000, max_ev=60, seed=1):
    H, D = profiles_to_tables(np.atleast_2d(H)); assign = np.ascontiguousarray(assign, dtype=np.int64); K = assign.shape[0]
    hmax = float(H.max()) + 1e-6; lip = float(np.abs(D).max())
    edges = np.linspace(0, 90, 91)
    BT = np.zeros((90, 90)); BR = np.zeros((90, 90)); lost = np.zeros(90)
    for i in range(90):
        oc, to = _batch(np.radians(i + 0.5), N, n, H, D, assign, K, hmax, lip, max_ev, seed * 7919 + i * 100003)
        d = np.degrees(to)
        BT[i] = np.histogram(d[oc == 0], edges)[0] / N
        BR[i] = np.histogram(d[oc == 1], edges)[0] / N
        lost[i] = np.mean(oc == 2)
    return BT, BR, lost


def spline_profile(xv):
    """Profile of the user's spline parameterisation: x2..x6 (0-1, sorted), y2..y6 (0-1.5), stretch_Z."""
    from scipy.interpolate import PchipInterpolator
    px = np.concatenate(([0.0], np.sort(xv[0:5]), [1.0])); py = np.concatenate(([1.0], xv[5:10], [0.0]))
    px, u = np.unique(px, return_index=True); py = py[u]
    r = np.linspace(0, 1, NT)
    return PchipInterpolator(px, py)(r) * xv[10]
