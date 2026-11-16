"""Ray tracer for a hexagonal close-packed array of axisymmetric lenslets whose cross-section is an arbitrary polyline
(undercuts allowed), index-matched to the substrate, with flat glass/air between the lenslets.
Units: cell radius 1 (pitch 2). The lens solid occupies r <= 1 between z = 0 and the revolved polyline; the polyline
segments (r1, z1) -> (r2, z2) are cones, cylinders (r1 = r2) or annuli (z1 = z2). Interfaces: the polyline surfaces
(lens/air) and the plane z = 0 outside the circles (glass/air). Unpolarized Fresnel with a Monte-Carlo R/T choice.
kernel() returns BT[in, air], K[in, sub, dphi], lost, KT[in, air, dphi] like hftrace_az.kernel(..., transmitted=True)."""
import numpy as np
from numba import njit, prange
SQ3 = np.sqrt(3.0)


@njit(cache=True)
def _nearest(x, y):
    j = 2.0 * y / (2.0 * SQ3); i = x / 2.0 - j / 2.0
    i0 = round(i); j0 = round(j); best = 1e18; bx = 0.0; by = 0.0
    for di in range(-1, 2):
        for dj in range(-1, 2):
            cx = (i0 + di) * 2.0 + (j0 + dj) * 1.0; cy = (j0 + dj) * SQ3
            d2 = (x - cx) ** 2 + (y - cy) ** 2
            if d2 < best:
                best = d2; bx = cx; by = cy
    return bx, by, np.sqrt(best)


@njit(cache=True)
def _lens_hit(ox, oy, oz, dx, dy, dz, seg, H, tmin):
    """Smallest t > tmin of the ray with one lenslet centred at the origin. Returns t, nx, ny, nz (t = 1e30 if none)."""
    dxy2 = dx * dx + dy * dy
    if dxy2 > 1e-20:
        tca = -(ox * dx + oy * dy) / dxy2
        cx = ox + tca * dx; cy = oy + tca * dy
        if cx * cx + cy * cy > 1.0 + 1e-9:
            return 1e30, 0.0, 0.0, 0.0
    elif ox * ox + oy * oy > 1.0 + 1e-9:
        return 1e30, 0.0, 0.0, 0.0
    best = 1e30; bnx = 0.0; bny = 0.0; bnz = 0.0
    for k in range(seg.shape[0]):
        r1 = seg[k, 0]; z1 = seg[k, 1]; r2 = seg[k, 2]; z2 = seg[k, 3]
        if abs(z2 - z1) < 1e-12:
            if abs(dz) < 1e-15:
                continue
            t = (z1 - oz) / dz
            if t > tmin and t < best:
                x = ox + t * dx; y = oy + t * dy; rr = np.sqrt(x * x + y * y)
                if rr >= min(r1, r2) and rr <= max(r1, r2):
                    best = t; bnx = 0.0; bny = 0.0; bnz = 1.0
            continue
        m = (r2 - r1) / (z2 - z1); a = r1 + m * (oz - z1); b = m * dz
        A = dxy2 - b * b; B = 2.0 * (ox * dx + oy * dy - a * b); C = ox * ox + oy * oy - a * a
        roots = np.empty(2); nr = 0
        if abs(A) < 1e-14:
            if abs(B) > 1e-15:
                roots[0] = -C / B; nr = 1
        else:
            disc = B * B - 4.0 * A * C
            if disc >= 0.0:
                sq = np.sqrt(disc); roots[0] = (-B - sq) / (2.0 * A); roots[1] = (-B + sq) / (2.0 * A); nr = 2
        for q in range(nr):
            t = roots[q]
            if t <= tmin or t >= best:
                continue
            z = oz + t * dz; s = (z - z1) / (z2 - z1)
            if s < 0.0 or s > 1.0:
                continue
            r = a + b * t
            if r < 0.0:
                continue
            x = ox + t * dx; y = oy + t * dy
            best = t; bnx = x; bny = y; bnz = -r * m
    nn = np.sqrt(bnx * bnx + bny * bny + bnz * bnz)
    if best < 1e29 and nn > 0:
        return best, bnx / nn, bny / nn, bnz / nn
    return best, 0.0, 0.0, 1.0


@njit(cache=True)
def _trace(th, phi, sx, sy, n, seg, H, offs, max_ev):
    dx = np.sin(th) * np.cos(phi); dy = np.sin(th) * np.sin(phi); dz = np.cos(th)
    px = sx; py = sy; pz = -1e-7
    glass = True
    for ev in range(max_ev):
        cx0, cy0, _ = _nearest(px, py)
        best = 1e30; nx = 0.0; ny = 0.0; nz = 1.0
        for k in range(offs.shape[0]):
            cx = cx0 + offs[k, 0]; cy = cy0 + offs[k, 1]
            t, a, b, c = _lens_hit(px - cx, py - cy, pz, dx, dy, dz, seg, H, 1e-9)
            if t < best:
                best = t; nx = a; ny = b; nz = c
        # plane z = 0 (glass/air only outside the circles)
        if abs(dz) > 1e-15:
            t0 = -pz / dz
            if t0 > 1e-9 and t0 < best:
                qx = px + t0 * dx; qy = py + t0 * dy; _, _, dist = _nearest(qx, qy)
                if dist > 1.0:
                    best = t0; nx = 0.0; ny = 0.0; nz = 1.0
                elif glass and dz < 0.0:
                    return 1, np.arccos(min(1.0, -dz)), np.arctan2(dy, dx)        # back into the substrate
        if best > 1e29:
            if (not glass) and dz > 0.0:
                return 0, np.arccos(min(1.0, dz)), np.arctan2(dy, dx)
            if glass and dz < 0.0:
                return 1, np.arccos(min(1.0, -dz)), np.arctan2(dy, dx)
            return 2, 0.0, 0.0
        px += best * dx; py += best * dy; pz += best * dz
        n1 = n if glass else 1.0; n2 = 1.0 if glass else n
        ci = -(dx * nx + dy * ny + dz * nz)
        if ci < 0.0:
            nx = -nx; ny = -ny; nz = -nz; ci = -ci
        eta = n1 / n2; k2 = 1.0 - eta * eta * (1.0 - ci * ci)
        if k2 < 0.0:
            Rf = 1.0; ct = 0.0
        else:
            ct = np.sqrt(k2)
            rs = (n1 * ci - n2 * ct) / (n1 * ci + n2 * ct); rp = (n1 * ct - n2 * ci) / (n1 * ct + n2 * ci)
            Rf = 0.5 * (rs * rs + rp * rp)
        if np.random.random() < Rf:
            dx += 2 * ci * nx; dy += 2 * ci * ny; dz += 2 * ci * nz
        else:
            dx = eta * dx + (eta * ci - ct) * nx; dy = eta * dy + (eta * ci - ct) * ny; dz = eta * dz + (eta * ci - ct) * nz
            glass = not glass
        mm = np.sqrt(dx * dx + dy * dy + dz * dz); dx /= mm; dy /= mm; dz /= mm
        px += 1e-9 * dx; py += 1e-9 * dy; pz += 1e-9 * dz
    return 2, 0.0, 0.0


@njit(parallel=True)
def _batch(th, N, n, seg, H, offs, max_ev, seed):
    oc = np.empty(N, np.int64); to = np.empty(N); dp = np.empty(N)
    for i in prange(N):
        np.random.seed(seed + i)
        u = np.random.random(); w = np.random.random()
        sx = u * 2.0 + w * 1.0; sy = w * SQ3
        phi = 2 * np.pi * np.random.random()
        a, b, po = _trace(th, phi, sx, sy, n, seg, H, offs, max_ev)
        oc[i] = a; to[i] = b; dp[i] = (po - phi) % (2 * np.pi)
    return oc, to, dp


def offsets(rings=2):
    o = []
    for i in range(-rings, rings + 1):
        for j in range(-rings, rings + 1):
            x = 2.0 * i + 1.0 * j; y = SQ3 * j
            if np.hypot(x, y) <= 2.0 * rings + 1e-9: o.append((x, y))
    return np.array(o)


def segments(r, z):
    return np.ascontiguousarray(np.c_[r[:-1], z[:-1], r[1:], z[1:]], dtype=np.float64)


def kernel(r, z, n=1.51, N=20000, max_ev=300, NA=36, seed=1, rings=2):
    seg = segments(np.asarray(r, float), np.asarray(z, float)); H = float(np.max(z)); offs = offsets(rings)
    BT = np.zeros((90, 90)); K = np.zeros((90, 90, NA)); KT = np.zeros((90, 90, NA)); lost = np.zeros(90)
    for i in range(90):
        oc, to, dp = _batch(np.radians(i + 0.5), N, n, seg, H, offs, max_ev, seed * 7919 + i * 100003)
        d = np.degrees(to); lost[i] = np.mean(oc == 2)
        BT[i] = np.histogram(d[oc == 0], np.linspace(0, 90, 91))[0] / N
        for o, A in ((1, K), (0, KT)):
            m = oc == o; rb = np.minimum(d[m].astype(int), 89); ab = np.minimum((dp[m] / (2 * np.pi) * NA).astype(int), NA - 1)
            np.add.at(A[i], (rb, ab), 1.0 / N)
    return BT, K, lost, KT


# ---------------------------------------------------------------------------------------------------------------
# Lattice-resolved kernels. The hexagonal lattice (nearest neighbour along 0 deg) has mirror lines every 30 deg, so an
# incident azimuth phi is reduced to a in [0, 30] deg with a sign s (+1 rotation, -1 mirror); the outgoing azimuth is
# phi_out = phi + s * psi. KTa[in, a_bin, air, psi_bin] and Ka[in, a_bin, sub, psi_bin] are fractions of the rays
# incident at (theta_in, a_bin).
def reduce_az(phi):
    r = np.degrees(phi) % 60.0; s = np.where(r <= 30.0, 1.0, -1.0); a = np.where(r <= 30.0, r, 60.0 - r)
    return a, s


@njit(parallel=True)
def _batch_abs(th, N, n, seg, H, offs, max_ev, seed):
    oc = np.empty(N, np.int64); to = np.empty(N); pi_ = np.empty(N); po_ = np.empty(N)
    for i in prange(N):
        np.random.seed(seed + i)
        u = np.random.random(); w = np.random.random()
        sx = u * 2.0 + w * 1.0; sy = w * SQ3
        phi = 2 * np.pi * np.random.random()
        a, b, po = _trace(th, phi, sx, sy, n, seg, H, offs, max_ev)
        oc[i] = a; to[i] = b; pi_[i] = phi; po_[i] = po
    return oc, to, pi_, po_


def kernel_abs(r, z, n=1.51, N=60000, max_ev=300, NA=72, NB=6, seed=1, rings=2):
    seg = segments(np.asarray(r, float), np.asarray(z, float)); H = float(np.max(z)); offs = offsets(rings)
    KTa = np.zeros((90, NB, 90, NA), np.float32); Ka = np.zeros((90, NB, 90, NA), np.float32); cnt = np.zeros((90, NB))
    for i in range(90):
        oc, to, pin, pout = _batch_abs(np.radians(i + 0.5), N, n, seg, H, offs, max_ev, seed * 7919 + i * 100003)
        a, s = reduce_az(pin); ab = np.minimum((a / (30.0 / NB)).astype(int), NB - 1)
        cnt[i] = np.bincount(ab, minlength=NB)
        psi = (s * (pout - pin)) % (2 * np.pi); pb = np.minimum((psi / (2 * np.pi) * NA).astype(int), NA - 1)
        d = np.minimum(np.degrees(to).astype(int), 89)
        for o, A in ((0, KTa), (1, Ka)):
            m = oc == o; np.add.at(A[i], (ab[m], d[m], pb[m]), 1.0)
        KTa[i] /= cnt[i][:, None, None]; Ka[i] /= cnt[i][:, None, None]
    return KTa, Ka
