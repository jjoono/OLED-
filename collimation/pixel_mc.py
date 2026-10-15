"""3-D Monte-Carlo of one micro-OLED pixel with the thin-film collimator:

    z = 0        device (specular reflector R(theta, lam); emits P(theta, lam))
    0 < z < s    spacer (thin-film encapsulation), clear
    s < z < t    scattering slab: mean free path ls, Henyey-Greenstein g, albedo w
    z = t        DBR: transmits T(theta, lam) into air, otherwise specular reflection
    x, y         square pixel of pitch P; side walls handled by unfolding

Unfolding: the path is propagated in an unbounded plane and every crossing of a pixel
boundary is counted as one wall hit.  A mirror wall of reflectance Rw then contributes the
weight Rw**hits, and with no walls the final unfolded cell gives the crosstalk.  Both come
from the same photon histories.  All media have n = N_MED (index-matched scatterer host).
"""
import numpy as np
from numba import njit, prange
import model as Mo

NB = Mo.NB


@njit(cache=True)
def _hg(g, u):
    if abs(g) < 1e-6:
        return 2 * u - 1
    s = (1 - g * g) / (1 - g + 2 * g * u)
    return (1 + g * g - s * s) / (2 * g)


@njit(cache=True)
def _rotate(ux, uy, uz, ct, phi):
    st = np.sqrt(max(0.0, 1 - ct * ct)); cp = np.cos(phi); sp = np.sin(phi)
    if abs(uz) > 0.99999:
        return st * cp, st * sp, ct * np.sign(uz)
    d = np.sqrt(1 - uz * uz)
    nx = st * (ux * uz * cp - uy * sp) / d + ux * ct
    ny = st * (uy * uz * cp + ux * sp) / d + uy * ct
    nz = -st * cp * d + uz * ct
    return nx, ny, nz


@njit(cache=True, parallel=True)
def run(n, seed, lamcdf, Pcdf, Rdev, Tdbr, P, s, t, ls, g, alb, max_ev, ngap1=1.0, ngap2=1.25, Rmet=0.95):
    """Returns per photon: status (0 escaped, 1 absorbed, 2 lost), air angle index (0.25 deg),
    wall hits, final cell (ix, iy), weight (scattering absorption)."""
    st = np.empty(n, np.int64); ang = np.empty(n, np.int64); hits = np.empty(n, np.int64)
    cx = np.empty(n, np.int64); cy = np.empty(n, np.int64); wt = np.empty(n)
    lw1 = np.empty(n); lw2 = np.empty(n)
    cth1 = np.sqrt(max(0.0, 1 - (ngap1 / Mo.N_MED) ** 2)); cth2 = np.sqrt(max(0.0, 1 - (ngap2 / Mo.N_MED) ** 2))
    lR = np.log(Rmet)
    nmed = Mo.N_MED
    for i in prange(n):
        np.random.seed(seed + i)
        k = np.searchsorted(lamcdf, np.random.random())
        b = np.searchsorted(Pcdf[k], np.random.random())
        th = (b + np.random.random()) * (np.pi / 2) / NB
        ph = 2 * np.pi * np.random.random()
        ux = np.sin(th) * np.cos(ph); uy = np.sin(th) * np.sin(ph); uz = np.cos(th)
        x = (np.random.random() - 0.5) * P; y = (np.random.random() - 0.5) * P; z = 0.0
        w = 1.0; status = 2; a_i = 0; nh = 0; l1 = 0.0; l2 = 0.0
        for ev in range(max_ev):
            # distance to the next horizontal boundary / scattering event
            if uz > 0:
                if z < s:
                    d = (s - z) / uz; zn = s; kind = 0
                else:
                    d = (t - z) / uz; zn = t; kind = 1
            else:
                if z > s:
                    d = (s - z) / uz; zn = s; kind = 0
                else:
                    d = (0.0 - z) / uz; zn = 0.0; kind = 2
            in_slab = (z >= s and uz > 0) or (z > s and uz < 0)
            if in_slab and ls > 0:
                ds = -ls * np.log(np.random.random() + 1e-300)
                if ds < d:
                    x0 = x; y0 = y
                    x += ds * ux; y += ds * uy; z += ds * uz
                    hx = abs(int(np.floor(x / P + 0.5)) - int(np.floor(x0 / P + 0.5)))
                    hy = abs(int(np.floor(y / P + 0.5)) - int(np.floor(y0 / P + 0.5)))
                    nh += hx + hy
                    if hx > 0:
                        if abs(ux) > cth1: l1 += hx * lR
                        if abs(ux) > cth2: l2 += hx * lR
                    if hy > 0:
                        if abs(uy) > cth1: l1 += hy * lR
                        if abs(uy) > cth2: l2 += hy * lR
                    w *= alb
                    ct = _hg(g, np.random.random())
                    ux, uy, uz = _rotate(ux, uy, uz, ct, 2 * np.pi * np.random.random())
                    if w < 1e-4:
                        status = 1; break
                    continue
            x0 = x; y0 = y
            x += d * ux; y += d * uy; z = zn
            hx = abs(int(np.floor(x / P + 0.5)) - int(np.floor(x0 / P + 0.5)))
            hy = abs(int(np.floor(y / P + 0.5)) - int(np.floor(y0 / P + 0.5)))
            nh += hx + hy
            if hx > 0:
                if abs(ux) > cth1: l1 += hx * lR
                if abs(ux) > cth2: l2 += hx * lR
            if hy > 0:
                if abs(uy) > cth1: l1 += hy * lR
                if abs(uy) > cth2: l2 += hy * lR
            if kind == 0:
                z = zn + 1e-9 * uz     # just cross the spacer/slab boundary
                continue
            bin_ = int(np.arccos(min(1.0, abs(uz))) / (np.pi / 2) * NB)
            if bin_ >= NB:
                bin_ = NB - 1
            if kind == 1:                       # DBR
                if np.random.random() < Tdbr[k, bin_]:
                    sa = nmed * np.sqrt(max(0.0, 1 - uz * uz))
                    a_i = int(np.degrees(np.arcsin(min(1.0, sa))) / 0.25)
                    status = 0; break
                uz = -uz; z = t - 1e-9
            else:                               # device
                if np.random.random() >= Rdev[k, bin_]:
                    status = 1; break
                uz = -uz; z = 1e-9
        st[i] = status; ang[i] = min(a_i, 359); wt[i] = w
        ix = int(np.floor(x / P + 0.5)); iy = int(np.floor(y / P + 0.5))
        cx[i] = ix; cy[i] = iy
        hits[i] = nh; lw1[i] = l1; lw2[i] = l2
    return st, ang, hits, cx, cy, wt, lw1, lw2
