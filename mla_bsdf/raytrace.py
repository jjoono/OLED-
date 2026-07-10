"""
Monte-Carlo geometric ray tracer for the BSDF of a hemispherical microlens array
(MLA) sitting on a glass substrate of the SAME refractive index.

Physics / model (see README for the full derivation of the simplifications):

  * Geometry  : full hemispheres, radius R, on a hexagonal close-packed lattice
                with centre-to-centre pitch = 2R (lenses touch).  Lengths are
                normalised to R = 1, pitch = 2.
  * Media     : lens = substrate = glass (n_glass).  Everything above the
                surface (spherical caps + flat inter-lens gaps) is air (n_air).
                Because the lens and the substrate share the same index, the
                ONLY optical interfaces are
                    (a) the spherical caps            glass <-> air
                    (b) the flat z = 0 gap regions    glass <-> air   (~9.3 %)
  * Source    : a plane wave inside the glass travelling upward at polar angle
                theta_i (measured from the +z substrate normal) and a random
                azimuth phi (=> azimuthal averaging).  Ray start points are
                sampled uniformly over one rhombic unit cell of the lattice,
                which reproduces uniform illumination of the periodic array and
                the correct 90.7 % / 9.3 % lens / gap area split.
  * Below z<0 : semi-infinite glass  -> a ray ending here (going down) is
                counted as REFLECTED (angle theta_r from the -z axis).
  * Above     : semi-infinite air    -> a ray leaving upward is TRANSMITTED
                (angle theta_t from the +z axis).
  * Interfaces: unpolarised Fresnel.  At each interface a single uniform deviate
                chooses reflection (prob = R_fresnel, incl. TIR where R=1) or
                refraction, so ray count stays constant and multiple bounces /
                neighbouring-lens re-entry are handled naturally.

No absorption anywhere => for every incident angle  T + R = 1  (up to the tiny
fraction of rays that exceed MAX_BOUNCE, reported separately as "lost").
"""

import numpy as np
from numba import njit, prange

SQRT3 = np.sqrt(3.0)

# outcome codes
TRANSMIT = 0
REFLECT = 1
LOST = 2


@njit(cache=True, inline="always")
def _nearest_center(x, y, pitch):
    """Return (cx, cy) of the hex-lattice centre nearest to (x, y).

    Basis: v1 = (pitch, 0), v2 = (pitch/2, pitch*sqrt3/2).
    """
    a = pitch
    j = 2.0 * y / (a * SQRT3)
    i = x / a - j / 2.0
    i0 = round(i)
    j0 = round(j)
    best = 1e18
    bcx = 0.0
    bcy = 0.0
    for di in range(-1, 2):
        for dj in range(-1, 2):
            ii = i0 + di
            jj = j0 + dj
            cx = ii * a + jj * a / 2.0
            cy = jj * a * SQRT3 / 2.0
            d2 = (x - cx) ** 2 + (y - cy) ** 2
            if d2 < best:
                best = d2
                bcx = cx
                bcy = cy
    return bcx, bcy


@njit(cache=True, inline="always")
def _in_footprint(x, y, R, pitch):
    cx, cy = _nearest_center(x, y, pitch)
    return (x - cx) ** 2 + (y - cy) ** 2 < R * R


@njit(cache=True)
def trace_ray(theta_i, phi, sx, sy, n_glass, n_air, R, pitch, max_bounce, concave):
    """Trace one ray. Returns (outcome_code, theta_out[rad]).

    concave = 0 : convex lenses (glass hemispheres bulge UP into air; the real
                  spherical surface is the UPPER hemisphere z >= 0).
    concave = 1 : concave pits (hemispherical AIR cavities carved DOWN into the
                  glass; the real spherical surface is the LOWER hemisphere
                  z <= 0, inside the sphere is the air cavity).
    In both cases the sphere centres sit on the z = 0 plane and the flat
    inter-lens gaps are identical glass/air interfaces.
    """
    st = np.sin(theta_i)
    ct = np.cos(theta_i)
    dx = st * np.cos(phi)
    dy = st * np.sin(phi)
    dz = ct                      # upward, dz > 0
    px = sx
    py = sy
    # start in the semi-infinite glass, below all geometry
    pz = -(1.0 + 1.0e-3) * R if concave else -1.0e-3 * R
    in_glass = True

    a = pitch
    for _ in range(max_bounce):
        best_t = 1.0e18
        hit = -1               # -1 none, 0 flat, 1 sphere
        nx = 0.0
        ny = 0.0
        nz = 0.0

        # (a) flat gap plane z = 0 (only a real interface outside lens footprints)
        if dz != 0.0:
            t = (0.0 - pz) / dz
            if 1.0e-7 < t < best_t:
                hx = px + t * dx
                hy = py + t * dy
                if not _in_footprint(hx, hy, R, pitch):
                    best_t = t
                    hit = 0
                    nx = 0.0
                    ny = 0.0
                    nz = 1.0

        # (b) spherical caps of lenses near the current (px, py)
        j = 2.0 * py / (a * SQRT3)
        i = px / a - j / 2.0
        i0 = int(round(i))
        j0 = int(round(j))
        for di in range(-4, 5):
            for dj in range(-4, 5):
                cx = (i0 + di) * a + (j0 + dj) * a / 2.0
                cy = (j0 + dj) * a * SQRT3 / 2.0
                ox = px - cx
                oy = py - cy
                oz = pz               # centres at z = 0
                bq = ox * dx + oy * dy + oz * dz
                cq = ox * ox + oy * oy + oz * oz - R * R
                disc = bq * bq - cq
                if disc <= 0.0:
                    continue
                sq = np.sqrt(disc)
                for t in (-bq - sq, -bq + sq):
                    if 1.0e-7 < t < best_t:
                        hz = pz + t * dz
                        # concave: lower hemisphere (z<=0); convex: upper (z>=0)
                        real = (hz <= 0.0) if concave else (hz >= 0.0)
                        if real:
                            best_t = t
                            hit = 1
                            hx = px + t * dx
                            hy = py + t * dy
                            nx = (hx - cx) / R
                            ny = (hy - cy) / R
                            nz = hz / R

        if hit == -1:
            # no interface ahead -> ray escapes to a semi-infinite half space
            if dz > 0.0:
                return TRANSMIT, np.arccos(min(1.0, dz))
            else:
                return REFLECT, np.arccos(min(1.0, -dz))

        # advance to the interface
        px += best_t * dx
        py += best_t * dy
        pz += best_t * dz

        # Fresnel (unpolarised) with the normal oriented against the incidence
        n1 = n_glass if in_glass else n_air
        n2 = n_air if in_glass else n_glass
        cosi = -(dx * nx + dy * ny + dz * nz)
        onx, ony, onz = nx, ny, nz
        if cosi < 0.0:
            onx, ony, onz = -nx, -ny, -nz
            cosi = -cosi
        eta = n1 / n2
        k = 1.0 - eta * eta * (1.0 - cosi * cosi)
        if k < 0.0:
            Rf = 1.0                      # total internal reflection
            cost = 0.0
        else:
            cost = np.sqrt(k)
            rs = (n1 * cosi - n2 * cost) / (n1 * cosi + n2 * cost)
            rp = (n1 * cost - n2 * cosi) / (n1 * cost + n2 * cosi)
            Rf = 0.5 * (rs * rs + rp * rp)

        if np.random.random() < Rf:
            # reflect: d + 2*cosi*on  (on points against d)
            dx = dx + 2.0 * cosi * onx
            dy = dy + 2.0 * cosi * ony
            dz = dz + 2.0 * cosi * onz
        else:
            # refract
            dx = eta * dx + (eta * cosi - cost) * onx
            dy = eta * dy + (eta * cosi - cost) * ony
            dz = eta * dz + (eta * cosi - cost) * onz
            nrm = np.sqrt(dx * dx + dy * dy + dz * dz)
            dx /= nrm
            dy /= nrm
            dz /= nrm
            in_glass = not in_glass

        # nudge off the surface to avoid immediate re-intersection
        px += 1.0e-7 * dx
        py += 1.0e-7 * dy
        pz += 1.0e-7 * dz

    return LOST, 0.0


@njit(cache=True, parallel=True)
def trace_batch(theta_i, N, n_glass, n_air, R, pitch, max_bounce, concave):
    """Trace N rays at incidence angle theta_i. Returns outcome[N], theta_out[N]."""
    outcome = np.empty(N, dtype=np.int64)
    tout = np.empty(N, dtype=np.float64)
    a = pitch
    # rhombic unit-cell vectors
    for n in prange(N):
        u = np.random.random()
        w = np.random.random()
        sx = u * a + w * (a / 2.0)
        sy = w * (a * SQRT3 / 2.0)
        phi = 2.0 * np.pi * np.random.random()
        oc, th = trace_ray(theta_i, phi, sx, sy, n_glass, n_air, R, pitch,
                           max_bounce, concave)
        outcome[n] = oc
        tout[n] = th
    return outcome, tout


def sweep(theta_i_deg, N=200_000, n_glass=1.5, n_air=1.0,
          R=1.0, pitch=2.0, max_bounce=200, nbins=90, seed=0, concave=False):
    """Run the full incident-angle sweep.

    Returns a dict with the 2-D BSDF maps (percent of incident power per 1 deg
    output bin) and the integrated reflectance / transmittance curves.
    """
    np.random.seed(seed)
    theta_i_deg = np.asarray(theta_i_deg, dtype=np.float64)
    ni = theta_i_deg.size
    bsdf_T = np.zeros((ni, nbins))
    bsdf_R = np.zeros((ni, nbins))
    Ttot = np.zeros(ni)
    Rtot = np.zeros(ni)
    lost = np.zeros(ni)
    edges = np.linspace(0.0, 90.0, nbins + 1)

    for idx, ti in enumerate(theta_i_deg):
        oc, th = trace_batch(np.deg2rad(ti), int(N), n_glass, n_air,
                             R, pitch, max_bounce, int(concave))
        thd = np.rad2deg(th)
        mT = oc == TRANSMIT
        mR = oc == REFLECT
        Ttot[idx] = mT.sum() / N
        Rtot[idx] = mR.sum() / N
        lost[idx] = (oc == LOST).sum() / N
        hT, _ = np.histogram(thd[mT], bins=edges)
        hR, _ = np.histogram(thd[mR], bins=edges)
        bsdf_T[idx] = hT / N * 100.0
        bsdf_R[idx] = hR / N * 100.0

    return {
        "theta_i": theta_i_deg,
        "theta_out_centers": 0.5 * (edges[:-1] + edges[1:]),
        "bsdf_T": bsdf_T,          # [%] power per 1 deg transmitted-angle bin
        "bsdf_R": bsdf_R,          # [%] power per 1 deg reflected-angle bin
        "T_total": Ttot,
        "R_total": Rtot,
        "lost": lost,
        "params": dict(N=N, n_glass=n_glass, n_air=n_air, R=R, pitch=pitch,
                       max_bounce=max_bounce, nbins=nbins, seed=seed,
                       concave=bool(concave)),
    }
