"""Single-pass escape probability of a microlens array, by Monte-Carlo ray tracing.

Hexagonally close-packed spherical caps of base radius r and height h on a
planar substrate, index-matched to it (n_MLA = n_sub), infinite in extent.
Aspect ratio AR = h/r; AR = 1 is a hemisphere.  The lens material and the
substrate are the same medium, so there is no interface at the base plane; the
9.3 % of the base plane that the close-packed circles do not cover is flat
glass/air.

p is defined as in the recycling law: the fraction of angularly randomised
light arriving at the extraction surface from inside that leaves into air in
one pass.  Light that comes back down through z = 0 is returned to the
substrate and is handled by the next term of the series, not here.

Rays carry no weight: reflection and refraction are sampled against the
unpolarised Fresnel transmittance, which keeps one ray per sample.
"""
import numpy as np

SQ3 = np.sqrt(3.0)


def _wrap(xy, a):
    """Fold a point back into one rhombic cell; the array is periodic, so this
    is exact and keeps the ray inside the 3x3 neighbourhood that is searched."""
    f1 = xy[:, 0] / a - xy[:, 1] / (a * SQ3)
    f2 = 2 * xy[:, 1] / (a * SQ3)
    f1 -= np.floor(f1)
    f2 -= np.floor(f2)
    out = np.empty_like(xy)
    out[:, 0] = a * (f1 + f2 / 2)
    out[:, 1] = a * SQ3 * f2 / 2
    return out


def _lattice_offsets(a):
    """Centres of the 3x3 neighbourhood of a hexagonal lattice of pitch a."""
    a1 = np.array([a, 0.0])
    a2 = np.array([a / 2, a * SQ3 / 2])
    rng = (-2, -1, 0, 1, 2)
    return np.array([i * a1 + j * a2 for i in rng for j in rng])


def _fresnel_T(ci, n1, n2):
    """Unpolarised transmittance for incidence cosine ci from n1 into n2."""
    s2 = (n1 / n2) ** 2 * (1.0 - ci ** 2)
    tir = s2 >= 1.0
    ct = np.sqrt(np.clip(1.0 - s2, 0.0, None))
    rs = (n1 * ci - n2 * ct) / (n1 * ci + n2 * ct + 1e-300)
    rp = (n1 * ct - n2 * ci) / (n1 * ct + n2 * ci + 1e-300)
    T = 1.0 - 0.5 * (rs ** 2 + rp ** 2)
    return np.where(tir, 0.0, np.clip(T, 0.0, 1.0)), ct


def escape_vs_angle(AR, n, th_deg, N=120000, max_events=40, seed=0):
    """B_T(theta): escape probability for light arriving at a fixed polar angle."""
    rng = np.random.default_rng(seed)
    return np.array([_trace_chunk(AR, n, N, max_events, rng, float(t))[0].mean()
                     for t in np.atleast_1d(th_deg)])


def escape_probability(AR, n, N=400000, max_events=40, seed=0, chunk=100000,
                       weight=None, th_deg=None):
    """Single-pass escape probability.

    weight = None   angularly randomised light (isotropic radiance), sampled directly
    weight = array  B_T(theta) is computed on th_deg and averaged with this weight,
                    which is what <P_sub . B_T> of the recycling law means."""
    if weight is not None:
        BT = escape_vs_angle(AR, n, th_deg, N=N, max_events=max_events, seed=seed)
        w = np.asarray(weight, float)
        return float(np.sum(w * BT) / np.sum(w))
    rng = np.random.default_rng(seed)
    esc = 0
    done = 0
    while done < N:
        m = min(chunk, N - done)
        esc += int(_trace_chunk(AR, n, m, max_events, rng)[0].sum())
        done += m
    return esc / N


def _trace_chunk(AR, n, N, max_events, rng, th_fixed=None):
    r = 1.0
    h = AR * r
    R = (r * r + h * h) / (2 * h)
    zc = h - R                      # sphere centre height
    a = 2 * r                       # pitch
    off = _lattice_offsets(a)       # (9,2)

    # start on the base plane, uniform over one rhombic cell
    u1, u2 = rng.random(N), rng.random(N)
    P = np.empty((N, 3))
    P[:, 0] = u1 * a + u2 * a / 2
    P[:, 1] = u2 * a * SQ3 / 2
    P[:, 2] = 0.0
    # direction: flux-weighted hemisphere (isotropic radiance), or a fixed polar angle
    if th_fixed is None:
        cz = np.sqrt(rng.random(N))
    else:
        cz = np.full(N, np.cos(np.radians(th_fixed)))
    sz = np.sqrt(np.clip(1 - cz ** 2, 0, None))
    ph = rng.random(N) * 2 * np.pi
    D = np.stack([sz * np.cos(ph), sz * np.sin(ph), cz], 1)

    inside = np.ones(N, bool)       # True: in the glass/lens; False: in air
    alive = np.ones(N, bool)
    escaped = np.zeros(N, bool)
    ret_cos = np.full(N, np.nan)    # |cos| of the direction on returning to the substrate

    for _ in range(max_events):
        if not alive.any():
            break
        idx = np.flatnonzero(alive)
        p, d, ins = P[idx], D[idx], inside[idx]

        # --- nearest cap-sphere intersection over the 3x3 neighbourhood
        best_t = np.full(len(idx), np.inf)
        best_n = np.zeros((len(idx), 3))
        for o in off:
            c = np.array([o[0], o[1], zc])
            L = p - c
            b = np.sum(L * d, 1)
            cq = np.sum(L * L, 1) - R * R
            disc = b * b - cq
            ok = disc > 0
            sq = np.sqrt(np.where(ok, disc, 0.0))
            for t in (-b - sq, -b + sq):
                hit = p + t[:, None] * d
                good = ok & (t > 1e-9) & (hit[:, 2] >= 0.0) & (t < best_t)
                if good.any():
                    best_t = np.where(good, t, best_t)
                    nrm = (hit - c) / R
                    best_n = np.where(good[:, None], nrm, best_n)

        # --- the base plane z = 0.  The close-packed circles leave 1 - pi/(2sqrt3)
        # --- = 9.3 % of it uncovered, and that part is flat glass/air.
        dz = d[:, 2]
        with np.errstate(divide='ignore', invalid='ignore'):
            t0 = -p[:, 2] / dz
        # downward crossing: a glass ray has returned, an air ray enters the glass
        t_down = np.where((dz < 0) & (t0 > 1e-9), t0, np.inf)
        # upward crossing while in the glass: an interface only where there is no cap
        up = ins & (dz > 0) & np.isfinite(t0)
        tu = np.where(up, np.maximum(t0, 0.0), np.inf)
        cross = p[:, :2] + np.where(np.isfinite(tu), tu, 0.0)[:, None] * d[:, :2]
        rho = np.full(len(idx), np.inf)
        for o in off:
            rho = np.minimum(rho, np.hypot(cross[:, 0] - o[0], cross[:, 1] - o[1]))
        t_up = np.where(up & (rho >= r), tu, np.inf)

        t_pl = np.minimum(t_down, t_up)
        take_sphere = best_t < t_pl
        t = np.where(take_sphere, best_t, t_pl)
        no_hit = ~np.isfinite(t)

        # a ray in air with nothing ahead of it has left; one in glass cannot happen
        if no_hit.any():
            g = idx[no_hit]
            esc_g = ~inside[g]
            escaped[g] = esc_g
            ret_cos[g[~esc_g]] = np.abs(d[no_hit][~esc_g, 2])   # defensive: back to the substrate
            alive[g] = False

        go = np.flatnonzero(~no_hit)
        if go.size == 0:
            continue
        gi = idx[go]
        hit = p[go] + t[go][:, None] * d[go]
        nrm = np.where(take_sphere[go][:, None], best_n[go],
                       np.tile([0.0, 0.0, 1.0], (go.size, 1)))
        ins_g = ins[go]
        dg = d[go]

        # a glass ray crossing the flat plane downwards has returned to the substrate;
        # one meeting it upwards is at the flat glass/air interface and carries on
        back = ins_g & ~take_sphere[go] & (dg[:, 2] < 0)
        if back.any():
            alive[gi[back]] = False
            ret_cos[gi[back]] = np.abs(dg[back][:, 2])

        act = np.flatnonzero(~back)
        if act.size == 0:
            continue
        ai = gi[act]
        nn = nrm[act]
        dd = dg[act]
        ig = ins_g[act]
        # outward normal relative to the direction of travel
        ci = -np.sum(dd * nn, 1)
        flip = ci < 0
        nn = np.where(flip[:, None], -nn, nn)
        ci = np.abs(ci)
        n1 = np.where(ig, n, 1.0)
        n2 = np.where(ig, 1.0, n)
        T, ct = _fresnel_T(ci, n1, n2)
        refract = rng.random(act.size) < T

        newd = np.empty((act.size, 3))
        # reflection
        newd[~refract] = (dd + 2 * ci[:, None] * nn)[~refract]
        # refraction
        eta = (n1 / n2)[:, None]
        tdir = eta * dd + (eta[:, 0] * ci - ct)[:, None] * nn
        nrm_len = np.linalg.norm(tdir, axis=1, keepdims=True)
        newd[refract] = (tdir / np.where(nrm_len > 0, nrm_len, 1.0))[refract]

        newp = hit[act] + 1e-9 * newd
        newp[:, :2] = _wrap(newp[:, :2], a)     # periodic array: fold back
        P[ai] = newp
        D[ai] = newd
        inside[ai] = np.where(refract, ~ig, ig)
        # a ray that has just entered the glass through the flat gap, heading down,
        # is back in the substrate
        home = inside[ai] & (newp[:, 2] <= 1e-8) & (newd[:, 2] < 0)
        if home.any():
            ret_cos[ai[home]] = np.abs(newd[home, 2])
            alive[ai[home]] = False

    return escaped, ret_cos


def bsdf(AR, n, N=40000, n_bin=90, max_events=40, seed=0):
    """Angle-resolved response of the array, as the recycling series needs it.

    B_T[i]    fraction of light arriving at polar angle bin i that escapes to air
    B_R[r, i] fraction that goes back into the substrate travelling at bin r

    Bins are 1 degree wide, bin i covering [i, i+1) degrees, so that they line up
    with the 0-89 degree grid the stack reflectance is computed on.  The two
    together sum to 1 for every column: the array neither absorbs nor traps.
    """
    rng = np.random.default_rng(seed)
    edges = np.linspace(0.0, 90.0, n_bin + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    BT = np.zeros(n_bin)
    BR = np.zeros((n_bin, n_bin))
    for i, th in enumerate(centres):
        esc, rc = _trace_chunk(AR, n, N, max_events, rng, float(th))
        BT[i] = esc.mean()
        good = np.isfinite(rc)
        if good.any():
            ang = np.degrees(np.arccos(np.clip(rc[good], 0.0, 1.0)))
            h, _ = np.histogram(ang, bins=edges)
            BR[:, i] = h / N
    return BT, BR
