"""
recycling.py -- multi-pass (angular photon recycling) models for a top-emitting OLED
whose bottom electrode is transparent and sits on a scattering / redistribution layer
(SL) backed by a mirror, and whose top (thin Ag + DBR) acts as an angle-selective filter.

Coordinate convention: the SL occupies 0 < z < d_SL; z = 0 is the mirror side,
z = d_SL is the interface to the coherent OLED stack (TCO / organics / thin Ag / DBR / air).
Angles theta_m are measured in the SL matrix (index n_m).

Two models:
  * analytic(): full angular randomisation at every visit of the SL (generalises Eq. (2)
    of the round-trip-loss manuscript; angle-resolved top transmission, polarisation-resolved).
  * monte_carlo(): photon Monte Carlo in an SL with Henyey-Greenstein scattering, absorption,
    a mirror with angle/polarisation-dependent reflectance and the coherent top stack as a
    boundary with R/T/A(theta, lambda, pol). Tracks lateral displacement (blur / PSF).
"""
from __future__ import annotations

import numpy as np

from oledsim import Medium, stack_RTA, radiant_intensity, KM, HC_E


def angle_grid(n=181):
    """theta grid (deg) in the SL, midpoints, for Lambertian averages."""
    edges = np.linspace(0, 90, n)
    th = 0.5 * (edges[1:] + edges[:-1])
    # Lambertian weight cos*sin dtheta, normalised
    w = np.cos(np.deg2rad(th)) * np.sin(np.deg2rad(th)) * np.deg2rad(edges[1] - edges[0])
    return th, w / w.sum(), edges


def boundary_tables(top_layers, n_m, wl, th_deg, bot_layers=(), mirror_nk=None):
    """R/T/A of the top stack and R of the bottom mirror for light inside the SL.

    Returns dict with arrays (2 pol, nwl, nth): Rt, Tt, At (and At_layers (2,nlay,nwl,nth)), Rb.
    pol index 0 = s, 1 = p."""
    wl = np.asarray(wl, float)
    nm = np.asarray(n_m, complex) * np.ones_like(wl)
    beta = np.real(nm)[:, None] * np.sin(np.deg2rad(th_deg))[None, :]
    air = Medium(np.ones_like(wl, dtype=complex))
    sl = Medium(nm)
    out = {k: [] for k in ("Rt", "Tt", "At", "Al", "Rb")}
    for pol in ("s", "p"):
        R, T, A = stack_RTA(top_layers, sl, air, beta, wl, pol)
        out["Rt"].append(R)
        out["Tt"].append(T)
        out["At"].append(A.sum(0))
        out["Al"].append(A)
        if mirror_nk is not None:
            Rb, _, _ = stack_RTA(list(bot_layers), sl, Medium(np.asarray(mirror_nk, complex)), beta, wl, pol)
        else:
            Rb = np.ones_like(R)
        out["Rb"].append(Rb)
    return {k: np.array(v) for k, v in out.items()}


def shift_tables(top_layers, n_m, wl, th_deg, bot_layers=(), mirror_nk=None, dbeta=1e-6, clip_um=10.0):
    """Lateral (Goos-Haenchen) displacement D = -(lambda/2pi) d arg(r)/d beta (in um) of a beam reflected
    by the coherent top stack and by the mirror(+spacer), as seen from the SL. Shape (2, nwl, nth).
    Light reflected obliquely by a ~0.6-um organic/TCO stack walks sideways inside it; ignoring this
    underestimates the point-spread function (blur)."""
    from oledsim import stack_r_tau
    wl = np.asarray(wl, float)
    nm = np.asarray(n_m, complex) * np.ones_like(wl)
    b = np.real(nm)[:, None] * np.sin(np.deg2rad(th_deg))[None, :]
    air = Medium(np.ones_like(wl, dtype=complex))
    sl = Medium(nm)
    Dt, Db = [], []
    for pol in ("s", "p"):
        r1, *_ = stack_r_tau(list(top_layers), sl, air, b, wl, pol)
        r2, *_ = stack_r_tau(list(top_layers), sl, air, b + dbeta, wl, pol)
        Dt.append(-(wl[:, None] / (2 * np.pi)) * np.angle(r2 / r1) / dbeta / 1000.0)
        if mirror_nk is not None:
            m = Medium(np.asarray(mirror_nk, complex))
            q1, *_ = stack_r_tau(list(bot_layers), sl, m, b, wl, pol)
            q2, *_ = stack_r_tau(list(bot_layers), sl, m, b + dbeta, wl, pol)
            Db.append(-(wl[:, None] / (2 * np.pi)) * np.angle(q2 / q1) / dbeta / 1000.0)
        else:
            Db.append(np.zeros_like(Dt[-1]))
    return np.clip(np.array(Dt), -clip_um, clip_um), np.clip(np.array(Db), -clip_um, clip_um)


def emission_into_sl(res, n_m, th_deg, edges_deg):
    """Photon power (per unit radiated power) emitted into the SL per angular bin and
    polarisation: shape (2, nwl, nth). Integrates the dipole density P_b(u)."""
    ne_ = res["n_eml"][:, None]
    nmr = np.real(np.asarray(n_m) * np.ones_like(res["n_eml"]))[:, None]
    u = res["u"]
    out = []
    for pol in ("s", "p"):
        P = res[f"Pb_{pol}"]
        # cumulative integral of P(u) 2u du on the u grid
        cum = np.cumsum(P * 2 * u[None, :] * res["du"], axis=1)
        ue = (nmr / ne_) * np.sin(np.deg2rad(edges_deg))[None, :]
        uu = np.r_[0.0, u + 0.5 * res["du"]]  # cumulative values live on the bin edges of the midpoint grid
        C = np.array([np.interp(ue[i], uu, np.r_[0.0, cum[i]]) for i in range(P.shape[0])])
        out.append(np.diff(C, axis=1))
    return np.array(out)


def analytic(res, tables, n_m, wl, th_deg, w_lamb, Pin, q_weights, a_sl_pass=0.0, R_b_extra=1.0):
    """Full-randomisation multi-pass model.

    res       : DipoleResult of the emitter in the stack with the SL as bottom medium
    tables    : boundary_tables(...)
    Pin       : emission_into_sl(...) (2, nwl, nth)
    q_weights : photons emitted at each wavelength per exciton (nwl,)
    a_sl_pass : extra loss per down-and-up visit of the SL (absorption in the SL)
    Returns dict of per-exciton quantities (spectrally resolved and totals).
    """
    Rt, Tt, At, Rb = tables["Rt"], tables["Tt"], tables["At"], tables["Rb"]
    w = w_lamb[None, None, :] * 0.5  # Lambertian, unpolarised (0.5 per pol)
    p = np.sum(Tt * w, axis=(0, 2))
    r = np.sum(Rt * w, axis=(0, 2))
    a = np.sum(At * w, axis=(0, 2))
    rb_refl = np.sum(Rt * Rb * w, axis=(0, 2)) / np.maximum(r, 1e-30)  # mirror reflectance seen by top-reflected light
    keep = (1 - a_sl_pass) * R_b_extra
    phi1 = np.sum(Pin * Rb, axis=(0, 2)) * keep  # first upward Lambertian flux (per unit radiated power)
    loop = r * rb_refl * keep
    phi_tot = phi1 / (1 - loop)
    T0 = 0.5 * (Tt[0, :, 0] + Tt[1, :, 0])
    nmr = np.real(np.asarray(n_m) * np.ones_like(wl))
    L_rec0 = T0 * phi_tot / (np.pi * nmr ** 2)  # on-axis photon radiance per unit radiated power
    Is, Ip = radiant_intensity(res, "top", 1.0, np.array([0.0]))
    L_dir0 = (Is + Ip)[:, 0]
    # spectra per exciton
    out = dict(p=p, r=r, a_top=a, rb=rb_refl, phi1=phi1, phi_tot=phi_tot,
               A_round=1 - (1 - a / np.maximum(1 - p, 1e-12)) * rb_refl * keep,
               eqe_dir=q_weights * res["top"], eqe_rec=q_weights * phi_tot * p,
               P_dn=q_weights * np.sum(Pin, axis=(0, 2)),
               L_dir0=q_weights * L_dir0, L_rec0=q_weights * L_rec0)
    # recycled light angular distribution in air (radiance ratio to on-axis): T_top(theta_m(theta_air))
    out["T0"] = T0
    return out


def summarize(out, wl, V):
    s = {}
    for k in ("eqe_dir", "eqe_rec", "P_dn"):
        s[k] = float(np.sum(out[k]))
    s["eqe"] = s["eqe_dir"] + s["eqe_rec"]
    s["CE_dir"] = float(KM * np.sum(V * (HC_E / wl) * out["L_dir0"]))
    s["CE_rec"] = float(KM * np.sum(V * (HC_E / wl) * out["L_rec0"]))
    s["CE"] = s["CE_dir"] + s["CE_rec"]
    wts = out["P_dn"] / max(np.sum(out["P_dn"]), 1e-30)
    s["p_avg"] = float(np.sum(wts * out["p"]))
    s["a_top_avg"] = float(np.sum(wts * out["a_top"]))
    s["rb_avg"] = float(np.sum(wts * out["rb"]))
    s["eta_rec"] = s["eqe_rec"] / max(s["P_dn"], 1e-30)
    return s


# ----------------------------------------------------------------------------
# Monte Carlo
# ----------------------------------------------------------------------------
def _hg_cos(g, xi):
    if abs(g) < 1e-6:
        return 2 * xi - 1
    t = (1 - g * g) / (1 - g + 2 * g * xi)
    return np.clip((1 + g * g - t * t) / (2 * g), -1, 1)


def monte_carlo(tables, th_deg, wl, Pin, q_weights, n_m, d_sl_um, mus_um, g, mua_um=0.0,
                n_photons=200000, max_steps=4000, seed=1, record_xy=True, shifts=None):
    """Photon Monte Carlo of the recycled light inside the SL.

    Photons are launched downward from the top boundary (z = d) with the joint
    (lambda, theta, pol) distribution of Pin weighted by q_weights.
    Returns per-exciton escape statistics (normalised by the launched power)."""
    rng = np.random.default_rng(seed)
    nwl, nth = len(wl), len(th_deg)
    Wjoint = (Pin * q_weights[None, :, None])  # (2, nwl, nth)
    total_in = Wjoint.sum()
    cdf = np.cumsum(Wjoint.ravel())
    cdf /= cdf[-1]
    idx = np.searchsorted(cdf, rng.random(n_photons))
    pol, il, it = np.unravel_index(idx, Wjoint.shape)
    dth = th_deg[1] - th_deg[0]
    theta = np.deg2rad(th_deg[it] + (rng.random(n_photons) - 0.5) * dth)
    phi = 2 * np.pi * rng.random(n_photons)
    ux, uy, uz = np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), -np.cos(theta)
    d = d_sl_um
    z = np.full(n_photons, d)
    x = np.zeros(n_photons)
    y = np.zeros(n_photons)
    alive = np.ones(n_photons, bool)
    mut = mus_um + mua_um
    th_edges = np.deg2rad(np.linspace(0, 90, nth + 1))
    Rt, Tt, Rb = tables["Rt"], tables["Tt"], tables["Rb"]
    esc = dict(il=[], th_air=[], x=[], y=[], pol=[], n_top_hits=[])
    fate = dict(top_abs=0, bot_abs=0, sl_abs=0, lost=0)
    hits = np.zeros(n_photons, int)
    nm = float(np.real(n_m))
    # radiance estimator: count top-boundary arrivals per (lambda, theta-bin, pol)
    arrivals = np.zeros((2, nwl, nth))
    for step in range(max_steps):
        ids = np.nonzero(alive)[0]
        if ids.size == 0:
            break
        uzi = uz[ids]
        with np.errstate(divide="ignore"):
            dist_b = np.where(uzi > 0, (d - z[ids]) / uzi, z[ids] / np.maximum(-uzi, 1e-300))
        if mut > 0:
            s = -np.log(rng.random(ids.size)) / mut
        else:
            s = np.full(ids.size, np.inf)
        hitb = s >= dist_b
        # --- boundary hits
        hb = ids[hitb]
        if hb.size:
            db = dist_b[hitb]
            x[hb] += ux[hb] * db
            y[hb] += uy[hb] * db
            up = uz[hb] > 0
            # top boundary
            ht = hb[up]
            if ht.size:
                z[ht] = d
                ct = np.clip(uz[ht], 0, 1)
                tbin = np.minimum(np.searchsorted(th_edges, np.arccos(ct)) - 1, nth - 1)
                np.add.at(arrivals, (pol[ht], il[ht], tbin), 1.0)
                hits[ht] += 1
                Rv = Rt[pol[ht], il[ht], tbin]
                Tv = Tt[pol[ht], il[ht], tbin]
                xi = rng.random(ht.size)
                refl = xi < Rv
                trans = (~refl) & (xi < Rv + Tv)
                if shifts is not None:  # lateral walk inside the coherent top stack
                    hr = ht[refl]
                    sh = shifts[0][pol[hr], il[hr], tbin[refl]]
                    upar = np.sqrt(ux[hr] ** 2 + uy[hr] ** 2) + 1e-300
                    x[hr] += sh * ux[hr] / upar
                    y[hr] += sh * uy[hr] / upar
                uz[ht[refl]] *= -1
                te = ht[trans]
                if te.size:
                    sin_air = np.clip(nm * np.sqrt(np.clip(1 - uz[te] ** 2, 0, 1)), 0, 1)
                    esc["il"].append(il[te]); esc["th_air"].append(np.arcsin(sin_air))
                    esc["pol"].append(pol[te]); esc["n_top_hits"].append(hits[te])
                    if record_xy:
                        esc["x"].append(x[te].copy()); esc["y"].append(y[te].copy())
                    alive[te] = False
                ab = ht[(~refl) & (~trans)]
                alive[ab] = False
                fate["top_abs"] += ab.size
            hbot = hb[~up]
            if hbot.size:
                z[hbot] = 0.0
                cb = np.clip(-uz[hbot], 0, 1)
                bbin = np.minimum(np.searchsorted(th_edges, np.arccos(cb)) - 1, nth - 1)
                Rv = Rb[pol[hbot], il[hbot], bbin]
                refl = rng.random(hbot.size) < Rv
                if shifts is not None:
                    hr = hbot[refl]
                    sh = shifts[1][pol[hr], il[hr], bbin[refl]]
                    upar = np.sqrt(ux[hr] ** 2 + uy[hr] ** 2) + 1e-300
                    x[hr] += sh * ux[hr] / upar
                    y[hr] += sh * uy[hr] / upar
                uz[hbot[refl]] *= -1
                ab = hbot[~refl]
                alive[ab] = False
                fate["bot_abs"] += ab.size
        # --- volume events
        hv = ids[~hitb]
        if hv.size:
            sv = s[~hitb]
            x[hv] += ux[hv] * sv
            y[hv] += uy[hv] * sv
            z[hv] += uz[hv] * sv
            scat = rng.random(hv.size) < (mus_um / mut)
            ab = hv[~scat]
            alive[ab] = False
            fate["sl_abs"] += ab.size
            sc = hv[scat]
            if sc.size:
                ct = _hg_cos(g, rng.random(sc.size))
                st = np.sqrt(np.clip(1 - ct * ct, 0, 1))
                ph = 2 * np.pi * rng.random(sc.size)
                cp, sp = np.cos(ph), np.sin(ph)
                ux0, uy0, uz0 = ux[sc], uy[sc], uz[sc]
                tmp = np.sqrt(np.clip(1 - uz0 * uz0, 1e-12, 1))
                big = np.abs(uz0) > 0.99999
                nux = np.where(big, st * cp, st * (ux0 * uz0 * cp - uy0 * sp) / tmp + ux0 * ct)
                nuy = np.where(big, st * sp, st * (uy0 * uz0 * cp + ux0 * sp) / tmp + uy0 * ct)
                nuz = np.where(big, np.sign(uz0) * ct, -st * cp * tmp + uz0 * ct)
                nrm = np.sqrt(nux * nux + nuy * nuy + nuz * nuz)
                ux[sc], uy[sc], uz[sc] = nux / nrm, nuy / nrm, nuz / nrm
                pol[sc] = (rng.random(sc.size) < 0.5).astype(pol.dtype)
    fate["lost"] = int(alive.sum())
    if fate["lost"] > 1e-4 * n_photons:
        print(f"WARNING monte_carlo: {fate['lost']} photons still alive after {max_steps} steps")
    for k in esc:
        esc[k] = np.concatenate(esc[k]) if len(esc[k]) else np.array([])
    scale = total_in / n_photons  # per exciton per photon
    return dict(esc=esc, fate=fate, scale=scale, n=n_photons, arrivals=arrivals * scale, steps=step)


def mc_on_axis_ce(mc, tables, wl, V, n_m, th_deg, cone_deg=10.0):
    """On-axis current efficiency of the recycled light from MC.

    Two estimators: (a) escaped photons within a small cone; (b) radiance theorem using
    the internal near-normal arrival flux x T_top(0) / n_m^2."""
    esc = mc["esc"]
    nwl = len(wl)
    th_air = esc["th_air"]
    il = esc["il"].astype(int)
    c = np.deg2rad(cone_deg)
    sel = th_air < c
    N = np.bincount(il[sel], minlength=nwl) * mc["scale"] / (np.pi * np.sin(c) ** 2)
    ce_a = KM * np.sum(V * (HC_E / wl) * N)
    # internal arrivals: radiance at normal from bins with theta_m < cone_m
    nm = float(np.real(n_m))
    cm = np.arcsin(np.sin(c) / nm)
    edges = np.deg2rad(np.linspace(0, 90, len(th_deg) + 1))
    nb = max(1, np.searchsorted(edges, cm) - 1)
    cm_eff = edges[nb]
    arr = mc["arrivals"][:, :, :nb].sum(axis=(0, 2))
    Lm = arr / (np.pi * np.sin(cm_eff) ** 2)
    T0 = 0.5 * (tables["Tt"][0, :, 0] + tables["Tt"][1, :, 0])
    ce_b = KM * np.sum(V * (HC_E / wl) * T0 * Lm / nm ** 2)
    return ce_a, ce_b


def psf_radii(mc, fractions=(0.5, 0.9)):
    esc = mc["esc"]
    if esc["x"].size == 0:
        return [np.nan for _ in fractions]
    r = np.sqrt(esc["x"] ** 2 + esc["y"] ** 2)
    return [float(np.quantile(r, f)) for f in fractions]
