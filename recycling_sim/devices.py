"""
devices.py -- stack builders and evaluators for
  (A) the proposed angular-recycling TEOLED (inverted):
      Ag mirror / [spacer] / SL (n_m) / TCO / EIL / ETL / EML / TCTA / HTL / HATCN / thin Ag / DBR / air
  (B) a conventional strong-microcavity TEOLED with the same organic stack and top anode:
      Ag / TCO(thin) / EIL / ETL / EML / TCTA / HTL / HATCN / thin Ag / capping (or DBR) / air
"""
from __future__ import annotations

import numpy as np

from oledsim import (MaterialDB, Layer, Medium, dipole_emission, radiant_intensity,
                     energy_to_photon_spectrum, gaussian_spectrum, KM, HC_E)
import recycling as rc

DEFAULT_MATS = dict(
    mirror="Ag_McPeak",          # thick Ag reflector (paper: McPeak, k=3.61 @550)
    ag_top="l_Ag_15_onHATCN",    # measured 15-nm Ag on HATCN
    tco="l_ITO",                 # ITO (Koenig), k=0.0032 @550
    eil="l_ZnOPEIE",
    etl=("l_B3_o_JO", "l_B3_e_JO"),
    eml=("l_TCTA_B3_o_JO", None),  # isotropic approximation for the emitter layer
    tcta=("l_TCTA_o_JO", "l_TCTA_e_JO"),
    htl=("l_TAPC_o_JO", "l_TAPC_e_JO"),
    hatcn="l_HATCN",
    H="l_ZnS",                   # evaporated ZnS (ellipsometry k~0; a k-floor of 1e-4 is applied)
    L="l_LiF",
    cpl="l_NPB",
    bH="l_TiO2_SJ_RTP",          # bottom (pre-OLED, sputter/e-beam) DBR for the hybrid cavity variant
    bL="l_SiO2_SJ_RTP",
)


class Ctx:
    """wavelength grid + materials + photometry"""

    def __init__(self, db: MaterialDB, wl=None, mats=None):
        self.db = db
        self.wl = np.arange(400.0, 701.0, 2.0) if wl is None else np.asarray(wl, float)
        self.mats = dict(DEFAULT_MATS)
        if mats:
            self.mats.update(mats)
        self.V = db.V(self.wl)
        self.ones = np.ones_like(self.wl)
        self.k_floor = {"H": 1e-4, "L": 1e-4, "bH": 1e-4, "bL": 1e-4}  # regularises loss-free DBR-guided modes
        korg = self.mats.pop("_korg", None)  # optional residual absorption (k floor) of all organic layers
        if korg:
            for k in ("eil", "etl", "eml", "tcta", "htl", "hatcn"):
                self.k_floor[k] = float(korg)

    def nk(self, spec):
        if spec is None:
            return None
        if isinstance(spec, tuple) and len(spec) == 2 and isinstance(spec[1], np.ndarray):
            arr = spec[1]  # custom array on 400:1:800 nm
            g = np.arange(400.0, 400.0 + arr.size)
            return np.interp(self.wl, g, arr.real) + 1j * np.interp(self.wl, g, arr.imag)
        if isinstance(spec, (int, float, complex)):
            return complex(spec) * self.ones
        if isinstance(spec, np.ndarray):
            return spec.astype(complex)
        return self.db.nk(spec, self.wl)

    def layer(self, key_or_spec, d, name=""):
        spec = self.mats.get(key_or_spec, key_or_spec) if isinstance(key_or_spec, str) else key_or_spec
        if isinstance(spec, tuple) and len(spec) == 2 and isinstance(spec[1], np.ndarray):
            no, ne = self.nk(spec), None
        elif isinstance(spec, tuple):
            no, ne = self.nk(spec[0]), self.nk(spec[1])
        else:
            no, ne = self.nk(spec), None
        kf = self.k_floor.get(key_or_spec) if isinstance(key_or_spec, str) else None
        if kf:
            no = no.real + 1j * np.maximum(no.imag, kf)
            if ne is not None:
                ne = ne.real + 1j * np.maximum(ne.imag, kf)
        return Layer(no, d, ne=ne, name=name or str(key_or_spec))

    def emitter(self, kind="narrow", **kw):
        """returns photon spectrum s0 (sum=1) on self.wl"""
        if kind == "narrow":
            sE = gaussian_spectrum(self.wl, kw.get("peak", 528), kw.get("fwhm", 25),
                                   shoulder=kw.get("shoulder", 0.15))
        else:
            sE = self.db.spectrum(kind, self.wl)
        return energy_to_photon_spectrum(self.wl, sE)


def _dbr_layers(ctx, dbr, hk="H", lk="L"):
    """dbr: list of thicknesses, alternating H, L starting with H (or dict with N, dH, dL, dH1)."""
    if dbr is None:
        return []
    if isinstance(dbr, dict):
        N = int(dbr["N"])
        th = [dbr.get("dH1", dbr["dH"])]
        for i in range(N):
            th += [dbr["dL"], dbr["dH"]]
        if dbr.get("end_L"):
            th += [dbr["dL"]]
    else:
        th = list(dbr)
    out = []
    for i, t in enumerate(th):
        key = hk if i % 2 == 0 else lk
        out.append(ctx.layer(key, t, name=key + str(i)))
    return out


def photon_weights(res, s0, plqy):
    """homogeneous broadening: photons emitted at each wavelength per exciton."""
    Fbar = np.sum(s0 * res["F"])
    qeff = plqy * Fbar / (1 - plqy + plqy * Fbar)
    return qeff * s0 * res["F"] / Fbar, qeff, Fbar


# ----------------------------------------------------------------------------
# (A) proposed structure
# ----------------------------------------------------------------------------
PROPOSED_DEFAULT = dict(
    n_m=1.85,        # SL matrix (planarisation) index
    d_tco=50.0, d_eil=10.0, d_etl=60.0, d_eml=20.0, z_frac=0.5,
    d_tcta=10.0, d_htl=60.0, d_hat=10.0, d_ag=12.0,
    dbr=dict(N=3, dH=60.0, dL=100.0, dH1=40.0),
    spacer=None,     # e.g. ("l_SiO2_SJ_RTP", 300.0) between SL and Ag mirror
    theta_h=0.8, plqy=1.0,
)


def proposed_layers(ctx, P):
    dn = [ctx.layer("etl", P["d_etl"]), ctx.layer("eil", P["d_eil"]), ctx.layer("tco", P["d_tco"])]
    # optional bottom partial mirror between the TCO and the SL (hybrid coherent cavity + recycling)
    if P.get("bdbr") is not None:
        dn += _dbr_layers(ctx, P["bdbr"], hk="bH", lk="bL")   # first layer (bH) touches the TCO
    if P.get("b_ag"):
        dn += [ctx.layer("ag_top", P["b_ag"], name="bAg")]     # thin Ag semi-mirror below the TCO
    dn = [l for l in dn if l.d > 0]
    up = [ctx.layer("tcta", P["d_tcta"]), ctx.layer("htl", P["d_htl"]), ctx.layer("hatcn", P["d_hat"]),
          ctx.layer("ag_top", P["d_ag"])] + _dbr_layers(ctx, P["dbr"])
    up = [l for l in up if l.d > 0]
    eml = ctx.layer("eml", P["d_eml"])
    top_stack = list(reversed(dn)) + [eml] + up  # from the SL upward
    bot = []
    if P.get("spacer"):
        bot = [ctx.layer(P["spacer"][0], P["spacer"][1])]
    return dn, up, eml, top_stack, bot


def run_proposed(ctx, P, s0, n_th=181, u_max=3.0, du=2.5e-4, a_sl_pass=0.0, mirror_scale=1.0):
    P = {**PROPOSED_DEFAULT, **P}
    dn, up, eml, top_stack, bot = proposed_layers(ctx, P)
    n_m = P["n_m"] * ctx.ones
    zdn = P["d_eml"] * P["z_frac"]
    zup = P["d_eml"] - zdn
    res = dipole_emission(ctx.wl, eml.no, zup, zdn, up, Medium(ctx.ones + 0j), dn, Medium(n_m + 0j),
                          theta_h=P["theta_h"], u_max=u_max, du=du)
    qw, qeff, Fbar = photon_weights(res, s0, P["plqy"])
    th, wL, edges = rc.angle_grid(n_th)
    tables = rc.boundary_tables(top_stack, n_m, ctx.wl, th, bot_layers=bot, mirror_nk=ctx.nk(ctx.mats["mirror"]))
    if mirror_scale != 1.0:  # optional: scale mirror absorption (1-R) for sensitivity studies
        tables["Rb"] = 1 - (1 - tables["Rb"]) * mirror_scale
    Pin = rc.emission_into_sl(res, n_m, th, edges)
    out = rc.analytic(res, tables, n_m, ctx.wl, th, wL, Pin, qw, a_sl_pass=a_sl_pass)
    summ = rc.summarize(out, ctx.wl, ctx.V)
    summ.update(qeff=float(qeff), Fbar=float(Fbar), top_frac=float(np.sum(qw * res["top"])),
                bot_frac=float(np.sum(qw * res["bot"])))
    summ["loss_frac"] = float(qeff) - summ["top_frac"] - summ["bot_frac"]
    return dict(summary=summ, res=res, tables=tables, Pin=Pin, qw=qw, out=out, th=th, wL=wL, edges=edges,
                P=P, top_stack=top_stack)


# ----------------------------------------------------------------------------
# (B) conventional strong-microcavity TEOLED (same top anode)
# ----------------------------------------------------------------------------
REF_DEFAULT = dict(
    d_tco=10.0, d_eil=10.0, d_etl=40.0, d_eml=20.0, z_frac=0.5,
    d_tcta=10.0, d_htl=40.0, d_hat=10.0, d_ag=12.0,
    cap=("cpl", 60.0), dbr=None, theta_h=0.8, plqy=1.0,
)


def run_reference(ctx, P, s0, u_max=3.0, du=2.5e-4):
    P = {**REF_DEFAULT, **P}
    dn = [ctx.layer("etl", P["d_etl"]), ctx.layer("eil", P["d_eil"]), ctx.layer("tco", P["d_tco"])]
    dn = [l for l in dn if l.d > 0]
    up = [ctx.layer("tcta", P["d_tcta"]), ctx.layer("htl", P["d_htl"]), ctx.layer("hatcn", P["d_hat"]),
          ctx.layer("ag_top", P["d_ag"])]
    if P.get("cap"):
        up.append(ctx.layer(P["cap"][0], P["cap"][1]))
    up += _dbr_layers(ctx, P.get("dbr"))
    up = [l for l in up if l.d > 0]
    eml = ctx.layer("eml", P["d_eml"])
    zdn = P["d_eml"] * P["z_frac"]
    zup = P["d_eml"] - zdn
    res = dipole_emission(ctx.wl, eml.no, zup, zdn, up, Medium(ctx.ones + 0j), dn,
                          Medium(ctx.nk(ctx.mats["mirror"])), theta_h=P["theta_h"], u_max=u_max, du=du)
    qw, qeff, Fbar = photon_weights(res, s0, P["plqy"])
    Is, Ip = radiant_intensity(res, "top", 1.0, np.array([0.0]))
    L0 = qw * (Is + Ip)[:, 0]
    CE = KM * np.sum(ctx.V * (HC_E / ctx.wl) * L0)
    eqe = float(np.sum(qw * res["top"]))
    return dict(summary=dict(eqe=eqe, CE=float(CE), qeff=float(qeff), Fbar=float(Fbar)), res=res, qw=qw, P=P)


def angular_profile_reference(ctx, run, theta_deg):
    Is, Ip = radiant_intensity(run["res"], "top", 1.0, theta_deg)
    I = run["qw"][:, None] * (Is + Ip)
    lum = KM * np.sum((ctx.V * HC_E / ctx.wl)[:, None] * I, axis=0)  # cd/A * cos(theta) (intensity)
    return I, lum


def angular_profile_proposed(ctx, run, theta_deg):
    """photon intensity per sr (per exciton) in air vs angle: direct + recycled."""
    Is, Ip = radiant_intensity(run["res"], "top", 1.0, theta_deg)
    I_dir = run["qw"][:, None] * (Is + Ip)
    out = run["out"]
    tables = run["tables"]
    nm = run["P"]["n_m"]
    th_m = np.rad2deg(np.arcsin(np.sin(np.deg2rad(theta_deg)) / nm))
    Tt = 0.5 * (tables["Tt"][0] + tables["Tt"][1])
    Tq = np.array([np.interp(th_m, run["th"], Tt[i]) for i in range(len(ctx.wl))])
    # radiance in air = T * phi_tot/(pi n^2); intensity = radiance * cos(theta_air)
    I_rec = (run["qw"] * out["phi_tot"])[:, None] * Tq / (np.pi * nm ** 2) * np.cos(np.deg2rad(theta_deg))[None, :]
    k = (ctx.V * HC_E / ctx.wl)[:, None] * KM
    return I_dir, I_rec, np.sum(k * I_dir, axis=0), np.sum(k * I_rec, axis=0)
