"""
analyze.py -- post-processing of optimised designs: angular profiles, A' budget,
cone-resolved efficiencies, Monte-Carlo blur, and generic design charts.

usage: python3 analyze.py   (reads results/opt_*.json, writes results/figs/*.png and results/summary.json)
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oledsim import KM, HC_E  # noqa: E402
from devices import (run_proposed, run_reference, angular_profile_proposed,  # noqa: E402
                     angular_profile_reference)
import recycling as rc  # noqa: E402
import optimize as opt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "results", "figs")
os.makedirs(FIG, exist_ok=True)
THETA = np.linspace(0, 89, 179)


def load(tag):
    p = os.path.join(HERE, "results", f"opt_{tag}.json")
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        d = json.load(fh)
    return d


def _restore(P):
    P = dict(P)
    if isinstance(P.get("cap"), list):
        P["cap"] = tuple(P["cap"])
    if isinstance(P.get("spacer"), list):
        P["spacer"] = tuple(P["spacer"])
    return P


def evaluate(d):
    """re-run an optimised design on the full grid and collect rich outputs"""
    pm, pf, _ = opt.preset_mats(d.get("preset"))
    ctx, s0 = opt.make_ctx(d["emitter"], mats=pm, full=True)
    P = _restore(d["params"])
    info = dict(tag=d.get("tag"), job=d["job"], emitter=d["emitter"], preset=d.get("preset"), params=d["params"])
    k = KM * ctx.V * HC_E / ctx.wl
    if d["job"] in ("ref", "refdbr"):
        run = run_reference(ctx, P, s0, du=2.5e-4)
        I, lum = angular_profile_reference(ctx, run, THETA)
        info.update(run["summary"])
        info["_run"] = run
        info["lum_theta"] = lum.tolist()
        I_tot = I
    else:
        run = run_proposed(ctx, P, s0, du=2.5e-4, a_sl_pass=d.get("a_sl", 0.0))
        I_dir, I_rec, lum_dir, lum_rec = angular_profile_proposed(ctx, run, THETA)
        info.update(run["summary"])
        info["lum_theta"] = (lum_dir + lum_rec).tolist()
        info["lum_theta_dir"] = lum_dir.tolist()
        I_tot = I_dir + I_rec
        # A' budget: spectrally (P_dn-) weighted Lambertian absorption per layer of the top stack
        wts = run["out"]["P_dn"] / run["out"]["P_dn"].sum()
        Al = run["tables"]["Al"]  # (2, nlay, nwl, nth)
        lam = 0.5 * np.sum(Al * run["wL"][None, None, None, :], axis=(0, 3))  # (nlay, nwl)
        names = [l.name for l in run["top_stack"]]
        budget = {}
        for n_, a in zip(names, lam @ wts):
            key = "DBR" if n_[0] in "HL" and n_[1:].isdigit() else n_
            budget[key] = budget.get(key, 0.0) + float(a)
        # absorption seen by top-reflected (recycled) light is what matters; give both
        info["Aprime_lambertian_layers"] = budget
        Rb = 0.5 * np.sum(run["tables"]["Rb"] * run["wL"][None, None, :], axis=(0, 2))
        info["mirror_loss_lambertian"] = float(1 - Rb @ wts)
        info["T0_peak"] = float(run["out"]["T0"][np.argmax(s0)])
        keep = 1.0 - d.get("a_sl", 0.0)  # SL absorption per visit belongs to the round-trip loss
        info["Aprime_eff"] = float(1 - (1 - info["a_top_avg"] / (1 - info["p_avg"])) * info["rb_avg"] * keep)
        nm = P.get("n_m", 1.85)
        info["G_radiance_gain"] = float(np.sum(wts * run["out"]["T0"]) / (nm ** 2 * (info["p_avg"] + (1 - info["p_avg"]) * info["Aprime_eff"])))
        info["_run"] = run
    # cone-resolved photon output (per exciton)
    th = np.deg2rad(THETA)
    dOm = 2 * np.pi * np.sin(th) * np.gradient(th)
    photons_theta = np.sum(I_tot, axis=0) * dOm
    for c in (10, 20, 30, 45):
        info[f"eqe_within_{c}deg"] = float(np.sum(photons_theta[THETA <= c]))
    info["eqe_total_from_profile"] = float(np.sum(photons_theta))
    info["CE_check"] = float(info["lum_theta"][0])
    info["_ctx"] = ctx
    info["_s0"] = s0
    return info


def figure_angular(infos, fname, title):
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for inf in infos:
        if inf is None:
            continue
        lum = np.array(inf["lum_theta"])
        ax[0].plot(THETA, lum / np.cos(np.deg2rad(THETA)).clip(1e-3), label=f"{inf['label']} ({inf['CE']:.0f} cd/A)")
        ax[1].plot(THETA, lum / lum[0], label=inf["label"])
    ax[1].plot(THETA, np.cos(np.deg2rad(THETA)), "k:", label="Lambertian")
    ax[0].set_xlabel("angle in air (deg)")
    ax[0].set_ylabel("luminance / current density (cd/A)")
    ax[1].set_xlabel("angle in air (deg)")
    ax[1].set_ylabel("normalised luminous intensity")
    ax[0].set_xlim(0, 80)
    ax[1].set_xlim(0, 80)
    ax[0].legend(fontsize=8)
    ax[1].legend(fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, fname), dpi=150)
    plt.close(fig)


def figure_gain_map(fname):
    """generic: on-axis radiance gain G = T0/(n^2 (p + (1-p)A')) and eta_rec vs p, A'"""
    p = np.linspace(0.005, 0.30, 200)
    A = np.linspace(0.0, 0.25, 200)
    PP, AA = np.meshgrid(p, A)
    nm, T0 = 1.85, 0.85
    G = T0 / (nm ** 2 * (PP + (1 - PP) * AA))
    eta = PP / (PP + (1 - PP) * AA)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    cs = ax[0].contourf(PP, AA, np.log2(G), levels=np.arange(-3, 4.01, 0.5), cmap="viridis")
    c2 = ax[0].contour(PP, AA, G, levels=[0.5, 1, 1.5, 2, 3, 4, 6], colors="w", linewidths=0.8)
    ax[0].clabel(c2, fmt="%.1fx")
    ax[0].set_xlabel("escape probability per pass p (Lambertian-averaged)")
    ax[0].set_ylabel("round-trip loss A'")
    ax[0].set_title("on-axis radiance gain G (T0=0.85, n_m=1.85)")
    fig.colorbar(cs, ax=ax[0], label="log2 G")
    cs = ax[1].contourf(PP, AA, eta, levels=np.linspace(0, 1, 11), cmap="magma")
    c2 = ax[1].contour(PP, AA, eta, levels=[0.2, 0.4, 0.6, 0.8], colors="w", linewidths=0.8)
    ax[1].clabel(c2, fmt="%.1f")
    ax[1].set_xlabel("p")
    ax[1].set_ylabel("A'")
    ax[1].set_title("recycled extraction efficiency p/(p+(1-p)A')")
    fig.colorbar(cs, ax=ax[1])
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, fname), dpi=150)
    plt.close(fig)


def mc_blur(info, settings, gh=True):
    run = info["_run"]
    ctx = info["_ctx"]
    rows = []
    shifts = None
    if gh:
        from devices import proposed_layers
        _, _, _, top_stack, bot = proposed_layers(ctx, run["P"])
        shifts = rc.shift_tables(top_stack, run["P"]["n_m"], ctx.wl, run["th"], bot_layers=bot,
                                 mirror_nk=ctx.nk(ctx.mats["mirror"]))
    for (d_um, mus, g, mua) in settings:
        mc = rc.monte_carlo(run["tables"], run["th"], ctx.wl, run["Pin"], run["qw"], run["P"]["n_m"],
                            d_um, mus, g, mua_um=mua, n_photons=120000, max_steps=20000, seed=7, shifts=shifts)
        eqe_rec = mc["esc"]["il"].size * mc["scale"]
        ce_a, ce_b = rc.mc_on_axis_ce(mc, run["tables"], ctx.wl, ctx.V, run["P"]["n_m"], run["th"], cone_deg=10)
        r50, r90 = rc.psf_radii(mc)
        rows.append(dict(d_um=d_um, mus=mus, g=g, mua=mua, S=mus * d_um, Sred=mus * d_um * (1 - g),
                         eqe_rec=float(eqe_rec), CE_rec_cone=float(ce_a), CE_rec_rad=float(ce_b),
                         r50_um=r50, r90_um=r90, fate={k: int(v) for k, v in mc["fate"].items()}))
    return rows


def clean(info):
    return {k: v for k, v in info.items() if not k.startswith("_")}


def cie_xy(ctx, spec):
    xb, yb, zb = ctx.db.xyz_bar(ctx.wl)
    e = spec * (HC_E / ctx.wl)  # photon -> energy weighting
    X, Y, Z = np.sum(xb * e), np.sum(yb * e), np.sum(zb * e)
    return X / (X + Y + Z), Y / (X + Y + Z)


def on_axis_spectrum(inf):
    """photon radiance spectrum on axis (per exciton per sr)"""
    run = inf["_run"]
    if inf["job"].startswith("ref"):
        from oledsim import radiant_intensity
        Is, Ip = radiant_intensity(run["res"], "top", 1.0, np.array([0.0]))
        return run["qw"] * (Is + Ip)[:, 0]
    o = run["out"]
    return o["L_dir0"] + o["L_rec0"]


def figure_budget(infos, fname):
    labs, rows = [], []
    keys = ["ag_top", "tco", "DBR", "mirror"]
    for inf in infos:
        if inf is None or "Aprime_lambertian_layers" not in inf:
            continue
        b = inf["Aprime_lambertian_layers"]
        rows.append([b.get("ag_top", 0), b.get("tco", 0), b.get("DBR", 0), inf["mirror_loss_lambertian"]])
        labs.append(inf["label"].replace("prop3_narrow528", "prop").replace("_", "\n"))
    if not rows:
        return
    rows = np.array(rows) * 100
    fig, ax = plt.subplots(figsize=(max(6, 1.1 * len(rows)), 4))
    bottom = np.zeros(len(rows))
    for j, k in enumerate(["top thin Ag/IZO", "bottom TCO", "DBR", "Ag mirror"]):
        ax.bar(range(len(rows)), rows[:, j], bottom=bottom, label=k)
        bottom += rows[:, j]
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(labs, fontsize=7)
    ax.set_ylabel("absorption per encounter, Lambertian (%)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, fname), dpi=150)
    plt.close(fig)


def figure_ce_eqe(allinfo, fname):
    fig, ax = plt.subplots(figsize=(7, 5))
    for t, inf in allinfo.items():
        mk = "o" if inf["job"].startswith("ref") else ("s" if "hyb" in inf["job"] else "^")
        col = "C0" if inf["job"].startswith("ref") else "C3"
        ax.scatter(inf["eqe"] * 100, inf["CE"], marker=mk, color=col)
        ax.annotate(t.replace("_narrow528", "").replace("prop3_", "p:").replace("ref_", "r:"), (inf["eqe"] * 100, inf["CE"]), fontsize=6)
    ax.set_xlabel("total EQE (%)")
    ax.set_ylabel("on-axis current efficiency (cd/A)")
    ax.set_title("optimised designs (blue o: conventional cavity, red ^: angular recycling)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, fname), dpi=150)
    plt.close(fig)


def loss_waterfall(tag):
    """Re-evaluate an optimised proposed design with individual loss channels switched off
    (geometry fixed): shows how much CE/EQE each A' component costs."""
    d = load(tag)
    pm, pf, _ = opt.preset_mats(d.get("preset"))
    P = _restore(d["params"])
    grid = np.arange(400.0, 801.0)
    base_mats = dict(pm)
    db = opt.DB
    agname = base_mats.get("ag_top", "l_Ag_15_onHATCN")
    tconame = base_mats.get("tco", "l_ITO")
    ag_arr = db.nk(agname, grid) if isinstance(agname, str) else agname[1]
    tco_arr = db.nk(tconame, grid) if isinstance(tconame, str) else tconame[1]
    steps = [
        ("as optimised", {}, 1.0),
        ("+ perfect mirror (R=1)", {}, 0.0),
        ("+ lossless TCO (k=0)", {"tco": ("tco_k0", tco_arr.real + 0j)}, 0.0),
        ("+ lossless top metal (n->0)", {"tco": ("tco_k0", tco_arr.real + 0j), "ag_top": ("ag_n0", 1e-4 + 1j * ag_arr.imag)}, 0.0),
    ]
    out = []
    for name, mats, ms in steps:
        ctx, s0 = opt.make_ctx(d["emitter"], mats={**base_mats, **mats}, full=True)
        ctx.k_floor = {"H": 1e-6, "L": 1e-6, "bH": 1e-6, "bL": 1e-6} if "lossless" in name else ctx.k_floor
        r = run_proposed(ctx, P, s0, du=2.5e-4, mirror_scale=ms, a_sl_pass=d.get("a_sl", 0.0))
        sm = r["summary"]
        out.append(dict(step=name, CE=sm["CE"], eqe=sm["eqe"], p=sm["p_avg"], a_top=sm["a_top_avg"], rb=sm["rb_avg"]))
    return out


def polarized_ce(info, t_good=0.95, r_bad=0.95, t_abs_pol=0.90):
    """On-axis CE in ONE linear polarisation for the recycling design with
    (a) an absorbing polariser on top (transmits t_abs_pol of the good polarisation), and
    (b) a reflective polariser (transmits t_good of the good polarisation, reflects r_bad of the
        other back into the depolarising recycling cavity).
    The polariser is treated as an incoherent sheet above the stack (normal-incidence values used
    for all angles)."""
    run = info["_run"]
    ctx = info["_ctx"]
    o = run["out"]
    nm = run["P"]["n_m"]
    k = KM * ctx.V * HC_E / ctx.wl
    p, r, rb, phi1, T0 = o["p"], o["r"], o["rb"], o["phi1"], o["T0"]
    qw = run["qw"]
    # (a) absorbing polariser: cavity unchanged, output halved x t_abs_pol
    phi_a = phi1 / (1 - r * rb)
    L_a = 0.5 * t_abs_pol * T0 * phi_a / (np.pi * nm ** 2)
    # (b) reflective polariser: half of the escaping light (bad pol) returns with r_bad
    #     (its re-entry through the stack has transmission ~ the same T; approximated by
    #      adding r_bad * p/2 to the reflected fraction)
    Tt = run["tables"]["Tt"]  # (2, nwl, nth); by reciprocity the air->SL transmission equals SL->air
    T2 = 0.5 * np.sum(Tt ** 2 * run["wL"][None, None, :], axis=(0, 2))
    r_eff = r + 0.5 * r_bad * T2
    phi_b = phi1 / (1 - r_eff * rb)
    L_b = 0.5 * t_good * T0 * phi_b / (np.pi * nm ** 2)
    ce_a = float(np.sum(k * qw * L_a))
    ce_b = float(np.sum(k * qw * L_b))
    return dict(CE_pol_absorbing=ce_a, CE_pol_reflective=ce_b, gain=ce_b / ce_a)


if __name__ == "__main__":
    figure_gain_map("gain_map.png")
    only = sys.argv[1:] if len(sys.argv) > 1 else None
    tags = sorted(os.path.basename(p)[4:-5] for p in glob.glob(os.path.join(HERE, "results", "opt_*.json")))
    if only:
        tags = [t for t in tags if t in only]
    allinfo = {}
    rows = []
    for t in tags:
        d = load(t)
        d["tag"] = t
        try:
            inf = evaluate(d)
        except Exception as e:  # keep going
            print("FAILED", t, e)
            continue
        inf["label"] = t
        if not inf["job"].startswith("ref"):
            inf.update(polarized_ce(inf))
        inf["CIE_onaxis"] = cie_xy(inf["_ctx"], on_axis_spectrum(inf))
        allinfo[t] = inf
        rows.append(t)
        print(t, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in clean(inf).items()
                  if k not in ("lum_theta", "lum_theta_dir", "params", "Aprime_lambertian_layers")})
    groups = [("narrow528", ["ref_narrow528_fine", "refdbr3_narrow528_fine", "prop3_narrow528_ext",
                             "prop3_narrow528_cons", "prop3_narrow528_real", "prop3_narrow528_best",
                             "prop3_narrow528_IZOtop"]),
              ("irppy2acac", ["ref_irppy2acac", "prop3_irppy2acac"]),
              ("blue", ["ref_vdabna_hf_DE", "prop3_vdabna_hf_DE", "prop3_vdabna_hf_best"]),
              ("red", ["ref_red_narrow_DE", "prop3_red_narrow_DE", "prop3_red_narrow_best"])]
    for name, members in groups:
        infos = [allinfo.get(m) for m in members if m in allinfo]
        if infos:
            figure_angular(infos, f"angular_{name}.png", f"optimised designs, emitter {name}")
    figure_budget([allinfo.get(t) for t in ["prop3_narrow528_ext", "prop3_narrow528_Ag8", "prop3_narrow528_Ag12",
                                           "prop3_narrow528_Ag12ODR", "prop3_narrow528_cons", "prop3_narrow528_real",
                                           "prop3_narrow528_bulkAg", "prop3_narrow528_ODR", "prop3_narrow528_lowkITO",
                                           "prop3_narrow528_lowHTL", "prop3_narrow528_best", "prop3_narrow528_IZOtop"]
                   if t in allinfo], "aprime_budget.png")
    figure_ce_eqe(allinfo, "ce_vs_eqe.png")
    wf = {}
    for t in ("prop3_narrow528_ext", "prop3_narrow528_cons"):
        if t in allinfo:
            wf[t] = loss_waterfall(t)
            for r in wf[t]:
                print("WATERFALL", t, r)
    blur = {}
    for t in ("prop3_narrow528_ext", "prop3_narrow528_cons"):
        if t in allinfo:
            blur[t] = mc_blur(allinfo[t], [(1.0, 2.5, 0.75, 0.0), (2.0, 2.5, 0.75, 0.0), (4.0, 2.5, 0.75, 0.0),
                                           (2.0, 5.0, 0.6, 0.0), (0.5, 10.0, 0.5, 0.0), (2.0, 2.5, 0.75, 0.002),
                                           (2.0, 0.5, 0.75, 0.0)])
            for r in blur[t]:
                print("MC", t, r)
    out = dict(designs={t: clean(i) for t, i in allinfo.items()}, mc_blur=blur, waterfall=wf)
    sp = os.path.join(HERE, "results", "summary.json")
    if only and os.path.exists(sp):  # merge into an existing summary
        old = json.load(open(sp))
        old["designs"].update(out["designs"])
        for k in ("mc_blur", "waterfall"):
            old.setdefault(k, {}).update(out[k])
        out = old
    with open(sp, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
