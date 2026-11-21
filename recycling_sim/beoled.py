"""
beoled.py -- bottom-emitting (BE) variant of the angle-selective recycling OLED:

    air / glass (thick, incoherent) / oxide DBR (sputtered TiO2/SiO2) / ITO cathode / EIL / ETL / EML /
    TCTA / HTL / HATCN / top anode (thin Ag or IZO) / SL (n_m, randomising) / SiO2 ODR spacer / Ag mirror

Light leaves through the DBR and the glass.  The SL + Ag mirror sit on top of the OLED.
The thick glass is handled incoherently.  Light that enters the glass outside the air cone (beta > 1) is
totally reflected at glass/air; two limits:
  'large'  -- it returns into the stack (large-area device, laterally infinite) and is recycled;
  'pixel'  -- it is lost (small pixel: it walks millimetres sideways in the glass before returning).
Same platform otherwise as the top-emitting comparison (compare_anodes.py): ODR, 15-nm low-k ITO,
organic k 1e-4.  DBR oxides get a k floor of 2e-5 (sputtered; evaporated ZnS/LiF in the TE design: 1e-4).

usage: python3 beoled.py run      (DE for 4 cases, writes results/be_*.json)
       python3 beoled.py table
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.optimize import differential_evolution

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import optimize as opt  # noqa: E402
import recycling as rc  # noqa: E402
from oledsim import Medium, dipole_emission, stack_RTA  # noqa: E402
from devices import _dbr_layers, photon_weights  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
N_GLASS = 1.52
K_OXIDE = 2e-5


def ctx_for(emitter, preset, full=False, k_oxide=K_OXIDE):
    """objective grid: 1 nm (the Ag/organic/DBR stack is a resonant microcavity; coarser grids alias);
    final grid (full=True): 0.5 nm over 400-700 nm."""
    from devices import Ctx
    pm, pf, ag_b = opt.preset_mats(preset)
    if full:
        ctx = Ctx(opt.DB, wl=np.arange(400.0, 700.01, 0.5), mats=pm)
        s0 = ctx.emitter(**opt.EMITTERS[emitter])
    else:
        ctx, s0 = opt.make_ctx(emitter, mats=pm, step=1.0)
    ctx.k_floor["bH"] = k_oxide; ctx.k_floor["bL"] = k_oxide
    return ctx, s0, pf, ag_b


def _fresnel_ga(beta, pol):
    """power reflectance glass -> air for in-plane wavevector beta (=n sin theta); 1 beyond the air cone."""
    ng = N_GLASS
    cg = np.sqrt(np.clip(1 - (beta / ng) ** 2, 0, 1)) + 0j
    ca = np.sqrt(1 - beta ** 2 + 0j)
    if pol == "s":
        r = (ng * cg - ca) / (ng * cg + ca)
    else:
        r = (cg - ng * ca) / (cg + ng * ca)
    R = np.abs(r) ** 2
    return np.where(beta >= 1.0, 1.0, R)


def run_be(ctx, s0, P, mode="pixel", n_th=91, du=1e-3):
    nm = P["n_m"] * ctx.ones
    up = [ctx.layer("etl", P["d_etl"]), ctx.layer("eil", 10.0), ctx.layer("tco", P.get("d_tco", 15.0))] + \
        _dbr_layers(ctx, P["dbr"], hk="bH", lk="bL")
    up = [l for l in up if l.d > 0]
    dn = [ctx.layer("tcta", 10.0), ctx.layer("htl", P["d_htl"]), ctx.layer("hatcn", P.get("d_hat", 10.0)),
          ctx.layer("ag_top", P["d_ag"])]
    eml = ctx.layer("eml", 20.0)
    glass = Medium(N_GLASS * ctx.ones + 0j)
    res = dipole_emission(ctx.wl, eml.no, 10.0, 10.0, up, glass, dn, Medium(nm + 0j), theta_h=0.8, u_max=3.0, du=du)
    qw, qeff, Fbar = photon_weights(res, s0, 1.0)
    th, wL, edges = rc.angle_grid(n_th)
    stack = list(reversed(dn)) + [eml] + up                    # from the SL towards the glass
    beta = np.real(nm)[:, None] * np.sin(np.deg2rad(th))[None, :]
    sl = Medium(nm + 0j)
    T = {k: [] for k in ("Rt", "Tt", "At", "Rb")}
    bot = [ctx.layer("l_SiO2_SJ_RTP", 800.0)]
    for pol in ("s", "p"):
        Rs, Ts, As = stack_RTA(stack, sl, glass, beta, ctx.wl, pol)
        As = As.sum(0)
        ok = beta < N_GLASS - 1e-9
        bg = np.where(ok, beta, 0.0)
        Rs2, _, _ = stack_RTA(list(reversed(stack)), glass, sl, bg, ctx.wl, pol)   # seen from the glass
        Rs2 = np.where(ok, Rs2, 0.0); Ts = np.where(ok, Ts, 0.0)
        Rga = _fresnel_ga(beta, pol)
        den = 1 - Rs2 * Rga
        Tt = Ts * (1 - Rga) / den
        if mode == "large":
            Rt = Rs + Ts * Ts * Rga / den
        else:   # glass-trapped light (beta >= 1) lost; inside the air cone keep the multiple reflections
            trap = beta >= 1.0
            Rt = np.where(trap, Rs, Rs + Ts * Ts * Rga / den)
        At = np.clip(1 - Rt - Tt, 0, 1)
        Rb, _, _ = stack_RTA(bot, sl, Medium(ctx.nk(ctx.mats["mirror"])), beta, ctx.wl, pol)
        for k, v in zip(("Rt", "Tt", "At", "Rb"), (Rt, Tt, At, Rb)):
            T[k].append(v)
    tables = {k: np.array(v) for k, v in T.items()}
    Pin = rc.emission_into_sl(res, nm, th, edges)
    out = rc.analytic(res, tables, nm, ctx.wl, th, wL, Pin, qw)
    # direct (first-pass) emission into the glass: only the air cone escapes, through glass/air
    Tga0 = 1 - 0.5 * (_fresnel_ga(np.array([0.0]), "s") + _fresnel_ga(np.array([0.0]), "p"))[0]
    out["L_dir0"] = out["L_dir0"] * Tga0
    u = res["u"]; ne = np.real(res["n_eml"])[:, None]; b = ne * u[None, :]
    esc = np.zeros(len(ctx.wl))
    for pol in ("s", "p"):
        Pt = res[f"Pt_{pol}"]
        esc += np.sum(Pt * (1 - _fresnel_ga(b, pol)) * (b < 1) * 2 * u[None, :] * res["du"], axis=1)
    out["eqe_dir"] = qw * esc
    s = rc.summarize(out, ctx.wl, ctx.V)
    s["Aprime"] = 1 - (1 - s["a_top_avg"] / (1 - s["p_avg"])) * s["rb_avg"]
    return s


def params(x, N):
    return dict(d_etl=x[0], d_htl=x[1], d_ag=x[2], dbr=dict(N=N, dH1=x[3], dH=x[4], dL=x[5]), n_m=x[6])


G = {}


def _obj(x):
    N = int(round(x[7]))
    return -run_be(G["ctx"], G["s0"], params(x, N), mode=G["mode"], n_th=121)["CE"]


CASES = {
    "be_n15_Ag8m": ("narrow532_15", "Ag8meas+ODR+lowkITO+korg1e-4"),
    "be_n15_IZO": ("narrow532_15", "IZOtop+ODR+lowkITO+korg1e-4"),
    "be_irppy_Ag8m": ("irppy2acac", "Ag8meas+ODR+lowkITO+korg1e-4"),
    "be_irppy_IZO": ("irppy2acac", "IZOtop+ODR+lowkITO+korg1e-4"),
}


def run(tag, mode):
    em, preset = CASES[tag]
    ctx, s0, pf, ag_b = ctx_for(em, preset)
    G.update(ctx=ctx, s0=s0, mode=mode)
    b = [(10, 500), (20, 600), ag_b or (6, 25), (0, 220), (20, 200), (20, 300), (1.6, 2.0), (2, 8.49)]
    r = differential_evolution(_obj, b, maxiter=30, popsize=10, seed=2, tol=1e-6, polish=False, workers=4,
                               updating="deferred")
    N = int(round(r.x[7])); P = params(r.x, N)
    out = dict(tag=tag, optimised_for=mode, emitter=em, preset=preset, N=N, x=list(map(float, r.x)), params=P)
    ctxf, s0f, _, _ = ctx_for(em, preset, full=True)
    for m in ("pixel", "large"):
        out[m] = run_be(ctxf, s0f, P, mode=m, n_th=181, du=2.5e-4)
    ctxk, s0k, _, _ = ctx_for(em, preset, full=True, k_oxide=1e-4)
    out["koxide1e-4"] = run_be(ctxk, s0k, P, mode=mode, n_th=181, du=2.5e-4)
    with open(os.path.join(HERE, "results", f"{tag}_{mode}.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(tag, "opt for", mode, "N", N, "| pixel CE", round(out["pixel"]["CE"]), "A'", round(100 * out["pixel"]["Aprime"], 2),
          "| large CE", round(out["large"]["CE"]), "A'", round(100 * out["large"]["Aprime"], 2),
          "| EQE pix/large", round(out["pixel"]["eqe"], 3), round(out["large"]["eqe"], 3), flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "run":
        for t in (sys.argv[2:] or CASES):
            for mode in ("pixel", "large"):
                run(t, mode)
