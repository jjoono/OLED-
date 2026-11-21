"""
arch.py -- unified architecture search for the angle-selective recycling OLED (maximum on-axis cd/A).

arch  'TE' : substrate / Ag mirror / SiO2 ODR / SL / ITO cathode / EIL / ETL / EML / TCTA / HTL / HATCN /
             top anode / evaporated ZnS-LiF DBR / exit                      (light leaves upward)
      'BE' : exit / sputtered TiO2-SiO2 DBR / ITO cathode / EIL / ETL / EML / TCTA / HTL / HATCN /
             top anode / SL / SiO2 ODR / Ag mirror                          (light leaves through the substrate)
exit  'air'      : the DBR faces air (TE lab device under a gas-filled lid)
      'glass'    : the DBR faces a thick n = 1.52 medium (BE substrate; TE with an OCA-bonded cover glass)
      'lowN'     : a low-index layer (n_low, thickness optimised) between the DBR and the thick glass
glass handling (incoherent): 'pixel' = light trapped in the glass (beta >= 1) is lost (it walks mm sideways);
                             'large' = it returns into the stack and is recycled (laterally infinite device)
Common platform: 800-nm SiO2 ODR, 15-nm low-k ITO cathode, organic k 1e-4, PLQY 1, horizontal dipoles 0.8,
SL index <= 2.0.  DBR k floor: evaporated ZnS/LiF 1e-4, sputtered TiO2/SiO2 2e-5.
Optimisation: differential evolution over thicknesses, SL index and DBR pair count, then bounded Nelder-Mead;
objective on a 1-nm (narrow emitter) or 2-nm (Ir(ppy)2acac) grid; final evaluation on a 0.5-nm grid.

usage: python3 arch.py run <case> [<case> ...]     (case names: see CASES; writes results/arch/<case>.json)
       python3 arch.py all                         (every case not yet done)
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
from scipy.optimize import differential_evolution, minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import optimize as opt  # noqa: E402
import recycling as rc  # noqa: E402
from oledsim import Medium, dipole_emission, stack_RTA  # noqa: E402
from devices import Ctx, _dbr_layers, photon_weights  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "arch")
N_GLASS = 1.52
K_EVAP, K_SPUT = 1e-4, 2e-5
EMITTERS = {"irppy": "irppy2acac", "n15": "narrow532_15"}
ANODES = {"Ag8m": "Ag8meas", "IZO": "IZOtop"}


def make(emitter, anode, final=False, k_dbr=None):
    pm, pf, ag_b = opt.preset_mats(f"{ANODES[anode]}+lowkITO+korg1e-4")
    em = EMITTERS[emitter]
    if final:
        ctx = Ctx(opt.DB, wl=np.arange(400.0, 700.01, 0.5), mats=pm)
        s0 = ctx.emitter(**opt.EMITTERS[em])
    else:
        step = 1.0 if emitter == "n15" else 2.0
        c0 = Ctx(opt.DB, wl=np.arange(400.0, 700.01, step), mats=pm)
        s_all = c0.emitter(**opt.EMITTERS[em])
        keep = s_all > 1e-2 * s_all.max()                  # spectral range for the objective
        ctx = Ctx(opt.DB, wl=c0.wl[keep.argmax(): len(keep) - keep[::-1].argmax()], mats=pm)
        s0 = ctx.emitter(**opt.EMITTERS[em])
    ctx.k_floor.update(H=K_EVAP, L=K_EVAP, bH=K_SPUT, bL=K_SPUT)
    if k_dbr is not None:
        ctx.k_floor.update(H=k_dbr, L=k_dbr, bH=k_dbr, bL=k_dbr)
    return ctx, s0, ag_b


def _r_glass_air(beta, pol):
    ng = N_GLASS
    cg = np.sqrt(np.clip(1 - (beta / ng) ** 2, 0, 1)) + 0j
    ca = np.sqrt(1 - beta ** 2 + 0j)
    r = (ng * cg - ca) / (ng * cg + ca) if pol == "s" else (cg - ng * ca) / (cg + ng * ca)
    return np.where(beta >= 1.0, 1.0, np.abs(r) ** 2)


def stacks(ctx, P, arch, exit_):
    cath = [ctx.layer("etl", P["d_etl"]), ctx.layer("eil", 10.0), ctx.layer("tco", 15.0)]
    anod = [ctx.layer("tcta", 10.0), ctx.layer("htl", P["d_htl"]), ctx.layer("hatcn", 10.0), ctx.layer("ag_top", P["d_ag"])]
    if arch == "TE":
        up = anod + _dbr_layers(ctx, P["dbr"], hk="H", lk="L"); dn = cath
    else:
        up = cath + _dbr_layers(ctx, P["dbr"], hk="bH", lk="bL"); dn = anod
    if exit_ == "lowN":
        up = up + [ctx.layer(complex(P["n_low"]), P["d_low"], name="lowN")]
    up = [l for l in up if l.d > 0]; dn = [l for l in dn if l.d > 0]
    return up, dn


def run(ctx, s0, P, arch, exit_, mode="pixel", n_th=121, du=1e-3):
    nm = P["n_m"] * ctx.ones
    up, dn = stacks(ctx, P, arch, exit_)
    eml = ctx.layer("eml", 20.0)
    ext = Medium(ctx.ones + 0j) if exit_ == "air" else Medium(N_GLASS * ctx.ones + 0j)
    res = dipole_emission(ctx.wl, eml.no, 10.0, 10.0, up, ext, dn, Medium(nm + 0j), theta_h=0.8, u_max=3.0, du=du)
    qw, qeff, Fbar = photon_weights(res, s0, 1.0)
    th, wL, edges = rc.angle_grid(n_th)
    stack = list(reversed(dn)) + [eml] + up
    beta = np.real(nm)[:, None] * np.sin(np.deg2rad(th))[None, :]
    sl = Medium(nm + 0j)
    bot = [ctx.layer("l_SiO2_SJ_RTP", 800.0)]
    T = {k: [] for k in ("Rt", "Tt", "At", "Rb")}
    for pol in ("s", "p"):
        Rs, Ts, As = stack_RTA(stack, sl, ext, beta, ctx.wl, pol)
        if exit_ == "air":
            Rt, Tt = Rs, Ts
        else:
            ok = beta < N_GLASS - 1e-9
            Rs2, _, _ = stack_RTA(list(reversed(stack)), ext, sl, np.where(ok, beta, 0.0), ctx.wl, pol)
            Rs2 = np.where(ok, Rs2, 0.0); Ts = np.where(ok, Ts, 0.0)
            Rga = _r_glass_air(beta, pol)
            den = 1 - Rs2 * Rga
            Tt = Ts * (1 - Rga) / den
            back = Rs + Ts * Ts * Rga / den
            Rt = back if mode == "large" else np.where(beta >= 1.0, Rs, back)
        At = np.clip(1 - Rt - Tt, 0, 1)
        Rb, _, _ = stack_RTA(bot, sl, Medium(ctx.nk(ctx.mats["mirror"])), beta, ctx.wl, pol)
        for k, v in zip(("Rt", "Tt", "At", "Rb"), (Rt, Tt, At, Rb)):
            T[k].append(v)
    tables = {k: np.array(v) for k, v in T.items()}
    Pin = rc.emission_into_sl(res, nm, th, edges)
    out = rc.analytic(res, tables, nm, ctx.wl, th, wL, Pin, qw)
    if exit_ != "air":   # first-pass light into the glass: only the air cone leaves (through glass/air)
        t0 = 1 - 0.5 * (_r_glass_air(np.array([0.0]), "s") + _r_glass_air(np.array([0.0]), "p"))[0]
        out["L_dir0"] = out["L_dir0"] * t0
        u = res["u"]; b = np.real(res["n_eml"])[:, None] * u[None, :]
        esc = np.zeros(len(ctx.wl))
        for pol in ("s", "p"):
            esc += np.sum(res[f"Pt_{pol}"] * (1 - _r_glass_air(b, pol)) * (b < 1) * 2 * u[None, :] * res["du"], axis=1)
        out["eqe_dir"] = qw * esc
    s = rc.summarize(out, ctx.wl, ctx.V)
    s["Aprime"] = 1 - (1 - s["a_top_avg"] / (1 - s["p_avg"])) * s["rb_avg"]
    return s


# ----------------------------------------------------------------------------------------------
def params(x, N, exit_, n_low):
    P = dict(d_etl=x[0], d_htl=x[1], d_ag=x[2], dbr=dict(N=N, dH1=x[3], dH=x[4], dL=x[5]), n_m=x[6])
    if exit_ == "lowN":
        P.update(n_low=n_low, d_low=x[8])
    return P


CASES = {}
for em in EMITTERS:
    for an in ANODES:
        CASES[f"TE_air_{an}_{em}"] = dict(arch="TE", exit="air", mode="pixel", anode=an, emitter=em)
        CASES[f"TE_glass_{an}_{em}"] = dict(arch="TE", exit="glass", mode="pixel", anode=an, emitter=em)
        CASES[f"BE_glass-pixel_{an}_{em}"] = dict(arch="BE", exit="glass", mode="pixel", anode=an, emitter=em)
        CASES[f"BE_glass-large_{an}_{em}"] = dict(arch="BE", exit="glass", mode="large", anode=an, emitter=em)
        for nl in (1.2, 1.1):
            CASES[f"BE_lowN{nl}_{an}_{em}"] = dict(arch="BE", exit="lowN", mode="pixel", anode=an, emitter=em, n_low=nl)
            CASES[f"TE_lowN{nl}_{an}_{em}"] = dict(arch="TE", exit="lowN", mode="pixel", anode=an, emitter=em, n_low=nl)

G = {}


def _obj(x):
    c = G["case"]
    P = params(x, int(round(x[7])), c["exit"], c.get("n_low"))
    return -run(G["ctx"], G["s0"], P, c["arch"], c["exit"], c["mode"], n_th=91)["CE"]


def optimise(name):
    path = os.path.join(OUT, f"{name}.json")
    if os.path.exists(path):
        return
    c = CASES[name]; t0 = time.time()
    ctx, s0, ag_b = make(c["emitter"], c["anode"])
    G.update(ctx=ctx, s0=s0, case=c)
    b = [(10, 500), (20, 600), ag_b, (0, 220), (20, 220), (20, 300), (1.6, 2.0), (2, 8.49)]
    if c["exit"] == "lowN":
        b.append((300, 2000))
    r = differential_evolution(_obj, b, maxiter=20, popsize=8, seed=3, tol=1e-6, polish=False,
                               workers=int(os.environ.get("ARCH_WORKERS", "4")), updating="deferred")
    N = int(round(r.x[7]))
    # bounded Nelder-Mead polish of the continuous variables at fixed N
    idx = [i for i in range(len(b)) if i != 7 and b[i][1] - b[i][0] > 0.1]
    lo = np.array([b[i][0] for i in idx]); hi = np.array([b[i][1] for i in idx])

    def f(z):
        x = np.array(r.x, float); x[idx] = z
        return _obj(x)
    z0 = np.array([r.x[i] for i in idx])
    sim = [z0] + [np.where(np.arange(len(z0)) == k, np.clip(z0 + 0.06 * np.maximum(abs(z0), 0.05 * (hi - lo)), lo, hi), z0)
                  for k in range(len(z0))]
    m = minimize(f, z0, method="Nelder-Mead", bounds=list(zip(lo, hi)),
                 options=dict(maxfev=200, xatol=0.3, fatol=0.05, initial_simplex=np.array(sim), adaptive=True))
    x = np.array(r.x, float)
    if m.fun < r.fun:
        x[idx] = m.x
    P = params(x, N, c["exit"], c.get("n_low"))
    ctxf, s0f, _ = make(c["emitter"], c["anode"], final=True)
    fin = run(ctxf, s0f, P, c["arch"], c["exit"], c["mode"], n_th=181, du=2.5e-4)
    other = None
    if c["exit"] != "air":
        other = run(ctxf, s0f, P, c["arch"], c["exit"], "large" if c["mode"] == "pixel" else "pixel", n_th=181, du=2.5e-4)
    out = dict(case=name, **c, N=N, x=list(map(float, x)), params=P, final=fin, other_glass_mode=other,
               objective_grid_CE=float(-min(m.fun, r.fun)), seconds=time.time() - t0)
    os.makedirs(OUT, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"{name:32s} N {N}  CE {fin['CE']:7.1f}  EQE {fin['eqe']:.3f}  A' {100 * fin['Aprime']:.2f}%  "
          f"p {fin['p_avg']:.3f}  ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    names = list(CASES) if sys.argv[1] == "all" else sys.argv[2:]
    for n in names:
        optimise(n)
