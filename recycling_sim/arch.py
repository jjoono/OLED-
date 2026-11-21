"""
arch.py -- unified architecture search for the angle-selective recycling OLED (maximum on-axis cd/A).

arch  'TE' : substrate / Ag mirror / SiO2 ODR / SL / ITO cathode / EIL / ETL / EML / TCTA / HTL / HATCN /
             top anode / DBR / exit                                         (light leaves upward)
      'BE' : exit / DBR / ITO cathode / EIL / ETL / EML / TCTA / HTL / HATCN /
             top anode / SL / SiO2 ODR / Ag mirror                          (light leaves through the substrate)
DBR   'ev' : evaporated ZnS/LiF, k floor 1e-4   (default for TE: deposited on the finished OLED)
      'ox' : sputtered TiO2/SiO2, k floor 2e-5  (default for BE: deposited on the bare glass first);
             case suffix '-ox' / '-ev' swaps the default (separates the DBR quality from the architecture)
exit  'air'      : the DBR faces air (TE lab device under a gas-filled lid)
      'glass'    : the DBR faces a thick n = 1.52 medium (BE substrate; TE with an OCA-bonded cover glass)
      'lowN'     : a low-index layer (n_low, thickness optimised) between the DBR and the thick glass
glass handling (incoherent): 'pixel' = light trapped in the glass (beta >= 1) is lost (it walks mm sideways);
                             'large' = it returns into the stack and is recycled (laterally infinite device)
                             first-pass emission reflected at glass/air re-enters the stack (beta < 1 always,
                             beta >= 1 only in 'large') and joins the recycled flux
Common platform: 800-nm SiO2 ODR, 15-nm low-k ITO cathode, organic k 1e-4, PLQY 1, horizontal dipoles 0.8,
SL index <= 2.0.  DBR k floor: evaporated ZnS/LiF 1e-4, sputtered TiO2/SiO2 2e-5.
Optimisation: differential evolution over thicknesses, SL index and DBR pair count (2..12), with the initial
population seeded by every finished design of the same arch/anode/emitter (and the validated compare_anodes TE
designs); then bounded Nelder-Mead at N-1, N, N+1.  Objective on a 1-nm (narrow emitter) or 3-nm (Ir(ppy)2acac)
grid; final evaluation on a 0.5-nm grid (dipole emission in wavelength chunks to bound memory).

usage: python3 arch.py run <case> [<case> ...]     (case names: see CASES; writes results/arch/<case>.json)
       python3 arch.py all                         (every case not yet done)
       python3 arch.py refinal <case ...|all>      (re-evaluate stored designs with the current model)
       python3 arch.py extend <case ...>           (continue the fixed-N polish while the best N is at the window edge)
"""
from __future__ import annotations

import glob
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.optimize import differential_evolution, minimize
from scipy.stats import qmc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import optimize as opt  # noqa: E402
import recycling as rc  # noqa: E402
from oledsim import DipoleResult, Layer, Medium, dipole_emission, stack_RTA  # noqa: E402
from devices import Ctx, _dbr_layers, photon_weights  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "arch")
OLD = os.path.join(HERE, "results", "arch_v1")          # first (unseeded, N <= 8) search, used as seeds
N_GLASS = 1.52
K_EVAP, K_SPUT = 1e-4, 2e-5
NMAX = 12.49
DBR_KEYS = {"ev": ("H", "L"), "ox": ("bH", "bL")}
EMITTERS = {"irppy": "irppy2acac", "n15": "narrow532_15"}
ANODES = {"Ag8m": "Ag8meas", "IZO": "IZOtop"}


def make(emitter, anode, final=False, k_dbr=None):
    pm, pf, ag_b = opt.preset_mats(f"{ANODES[anode]}+lowkITO+korg1e-4")
    em = EMITTERS[emitter]
    if final:
        ctx = Ctx(opt.DB, wl=np.arange(400.0, 700.01, 0.5), mats=pm)
        s0 = ctx.emitter(**opt.EMITTERS[em])
    else:
        step = 1.0 if emitter == "n15" else 3.0       # 3 nm: within 0.2 % of 2 nm for Ir(ppy)2acac
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


def _fresnel(n1, n2, beta, pol):
    """power reflectance n1 -> n2 (both real, n2 > n1) for in-plane wavevector beta < n1"""
    c1 = np.sqrt(np.clip(1 - (beta / n1) ** 2, 0, 1)); c2 = np.sqrt(np.clip(1 - (beta / n2) ** 2, 0, 1))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (n1 * c1 - n2 * c2) / (n1 * c1 + n2 * c2) if pol == "s" else (n2 * c1 - n1 * c2) / (n2 * c1 + n1 * c2)
    return np.where(beta < n1, np.nan_to_num(r * r, nan=1.0), 1.0)


def exit_index(exit_, P):
    return 1.0 if exit_ == "air" else (P["n_low"] if exit_ == "lowNi" else N_GLASS)


def exit_props(beta, pol, exit_, n_low, mode):
    """incoherent region beyond the coherent stack, for light leaving the stack into the exit medium at beta:
    Rx = fraction that comes back to the stack (multiple incoherent bounces; in 'pixel' mode light that enters the
    glass at beta >= 1 is lost), Tx = fraction that reaches air."""
    Rga = _r_glass_air(beta, pol)
    if exit_ == "lowNi":
        R1 = _fresnel(n_low, N_GLASS, beta, pol); T1 = 1 - R1
        d = 1 - R1 * Rga
        Rx = R1 + _div(T1 * T1 * Rga, d); Tx = _div(T1 * (1 - Rga), d); Rpix = R1
    else:
        Rx = Rga; Tx = 1 - Rga; Rpix = 0.0
    if mode == "pixel":
        Rx = np.where(beta >= 1.0, Rpix, Rx)
    return Rx, Tx


ODR = [["l_SiO2_SJ_RTP", 800.0]]
_QH, _QL = ["l_TiO2_SJ_RTP", 530.0 / (4 * 2.35)], ["l_SiO2_SJ_RTP", 530.0 / (4 * 1.45)]
MIRRORS = {"odr": ODR,                                              # 800-nm SiO2 ODR on Ag (default)
           "edr": ODR + [_QH, _QL] * 2,                             # dense: SiO2 800 + 2 quarter-wave TiO2/SiO2 pairs on Ag
           "epr": [[1.2, 1000.0]] + [_QH, _QL] * 2}                 # porous n 1.2 spacer (sealed) + 2 pairs on Ag


def mirror_layers(ctx, P):
    """bottom reflector below the SL (from the SL towards the Ag mirror); default 800-nm SiO2 ODR.
    P['bot'] = [[material key | real index, thickness nm], ...]"""
    out = []
    for m, d in P.get("bot") or ODR:
        out.append(ctx.layer(complex(m) if isinstance(m, (int, float)) else m, float(d), name=f"bot_{m}"))
    return [l for l in out if l.d > 0]


def _div(a, b):
    a, b = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float))
    return np.divide(a, b, out=np.zeros(a.shape), where=b > 0)


def stacks(ctx, P, arch, exit_, dbr=None):
    hk, lk = DBR_KEYS[dbr or ("ev" if arch == "TE" else "ox")]
    cath = [ctx.layer("etl", P["d_etl"]), ctx.layer("eil", 10.0), ctx.layer("tco", 15.0)]
    anod = [ctx.layer("tcta", 10.0), ctx.layer("htl", P["d_htl"]), ctx.layer("hatcn", 10.0), ctx.layer("ag_top", P["d_ag"])]
    D = _dbr_layers(ctx, P["dbr"], hk=hk, lk=lk)
    if arch == "TE":
        up = anod + D; dn = cath
    else:
        up = cath + D; dn = anod
    if exit_ == "lowN":
        up = up + [ctx.layer(complex(P["n_low"]), P["d_low"], name="lowN")]
    up = [l for l in up if l.d > 0]; dn = [l for l in dn if l.d > 0]
    return up, dn


def _dipole(wl, n_eml, up, top, dn, bot, du, chunk=None):
    """dipole emission, optionally in wavelength chunks (same result, bounded memory)."""
    kw = dict(theta_h=0.8, u_max=3.0, du=du)
    if not chunk or len(wl) <= chunk:
        return dipole_emission(wl, n_eml, 10.0, 10.0, up, top, dn, bot, **kw)
    parts = []
    for i in range(0, len(wl), chunk):
        s = slice(i, i + chunk)
        cut = lambda ls: [Layer(l.no[s], l.d, l.ne[s], l.name) for l in ls]  # noqa: E731
        parts.append(dipole_emission(wl[s], n_eml[s], 10.0, 10.0, cut(up), Medium(top.no[s], top.ne[s]),
                                     cut(dn), Medium(bot.no[s], bot.ne[s]), **kw))
    res = DipoleResult()
    for k, v in parts[0].items():
        res[k] = v if k in ("u", "du") else np.concatenate([p[k] for p in parts], axis=0)
    return res


def _glass_first_pass(res, nm, beta, edges, gt, mode, exit_="glass", n_low=None):
    """fate of the first-pass emission into the exit medium (incoherent region beyond the stack).
    Returns Pin_x (2, nwl, nth): power that re-enters the SL after coming back from the exit side, by SL angle;
            esc (nwl,): first-pass power that reaches air, including exit/stack multiple bounces;
            f0 (nwl,): on-axis factor Tx(0) / (1 - Rs2(0) Rx(0)) for the direct radiance."""
    u = res["u"]; du = res["du"]; ne = np.real(res["n_eml"]); nmr = np.real(nm)
    n_x = n_low if exit_ == "lowNi" else N_GLASS
    nwl, nth = len(ne), len(edges) - 1
    Pin_x = np.zeros((2, nwl, nth)); esc = np.zeros(nwl)
    sin_e = np.sin(np.deg2rad(edges)); uu = np.r_[0.0, u + 0.5 * du]
    for ip, pol in enumerate(("s", "p")):
        Rs2, Ts2, L2 = gt[pol]
        Pt = res[f"Pt_{pol}"] * (2 * u * du)[None, :]
        for i in range(nwl):
            b = ne[i] * u
            Rx, Tx = exit_props(b, pol, exit_, n_low, mode)
            ts2 = np.interp(b, beta[i], Ts2[i]); l2 = np.interp(b, beta[i], L2[i])
            den = (1.0 - Rx) + Rx * l2                                 # 1 - Rs2 Rx
            P = np.where(b < n_x, Pt[i], 0.0)
            esc[i] += np.sum(_div(P * Tx, den))
            cum = np.r_[0.0, np.cumsum(_div(P * Rx * ts2, den))]
            Pin_x[ip, i] = np.diff(np.interp((nmr[i] / ne[i]) * sin_e, uu, cum))
    Rx0, Tx0 = exit_props(np.array([0.0]), "s", exit_, n_low, mode)
    rs2_0 = 0.5 * (gt["s"][0][:, 0] + gt["p"][0][:, 0])
    return Pin_x, esc, Tx0[0] / (1.0 - rs2_0 * Rx0[0])


def run(ctx, s0, P, arch, exit_, mode="pixel", n_th=121, du=1e-3, dbr=None, chunk=None):
    nm = P["n_m"] * ctx.ones
    up, dn = stacks(ctx, P, arch, exit_, dbr)
    eml = ctx.layer("eml", 20.0)
    n_x = exit_index(exit_, P)
    ext = Medium(n_x * ctx.ones + 0j)
    res = _dipole(ctx.wl, eml.no, up, ext, dn, Medium(nm + 0j), du, chunk)
    qw, qeff, Fbar = photon_weights(res, s0, 1.0)
    th, wL, edges = rc.angle_grid(n_th)
    stack = list(reversed(dn)) + [eml] + up
    beta = np.real(nm)[:, None] * np.sin(np.deg2rad(th))[None, :]
    sl = Medium(nm + 0j)
    bot = mirror_layers(ctx, P)
    mirror = Medium(ctx.nk(ctx.mats["mirror"]))
    T = {k: [] for k in ("Rt", "Tt", "At", "Rb")}
    gt = {}
    for pol in ("s", "p"):
        Rs, Ts, _ = stack_RTA(stack, sl, ext, beta, ctx.wl, pol)
        if exit_ == "air":
            Rt, Tt = Rs, Ts
        else:
            ok = beta < n_x - 1e-9
            Rs2, Ts2, As2 = stack_RTA(list(reversed(stack)), ext, sl, np.where(ok, beta, 0.0), ctx.wl, pol)
            Ts = np.where(ok, Ts, 0.0); Rs2 = np.where(ok, Rs2, 0.0); Ts2 = np.where(ok, Ts2, 0.0)
            L2 = np.where(ok, np.clip(Ts2 + As2.sum(0), 0.0, 1.0), 1.0)    # 1 - Rs2 without cancellation
            Rx, Tx = exit_props(beta, pol, exit_, P.get("n_low"), mode)
            den = (1.0 - Rx) + Rx * L2                                        # 1 - Rs2 Rx
            Tt = _div(Ts * Tx, den)
            Rt = Rs + _div(Ts * Ts * Rx, den)
            gt[pol] = (Rs2, Ts2, L2)
        At = np.clip(1 - Rt - Tt, 0, 1)
        Rb, _, _ = stack_RTA(bot, sl, mirror, beta, ctx.wl, pol)
        for k, v in zip(("Rt", "Tt", "At", "Rb"), (Rt, Tt, At, Rb)):
            T[k].append(v)
    tables = {k: np.array(v) for k, v in T.items()}
    Pin = rc.emission_into_sl(res, nm, th, edges)
    if exit_ != "air":
        Pin_x, esc, f0 = _glass_first_pass(res, nm, beta, edges, gt, mode, exit_, P.get("n_low"))
        Pin = Pin + Pin_x
    out = rc.analytic(res, tables, nm, ctx.wl, th, wL, Pin, qw)
    if exit_ != "air":
        out["L_dir0"] = out["L_dir0"] * f0
        out["eqe_dir"] = qw * esc
    s = rc.summarize(out, ctx.wl, ctx.V)
    s["Aprime"] = 1 - (1 - s["a_top_avg"] / (1 - s["p_avg"])) * s["rb_avg"]
    return s


# ----------------------------------------------------------------------------------------------
def params(x, N, exit_, n_low, mir=None):
    P = dict(d_etl=x[0], d_htl=x[1], d_ag=x[2], dbr=dict(N=N, dH1=x[3], dH=x[4], dL=x[5]), n_m=x[6])
    if mir and mir != "odr":
        P["bot"] = MIRRORS[mir]
    if exit_ == "lowN":
        P.update(n_low=n_low, d_low=x[8])
    elif exit_ == "lowNi":
        P.update(n_low=n_low)
    return P


CASES = {}
for em in EMITTERS:
    for an in ANODES:
        CASES[f"TE_air_{an}_{em}"] = dict(arch="TE", exit="air", mode="pixel", anode=an, emitter=em)
        CASES[f"TE_air-ox_{an}_{em}"] = dict(arch="TE", exit="air", mode="pixel", anode=an, emitter=em, dbr="ox")
        CASES[f"TE_glass_{an}_{em}"] = dict(arch="TE", exit="glass", mode="pixel", anode=an, emitter=em)
        CASES[f"BE_glass-pixel_{an}_{em}"] = dict(arch="BE", exit="glass", mode="pixel", anode=an, emitter=em)
        CASES[f"BE_glass-large_{an}_{em}"] = dict(arch="BE", exit="glass", mode="large", anode=an, emitter=em)
        for nl in (1.2, 1.1):
            CASES[f"BE_lowN{nl}_{an}_{em}"] = dict(arch="BE", exit="lowN", mode="pixel", anode=an, emitter=em, n_low=nl)
            CASES[f"TE_lowN{nl}_{an}_{em}"] = dict(arch="TE", exit="lowN", mode="pixel", anode=an, emitter=em, n_low=nl)
        for mir in ("edr", "epr"):
            CASES[f"TE_air-ox-{mir}_{an}_{em}"] = dict(arch="TE", exit="air", mode="pixel", anode=an, emitter=em, dbr="ox", mir=mir)
            CASES[f"TE_air-{mir}_{an}_{em}"] = dict(arch="TE", exit="air", mode="pixel", anode=an, emitter=em, mir=mir)
            CASES[f"BE_lowN1.05-inc-{mir}_{an}_{em}"] = dict(arch="BE", exit="lowNi", mode="pixel", anode=an, emitter=em,
                                                            n_low=1.05, mir=mir)
        CASES[f"BE_lowN1.05_{an}_{em}"] = dict(arch="BE", exit="lowN", mode="pixel", anode=an, emitter=em, n_low=1.05)
        for nl in (1.2, 1.1, 1.05):
            CASES[f"BE_lowN{nl}-inc_{an}_{em}"] = dict(arch="BE", exit="lowNi", mode="pixel", anode=an, emitter=em, n_low=nl)
            CASES[f"TE_lowN{nl}-inc_{an}_{em}"] = dict(arch="TE", exit="lowNi", mode="pixel", anode=an, emitter=em, n_low=nl)
        CASES[f"BE_lowN1.1-ev_{an}_{em}"] = dict(arch="BE", exit="lowN", mode="pixel", anode=an, emitter=em, n_low=1.1,
                                                 dbr="ev")

G = {}


def _obj(x):
    c = G["case"]
    P = params(x, int(round(x[7])), c["exit"], c.get("n_low"), c.get("mir"))
    return -run(G["ctx"], G["s0"], P, c["arch"], c["exit"], c["mode"], n_th=91, dbr=c.get("dbr"))["CE"]


def bounds(c, ag_b):
    b = [(10, 500), (20, 600), ag_b, (0, 220), (20, 220), (20, 300), (1.6, 2.0), (2, NMAX)]
    if c["exit"] == "lowN":
        b.append((300, 2000))
    return b


def _vec(d, c):
    """parameter vector of a stored design, mapped onto case c"""
    P = d["params"]; D = P["dbr"]
    x = [P["d_etl"], P["d_htl"], P["d_ag"], D["dH1"], D["dH"], D["dL"], P["n_m"], D["N"]]
    if c["exit"] == "lowN":
        x.append(P.get("d_low") or 1000.0)
    return x


def seeds(name, c):
    """finished designs with the same arch, anode and emitter (any exit or DBR), plus the compare_anodes TE designs"""
    xs = []
    for p in sorted(glob.glob(os.path.join(OUT, "*.json")) + glob.glob(os.path.join(OLD, "*.json"))):
        try:
            d = json.load(open(p))
            if (d["arch"], d["anode"], d["emitter"]) == (c["arch"], c["anode"], c["emitter"]):
                xs.append(_vec(d, c))
        except (KeyError, ValueError):
            pass
    if c["arch"] == "TE":
        for p in sorted(glob.glob(os.path.join(HERE, "results", f"opt_cmp_{c['emitter']}_{c['anode']}_N*.json"))):
            try:
                xs.append(_vec(json.load(open(p)), c))
            except (KeyError, ValueError):
                pass
    return xs


def init_pop(b, xs, n, seed=3):
    lo = np.array([v[0] for v in b], float); hi = np.array([v[1] for v in b], float)
    pop = lo + qmc.LatinHypercube(d=len(b), seed=seed).random(n) * (hi - lo)
    for k, x in enumerate(xs[: n // 3]):
        pop[k] = np.clip(np.array(x, float), lo, hi)
    return pop


def _polish(x0, N, b, maxfev=300):
    """bounded Nelder-Mead of the continuous variables at a fixed pair count N"""
    idx = [i for i in range(len(b)) if i != 7 and b[i][1] - b[i][0] > 0.1]
    lo = np.array([b[i][0] for i in idx]); hi = np.array([b[i][1] for i in idx])
    xb = np.array(x0, float); xb[7] = N

    def f(z):
        x = xb.copy(); x[idx] = z
        return _obj(x)
    z0 = xb[idx]
    sim = [z0] + [np.where(np.arange(len(z0)) == k, np.clip(z0 + 0.06 * np.maximum(abs(z0), 0.05 * (hi - lo)), lo, hi), z0)
                  for k in range(len(z0))]
    m = minimize(f, z0, method="Nelder-Mead", bounds=list(zip(lo, hi)),
                 options=dict(maxfev=maxfev, xatol=0.3, fatol=0.05, initial_simplex=np.array(sim), adaptive=True))
    x = xb.copy(); x[idx] = m.x
    f0 = f(z0)
    return (m.fun, x) if m.fun < f0 else (f0, xb)


def final(c, P):
    ctxf, s0f, _ = make(c["emitter"], c["anode"], final=True)
    kw = dict(n_th=181, du=2.5e-4, dbr=c.get("dbr"), chunk=40)
    fin = run(ctxf, s0f, P, c["arch"], c["exit"], c["mode"], **kw)
    other = None
    if c["exit"] != "air":
        other = run(ctxf, s0f, P, c["arch"], c["exit"], "large" if c["mode"] == "pixel" else "pixel", **kw)
    return fin, other


def optimise(name, popsize=10, maxiter=25):
    path = os.path.join(OUT, f"{name}.json")
    if os.path.exists(path):
        return
    c = CASES[name]; t0 = time.time()
    ctx, s0, ag_b = make(c["emitter"], c["anode"])
    G.update(ctx=ctx, s0=s0, case=c)
    b = bounds(c, ag_b)
    xs = seeds(name, c)
    W = int(os.environ.get("ARCH_WORKERS", "4"))
    r = differential_evolution(_obj, b, maxiter=maxiter, init=init_pop(b, xs, popsize * len(b)), seed=3, tol=1e-6,
                               polish=False, workers=W, updating="deferred")
    N0 = int(round(r.x[7]))
    Ns = [n for n in (N0 - 1, N0, N0 + 1) if 2 <= n <= int(NMAX)]
    with Pool(len(Ns)) as pool:
        pol = pool.starmap(_polish, [(r.x, n, b) for n in Ns])
    fun, x = min([(r.fun, np.array(r.x, float))] + list(pol), key=lambda t: t[0])
    N = int(round(x[7]))
    P = params(x, N, c["exit"], c.get("n_low"), c.get("mir"))
    fin, other = final(c, P)
    out = dict(case=name, **c, N=N, x=list(map(float, x)), params=P, final=fin, other_glass_mode=other,
               objective_grid_CE=float(-fun), n_seeds=len(xs), de_CE=float(-r.fun), de_N=N0,
               polish_CE={int(round(xp[7])): float(-fp) for fp, xp in pol}, seconds=time.time() - t0)
    os.makedirs(OUT, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"{name:32s} N {N}  CE {fin['CE']:7.1f}  EQE {fin['eqe']:.3f}  A' {100 * fin['Aprime']:.2f}%  "
          f"p {fin['p_avg']:.3f}  (DE {-r.fun:.0f} @N{N0}, polish {max(-fp for fp, _ in pol):.0f}, "
          f"{len(xs)} seeds, {time.time() - t0:.0f}s)", flush=True)


def extend(name):
    """continue the fixed-N polish beyond the N-1..N+1 window while the best N sits at its edge"""
    path = os.path.join(OUT, f"{name}.json")
    d = json.load(open(path)); c = CASES[name]
    pol = {int(k): v for k, v in (d.get("polish_CE") or {}).items()}
    if not pol:
        return
    ctx, s0, ag_b = make(c["emitter"], c["anode"])
    G.update(ctx=ctx, s0=s0, case=c)
    b = bounds(c, ag_b)
    x = np.array(d["x"], float); N = int(d["N"]); best = d["objective_grid_CE"]; changed = False
    for step in (+1, -1):
        if not ((step > 0 and N == max(pol)) or (step < 0 and N == min(pol))):
            continue
        while 2 <= N + step <= int(NMAX):
            fun, xn = _polish(x, N + step, b)
            pol[N + step] = -fun
            if -fun <= best + 0.5:
                break
            best, x, N, changed = -fun, xn, N + step, True
    d["polish_CE"] = pol
    if changed:
        P = params(x, N, c["exit"], c.get("n_low"), c.get("mir"))
        fin, other = final(c, P)
        d["pre_extend"] = dict(N=d["N"], x=d["x"], final_CE=d["final"]["CE"], objective_grid_CE=d["objective_grid_CE"])
        d.update(N=N, x=list(map(float, x)), params=P, final=fin, other_glass_mode=other, objective_grid_CE=float(best))
    with open(path, "w") as fh:
        json.dump(d, fh, indent=1, default=float)
    print(f"{name:32s} extend: N {d.get('pre_extend', {}).get('N', N)} -> {N}  CE {d['final']['CE']:7.1f}  "
          f"polish {dict(sorted((k, round(v)) for k, v in pol.items()))}", flush=True)


def refinal(name):
    """re-evaluate a stored design with the current model (keeps the previous evaluation as final_v1)"""
    path = os.path.join(OUT, f"{name}.json")
    d = json.load(open(path))
    c = CASES[name]
    fin, other = final(c, d["params"])
    d.setdefault("final_v1", d["final"]); d.setdefault("other_glass_mode_v1", d["other_glass_mode"])
    d["final"], d["other_glass_mode"] = fin, other
    with open(path, "w") as fh:
        json.dump(d, fh, indent=1, default=float)
    print(f"{name:32s} refinal CE {d['final_v1']['CE']:7.1f} -> {fin['CE']:7.1f}  EQE {fin['eqe']:.3f}  "
          f"A' {100 * fin['Aprime']:.2f}%  other {None if other is None else round(other['CE'], 1)}", flush=True)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "extend":
        for n in sys.argv[2:]:
            extend(n)
    elif cmd == "refinal":
        names = sys.argv[2:]
        if names == ["all"]:
            names = [os.path.basename(p)[:-5] for p in sorted(glob.glob(os.path.join(OUT, "*.json")))]
        for n in names:
            refinal(n)
    else:
        for n in (list(CASES) if cmd == "all" else sys.argv[2:]):
            optimise(n)
