"""
arch_budget.py -- where the recycled light goes in one round trip, for the optimised architectures (arch.py).

A Lambertian upward flux in the SL (spectrally weighted like the light that enters the SL, P_dn) meets the top
stack.  Per round trip it is split into
    escape into air (useful, = p)  |  absorption in the DBR, anode, ITO cathode, organics (+EIL), low-n layer
    |  glass-trapped (pixel mode: beta >= 1 light that walks out of the pixel)  |  mirror loss (1 - Rb of the
    reflected part)  |  recirculated (reflected and returned by the mirror).
Glass exits are incoherent: light reflected at glass/air re-enters the stack (multiple bounces), as in arch.run.
Also reported: on-axis transmission T0 (photopically weighted) and n_m (CE ~ T0 / (n_m^2 (p + losses))).

usage: python3 arch_budget.py [case ...]     (default: every results/arch/*.json; writes results/arch_budget.json)
       python3 arch_budget.py plot          (stacked bars, results/figs/arch_budget_IZO_<emitter>.*)
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import arch as A  # noqa: E402
import recycling as rc  # noqa: E402
from oledsim import Medium, stack_RTA  # noqa: E402
from devices import photon_weights  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
GROUPS = ("DBR", "anode", "ITO", "organics", "lowN")


def group(name):
    if name == "ag_top":
        return "anode"
    if name == "tco":
        return "ITO"
    if name == "lowN":
        return "lowN"
    if name[:1] in ("H", "L") or name[:2] in ("bH", "bL"):
        return "DBR"
    return "organics"


def budget(name, n_th=121):
    d = json.load(open(os.path.join(A.OUT, f"{name}.json")))
    c = A.CASES[name]; P = d["params"]
    ctx, s0, _ = A.make(c["emitter"], c["anode"], final=True)
    nm = P["n_m"] * ctx.ones
    up, dn = A.stacks(ctx, P, c["arch"], c["exit"], c.get("dbr"))
    eml = ctx.layer("eml", 20.0)
    ext = Medium(ctx.ones + 0j) if c["exit"] == "air" else Medium(A.N_GLASS * ctx.ones + 0j)
    res = A._dipole(ctx.wl, eml.no, up, ext, dn, Medium(nm + 0j), 1e-3, chunk=60)   # spectral weights only
    qw, _, _ = photon_weights(res, s0, 1.0)
    th, wL, edges = rc.angle_grid(n_th)
    Pin = rc.emission_into_sl(res, nm, th, edges)
    stack = list(reversed(dn)) + [eml] + up
    gidx = {g: [i for i, l in enumerate(stack) if group(l.name) == g] for g in GROUPS}
    beta = np.real(nm)[:, None] * np.sin(np.deg2rad(th))[None, :]
    sl = Medium(nm + 0j)
    mirror = Medium(ctx.nk(ctx.mats["mirror"]))
    w = 0.5 * wL[None, :]
    acc = {k: np.zeros(len(ctx.wl)) for k in ("escape", "trap", "mirror", "recirc") + GROUPS}
    T0 = np.zeros(len(ctx.wl)); gt = {}
    for pol in ("s", "p"):
        R, T, Al = stack_RTA(stack, sl, ext, beta, ctx.wl, pol)
        if c["exit"] == "air":
            esc, back, trap, Aeff = T, R, 0 * T, Al
        else:
            ok = beta < A.N_GLASS - 1e-9
            R2, T2, A2 = stack_RTA(list(reversed(stack)), ext, sl, np.where(ok, beta, 0.0), ctx.wl, pol)
            T = np.where(ok, T, 0.0); R2 = np.where(ok, R2, 0.0); T2 = np.where(ok, T2, 0.0)
            A2 = np.where(ok[None], A2[::-1], 0.0)                     # back to the stack's layer order
            L2 = np.where(ok, np.clip(T2 + A2.sum(0), 0, 1), 1.0)
            Rga = A._r_glass_air(beta, pol)
            if c["mode"] == "pixel":
                ret = np.where(beta < 1.0, Rga, 0.0); trap = np.where(beta >= 1.0, T, 0.0)
            else:
                ret = Rga; trap = 0 * T
            den = (1 - Rga) + Rga * L2
            esc = A._div(T * (1 - Rga), den)
            back = R + A._div(T * ret * T2, den)
            Aeff = Al + A._div(T * ret, den)[None] * A2
            gt[pol] = (R2, T2, L2)
        Rb, _, _ = stack_RTA([ctx.layer("l_SiO2_SJ_RTP", 800.0)], sl, mirror, beta, ctx.wl, pol)
        acc["escape"] += np.sum(esc * w, 1); acc["trap"] += np.sum(trap * w, 1)
        acc["mirror"] += np.sum(back * (1 - Rb) * w, 1); acc["recirc"] += np.sum(back * Rb * w, 1)
        for g in GROUPS:
            if gidx[g]:
                acc[g] += np.sum(Aeff[gidx[g]].sum(0) * w, 1)
        T0 += 0.5 * esc[:, 0]
    if c["exit"] != "air":
        Pin = Pin + A._glass_first_pass(res, nm, beta, edges, gt, c["mode"])[0]
    wts = qw * Pin.sum(axis=(0, 2)); wts = wts / wts.sum()
    out = {k: float(np.sum(wts * v)) for k, v in acc.items()}
    out["closure"] = float(sum(out[k] for k in ("escape", "trap", "mirror", "recirc") + GROUPS))
    vw = wts * ctx.V; out["T0_phot"] = float(np.sum(vw * T0) / np.sum(vw))
    out.update(n_m=float(P["n_m"]), CE=float(d["final"]["CE"]), Aprime=float(d["final"]["Aprime"]),
               p_avg=float(d["final"]["p_avg"]), N=int(d["N"]))
    return out


SEG = (("escape", "Escape (useful)", "#2e7d32"), ("DBR", "DBR absorption", "#1f5fa8"),
       ("anode", "Anode absorption", "#7f8c8d"), ("ITO", "ITO absorption", "#b9770e"),
       ("organics", "Organic/EIL absorption", "#8e44ad"), ("lowN", "Low-n layer", "#76d7c4"),
       ("trap", "Trapped in glass", "#c0392b"), ("mirror", "Mirror (ODR/Ag) loss", "#f1c40f"))
SHORT = {"TE_air": "TE, air", "TE_air-ox": "TE, air, oxide DBR", "TE_glass": "TE + OCA glass",
         "TE_lowN1.2": "TE + n1.2 + glass", "BE_glass-pixel": "BE glass (pixel)", "BE_glass-large": "BE glass (large)",
         "BE_lowN1.2": "BE + n1.2", "BE_lowN1.1": "BE + n1.1", "BE_lowN1.05": "BE + n1.05", "BE_lowN1.1-ev": "BE + n1.1, ZnS/LiF DBR"}


def plot(allb, emitter="n15", anode="IZO"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import plot_phase_origin  # noqa: F401  (Origin-style rcParams)
    rows = [(k.rsplit("_", 2)[0], v) for k, v in allb.items() if k.endswith(f"_{anode}_{emitter}")]
    rows.sort(key=lambda r: r[1]["CE"])
    fig, ax = plt.subplots(figsize=(12, 0.62 * len(rows) + 2.2))
    y = np.arange(len(rows))
    left = np.zeros(len(rows))
    for key, lab, col in SEG:
        v = np.array([100 * r[1][key] for r in rows])
        ax.barh(y, v, left=left, color=col, edgecolor="k", lw=0.8, label=lab, height=0.68)
        left += v
    for yi, (k, b), x in zip(y, rows, left):
        ax.text(x + 0.3, yi, f"{b['CE']:.0f} cd/A", va="center", fontsize=14)
    ax.set_yticks(y); ax.set_yticklabels([SHORT.get(k, k) for k, _ in rows], fontsize=15)
    ax.set_xlabel("Share of the upward flux per round trip (%)")
    ax.set_xlim(0, max(left) * 1.28)
    ax.legend(loc="lower right", frameon=False, fontsize=12, ncol=2)
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(HERE, "results", "figs", f"arch_budget_{anode}_{emitter}.{ext}"),
                    dpi=300 if ext == "png" else None)
    plt.close(fig)


if __name__ == "__main__":
    if sys.argv[1:2] == ["plot"]:
        allb = json.load(open(os.path.join(HERE, "results", "arch_budget.json")))
        for em in ("n15", "irppy"):
            plot(allb, em, "IZO")
        sys.exit()
    names = sys.argv[1:] or [os.path.basename(p)[:-5] for p in sorted(glob.glob(os.path.join(A.OUT, "*.json")))]
    path = os.path.join(HERE, "results", "arch_budget.json")
    allb = json.load(open(path)) if os.path.exists(path) else {}
    for n in names:
        b = budget(n); allb[n] = b
        print(f"{n:30s} CE {b['CE']:7.1f} | esc {100 * b['escape']:.2f}  DBR {100 * b['DBR']:.2f}  anode {100 * b['anode']:.2f}"
              f"  ITO {100 * b['ITO']:.2f}  org {100 * b['organics']:.2f}  lowN {100 * b['lowN']:.3f}  trap {100 * b['trap']:.2f}"
              f"  mirror {100 * b['mirror']:.2f}  recirc {100 * b['recirc']:.1f}  (sum {b['closure']:.4f})  T0 {b['T0_phot']:.3f}",
              flush=True)
        json.dump(allb, open(path, "w"), indent=1)
