"""
amap.py -- generalised map of the maximum on-axis current efficiency versus round-trip loss A' and
emitter linewidth (FWHM), for
  (R) the angle-selective recycling TEOLED, and
  (S) the strong-microcavity reference TEOLED (measured-Ag and ideal-Ag top anode).

Recycling: the metal-free platform (sputtered IZO top anode on HATCN, ZnS/LiF DBR, 800-nm SiO2 ODR,
15-nm low-k bottom ITO, organic k 1e-4) has the lowest A' (~3 %).  A' is raised continuously by an
angle-independent extra loss per SL round trip (a_extra = a_sl_pass), and for every (FWHM, a_extra, N)
the stack is re-optimised (bounded Nelder-Mead, continuation in a_extra).  The resulting effective A'
(1 - (1 - a_top/(1 - p)) * R_b * (1 - a_extra), spectrally weighted) is the abscissa, so the map does not
depend on what causes the loss.  Measured/ideal thin-Ag and IZO devices (compare_anodes.py) are overlaid
as a check of this generalisation.

Emitter family: Gaussian in photon energy, peak 530 nm, no vibronic shoulder, FWHM 10-80 nm.

usage: python3 amap.py run      (optimisation; writes results/amap/*.json, resumable)
       python3 amap.py plot     (figures in results/figs/amap_*.png)
"""
from __future__ import annotations

import glob
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import optimize as opt  # noqa: E402
import refine as rf  # noqa: E402
from devices import Ctx, run_proposed, run_reference  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "amap")
FIG = os.path.join(HERE, "results", "figs")
PEAK = 530.0
FWHMS = [10, 15, 20, 25, 35, 50, 65, 80]
AEXTRA = [0.0, 0.015, 0.03, 0.05, 0.08, 0.12, 0.18]
NS = [3, 5, 8]
PRESET_R = "IZOtop+ODR+lowkITO+korg1e-4"
REF_PRESETS = {"realAg": "korg1e-4", "idealAg": "bulkAg+korg1e-4"}
SEEDS_R = ["prop6_narrow532_15_IZOtop_lowk_korg_tco15", "cmp_n15_IZO_N8", "cmp_irppy_IZO_N3",
           "cmp_irppy_IZO_N6p", "cmp_irppy_IZO_N8"]
SEEDS_S = ["ref_narrow532_15_korg", "ref_irppy2acac_korg", "ref_narrow532_15_bulkAg_korg_s1",
           "ref_irppy2acac_bulkAg_korg_s2"]


def emitter_key(fwhm):
    key = f"g{int(PEAK)}_{fwhm}"
    opt.EMITTERS[key] = dict(kind="narrow", peak=PEAK, fwhm=float(fwhm), shoulder=0.0)
    return key


def _x(tag):
    with open(os.path.join(HERE, "results", f"opt_{tag}.json")) as fh:
        return np.array(json.load(fh)["x"], float)


def _nm(f, x0, b, maxfev):
    lo = np.array([l for l, h in b]); hi = np.array([h for l, h in b])
    x0 = np.clip(x0, lo, hi)
    sim = [x0]
    for i in range(len(x0)):
        xi = x0.copy(); st = 0.08 * max(abs(x0[i]), 0.05 * (hi[i] - lo[i]))
        xi[i] = x0[i] + st if x0[i] + st <= hi[i] else x0[i] - st
        sim.append(xi)
    f0 = f(x0)
    r = minimize(f, x0, method="Nelder-Mead", bounds=b,
                 options=dict(maxfev=maxfev, xatol=0.3, fatol=0.05, initial_simplex=np.array(sim), adaptive=True))
    return (r.x, r.fun) if r.fun < f0 else (x0, f0)


def chain_recycling(task):
    fwhm, N = task
    path = os.path.join(OUT, f"prop_f{fwhm}_N{N}.json")
    if os.path.exists(path):
        return path
    em = emitter_key(fwhm)
    pm, pf, ag_b = opt.preset_mats(PRESET_R)
    fixed = {**pf, "d_tco": 15.0}
    ctx, s0 = opt.make_ctx(em, mats=pm, step=min(4.0, fwhm / 5.0))
    ctxf, s0f = opt.make_ctx(em, mats=pm, full=True)
    b = rf.bounds_for("prop", ag_b)
    rows = []
    x = None
    for a in AEXTRA:
        def f(xx):
            P = {**fixed, **rf.params_for("prop", xx, N, 2)}
            return -run_proposed(ctx, P, s0, n_th=91, du=1e-3, a_sl_pass=a)["summary"]["CE"]
        cands = [_x(t) for t in SEEDS_R] + ([x] if x is not None else [])
        x0 = min(cands, key=f)
        x, fun = _nm(f, x0, b, 100)
        P = {**fixed, **rf.params_for("prop", x, N, 2)}
        s = run_proposed(ctxf, P, s0f, du=2.5e-4, a_sl_pass=a)["summary"]
        Ap = 1 - (1 - s["a_top_avg"] / (1 - s["p_avg"])) * s["rb_avg"] * (1 - a)
        rows.append(dict(a_extra=a, x=list(map(float, x)), CE=s["CE"], eqe=s["eqe"], p=s["p_avg"],
                         a_top=s["a_top_avg"], rb=s["rb_avg"], Aprime=Ap))
        print(f"R fwhm {fwhm} N {N} a+ {a:.3f}: A' {100 * Ap:.2f}%  CE {s['CE']:.0f}  EQE {s['eqe']:.3f}", flush=True)
    with open(path, "w") as fh:
        json.dump(dict(kind="recycling", fwhm=fwhm, N=N, preset=PRESET_R, rows=rows), fh, indent=1)
    return path


def chain_reference(task):
    fwhm, mode = task
    path = os.path.join(OUT, f"ref_f{fwhm}_{mode}.json")
    if os.path.exists(path):
        return path
    em = emitter_key(fwhm)
    pm, pf, _ = opt.preset_mats(REF_PRESETS[mode])
    ctx, s0 = opt.make_ctx(em, mats=pm, step=1.0 if fwhm <= 25 else 2.0)
    b = rf.bounds_for("ref", None, d_ag_min=10.0)

    def f(xx):
        return -run_reference(ctx, {**pf, **rf.params_for("ref", xx, 3, 2)}, s0, du=1e-3)["summary"]["CE"]
    x0 = min([_x(t) for t in SEEDS_S], key=f)
    x, fun = _nm(f, x0, b, 150)
    c0 = Ctx(opt.DB, wl=np.arange(400.0, 700.01, 0.5 if fwhm <= 25 else 2.0), mats=pm)
    s0f = c0.emitter(**opt.EMITTERS[em])
    s = run_reference(c0, {**pf, **rf.params_for("ref", x, 3, 2)}, s0f, du=2.5e-4)["summary"]
    print(f"S fwhm {fwhm} {mode}: CE {s['CE']:.0f}  EQE {s['eqe']:.3f}  x {np.round(x, 1)}", flush=True)
    with open(path, "w") as fh:
        json.dump(dict(kind="reference", fwhm=fwhm, mode=mode, x=list(map(float, x)), CE=s["CE"], eqe=s["eqe"]),
                  fh, indent=1)
    return path


def _task(t):
    return chain_recycling(t[1:]) if t[0] == "R" else chain_reference(t[1:])


def run():
    os.makedirs(OUT, exist_ok=True)
    tasks = [("R", f, N) for f in sorted(FWHMS, reverse=True) for N in NS] + \
            [("S", f, m) for f in sorted(FWHMS, reverse=True) for m in REF_PRESETS]
    with Pool(4) as pool:
        for p in pool.imap_unordered(_task, tasks):
            print("done", os.path.basename(p), flush=True)


# ----------------------------------------------------------------------------------------------
def load():
    R, S = {}, {}
    for p in glob.glob(os.path.join(OUT, "*.json")):
        d = json.load(open(p))
        if d["kind"] == "recycling":
            R.setdefault(d["fwhm"], {})[d["N"]] = d["rows"]
        else:
            S.setdefault(d["mode"], {})[d["fwhm"]] = d["CE"]
    return R, S


def envelope(R, fwhm, Agrid):
    """best-over-N recycling CE on a common A' grid (monotone interpolation per N, then max)."""
    best = np.full(Agrid.shape, np.nan)
    for N, rows in R[fwhm].items():
        A = np.array([r["Aprime"] for r in rows]); C = np.array([r["CE"] for r in rows])
        o = np.argsort(A); A, C = A[o], C[o]
        v = np.interp(Agrid, A, C, left=np.nan, right=np.nan)
        best = np.fmax(best, v)
    return best


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    os.makedirs(FIG, exist_ok=True)
    R, S = load()
    fw = np.array(sorted(R))
    Ag = np.linspace(0.03, 0.20, 69)
    M = np.array([envelope(R, f, Ag) for f in fw])
    sR = np.array([S.get("realAg", {}).get(f, np.nan) for f in fw])
    sI = np.array([S.get("idealAg", {}).get(f, np.nan) for f in fw])
    cmp = json.load(open(os.path.join(HERE, "results", "compare_anodes.json")))
    dev = []  # (label, A', fwhm-equivalent, CE, marker)
    for key, lab, mk in (("Ag5", "meas. Ag 5 nm", "v"), ("Ag7", "meas. Ag 7 nm", "<"), ("Ag8", "meas. Ag 8 nm", "s"),
                         ("Ag_ideal", "ideal thin Ag", "D"), ("IZO", "IZO", "o")):
        for em, fwe in (("n15", 15), ("irppy", 64)):
            r = cmp[key][em]
            if r and r.get("Aprime") is not None:
                dev.append((lab, r["Aprime"], fwe, r["CE"], mk, em))

    # Figure 1: map
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    lv = [200, 300, 400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800]
    cs = ax.contourf(100 * Ag, fw, M, levels=40, cmap="viridis")
    cl = ax.contour(100 * Ag, fw, M, levels=lv, colors="w", linewidths=0.6)
    ax.clabel(cl, fmt="%d", fontsize=7)
    for sv, ls, lab in ((sR, "-", "= strong cavity (measured Ag)"), (sI, "--", "= strong cavity (ideal Ag)")):
        ratio = M / sv[:, None]
        if np.nanmin(ratio) < 1 < np.nanmax(ratio):
            c = ax.contour(100 * Ag, fw, ratio, levels=[1.0], colors="r", linestyles=ls, linewidths=1.6)
            c.collections[0].set_label("recycling " + lab) if hasattr(c, "collections") else None
    for lab, A, f, ce, mk, em in dev:
        ax.plot(100 * A, f, mk, mfc="none", mec="w" if em == "n15" else "orange", ms=7)
    ax.set_xlabel("round-trip loss A' (%)")
    ax.set_ylabel("emitter FWHM (nm), Gaussian @530 nm")
    ax.set_title("Recycling TEOLED: max on-axis cd/A (organic k 1e-4)")
    fig.colorbar(cs, ax=ax, label="on-axis cd/A")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color="r", ls="-", label="= strong cavity, measured Ag"),
         Line2D([], [], color="r", ls="--", label="= strong cavity, ideal Ag")]
    h += [Line2D([], [], marker=mk, ls="", mfc="none", mec="k", label=lab) for lab, mk in
          (("meas. Ag 5/7/8 nm", "s"), ("ideal thin Ag", "D"), ("IZO", "o"))]
    ax.legend(handles=h, fontsize=7, loc="upper right")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "amap_map.png"), dpi=160); plt.close(fig)

    # Figure 2: cuts cd/A vs A' for several FWHM, with strong-cavity levels and device points
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    cols = plt.cm.plasma(np.linspace(0.05, 0.85, len(fw)))
    for i, f in enumerate(fw):
        ax.plot(100 * Ag, M[i], color=cols[i], lw=1.8, label=f"{f} nm")
        ax.plot([3, 20], [sR[i]] * 2, color=cols[i], ls=":", lw=1.0)
    for lab, A, f, ce, mk, em in dev:
        i = int(np.argmin(abs(fw - f)))
        ax.plot(100 * A, ce, mk, color="k" if em == "n15" else "tab:orange", ms=6, mfc="none")
    ax.set_xlabel("round-trip loss A' (%)"); ax.set_ylabel("on-axis cd/A")
    ax.set_title("solid: recycling (best DBR pairs)   dotted: strong cavity, measured Ag")
    ax.legend(title="FWHM", fontsize=7, ncol=2)
    ax.set_xlim(3, 20); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "amap_cuts.png"), dpi=160); plt.close(fig)

    # Figure 3: cd/A vs FWHM for strong cavity and recycling at fixed A'
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(fw, sR, "k-o", label="strong cavity, measured Ag")
    ax.plot(fw, sI, "k--s", label="strong cavity, ideal Ag")
    for Av, c in ((0.035, "tab:green"), (0.05, "tab:blue"), (0.08, "tab:purple"), (0.12, "tab:red")):
        j = int(np.argmin(abs(Ag - Av)))
        ax.plot(fw, M[:, j], "-", color=c, lw=2, label=f"recycling, A' = {100 * Av:.1f}%")
    ax.set_xlabel("emitter FWHM (nm)"); ax.set_ylabel("on-axis cd/A"); ax.grid(alpha=0.3)
    ax.legend(fontsize=8); ax.set_title("Gaussian emitters @530 nm, organic k 1e-4")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "amap_vs_fwhm.png"), dpi=160); plt.close(fig)

    # crossover A'* where recycling = strong cavity
    cross = {}
    for i, f in enumerate(fw):
        for name, sv in (("realAg", sR[i]), ("idealAg", sI[i])):
            v = M[i]; ok = ~np.isnan(v)
            above = v[ok] > sv
            cross.setdefault(name, {})[int(f)] = (float(100 * np.interp(0.0, -(v[ok] - sv), Ag[ok]))
                                                  if above.any() and not above.all() else
                                                  (">20" if above.all() else "<3"))
    summ = dict(fwhm=fw.tolist(), strong_real=sR.tolist(), strong_ideal=sI.tolist(),
                recyc_at={f"{100 * Av:.1f}": M[:, int(np.argmin(abs(Ag - Av)))].tolist() for Av in (0.035, 0.05, 0.08, 0.12)},
                crossover_Aprime_percent=cross)
    json.dump(summ, open(os.path.join(HERE, "results", "amap_summary.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    run() if sys.argv[1] == "run" else plot()
