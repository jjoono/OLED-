"""
plot_phase_origin.py -- Origin-style redraw of the recycling vs strong-cavity win/lose diagram (amap_phase).
Calibri (Carlito, metric-compatible, used for the PNG); SVG/PDF keep text as 'Calibri'.
usage: python3 plot_phase_origin.py   (needs results/amap/*.json and results/compare_anodes.json)
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy.interpolate import PchipInterpolator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import amap  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "results", "figs")

plt.rcParams.update({
    "font.family": ["Calibri", "Carlito", "DejaVu Sans"],
    "font.size": 17, "axes.labelsize": 21, "xtick.labelsize": 18, "ytick.labelsize": 18, "legend.fontsize": 13,
    "axes.linewidth": 1.6, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.major.size": 7, "ytick.major.size": 7, "xtick.minor.size": 3.5, "ytick.minor.size": 3.5,
    "xtick.major.width": 1.4, "ytick.major.width": 1.4, "xtick.minor.width": 1.1, "ytick.minor.width": 1.1,
    "xtick.top": True, "ytick.right": True, "svg.fonttype": "none", "pdf.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Carlito", "mathtext.it": "Carlito:italic",
})

TOP = 20.0


def crossing(M, Ag, ref, g):
    out = []
    for i in range(M.shape[0]):
        d = M[i] - g * ref[i]; ok = ~np.isnan(d); dd, aa = d[ok], 100 * Ag[ok]
        if (dd > 0).all():
            out.append(np.inf)
        elif (dd < 0).all():
            out.append(np.nan)
        else:
            j = int(np.argmax(dd < 0))
            out.append(float(np.interp(0.0, [-dd[j - 1], -dd[j]], [aa[j - 1], aa[j]])))
    return np.array(out)


def smooth(fw, y, x0=10, x1=80):
    """PCHIP through finite points; a line leaving the frame (inf) is continued linearly to y = TOP."""
    xs, ys = [], []
    for k, (x, v) in enumerate(zip(fw, y)):
        if np.isnan(v):
            continue
        if np.isinf(v):
            if xs and ys[-1] < TOP and k > 0:
                xe = xs[-1] + (x - xs[-1]) * (TOP + 0.6 - ys[-1]) / (25.0 - ys[-1])
                xs.append(xe); ys.append(TOP + 0.6)
            break
        xs.append(x); ys.append(v)
    xs, ys = np.array(xs), np.array(ys)
    xx = np.linspace(xs[0], xs[-1], 400)
    return xx, PchipInterpolator(xs, ys)(xx)


def main():
    R, S = amap.load()
    fw = np.array(sorted(R), float)
    Ag = np.linspace(0.035, 0.20, 67)
    M = np.array([amap.envelope(R, int(f), Ag) for f in fw])
    sR = np.array([S["realAg"][int(f)] for f in fw]); sI = np.array([S["idealAg"][int(f)] for f in fw])
    cmp = json.load(open(os.path.join(HERE, "results", "compare_anodes.json")))

    b1 = crossing(M, Ag, sR, 1.0); bI = crossing(M, Ag, sI, 1.0)
    b15 = crossing(M, Ag, sR, 1.5); b2 = crossing(M, Ag, sR, 2.0)

    fig, ax = plt.subplots(figsize=(9.0, 6.9))
    x1, y1 = smooth(fw, b1)
    xf = np.linspace(8, 82, 600)
    yb = np.interp(xf, x1, y1, left=y1[0], right=TOP + 5)
    ax.fill_between(xf, 2.5, np.minimum(yb, TOP + 5), color="#dcefdc", lw=0, zorder=0)
    ax.fill_between(xf, np.minimum(yb, TOP + 5), TOP + 5, color="#f8dcd8", lw=0, zorder=0)

    for b, lab in ((b15, "×1.5"), (b2, "×2")):
        xs, ys = smooth(fw, b)
        ax.plot(xs, ys, ":", color="0.35", lw=1.8, zorder=2)
        ax.text(xs[-1] - 0.5, ys[-1] - 0.35, lab, fontsize=15, color="0.25", ha="right", va="bottom")
    ax.plot([], [], ":", color="0.35", lw=1.8, label="Recycling = ×1.5 / ×2 of strong cavity")

    ax.text(60, 5.9, "Recycling wins", fontsize=22, weight="bold", color="#2e7d32", ha="center")
    ax.text(23, 17.6, "Strong cavity wins", fontsize=22, weight="bold", color="#b23b2e", ha="center")

    for xv in (15, 64):
        ax.axvline(xv, color="0.55", lw=1.0, zorder=1)
    pts = [("IZO", "o", "IZO"), ("Ag8", "s", "Ag 8 nm")]
    off = {("n15", "IZO"): (9, 0), ("n15", "Ag_ideal"): (9, 0), ("n15", "Ag8"): (12, 0), ("n15", "Ag7"): (-12, -14),
           ("n15", "Ag5"): (9, 0), ("irppy", "IZO"): (9, 0), ("irppy", "Ag_ideal"): (9, -6),
           ("irppy", "Ag8"): (-11, -2), ("irppy", "Ag7"): (9, -2), ("irppy", "Ag5"): (9, 0)}
    for em, x in (("n15", 15.0), ("irppy", 64.0)):
        ref = cmp["ref_real"][em]["CE"]
        for key, mk, name in pts:
            r = cmp[key][em]; y = 100 * r["Aprime"]; g = r["CE"] / ref
            ax.plot(x, y, mk, ms=11, mfc="#2e7d32" if g > 1 else "#c0392b", mec="k", mew=1.3, zorder=5)
            dx, dy = off[(em, key)]
            if (em, key) == ("n15", "Ag7"):
                continue                     # overlaps Ag 8 nm (A' 9.5 vs 9.6 %): labelled together below
            gl = f"×{g:.2f}"
            ax.annotate(f"{name}  {gl}", (x, y), xytext=(dx, dy), textcoords="offset points", fontsize=14,
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75),
                        va="center", ha="right" if dx < 0 else "left", zorder=6)

    ax.set_xlim(8, 82); ax.set_ylim(TOP, 2.5)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2)); ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_xlabel("Emitter FWHM (nm)")
    ax.set_ylabel("Round-trip loss $A'$ (%)")
    sec = ax.secondary_xaxis("top")
    sec.set_xticks([15, 64]); sec.set_xticklabels(["532 nm / 15 nm", "Ir(ppy)$_2$acac"], fontsize=15)
    sec.tick_params(direction="in", length=0)
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(FIG, f"amap_phase_origin.{ext}"), dpi=300 if ext == "png" else None)
    plt.close(fig)


def panels():
    """Origin-style per-FWHM cuts, stacked vertically (4 x 1, shared A' axis)."""
    R, S = amap.load()
    fw = np.array(sorted(R), float)
    Ag = np.linspace(0.035, 0.20, 67)
    M = np.array([amap.envelope(R, int(f), Ag) for f in fw])
    sR = np.array([S["realAg"][int(f)] for f in fw]); sI = np.array([S["idealAg"][int(f)] for f in fw])
    b1 = crossing(M, Ag, sR, 1.0)
    cmp = json.load(open(os.path.join(HERE, "results", "compare_anodes.json")))
    pts = [("IZO", "o"), ("Ag_ideal", "D"), ("Ag8", "s"), ("Ag7", "^"), ("Ag5", "v")]
    sel = [15, 25, 50, 65]
    fig, axs = plt.subplots(4, 1, figsize=(7.2, 15.5), sharex=True)
    A = 100 * Ag
    for ax, f in zip(axs, sel):
        i = int(np.nonzero(fw == f)[0][0]); c = M[i]
        ax.fill_between(A, sR[i], c, where=c >= sR[i], color="#2e7d32", alpha=0.15, lw=0, interpolate=True)
        ax.fill_between(A, c, sR[i], where=c < sR[i], color="#b23b2e", alpha=0.15, lw=0, interpolate=True)
        ax.plot(A, c, color="#2e7d32", lw=3.0, label="Recycling")
        ax.axhline(sR[i], color="#b23b2e", lw=2.4, label="Strong cavity (measured Ag)")
        ax.axhline(sI[i], color="#b23b2e", lw=1.8, ls="--", label="Strong cavity (ideal Ag)")
        if np.isfinite(b1[i]) and b1[i] < 20:
            ax.axvline(b1[i], color="k", ls=":", lw=1.6)
            right = b1[i] < 15
            ax.text(b1[i] + (0.4 if right else -0.4), 0.42, f"$A'^*$ = {b1[i]:.1f}%", transform=ax.get_xaxis_transform(),
                    fontsize=17, va="center", ha="left" if right else "right",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))
        em = "n15" if f == 15 else ("irppy" if f == 65 else None)
        if em:
            for key, mk in pts:
                r = cmp[key][em]
                ax.plot(100 * r["Aprime"], r["CE"], mk, ms=11, mfc="none", mec="k", mew=1.6, zorder=5)
        lab = {15: "FWHM 15 nm  (532 nm / 15 nm devices)", 65: "FWHM 65 nm  (Ir(ppy)$_2$acac devices)"}.get(f, f"FWHM {f} nm")
        ax.text(0.97, 0.93, lab, transform=ax.transAxes, ha="right", va="top", fontsize=17)
        ax.set_ylabel("On-axis cd/A")
        from matplotlib.ticker import MultipleLocator
        ax.xaxis.set_major_locator(MultipleLocator(5)); ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.set_xlim(3, 20)
        lo = min(np.nanmin(c), sR[i]); hi = max(np.nanmax(c), sI[i])
        if em:
            hi = max(hi, max(cmp[k][em]["CE"] for k, _ in pts))
        ax.set_ylim(lo - 0.08 * (hi - lo), hi + 0.22 * (hi - lo))
    axs[-1].set_xlabel("Round-trip loss $A'$ (%)")
    axs[1].legend(loc="upper right", frameon=False, fontsize=14, bbox_to_anchor=(1.0, 0.86))
    from matplotlib.lines import Line2D
    mh = [Line2D([], [], marker=m, ls="", ms=11, mfc="none", mec="k", mew=1.6, label=t) for m, t in
          (("o", "IZO"), ("D", "Ideal thin Ag"), ("s", "Measured Ag 8 nm"), ("^", "Measured Ag 7 nm"), ("v", "Measured Ag 5 nm"))]
    axs[-1].legend(handles=mh, loc="upper right", bbox_to_anchor=(1.0, 0.86), frameon=False, fontsize=14,
                   handletextpad=0.4, labelspacing=0.35)
    fig.tight_layout(h_pad=0.4)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(FIG, f"amap_panels_origin.{ext}"), dpi=300 if ext == "png" else None)
    plt.close(fig)


def ratio_map():
    """filled contour of on-axis cd/A ratio recycling / strong cavity (measured Ag) vs FWHM and A'."""
    from matplotlib.colors import TwoSlopeNorm
    R, S = amap.load()
    fw = np.array(sorted(R), float)
    Ag = np.linspace(0.035, 0.20, 67)
    M = np.array([amap.envelope(R, int(f), Ag) for f in fw])
    sR = np.array([S["realAg"][int(f)] for f in fw])
    ratio = M / sR[:, None]
    ff = np.linspace(10, 80, 281)
    Z = np.array([PchipInterpolator(fw, ratio[:, j])(ff) for j in range(len(Ag))])  # (nA, nF)
    # extend the top edge (A' 3.0-3.5 %, below the 8-pair floor) with the 3.5 % values so the IZO points sit on the map
    Ag = np.r_[0.030, Ag]; Z = np.vstack([Z[:1], Z])
    cmp = json.load(open(os.path.join(HERE, "results", "compare_anodes.json")))

    fig, ax = plt.subplots(figsize=(9.6, 6.9))
    norm = TwoSlopeNorm(vmin=0.4, vcenter=1.0, vmax=3.0)
    levels = np.concatenate([np.linspace(0.4, 1.0, 13), np.linspace(1.0, 3.0, 21)[1:]])
    cf = ax.contourf(ff, 100 * Ag, np.clip(Z, 0.4, 3.0), levels=levels, cmap="RdYlGn", norm=norm, extend="both")
    cl = ax.contour(ff, 100 * Ag, Z, levels=[0.6, 0.8, 1.25, 1.5, 2.0, 2.5], colors="0.25", linewidths=1.0)
    ax.clabel(cl, fmt=lambda v: f"×{v:g}", fontsize=14, inline_spacing=6)
    c1 = ax.contour(ff, 100 * Ag, Z, levels=[1.0], colors="k", linewidths=2.6)
    ax.clabel(c1, fmt=lambda v: "×1", fontsize=15, inline_spacing=6)

    ax.text(23, 17.6, "Strong cavity wins\n(ratio < 1)", fontsize=20, weight="bold", color="#8b1a1a", ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.55))
    ax.text(60, 11.0, "Recycling wins\n(ratio > 1)", fontsize=20, weight="bold", color="#1b5e20", ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.55))

    for em, x in (("n15", 15.0), ("irppy", 64.0)):
        ref = cmp["ref_real"][em]["CE"]
        for key, mk, name, dx in (("IZO", "o", "IZO", 9), ("Ag8", "s", "Ag 8 nm", 9 if em == "n15" else -11)):
            r = cmp[key][em]; y = 100 * r["Aprime"]; g = r["CE"] / ref
            ax.plot(x, y, mk, ms=11, mfc="white", mec="k", mew=1.6, zorder=5)
            ax.annotate(f"{name}  ×{g:.2f}", (x, y), xytext=(dx, 0), textcoords="offset points", fontsize=14,
                        va="center", ha="right" if dx < 0 else "left", zorder=6,
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))
    ax.set_xlim(10, 80); ax.set_ylim(20, 3.0)
    ax.tick_params(pad=7)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2)); ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_xlabel("Emitter FWHM (nm)")
    ax.set_ylabel("Round-trip loss $A'$ (%)")
    sec = ax.secondary_xaxis("top")
    sec.set_xticks([15, 64]); sec.set_xticklabels(["532 nm / 15 nm", "Ir(ppy)$_2$acac"], fontsize=15)
    sec.tick_params(direction="in", length=0)
    cb = fig.colorbar(cf, ax=ax, ticks=[0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0], pad=0.02)
    cb.set_label("On-axis cd/A ratio (recycling / strong cavity)", fontsize=17)
    cb.ax.tick_params(labelsize=15, direction="in")
    cb.outline.set_linewidth(1.4)
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(FIG, f"amap_ratio_origin.{ext}"), dpi=300 if ext == "png" else None)
    plt.close(fig)


if __name__ == "__main__":
    main()
    panels()
    ratio_map()
