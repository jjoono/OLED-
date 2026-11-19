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
    "font.size": 16, "axes.labelsize": 20, "xtick.labelsize": 16, "ytick.labelsize": 16, "legend.fontsize": 13,
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

    fig, ax = plt.subplots(figsize=(9.0, 7.6))
    x1, y1 = smooth(fw, b1)
    xf = np.linspace(8, 82, 600)
    yb = np.interp(xf, x1, y1, left=y1[0], right=TOP + 5)
    ax.fill_between(xf, 2.5, np.minimum(yb, TOP + 5), color="#dcefdc", lw=0, zorder=0)
    ax.fill_between(xf, np.minimum(yb, TOP + 5), TOP + 5, color="#f8dcd8", lw=0, zorder=0)

    ax.plot(x1, y1, "-", color="k", lw=2.8, label="Break-even (strong cavity, measured Ag)", zorder=3)
    xi, yi = smooth(fw, bI)
    ax.plot(xi, yi, "--", color="k", lw=2.0, label="Break-even (strong cavity, ideal Ag)", zorder=3)
    for b, lab in ((b15, "×1.5"), (b2, "×2")):
        xs, ys = smooth(fw, b)
        ax.plot(xs, ys, ":", color="0.35", lw=1.8, zorder=2)
        ax.text(xs[-1] - 0.5, ys[-1] - 0.35, lab, fontsize=15, color="0.25", ha="right", va="bottom")
    ax.plot([], [], ":", color="0.35", lw=1.8, label="Recycling = ×1.5 / ×2 of strong cavity")

    ax.text(60, 5.9, "Recycling wins", fontsize=22, weight="bold", color="#2e7d32", ha="center")
    ax.text(23, 17.6, "Strong cavity wins", fontsize=22, weight="bold", color="#b23b2e", ha="center")

    for xv in (15, 64):
        ax.axvline(xv, color="0.55", lw=1.0, zorder=1)
    pts = [("IZO", "o", "IZO"), ("Ag_ideal", "D", "Ideal Ag"), ("Ag8", "s", "Ag 8 nm"),
           ("Ag7", "^", "Ag 7 nm"), ("Ag5", "v", "Ag 5 nm")]
    off = {("n15", "IZO"): (9, 0), ("n15", "Ag_ideal"): (9, 0), ("n15", "Ag8"): (12, 2), ("n15", "Ag7"): (-12, -14),
           ("n15", "Ag5"): (9, 0), ("irppy", "IZO"): (9, 6), ("irppy", "Ag_ideal"): (9, -6),
           ("irppy", "Ag8"): (-11, -2), ("irppy", "Ag7"): (9, -2), ("irppy", "Ag5"): (9, 0)}
    for em, x in (("n15", 15.0), ("irppy", 64.0)):
        ref = cmp["ref_real"][em]["CE"]
        for key, mk, name in pts:
            r = cmp[key][em]; y = 100 * r["Aprime"]; g = r["CE"] / ref
            ax.plot(x, y, mk, ms=11, mfc="#2e7d32" if g > 1 else "#c0392b", mec="k", mew=1.3, zorder=5)
            dx, dy = off[(em, key)]
            if (em, key) == ("n15", "Ag7"):
                continue                     # overlaps Ag 8 nm (A' 9.5 vs 9.6 %): labelled together below
            if (em, key) == ("n15", "Ag8"):
                g7 = cmp["Ag7"][em]["CE"] / ref
                name, gl = "Ag 7 / 8 nm", f"×{g7:.2f} / {g:.2f}"
            else:
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
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False, handlelength=2.6, fontsize=13,
              columnspacing=1.5)
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(FIG, f"amap_phase_origin.{ext}"), dpi=300 if ext == "png" else None)
    plt.close(fig)


if __name__ == "__main__":
    main()
