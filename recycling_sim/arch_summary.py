"""
arch_summary.py -- ranking of all architectures from arch.py (results/arch/*.json):
markdown table, Excel sheet (results/arch_summary.xlsx) and an Origin-style bar chart.
TE-air falls back to the validated compare_anodes designs when the arch.py run is missing.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "results", "figs")

LABEL = {
    "TE_air": "TE, ZnS/LiF DBR facing air (lab device under a gas-filled lid)",
    "TE_air-ox": "TE, low-loss TiO2/SiO2 DBR (k 2e-5) facing air",
    "TE_glass": "TE, cover glass bonded with OCA (n 1.5) on the DBR",
    "TE_lowN1.2": "TE, n 1.2 layer between DBR and cover glass",
    "BE_glass-pixel": "BE, TiO2/SiO2 DBR on glass, small pixel (glass-trapped light lost)",
    "BE_glass-large": "BE, TiO2/SiO2 DBR on glass, large area (glass-trapped light recycled)",
    "BE_lowN1.2": "BE, n 1.2 layer between TiO2/SiO2 DBR and glass (pixel)",
    "BE_lowN1.1": "BE, n 1.1 layer between TiO2/SiO2 DBR and glass (pixel)",
    "BE_lowN1.1-ev": "BE, n 1.1 layer, ZnS/LiF-quality DBR (k 1e-4) (pixel)",
}
ORDER = list(LABEL)
REF = {"irppy": {"real": 263.5, "ideal": 302.5}, "n15": {"real": 775.3, "ideal": 1030.1}}
TE_AIR_FALLBACK = {("IZO", "n15"): ("cmp_n15_IZO_N8", 1474.6), ("Ag8m", "n15"): ("cmp_n15_Ag8m_N4", 721.1),
                   ("IZO", "irppy"): ("cmp_irppy_IZO_N8", None), ("Ag8m", "irppy"): ("cmp_irppy_Ag8m_N6", 425.0)}


def load():
    rows = []
    for p in glob.glob(os.path.join(HERE, "results", "arch", "*.json")):
        d = json.load(open(p)); f = d["final"]
        key = d["case"].rsplit("_", 2)[0]
        rows.append(dict(variant=key, anode=d["anode"], emitter=d["emitter"], CE=f["CE"], EQE=f["eqe"],
                         Aprime=100 * f["Aprime"], p=f["p_avg"], N=d["N"], n_m=d["params"]["n_m"],
                         d_etl=d["params"]["d_etl"], d_htl=d["params"]["d_htl"], d_anode=d["params"]["d_ag"],
                         dbr=d["params"]["dbr"], d_low=d["params"].get("d_low"),
                         CE_other_glass_mode=(d["other_glass_mode"] or {}).get("CE")))
    have = {(r["variant"], r["anode"], r["emitter"]) for r in rows}
    cmp = json.load(open(os.path.join(HERE, "results", "compare_anodes.json")))
    for (an, em), (tag, ce05) in TE_AIR_FALLBACK.items():
        if ("TE_air", an, em) in have:
            continue
        r = cmp["IZO" if an == "IZO" else "Ag8"][em]
        d = json.load(open(os.path.join(HERE, "results", f"opt_{tag}.json")))
        rows.append(dict(variant="TE_air", anode=an, emitter=em, CE=ce05 or r["CE"], EQE=r["eqe"], Aprime=100 * r["Aprime"],
                         p=r["p"], N=r["N"], n_m=d["params"]["n_m"], d_etl=d["params"]["d_etl"],
                         d_htl=d["params"]["d_htl"], d_anode=d["params"]["d_ag"], dbr=d["params"]["dbr"], d_low=None,
                         CE_other_glass_mode=None))
    df = pd.DataFrame(rows)
    df["vs strong cavity (measured Ag)"] = [r.CE / REF[r.emitter]["real"] for r in df.itertuples()]
    df["vs strong cavity (ideal Ag)"] = [r.CE / REF[r.emitter]["ideal"] for r in df.itertuples()]
    df["order"] = [ORDER.index(v) if v in ORDER else 99 for v in df.variant]
    return df.sort_values(["emitter", "anode", "order"])


def table(df):
    out = []
    for em, emn in (("n15", "532 nm / 15 nm"), ("irppy", "Ir(ppy)2acac")):
        out.append(f"\n### {emn}  (strong cavity: measured Ag {REF[em]['real']:.0f}, ideal Ag {REF[em]['ideal']:.0f} cd/A)\n")
        out.append("| structure | IZO anode: cd/A (x vs SC) EQE, A' | measured Ag 8 nm: cd/A (x vs SC) EQE, A' |")
        out.append("|---|---|---|")
        for v in ORDER:
            cells = []
            for an in ("IZO", "Ag8m"):
                s = df[(df.variant == v) & (df.anode == an) & (df.emitter == em)]
                if len(s):
                    r = s.iloc[0]
                    cells.append(f"**{r.CE:.0f}** (x{r['vs strong cavity (measured Ag)']:.2f}), {100 * r.EQE:.0f}%, "
                                 f"{r.Aprime:.1f}%, N={int(r.N)}")
                else:
                    cells.append("-")
            out.append(f"| {LABEL[v]} | {cells[0]} | {cells[1]} |")
    return "\n".join(out)


def figure(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import plot_phase_origin  # noqa: F401  (Origin-style rcParams, Calibri/Carlito)
    short = {"TE_air": "TE\nair", "TE_air-ox": "TE air\noxide DBR", "TE_glass": "TE\n+OCA glass",
             "TE_lowN1.2": "TE\n+n1.2+glass", "BE_glass-pixel": "BE glass\npixel",
             "BE_glass-large": "BE glass\nlarge area", "BE_lowN1.2": "BE\n+n1.2", "BE_lowN1.1": "BE\n+n1.1",
             "BE_lowN1.1-ev": "BE +n1.1\nZnS/LiF DBR"}
    vs = [v for v in ORDER if v in set(df.variant)]
    fig, axs = plt.subplots(2, 1, figsize=(13, 10), sharex=True)
    for ax, em, emn in zip(axs, ("n15", "irppy"), ("532 nm / 15 nm", "Ir(ppy)$_2$acac")):
        x = np.arange(len(vs)); w = 0.38
        for k, (an, col, lab) in enumerate((("IZO", "#2e7d32", "IZO anode"), ("Ag8m", "#7f8c8d", "Measured Ag 8 nm"))):
            y = [df[(df.variant == v) & (df.anode == an) & (df.emitter == em)].CE.max() for v in vs]
            ax.bar(x + (k - 0.5) * w, y, w, color=col, edgecolor="k", lw=1.2, label=lab)
            for xi, yi in zip(x + (k - 0.5) * w, y):
                if np.isfinite(yi):
                    ax.text(xi, yi + 15, f"{yi:.0f}", ha="center", va="bottom", fontsize=12)
        ax.axhline(REF[em]["real"], color="#b23b2e", lw=2.2, label="Strong cavity (measured Ag)")
        ax.axhline(REF[em]["ideal"], color="#b23b2e", lw=1.6, ls="--", label="Strong cavity (ideal Ag)")
        ax.set_ylabel("On-axis cd/A")
        ax.text(0.01, 0.95, emn, transform=ax.transAxes, fontsize=18, va="top")
        ax.set_ylim(0, max(df[df.emitter == em].CE.max(), REF[em]["ideal"]) * 1.18)
    axs[0].legend(loc="upper right", frameon=False, fontsize=13, ncol=2)
    axs[1].set_xticks(np.arange(len(vs))); axs[1].set_xticklabels([short[v] for v in vs], fontsize=13)
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(FIG, f"arch_ranking.{ext}"), dpi=300 if ext == "png" else None)
    plt.close(fig)


if __name__ == "__main__":
    df = load()
    print(table(df))
    cols = ["variant", "emitter", "anode", "CE", "vs strong cavity (measured Ag)", "vs strong cavity (ideal Ag)", "EQE",
            "Aprime", "p", "N", "n_m", "d_etl", "d_htl", "d_anode", "dbr", "d_low", "CE_other_glass_mode"]
    x = df[cols].copy(); x["dbr"] = x["dbr"].astype(str); x["variant"] = [LABEL.get(v, v) for v in x.variant]
    x.to_excel(os.path.join(HERE, "results", "arch_summary.xlsx"), index=False)
    figure(df)
