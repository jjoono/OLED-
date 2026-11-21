"""
arch_summary.py -- ranking of all architectures from arch.py (results/arch/*.json):
markdown tables, Excel workbook (results/arch_summary.xlsx) and Origin-style figures
  results/figs/arch_ranking.*   (A) fair architecture comparison (periodic DBR, 800-nm SiO2 ODR mirror,
                                    thick incoherent low-n films)
  results/figs/arch_upgrades.*  (B) component upgrades on the best architecture (waterfall)
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

LABEL = {   # (A) architecture comparison: same periodic DBR parametrisation, 800-nm SiO2 ODR mirror
    "TE_air": "TE, ZnS/LiF DBR facing air (N2-gap glass lid)",
    "TE_air-ox": "TE, low-loss TiO2/SiO2 DBR (k 2e-5) facing air",
    "TE_lowN1.1-inc": "TE, n 1.1 film + cover glass",
    "TE_lowN1.2-inc": "TE, n 1.2 film + cover glass",
    "TE_glass": "TE, cover glass bonded with OCA (n 1.5)",
    "BE_glass-pixel": "BE, TiO2/SiO2 DBR on glass, small pixel (glass-trapped light lost)",
    "BE_glass-large": "BE, TiO2/SiO2 DBR on glass, large area (glass-trapped light recycled)",
    "BE_lowN1.2-inc": "BE, n 1.2 film between DBR and glass",
    "BE_lowN1.1-inc": "BE, n 1.1 film between DBR and glass",
    "BE_lowN1.05-inc": "BE, n 1.05 film (SiO2 nanorod limit) between DBR and glass",
}
ORDER = list(LABEL)
SHORT = {"TE_air": "TE air\nZnS/LiF", "TE_air-ox": "TE air\noxide DBR", "TE_lowN1.1-inc": "TE\n+n1.1+glass",
         "TE_lowN1.2-inc": "TE\n+n1.2+glass", "TE_glass": "TE\n+OCA glass", "BE_glass-pixel": "BE glass\npixel",
         "BE_glass-large": "BE glass\nlarge area", "BE_lowN1.2-inc": "BE\n+n1.2", "BE_lowN1.1-inc": "BE\n+n1.1",
         "BE_lowN1.05-inc": "BE\n+n1.05"}
EXTRA = {   # (B) beyond the fair comparison
    "TE_air-ox-edr": "TE air, oxide DBR + enhanced dense mirror (SiO2 800 + 2 TiO2/SiO2 pairs on Ag)",
    "TE_air-ox-epr": "TE air, oxide DBR + enhanced porous mirror (n 1.2 1 um + 2 pairs on Ag)",
    "TE_air-epr": "TE air, ZnS/LiF DBR + enhanced porous mirror",
    "BE_lowN1.05-inc-epr": "BE n 1.05 film + enhanced porous mirror",
    "BE_lowN1.1": "BE, n 1.1 film, coherent (uniform thickness: its far interface is one more filter element)",
    "BE_lowN1.05": "BE, n 1.05 film, coherent (uniform thickness)",
    "BE_lowN1.2": "BE, n 1.2 film, coherent",
    "TE_lowN1.2": "TE, n 1.2 film + glass, coherent",
    "BE_lowN1.1-ev": "BE, n 1.1 film coherent, ZnS/LiF-quality DBR",
}
REF = {"irppy": {"real": 263.5, "ideal": 302.5}, "n15": {"real": 775.3, "ideal": 1030.1}}
TE_AIR_FALLBACK = {("IZO", "n15"): ("cmp_n15_IZO_N8", 1474.6), ("Ag8m", "n15"): ("cmp_n15_Ag8m_N4", 721.1),
                   ("IZO", "irppy"): ("cmp_irppy_IZO_N8", None), ("Ag8m", "irppy"): ("cmp_irppy_Ag8m_N6", 425.0)}
EMN = (("n15", "532 nm / 15 nm"), ("irppy", "Ir(ppy)2acac"))


def load():
    rows = []
    for p in glob.glob(os.path.join(HERE, "results", "arch", "*.json")):
        d = json.load(open(p)); f = d["final"]
        P = d["params"]
        rows.append(dict(variant=d["case"].rsplit("_", 2)[0], anode=d["anode"], emitter=d["emitter"], CE=f["CE"],
                         EQE=f["eqe"], Aprime=100 * f["Aprime"], p=100 * f["p_avg"], N=d["N"], n_m=P["n_m"],
                         d_etl=P["d_etl"], d_htl=P["d_htl"], d_anode=P["d_ag"], dH1=P["dbr"]["dH1"], dH=P["dbr"]["dH"],
                         dL=P["dbr"]["dL"], n_low=P.get("n_low"), d_low=P.get("d_low"),
                         mirror=json.dumps(P.get("bot")) if P.get("bot") else "SiO2 800 nm ODR",
                         CE_other_glass_mode=(d["other_glass_mode"] or {}).get("CE"),
                         source=f"adopted from {d['adopted_from']} (not re-optimised)" if d.get("adopted_from") else "arch.py"))
    have = {(r["variant"], r["anode"], r["emitter"]) for r in rows}
    cmp = json.load(open(os.path.join(HERE, "results", "compare_anodes.json")))
    for (an, em), (tag, ce05) in TE_AIR_FALLBACK.items():
        if ("TE_air", an, em) in have:
            continue
        r = cmp["IZO" if an == "IZO" else "Ag8"][em]
        P = json.load(open(os.path.join(HERE, "results", f"opt_{tag}.json")))["params"]
        rows.append(dict(variant="TE_air", anode=an, emitter=em, CE=ce05 or r["CE"], EQE=r["eqe"], Aprime=100 * r["Aprime"],
                         p=100 * r["p"], N=r["N"], n_m=P["n_m"], d_etl=P["d_etl"], d_htl=P["d_htl"], d_anode=P["d_ag"],
                         dH1=P["dbr"]["dH1"], dH=P["dbr"]["dH"], dL=P["dbr"]["dL"], n_low=None, d_low=None,
                         mirror="SiO2 800 nm ODR", CE_other_glass_mode=None, source="compare_anodes.py"))
    df = pd.DataFrame(rows)
    df["x vs strong cavity (measured Ag)"] = [r.CE / REF[r.emitter]["real"] for r in df.itertuples()]
    df["x vs strong cavity (ideal Ag)"] = [r.CE / REF[r.emitter]["ideal"] for r in df.itertuples()]
    df["group"] = ["A" if v in LABEL else "B" for v in df.variant]
    df["order"] = [ORDER.index(v) if v in ORDER else 50 + list(EXTRA).index(v) if v in EXTRA else 99 for v in df.variant]
    df["structure"] = [LABEL.get(v) or EXTRA.get(v) or v for v in df.variant]
    return df.sort_values(["emitter", "group", "order", "anode"])


def cell(df, v, an, em):
    s = df[(df.variant == v) & (df.anode == an) & (df.emitter == em)]
    if not len(s):
        return "-"
    r = s.iloc[0]
    return (f"**{r.CE:.0f}** (x{r['x vs strong cavity (measured Ag)']:.2f}), {100 * r.EQE:.0f}%, A' {r.Aprime:.1f}%, "
            f"N={int(r.N)}")


def table(df):
    out = []
    for em, emn in EMN:
        out.append(f"\n### {emn}  (strong cavity: measured Ag {REF[em]['real']:.0f}, ideal Ag {REF[em]['ideal']:.0f} cd/A)\n")
        out.append("| (A) structure | IZO anode: cd/A (x vs SC), EQE, A', N | measured Ag 8 nm |")
        out.append("|---|---|---|")
        for v in ORDER:
            out.append(f"| {LABEL[v]} | {cell(df, v, 'IZO', em)} | {cell(df, v, 'Ag8m', em)} |")
        out.append("\n| (B) beyond the fair comparison | IZO anode | measured Ag 8 nm |")
        out.append("|---|---|---|")
        for v in EXTRA:
            if len(df[(df.variant == v) & (df.emitter == em)]):
                out.append(f"| {EXTRA[v]} | {cell(df, v, 'IZO', em)} | {cell(df, v, 'Ag8m', em)} |")
    return "\n".join(out)


def figure(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import plot_phase_origin  # noqa: F401  (Origin-style rcParams, Calibri/Carlito)
    vs = [v for v in ORDER if v in set(df.variant)]
    fig, axs = plt.subplots(2, 1, figsize=(15, 10.5), sharex=True)
    for ax, (em, emn) in zip(axs, (("n15", "532 nm / 15 nm"), ("irppy", "Ir(ppy)$_2$acac"))):
        x = np.arange(len(vs)); w = 0.38
        for k, (an, col, lab) in enumerate((("IZO", "#2e7d32", "IZO anode"), ("Ag8m", "#7f8c8d", "Measured Ag 8 nm"))):
            y = np.array([df[(df.variant == v) & (df.anode == an) & (df.emitter == em)].CE.max() for v in vs], float)
            ax.bar(x + (k - 0.5) * w, np.nan_to_num(y), w, color=col, edgecolor="k", lw=1.2, label=lab)
            for xi, yi in zip(x + (k - 0.5) * w, y):
                if np.isfinite(yi):
                    ax.text(xi, yi + 18, f"{yi:.0f}", ha="center", va="bottom", fontsize=12)
        ax.axhline(REF[em]["real"], color="#b23b2e", lw=2.2, label="Strong cavity (measured Ag)")
        ax.axhline(REF[em]["ideal"], color="#b23b2e", lw=1.6, ls="--", label="Strong cavity (ideal Ag)")
        ax.set_ylabel("On-axis cd/A")
        ax.text(0.01, 0.95, emn, transform=ax.transAxes, fontsize=19, va="top")
        ax.set_ylim(0, np.nanmax([df[(df.emitter == em) & (df.group == "A")].CE.max(), REF[em]["ideal"]]) * 1.2)
    axs[0].legend(loc="upper right", frameon=False, fontsize=13, ncol=2)
    axs[1].set_xticks(np.arange(len(vs))); axs[1].set_xticklabels([SHORT[v] for v in vs], fontsize=13)
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(FIG, f"arch_ranking.{ext}"), dpi=300 if ext == "png" else None)
    plt.close(fig)


STEPS = (("TE_air", "TE air, ZnS/LiF DBR\n(reference)"), ("TE_air-ox", "+ low-loss oxide DBR"),
         ("TE_air-ox-edr", "+ dense enhanced mirror"), ("TE_air-ox-epr", "+ porous enhanced mirror"))


def upgrades(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import plot_phase_origin  # noqa: F401
    fig, axs = plt.subplots(1, 2, figsize=(15, 6.2))
    for ax, (em, emn) in zip(axs, (("n15", "532 nm / 15 nm"), ("irppy", "Ir(ppy)$_2$acac"))):
        st = [(v, lab, df[(df.variant == v) & (df.anode == "IZO") & (df.emitter == em)].CE.max()) for v, lab in STEPS]
        ff = os.path.join(HERE, "results", f"freeform_TE_air-ox-epr_IZO_{em}.json")
        if os.path.exists(ff):
            st.append(("freeform", "+ aperiodic filter\n(thick top spacer)", json.load(open(ff))["final"]["CE"]))
        st = [s for s in st if np.isfinite(s[2])]
        y = np.arange(len(st))[::-1]
        for yi, (v, lab, ce) in zip(y, st):
            ax.barh(yi, ce, color="#2e7d32" if v != "TE_air" else "#7f8c8d", edgecolor="k", lw=1.2, height=0.62)
            ax.text(ce + 25, yi, f"{ce:.0f}", va="center", fontsize=15)
        ax.axvline(REF[em]["real"], color="#b23b2e", lw=2.2)
        ax.axvline(REF[em]["ideal"], color="#b23b2e", lw=1.6, ls="--")
        ax.set_yticks(y); ax.set_yticklabels([s[1] for s in st], fontsize=14)
        ax.set_xlabel("On-axis cd/A (IZO anode)")
        ax.set_title(emn, fontsize=19)
        ax.set_xlim(0, max(s[2] for s in st) * 1.18)
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(FIG, f"arch_upgrades.{ext}"), dpi=300 if ext == "png" else None)
    plt.close(fig)


if __name__ == "__main__":
    df = load()
    print(table(df))
    cols = ["group", "structure", "variant", "emitter", "anode", "CE", "x vs strong cavity (measured Ag)",
            "x vs strong cavity (ideal Ag)", "EQE", "Aprime", "p", "N", "n_m", "d_etl", "d_htl", "d_anode", "dH1", "dH",
            "dL", "n_low", "d_low", "mirror", "CE_other_glass_mode", "source"]
    with pd.ExcelWriter(os.path.join(HERE, "results", "arch_summary.xlsx")) as xw:
        df[cols].to_excel(xw, sheet_name="designs", index=False)
        bp = os.path.join(HERE, "results", "arch_budget.json")
        if os.path.exists(bp):
            b = pd.DataFrame(json.load(open(bp))).T
            b.index.name = "case"
            b.to_excel(xw, sheet_name="loss budget per round trip")
    figure(df)
    upgrades(df)
