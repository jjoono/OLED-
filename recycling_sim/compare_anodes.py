"""
compare_anodes.py -- summary of the maximum on-axis current efficiency for
  (i)  strong-microcavity reference TEOLED,
  (ii) angle-selective recycling TEOLED with the top anode = measured thin Ag on HATCN (5 / 7 / 8 nm),
       ideal thin Ag (bulk McPeak optics, thickness optimised) or sputtered IZO,
for Ir(ppy)2acac and an ultra-narrow 532/15-nm emitter.  Common platform for (ii): 800-nm SiO2 ODR under
the SL, 15-nm low-k bottom ITO, organic k floor 1e-4.  For each anode the best DBR pair count is kept.

Writes results/compare_anodes.json and prints a markdown table.
usage: python3 compare_anodes.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze as an  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

ROWS = [  # (label, {emitter_key: [candidate tags]})
    ("ref_real", "강공진 기준(실측 15 nm Ag nk, 두께 최적)",
     {"irppy": ["ref_irppy2acac_korg"], "n15": ["ref_narrow532_15_korg"]}),
    ("ref_ideal", "강공진 기준(이상적 Ag: 벌크 nk, 두께 최적)",
     {"irppy": ["ref_irppy2acac_bulkAg_korg", "ref_irppy2acac_bulkAg_korg_s2"],
      "n15": ["ref_narrow532_15_bulkAg_korg", "ref_narrow532_15_bulkAg_korg_s1"]}),
    ("Ag5", "재순환 + 실측 HATCN/Ag 5 nm",
     {"irppy": ["cmp_irppy_Ag5m_N3", "cmp_irppy_Ag5m_N6"], "n15": ["cmp_n15_Ag5m_N3", "cmp_n15_Ag5m_N6"]}),
    ("Ag7", "재순환 + 실측 HATCN/Ag 7 nm",
     {"irppy": ["cmp_irppy_Ag7m_N3", "cmp_irppy_Ag7m_N6"], "n15": ["cmp_n15_Ag7m_N3", "cmp_n15_Ag7m_N6"]}),
    ("Ag8", "재순환 + 실측 HATCN/Ag 8 nm",
     {"irppy": ["cmp_irppy_Ag8m_N3", "cmp_irppy_Ag8m_N6"], "n15": ["cmp_n15_Ag8m_N3", "cmp_n15_Ag8m_N6"]}),
    ("Ag_ideal", "재순환 + 이상적 thin Ag(벌크 nk, 두께 최적)",
     {"irppy": ["cmp_irppy_Ag_N3", "cmp_irppy_Ag_N6", "cmp_irppy_Ag_N6p", "cmp_irppy_Ag_N8"],
      "n15": ["cmp_n15_Ag_N3", "cmp_n15_Ag_N6", "cmp_n15_Ag_N6p", "cmp_n15_Ag_N8"]}),
    ("IZO", "재순환 + 스퍼터 IZO(HATCN 버퍼)",
     {"irppy": ["cmp_irppy_IZO_N3", "cmp_irppy_IZO_N6", "cmp_irppy_IZO_N6p", "cmp_irppy_IZO_N8"],
      "n15": ["prop6_narrow532_15_IZOtop_lowk_korg_tco15", "cmp_n15_IZO_N8"]}),
]


def best_of(tags):
    best = None
    for t in tags:
        d = an.load(t)
        if d is None:
            continue
        if best is None or d["final"]["CE"] > best[1]["final"]["CE"]:
            best = (t, d)
    return best


def detail(tag, d):
    d = dict(d); d["tag"] = tag
    inf = an.evaluate(d)
    x, y = an.cie_xy(inf["_ctx"], an.on_axis_spectrum(inf))
    out = dict(tag=tag, CE=inf["CE"], eqe=inf["eqe"], eqe20=inf["eqe_within_20deg"], cie=(x, y),
               N=d.get("N"), d_ag=d["params"].get("d_ag"), preset=d.get("preset"))
    if "Aprime_eff" in inf:
        out.update(Aprime=inf["Aprime_eff"], p=inf["p_avg"], top_anode_abs=inf["Aprime_lambertian_layers"].get("ag_top"))
    return out


if __name__ == "__main__":
    res = {}
    for key, lab, cand in ROWS:
        res[key] = dict(label=lab)
        for em, tags in cand.items():
            b = best_of(tags)
            res[key][em] = detail(*b) if b else None
            if b:
                r = res[key][em]
                print(key, em, r["tag"], round(r["CE"], 1), flush=True)
    with open(os.path.join(HERE, "results", "compare_anodes.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    ref = {em: res["ref_real"][em]["CE"] for em in ("irppy", "n15")}

    def cell(r, em):
        if r is None:
            return "–"
        s = f"**{r['CE']:.0f}** (×{r['CE'] / ref[em]:.2f}), EQE {100 * r['eqe']:.0f}%, ±20° {100 * r['eqe20']:.1f}%"
        s += f", CIE ({r['cie'][0]:.3f}, {r['cie'][1]:.3f})"
        if r.get("N") and not r["tag"].startswith("ref"):
            s += f", DBR {r['N']}쌍"
        if r["tag"].startswith(("cmp_irppy_Ag_", "cmp_n15_Ag_", "ref_")) and r.get("d_ag"):
            s += f", Ag {r['d_ag']:.1f} nm"
        if r.get("Aprime") is not None:
            s += f", A′ {100 * r['Aprime']:.1f}% (양극 {100 * (r['top_anode_abs'] or 0):.1f}%p)"
        return s
    print("\n| 구조 | Ir(ppy)₂acac | 532/15 nm |\n|---|---|---|")
    for key, lab, _ in ROWS:
        print(f"| {lab} | {cell(res[key]['irppy'], 'irppy')} | {cell(res[key]['n15'], 'n15')} |")
