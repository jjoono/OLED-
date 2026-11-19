"""export_amap_xlsx.py -- raw data of the A' x FWHM study (recycling vs strong cavity) to results/amap_rawdata.xlsx"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import amap  # noqa: E402
from plot_phase_origin import crossing  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
R, S = amap.load()
fw = np.array(sorted(R), float)
Ag = np.round(np.arange(0.035, 0.2001, 0.0025), 4)
M = np.array([amap.envelope(R, int(f), Ag) for f in fw])
sR = np.array([S["realAg"][int(f)] for f in fw]); sI = np.array([S["idealAg"][int(f)] for f in fw])
ratio = M / sR[:, None]
cols = [f"FWHM {int(f)} nm" for f in fw]
idx = pd.Index(np.round(100 * Ag, 2), name="A' (%)")

ff = np.arange(10, 80.01, 1.0)
Zf = np.array([PchipInterpolator(fw, ratio[:, j])(ff) for j in range(len(Ag))])
cmp = json.load(open(os.path.join(HERE, "results", "compare_anodes.json")))

raw = []
for f in sorted(R):
    for N, rows in sorted(R[f].items()):
        for r in rows:
            raw.append({"FWHM (nm)": f, "DBR pairs N": N, "extra loss a_extra": r["a_extra"], "A' (%)": 100 * r["Aprime"],
                        "on-axis cd/A": r["CE"], "total EQE": r["eqe"], "p": r["p"], "a_top": r["a_top"], "R_b": r["rb"]})

dev = []
for key, name in (("IZO", "IZO top anode"), ("Ag8", "measured Ag 8 nm"), ("Ag7", "measured Ag 7 nm"),
                  ("Ag5", "measured Ag 5 nm"), ("Ag_ideal", "ideal thin Ag (bulk nk)")):
    for em, emn, fwe in (("n15", "532 nm / 15 nm (Gaussian)", 15), ("irppy", "Ir(ppy)2acac (measured spectrum)", 64)):
        r = cmp[key][em]; ref = cmp["ref_real"][em]["CE"]; refi = cmp["ref_ideal"][em]["CE"]
        dev.append({"top anode": name, "emitter": emn, "FWHM (nm)": fwe, "A' (%)": 100 * r["Aprime"],
                    "recycling cd/A": r["CE"], "total EQE": r["eqe"], "strong cavity (measured Ag) cd/A": ref,
                    "strong cavity (ideal Ag) cd/A": refi, "ratio vs measured-Ag cavity": r["CE"] / ref,
                    "ratio vs ideal-Ag cavity": r["CE"] / refi, "DBR pairs": r["N"], "CIE x": r["cie"][0], "CIE y": r["cie"][1]})

lines = pd.DataFrame({"FWHM (nm)": fw.astype(int),
                      "A' at ratio 1 (vs measured Ag) %": crossing(M, Ag, sR, 1.0),
                      "A' at ratio 1 (vs ideal Ag) %": crossing(M, Ag, sI, 1.0),
                      "A' at ratio 1.5 (vs measured Ag) %": crossing(M, Ag, sR, 1.5),
                      "A' at ratio 2 (vs measured Ag) %": crossing(M, Ag, sR, 2.0)}).replace(np.inf, ">20 (beyond range)")

notes = pd.DataFrame({"note": [
    "Model values (optical simulation), not measurements. Organic k = 1e-4, PLQY = 1, horizontal dipole ratio 0.8.",
    "Recycling platform: IZO top anode on HATCN, ZnS/LiF DBR (3/5/8 pairs, best kept), 800-nm SiO2 ODR, 15-nm low-k ITO cathode.",
    "A' raised by an angle-independent extra loss per round trip and the stack re-optimised at every point.",
    "Emitters for the map: Gaussian in photon energy, peak 530 nm, no vibronic shoulder.",
    "Strong cavity: Ag thickness, organic thicknesses and capping layer optimised per FWHM (measured 15-nm Ag nk or bulk McPeak nk).",
    "Sheets ratio_grid / recycling_cdA_grid: A' grid every 0.25 %, FWHM at the 8 computed values (no interpolation in FWHM).",
    "Sheet ratio_grid_fine_FWHM: same, PCHIP-interpolated in FWHM every 1 nm (what the contour figure uses).",
    "Below A' ~3.3-3.5 % the 8-pair design cannot reach; the grid starts at 3.5 %. The IZO devices sit at A' 3.2-3.4 %.",
    "Ir(ppy)2acac devices use the measured spectrum, so they do not lie exactly on the Gaussian map (map overestimates by up to ~20 %).",
]})

out = os.path.join(HERE, "results", "amap_rawdata.xlsx")
with pd.ExcelWriter(out, engine="openpyxl") as xw:
    notes.to_excel(xw, sheet_name="README", index=False)
    pd.DataFrame(ratio.T, index=idx, columns=cols).to_excel(xw, sheet_name="ratio_grid")
    pd.DataFrame(Zf, index=idx, columns=[f"{int(f)}" for f in ff]).to_excel(xw, sheet_name="ratio_grid_fine_FWHM")
    pd.DataFrame(M.T, index=idx, columns=cols).to_excel(xw, sheet_name="recycling_cdA_grid")
    pd.DataFrame({"FWHM (nm)": fw.astype(int), "strong cavity measured Ag cd/A": sR, "strong cavity ideal Ag cd/A": sI}
                 ).to_excel(xw, sheet_name="strong_cavity", index=False)
    lines.to_excel(xw, sheet_name="ratio_lines", index=False)
    pd.DataFrame(dev).to_excel(xw, sheet_name="devices", index=False)
    pd.DataFrame(raw).to_excel(xw, sheet_name="raw_optimisation_points", index=False)
    for ws in xw.book.worksheets:
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = max(10, min(45, max(len(str(c.value or "")) for c in col) + 2))
print(out)
