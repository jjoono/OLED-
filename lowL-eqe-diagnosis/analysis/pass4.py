# -*- coding: utf-8 -*-
"""Pass 3: PD-voltage-based noise floor + point-level EQE error, full reliable set."""
import os
import numpy as np
import pandas as pd
from parse_all import parse_file, col

# 원본 eval 데이터 폴더. OLED_EVAL_DIR 환경변수로 덮어쓸 수 있다.
DATA = os.environ.get(
    "OLED_EVAL_DIR", r"C:\Users\Junho\Dropbox\개인자료\HCNB\Experiments\eval"
)
d = os.path.dirname(os.path.abspath(__file__))

mdf = pd.read_csv(os.path.join(d, "metrics.csv"))
rel = mdf[(mdf.Lmax >= 500) & (mdf.EQE_100.notna()) & (mdf.EQE_100 > 0.05) & (mdf.EQE_100 < 40)]
files = rel.file.tolist()
print(f"analyzing {len(files)} pixels")

def get_gain(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            head = f.read(600)
        import re
        m = re.search(r"PD gain: *([0-9.]+) *dB", head)
        return float(m.group(1)) if m else float("nan")
    except Exception:
        return float("nan")

pts, pix = [], []
for k, fn in enumerate(files):
    res = parse_file(os.path.join(DATA, fn))
    if res is None:
        continue
    arr, cols = res
    V = col(arr, cols, "Voltage"); PD = col(arr, cols, "PD Voltage")
    J = col(arr, cols, "Current Density")
    if J is None: J = col(arr, cols, "Abs. Current Density")
    L = col(arr, cols, "Luminance")
    E = col(arr, cols, "EQE")
    if any(x is None for x in (V, PD, J, L, E)):
        continue
    with np.errstate(all="ignore"):
        imax = int(np.nanargmax(L))
        V, PD, J, L, E = V[:imax+1], PD[:imax+1], J[:imax+1], L[:imax+1], E[:imax+1]
        # L per PD volt scale from bright points
        bright = (PD > 0.002) & np.isfinite(L) & (L > 0)
        if bright.sum() < 2:
            continue
        LperPD = float(np.nanmedian(L[bright] / PD[bright]))
        # PD noise from off region
        off = np.abs(J) < 1e-3
        if off.sum() < 4:
            off = V < 1.0
        if off.sum() < 4:
            continue
        pdnoise = float(np.nanstd(PD[off]))
        pdoffset = float(np.nanmedian(PD[off]))
        noiseL = pdnoise * LperPD
        offL = pdoffset * LperPD
        # plateau EQE
        m = (L >= 200) & (L <= 2000) & np.isfinite(E) & (E > 0)
        if m.sum() < 2:
            continue
        Ep = float(np.nanmedian(E[m]))
        if not (0.05 < Ep < 40):
            continue
        n_1_10 = int(((L >= 1) & (L < 10)).sum())
        n_10_100 = int(((L >= 10) & (L < 100)).sum())
        sel = (L > 0.3) & (L <= 300) & np.isfinite(E) & (np.abs(J) > 1e-4)
        for i in np.where(sel)[0]:
            pts.append((fn, L[i], E[i]/Ep, noiseL, J[i], V[i]))
        pix.append((fn, noiseL, offL, pdnoise, LperPD, Ep, n_1_10, n_10_100, get_gain(os.path.join(DATA, fn))))
    if (k+1) % 3000 == 0:
        print(f"{k+1}", flush=True)

pd.DataFrame(pts, columns=["file","L","Erel","noiseL","J","V"]).to_csv(os.path.join(d,"points4.csv"), index=False)
pd.DataFrame(pix, columns=["file","noiseL","offL","pdnoise","LperPD","Eplateau","n_1_10","n_10_100","gain"]).to_csv(os.path.join(d,"pixels4.csv"), index=False)
print(f"done: {len(pts)} points, {len(pix)} pixels")
