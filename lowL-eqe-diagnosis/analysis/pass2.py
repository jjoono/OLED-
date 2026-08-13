# -*- coding: utf-8 -*-
"""Pass 2: point-level low-luminance EQE error analysis over reliable pixels."""
import os, csv
import numpy as np
from parse_all import parse_file, col

# 원본 eval 데이터 폴더. OLED_EVAL_DIR 환경변수로 덮어쓸 수 있다.
DATA = os.environ.get(
    "OLED_EVAL_DIR", r"C:\Users\Junho\Dropbox\개인자료\HCNB\Experiments\eval"
)
d = os.path.dirname(os.path.abspath(__file__))

import pandas as pd
mdf = pd.read_csv(os.path.join(d, "metrics.csv"))
rel = mdf[(mdf.Lmax >= 500) & (mdf.EQE_100.notna()) & (mdf.EQE_100 > 0.05) & (mdf.EQE_100 < 40)]
files = rel.file.tolist()
print(f"analyzing {len(files)} reliable pixels")

pts = []       # per-point records
pix = []       # per-pixel records
for k, fn in enumerate(files):
    res = parse_file(os.path.join(DATA, fn))
    if res is None:
        continue
    arr, cols = res
    V = col(arr, cols, "Voltage"); PD = col(arr, cols, "PD Voltage")
    J = col(arr, cols, "Current Density"); L = col(arr, cols, "Luminance")
    E = col(arr, cols, "EQE")
    if V is None or L is None or E is None or J is None:
        continue
    with np.errstate(all="ignore"):
        # forward sweep only: stop at index of Lmax
        imax = int(np.nanargmax(L))
        V, PD, J, L, E = V[:imax+1], PD[:imax+1], J[:imax+1], L[:imax+1], E[:imax+1]
        # noise floor: std of L where luminance signal absent (|L|<5 and J<1e-3)
        off = (np.abs(J) < 1e-3)
        if off.sum() >= 4:
            noiseL = float(np.nanstd(L[off]))
            offsetL = float(np.nanmedian(L[off]))
        else:
            noiseL = np.nan; offsetL = np.nan
        # plateau EQE: median EQE where L in [200, 2000]
        m = (L >= 200) & (L <= 2000) & np.isfinite(E) & (E > 0)
        if m.sum() < 2:
            continue
        Ep = float(np.nanmedian(E[m]))
        if not (0.05 < Ep < 40):
            continue
        # count points per luminance decade (real signal region)
        n_1_10 = int(((L >= 1) & (L < 10)).sum())
        n_10_100 = int(((L >= 10) & (L < 100)).sum())
        # low-L points: 0.5 <= L <= 150, device conducting (J>1e-4)
        sel = (L > 0.3) & (L <= 150) & np.isfinite(E) & (np.abs(J) > 1e-4)
        for i in np.where(sel)[0]:
            pts.append((L[i], E[i]/Ep, noiseL, J[i], V[i]))
        # negative-L points while conducting (pure noise evidence)
        n_negL_conducting = int(((L < 0) & (np.abs(J) > 1e-3)).sum())
        pix.append((fn, noiseL, offsetL, Ep, n_1_10, n_10_100, n_negL_conducting))
    if (k+1) % 2000 == 0:
        print(f"{k+1}/{len(files)}", flush=True)

pd.DataFrame(pts, columns=["L","Erel","noiseL","J","V"]).to_csv(os.path.join(d,"points.csv"), index=False)
pd.DataFrame(pix, columns=["file","noiseL","offsetL","Eplateau","n_1_10","n_10_100","n_negL"]).to_csv(os.path.join(d,"pixels.csv"), index=False)
print(f"done: {len(pts)} low-L points from {len(pix)} pixels")
