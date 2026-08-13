# -*- coding: utf-8 -*-
"""Parse all OLED jvl eval files, compute per-pixel metrics for classification
and low-luminance EQE anomaly analysis."""
import os, sys, math, csv, re
import numpy as np

# 원본 eval 데이터 폴더. OLED_EVAL_DIR 환경변수로 덮어쓸 수 있다.
DATA = os.environ.get(
    "OLED_EVAL_DIR", r"C:\Users\Junho\Dropbox\개인자료\HCNB\Experiments\eval"
)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metrics.csv")

def parse_file(path):
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()
    except Exception:
        return None
    # find measurement data marker
    start = None
    for i, ln in enumerate(lines[:10]):
        if "### Measurement data ###" in ln:
            start = i
            break
    if start is None:
        return None
    header = [h.strip() for h in lines[start+1].split("\t")]
    data_lines = lines[start+3:]
    rows = []
    for ln in data_lines:
        parts = ln.split("\t")
        if len(parts) < 5:
            continue
        try:
            vals = [float(p) for p in parts[:len(header)] if p.strip() != ""]
        except ValueError:
            continue
        rows.append(vals)
    if len(rows) < 5:
        return None
    n = min(len(r) for r in rows)
    arr = np.array([r[:n] for r in rows], dtype=float)
    cols = {h: j for j, h in enumerate(header[:n])}
    return arr, cols

def col(arr, cols, name):
    j = cols.get(name)
    if j is None:
        return None
    return arr[:, j]

def interp_eqe_at_L(L, EQE, target):
    """EQE at target luminance via interpolation in log(L). Only forward sweep,
    increasing region."""
    ok = np.isfinite(L) & np.isfinite(EQE) & (L > 0) & (EQE > 0)
    if ok.sum() < 3:
        return np.nan
    L2, E2 = L[ok], EQE[ok]
    # ensure roughly monotonic: sort by L
    idx = np.argsort(L2)
    L2, E2 = L2[idx], E2[idx]
    if target < L2[0] or target > L2[-1]:
        return np.nan
    return float(np.interp(np.log10(target), np.log10(L2), E2))

def main():
    files = [f for f in os.listdir(DATA)
             if f.endswith("_eval.csv") and "spec" not in f]
    out = open(OUT, "w", newline="", encoding="utf-8")
    w = csv.writer(out)
    w.writerow(["file", "npts", "vmin", "vmax", "vstep_med",
                "Lmax", "Jmax", "PDmax",
                "pd_noise_off",  # std of PD voltage where device off
                "L_per_PDV",     # luminance per PD volt scale
                "noise_L_equiv", # pd_noise_off * L_per_PDV
                "J_leak",        # |J| at ~1V (or lowest fwd V >0.5)
                "J_rev",         # |J| at most negative V
                "EQEmax", "L_at_EQEmax",
                "EQE_1", "EQE_10", "EQE_100", "EQE_1000",
                "lum_per_mW",    # Lmax / powerdensity at Lmax (cd/m2 per mW/cm2) - NIR flag
                ])
    nerr = 0
    for k, fn in enumerate(files):
        res = parse_file(os.path.join(DATA, fn))
        if res is None:
            nerr += 1
            continue
        arr, cols = res
        V = col(arr, cols, "Voltage")
        PD = col(arr, cols, "PD Voltage")
        J = col(arr, cols, "Current Density")
        L = col(arr, cols, "Luminance")
        EQE = col(arr, cols, "EQE")
        PWR = col(arr, cols, "Power Density")
        if V is None or L is None or EQE is None:
            nerr += 1
            continue
        with np.errstate(all="ignore"):
            npts = len(V)
            vmin, vmax = float(np.nanmin(V)), float(np.nanmax(V))
            dv = np.diff(V)
            vstep = float(np.nanmedian(np.abs(dv[dv != 0]))) if len(dv) else np.nan
            Lmax = float(np.nanmax(L)) if L is not None else np.nan
            Jmax = float(np.nanmax(np.abs(J))) if J is not None else np.nan
            PDmax = float(np.nanmax(np.abs(PD))) if PD is not None else np.nan
            # off region: |J| < 1e-4 mA/cm2 or V < 0
            pd_noise = np.nan
            LperPD = np.nan
            noiseL = np.nan
            if PD is not None and J is not None:
                off = np.abs(J) < 1e-4
                if off.sum() >= 4:
                    pd_noise = float(np.nanstd(PD[off]))
                on = (np.abs(PD) > 0.001) & np.isfinite(L)
                if on.sum() >= 3:
                    ratio = L[on] / PD[on]
                    LperPD = float(np.nanmedian(ratio))
                if np.isfinite(pd_noise) and np.isfinite(LperPD):
                    noiseL = pd_noise * LperPD
            # leakage at ~1V forward
            J_leak = np.nan
            if J is not None:
                m = (V >= 0.9) & (V <= 1.6)
                if m.sum():
                    J_leak = float(np.nanmax(np.abs(J[m])))
            J_rev = np.nan
            if J is not None:
                m = V <= -0.5
                if m.sum():
                    J_rev = float(np.nanmax(np.abs(J[m])))
            # EQE metrics: restrict to L > noise*5 region? Use raw for now
            okE = np.isfinite(EQE) & (EQE > 0) & np.isfinite(L) & (L > 0)
            EQEmax = np.nan; L_at = np.nan
            if okE.sum():
                # EQEmax only where L > 0.5 to avoid pure-noise picks... keep raw, filter later
                i = int(np.nanargmax(np.where(okE, EQE, -np.inf)))
                EQEmax = float(EQE[i]); L_at = float(L[i])
            e1 = interp_eqe_at_L(L, EQE, 1.0)
            e10 = interp_eqe_at_L(L, EQE, 10.0)
            e100 = interp_eqe_at_L(L, EQE, 100.0)
            e1000 = interp_eqe_at_L(L, EQE, 1000.0)
            lum_per_mW = np.nan
            if PWR is not None and np.isfinite(Lmax) and Lmax > 0:
                i = int(np.nanargmax(L))
                if np.isfinite(PWR[i]) and PWR[i] > 0:
                    lum_per_mW = float(L[i] / PWR[i])
        w.writerow([fn, npts, vmin, vmax, round(vstep,4) if np.isfinite(vstep) else "",
                    Lmax, Jmax, PDmax, pd_noise, LperPD, noiseL,
                    J_leak, J_rev, EQEmax, L_at, e1, e10, e100, e1000, lum_per_mW])
        if (k+1) % 2000 == 0:
            print(f"{k+1}/{len(files)}", flush=True)
    out.close()
    print(f"done, {len(files)} files, {nerr} parse errors")

if __name__ == "__main__":
    main()
