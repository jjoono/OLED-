"""Consolidate the 2026-08-20 T/R campaign into one wavelength-resolved table.

Writes data/TR_20260820/ALL_SAMPLES_TRA.csv with, for every sample and every
wavelength: absolute T as measured, R as the instrument reported it, R after
the back-surface correction, and the resulting absorptance.

The reflection accessory collects only 84.3 % of the substrate's back-surface
beam -- established on bare glass, where assuming that value drives the fitted
absorptance across 450-700 nm to +0.01 +/- 0.15 %p, i.e. to zero, while
leaving physical Fe3+ absorption below 400 nm and an Fe2+ tail beyond 750 nm.
"""
import csv, os
import numpy as np

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RAW = os.path.join(BASE, "data", "TR_20260820", "raw")
OUT = os.path.join(BASE, "data", "TR_20260820", "ALL_SAMPLES_TRA.csv")
KEEP = 0.843

SAMPLES = [("1-2","HATCN",4), ("1-3","HATCN",5), ("1-4","HATCN",6),
           ("2-1","HATCN",7), ("2-2","HATCN",8), ("2-3","HATCN",10), ("2-4","HATCN",12),
           ("1-9","MoOx",0),  ("1-10","MoOx",4), ("1-11","MoOx",5), ("1-12","MoOx",6),
           ("2-9","MoOx",7),  ("2-10","MoOx",8), ("2-11","MoOx",10), ("2-12","MoOx",12)]


def read(p):
    out = {}
    if not os.path.exists(p): return out
    for row in csv.reader(open(p)):
        try: out[float(row[0])] = float(row[1])
        except (ValueError, IndexError): pass
    return out


def n_glass(l): return 1.5220 + 3900./l**2


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    data, lams = {}, set()
    for sid, seed, d in SAMPLES:
        T = read(os.path.join(RAW, f"{sid}T.csv"))
        R = read(os.path.join(RAW, f"{sid}R.csv"))
        data[sid] = (T, R)
        lams |= set(T) | set(R)
    lam = sorted((l for l in lams if 350 <= l <= 850), reverse=True)

    cols = ["wavelength_nm"]
    for sid, seed, d in SAMPLES:
        tag = f"{seed}5_Ag{d}" if d else f"{seed}5_bare"
        cols += [f"{tag}_T_pct", f"{tag}_Rmeas_pct", f"{tag}_Rcorr_pct", f"{tag}_A_pct"]

    rows = []
    for l in lam:
        row = [f"{l:.0f}"]
        Rb = ((n_glass(l)-1)/(n_glass(l)+1))**2
        for sid, seed, d in SAMPLES:
            T, R = data[sid]
            t, r = T.get(l), R.get(l)
            if t is None or r is None:
                row += ["", "", "", ""]; continue
            back = (t/100.0)**2 * Rb              # back-surface beam, two film passes
            rc = r + 100*(1-KEEP)*back
            row += [f"{t:.4f}", f"{r:.4f}", f"{rc:.4f}", f"{100-t-rc:.4f}"]
        rows.append(row)

    with open(OUT, "w", newline="") as f:
        w = csv.writer(f); w.writerow(cols); w.writerows(rows)

    have = sum(1 for sid, _, _ in SAMPLES if data[sid][0] and data[sid][1])
    print(f"wrote {os.path.relpath(OUT, BASE)}")
    print(f"  {len(rows)} wavelengths, {lam[-1]:.0f}-{lam[0]:.0f} nm, {len(cols)} columns")
    print(f"  {have}/{len(SAMPLES)} samples complete "
          f"(1-9 has no T, 2-12 has no R)")


if __name__ == "__main__":
    main()
