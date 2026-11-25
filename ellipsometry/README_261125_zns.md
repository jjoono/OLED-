# Session 2026-11-25 — ZnS on Si, two pieces

`scripts/zns_261125/` fits `ZnS_fitting_Sanghee.xlsx`: two Si pieces with ZnS,
stated to be the same thickness (expected 100–250 nm), 65 / 70 / 75°,
192.4–1688.3 nm. Fitted 205–1000 nm (beyond that the wafer is transparent).

```bash
export ZNS_XLSX=/path/to/ZnS_fitting_Sanghee.xlsx
export ELLIPS_DATA=/path/to/data      # the lab .mod with SI_JAW / NTVE_JAW
export ELLIPS_OUT=/path/to/output
cd scripts/zns_261125
python zn_load.py       # inventory, depolarization, difference between the pieces
python zn_fit.py        # thickness scan, then full-range fits and variants -> zns_fits.json
python zn_report.py     # -> $ELLIPS_OUT/zns/ZnS_nk.csv, zns_summary.png
```

Model: air / roughness (ZnS + void 50:50) / ZnS / native oxide 1.5 nm
(NTVE_JAW) / Si (SI_JAW). ZnS is a Kramers-Kronig-consistent Gaussian basis
from 3.0 eV up with a UV pole and ε∞ (`gp_model.KKBasis`).

## Result

| | Sheet1 | Sheet2 |
|---|---|---|
| ZnS thickness, homogeneous film | 95.96 ± 0.04 nm | 93.76 ± 0.04 nm |
| ZnS thickness, graded film | 97.02 ± 0.04 nm | 94.81 ± 0.04 nm |
| roughness layer | 1.0–1.7 nm | 1.1–1.8 nm |
| MSE, homogeneous / graded | 10.8 / 7.5 | 10.9 / 7.4 |

n,k (depth average, both pieces with one dispersion): 2.539 + 0.007i at 400 nm,
2.445 at 450, 2.363 at 550, 2.330 at 633, 2.296 at 800, 2.278 at 1000 nm;
absorption edge near 340–360 nm (k = 0.17 at 340, 0.06 at 360 nm).

- **The thickness is about 95 nm, just under the expected range.** In the
  transparent range every thickness from 40 to 400 nm was imposed in turn with
  the index refitted: one minimum per piece, at 95 and 93 nm (MSE 5); the next
  best, near 265 and 308 nm, are at MSE 370–390 and need n(550) = 2.5–2.6.
  Nothing inside 100–250 nm fits better than MSE 30.
- **The two pieces are not the same thickness**: 2.2 nm apart (2.3 %), fifty
  times the statistical interval, with either depth profile. They do share one
  dispersion (fitting them separately gives the same n,k and MSE).
- **The film is slightly denser at the bottom.** A linear gradient of 3 % in
  (ε − 1), top below bottom, lowers the MSE from 10.8 to 7.5; a denser 25 nm
  at the substrate (+9 %) does the same (5.9). The total thickness stays at
  96–97 nm (Sheet1) for every profile tried.
- The absolute thickness depends on the model at the ±1.5 nm level: 96.0 nm
  homogeneous, 97.0 graded, and if the native oxide is left free it trades
  with the film (oxide 7 nm, ZnS 91.5 nm, sum unchanged). 1.5 nm of oxide is
  assumed.
- n(550) = 2.36 is 1 % below bulk cubic ZnS (2.38–2.39): a dense film.
- k below ~0.005 (above ~420 nm) is not determined: the basis has no
  absorption below 3.0 eV, so k = 0 there by construction.
- Depolarization reaches 12–15 % at 300–400 nm, where Ψ and Δ swing fastest.
  A lateral thickness spread of 2.8 % reproduces it without changing the
  thickness or the index.
