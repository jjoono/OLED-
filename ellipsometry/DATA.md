# Data the ellipsometry scripts need

Two categories, handled differently.

## 1. Tracked in this repository

The 2026-08-20 absolute T/R campaign and the optical-constant library are
committed under `dft_seedlayer_screen/data/`. The scripts locate them
automatically, relative to their own file; no environment variable is required.

```
dft_seedlayer_screen/data/
  TR_20260820/
    raw/{id}T.csv, {id}R.csv     28 untouched Cary 6000i / UMA exports
    raw/glass.csv, raw/glassR.csv  bare-substrate reference, same campaign
    ALL_SAMPLES_TRA.csv          consolidated T, Rmeas, Rcorr, A  (see below)
    HATCN5_Ag5_TR.txt, MoOx5_Ag5_TR.txt
  nk/
    Ag_McPeak.csv, Ag_Palik.csv  literature silver
    Ag{5,7,8}nm_on_HATCN5_measured.csv   n,k inverted from this campaign's T/R
    l_HATCN.csv, MoO3.csv, ITO.csv, IZO.csv, NPB.csv, Alq3.csv, TPBi.csv, ...
```

Overrides, if the data lives elsewhere:

```bash
export TR20260820_DIR=/path/to/TR_20260820
export NK_DIR=/path/to/nk
```

### Rebuilding `ALL_SAMPLES_TRA.csv`

`scripts/tr_crosscheck/build_tra.py` regenerates it from `raw/`. It applies the
back-surface reflection correction described in the handoff (the UMA collects
only ~84 % of the substrate's back-surface beam at 6°, so measured R is low by
the uncollected part, weighted by the sample's own transmittance squared).

```bash
python scripts/tr_crosscheck/build_tra.py --verify   # rebuild to a scratch file and compare
python scripts/tr_crosscheck/build_tra.py            # overwrite the consolidated CSV
```

The original consolidation script was never committed, so the correction
constant was recovered by least squares against the delivered file. The rebuild
matches it to **max 0.012 %p across every column** — far inside the campaign's
σ(A) = 0.15 %p noise — and the implied collection factor, 0.837, sits inside the
handoff's 0.843 ± 0.024. `--verify` re-checks this and never touches the
delivered file.

The bare-substrate exports `raw/glass.csv` / `raw/glassR.csv` check that
constant directly, without reference to the delivered file: a transparent slab
must show zero absorptance, and over 450–700 nm the bare glass gives
**+0.581 %p with the raw R and −0.007 %p with the corrected R**
(`scripts/joint_nk_260910/jnk_glass.py` prints it). So f is right and the
substrate is lossless over the working window.

Two exports are missing at source and stay blank: **1-9 has no T**,
**2-12 has no R**.

### Which `HATCN` entry to use

`nk/HATCN.csv` and `HATCN_NIR` both have k rising toward the red, which a
wide-gap organic cannot do. Use **`l_HATCN.csv`** (k ≡ 0 over 430–830 nm,
n = 1.849 at 550 nm — which the ellipsometry independently reproduces as 1.867).

## 2. Not in this repository

Kept out deliberately: raw ellipsometer exports, fitted results, figures.
`.gitignore` blocks `ellipsometry/**/*.{xlsx,csv,mod,png,json,npz,SE}`.

```bash
export ELLIPS_DATA=/path/to/measurements   # .SE xlsx exports, CompleteEASE .mod templates
export ELLIPS_OUT=/path/to/output          # fits (.json/.npz), n,k CSV, figures, generated .mod
```

| What | Used by | Notes |
|---|---|---|
| `summary.xlsx` (260819 ThinAg) | `thin_ag_260819/*`, `joint_nk_260910/*` | 16 sheets, 5 angles, 245.8–1688.1 nm |
| `se추출.xlsx` (260813 ITO/IZO) | `ito_izo_v2/*` | 5 sheets, 3 angles, 192.4–1688.3 nm |
| `*_GenOSC_*.mod` | `completeease_mod/*` | CompleteEASE model file used as a structural template |
| `ag_seed_result.json`, `ag_v3_result.json`, `ce_v13_result.json` | most fit scripts | written by the pipeline; regenerate by running it |
| `k_SI_digitized.csv` | `k_anchored.py`, `final_window.py`, `err_budget.py` | perovskite absorption digitised from the source paper's SI; not redistributed |

### Order to regenerate the fits

```bash
# thin Ag, 260819
python scripts/thin_ag_260819/ag_load.py      # builds the npz cache, prints depolarization
python scripts/thin_ag_260819/ag_seed.py      # bare-seed dispersion + thickness  -> ag_seed_result.json
python scripts/thin_ag_260819/ag_fit.py       # first pass, 3 oscillators          -> ag_fit_result.json
python scripts/thin_ag_260819/ag_polish.py    # cross-seeded local-minimum check    -> ag_polish_result.json
python scripts/thin_ag_260819/ag_final.py     # 4 oscillators, 260-1080 nm          -> ag_final_result.json
python scripts/thin_ag_260819/ag_v3.py        # seed scan + final refit             -> ag_v3_result.json
python scripts/completeease_mod/make_mod_ag.py
python scripts/thin_ag_260819/ag_roundtrip.py # verify the .mod reproduces the fit

# sputtered ITO / IZO, 260813
python scripts/lib/ce_mat.py                  # decode SI_JAW / NTVE_JAW from a .mod -> ce_mat.npz
python scripts/ito_izo_v2/ce_v12.py           # geometry scan + fit
python scripts/ito_izo_v2/ce_v13.py           # lock unidentifiable parameters
python scripts/completeease_mod/make_mod_v13.py
```

## 3. A note on seed thickness

Both seeds were deposited to a **quartz-monitor reading of 5 nm**. The bare-seed
samples fit to 6.4 nm (HATCN) and 6.9–7.0 nm (MoOx); the Ag-covered samples need
about 1 nm more still (7.5 / 8.0 nm). Forcing 5 nm makes the fit 5.5× (HATCN)
and 12× (MoOx) worse and demands n = 2.16–2.18, impossible for HATCN. The
`SPEC` table in `thin_ag_260819/ag_load.py` records the **QCM** value, 5 nm;
the fitted values are what the pipeline uses. This gap is the reason to check
the evaporator's tooling factor.
