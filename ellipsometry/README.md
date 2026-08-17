# Ellipsometry analysis — n,k extraction

Spectroscopic-ellipsometry (VASE) analysis for thin films measured in 2026:
perovskite/GraHIL, IZO, Ag-on-HATCN, and sputtered ITO / IZO.

## Method

All samples are analysed with the same validated pipeline:

1. **NCS-space residuals** (`N=cos2Ψ, C=sin2Ψ·cosΔ, S=sin2Ψ·sinΔ`) — bounded, immune
   to the Ψ→90° / Δ-wrap pathologies that break Ψ/Δ-space fitting.
2. **Depolarization weighting** — measured `% Depolarization` down-weights unreliable
   points (Si turns transparent >1107 nm → backside reflection).
3. **Geometry scan** — chain-consistency residual over a (thickness, roughness) grid
   locates the true geometry and exposes n·d alias branches.
4. **Model-free chain extraction** — per-wavelength (n,k) from the multi-angle data with
   a continuity penalty; no dispersion-model assumption.
5. **KK-consistent Gen-Osc** — Drude + Tauc-Lorentz/Gaussian with a singularity-subtracted
   Kramers-Kronig integral (validated against the analytic Lorentz to 2e-3).

Independent anchors are used whenever available: SEM thickness, UV-Vis absorption,
and glass-baseline transmittance.

## Usage

Measurement data and fitted results are not in this repository. Point the scripts at
your own copies:

```bash
export ELLIPS_DATA=/path/to/measurements   # xlsx exports, .mod/.mat files
export ELLIPS_OUT=/path/to/output          # figures, CSV, intermediate .npz
python scripts/ito_all5.py
```

Requires `numpy`, `scipy`, `matplotlib`, `openpyxl`.

## Layout

```
scripts/   analysis code (see table below)
```

### Core

| File | Purpose |
|---|---|
| `ellipsometry_fit.py` | Shared library: TMM, Bruggeman EMA, KK operators, Gen-Osc, loaders |
| `kkcheck.py` | Validates the KK integral against the analytic Lorentz oscillator |

### Per-material

| Material | Scripts | Result |
|---|---|---|
| Perovskite on GraHIL | `perovfit*.py`, `k_anchored.py`, `final_window.py`, `err_budget.py` | `Perov_FINAL_nk.csv` |
| GraHIL (PEDOT:PSS) | (in `ellipsometry_fit.py` presets) | `GraHIL_genosc_nk.csv` |
| IZO (O2-annealed) | `izo_*.py` | `IZO_O2_nk.csv` |
| Ag on HATCN | `hatcn_*.py`, `ag_*.py` | `AgHATCN_nk.csv` |
| sputtered ITO / IZO | `ito_*.py`, `substrate_sens.py` | `ITO_IZO_all5_nk.csv` |

## Data-file formats

`.SE` files (J.A. Woollam CompleteEASE measurement files) are **proprietary-obfuscated**:
a fixed-position keystream, no readable tokens, no deflate blocks. They cannot be parsed
without CompleteEASE — export Ψ/Δ (+ depolarization) to xlsx instead.

What *is* parseable:

- **`.mod` / `.mat`** (CompleteEASE model and material files) — plain text with
  base64+gzip arrays. Decoders live in `izo_extract2.py` (JAW substrate tables:
  big-endian float32, wavelength in Ångström) and `hatcn_mat.py` (uniaxial B-spline
  material file: `spline_e2(E)` nodes + `E Inf`).
- **xlsx exports** — `load_data_xlsx()` in `ellipsometry_fit.py` auto-detects the header
  row and both column layouts (interleaved `wl|Ψ Δ|Ψ Δ` and grouped `wl|Ψ Ψ Ψ|Δ Δ Δ`);
  `depol_weights()` builds fit weights from the depolarization block.

## Notes on reliability

Each CSV carries the residual (and, where available, an uncertainty band). Values quoted
as `0.0000` for k mean *below the detection limit*, not literally zero — the limit is
~0.01 for these films, established by forcing k to fixed values and watching the residual.
