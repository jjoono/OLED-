# Session 2026-08-20/21 — TMM fix, CompleteEASE `.mod` generation, thin Ag on seed layers

This document records the analysis process added on branch
`claude/ellipsometry-tmm-fix-ag-seed-260820`. Measurement data and fitted
results are **not** in the repository; every script reads its inputs from
`ELLIPS_DATA` and writes to `ELLIPS_OUT`.

```bash
export ELLIPS_DATA=/path/to/measurements     # xlsx exports of .SE files, CompleteEASE .mod files
export ELLIPS_OUT=/path/to/output            # fits (.json/.npz), n,k CSV, figures, generated .mod
```

The 2026-08-20 absolute T/R campaign and the n,k library ARE tracked, under
`dft_seedlayer_screen/data/`; the T/R and device-optics scripts find them
automatically. **See [DATA.md](DATA.md)** for what every script needs, how to
rebuild `ALL_SAMPLES_TRA.csv` from the raw Cary exports, and the order to
regenerate the fits.

---

## 1. The transfer-matrix bug (read this first)

`ellipsometry_fit._tmm` (and `_tmm_uniaxial`) had the two propagation-phase
exponents swapped. For the convention `(E0+,E0-) = M (EN+,EN-)` the forward
wave must pick up `exp(-i·k0·q·d)` travelling back up; the old code amplified
the *backward* wave inside absorbing layers.

Why it hid: ellipsometry only uses the **ratio** `rp/rs`, and the fit simply
absorbed the error into n and k. Ψ/Δ fits looked excellent (chain residual
0.003) while the extracted constants were biased. It was exposed by computing an
**absolute** quantity — the modelled reflectance exceeded 1 in the UV.

Regression tests (`scripts/lib/tmm_fix.py`, run as `__main__`):

| test | buggy | fixed |
|---|---|---|
| 200 nm Ag on glass, 0°, must equal bulk Fresnel R = 0.9685 | 1.0326 | 0.9685 |
| bare Si vs analytic Fresnel | passes either way (no finite layer) | passes |
| single film vs Airy formula, `Im(q) > 0` | matched the *wrong* sign variant | matches |

The fix swaps the exponents **and** conjugates the returned `rp, rs` so Δ keeps
the instrument's sign convention. Impact on the sputtered-ITO batch: n moved
< 1 %, **k was inflated ~20× (0.022 → 0.001 at 550 nm), thickness biased −10 %,
MSE 22 → 2.3.** Everything computed before 2026-08-20 with this library needs
re-checking, absorbing films most of all. The pre-fix library is kept as
`scripts/archive/ellipsometry_fit_prebugfix.py` for reproducing old numbers.

## 2. Writing CompleteEASE `.mod` files from Python

`scripts/completeease_mod/` generates Gen-Osc model files that CompleteEASE
loads directly (`make_mod_v13.py` — single Gen-Osc layer; `make_mod_ag.py` —
four-layer substrate/oxide/seed/metal stack). Format facts that each cost a
failed load, so they are documented in the code:

* UTF-8 **with BOM**, CRLF, tab-delimited; leading tabs are indentation, not
  fields. Thicknesses in **Å**.
* `start_Gen-Osc Fit Parms` line 0 = `<oscillator count>` TAB `<Einf>` … .
* **Tauc-Lorentz** = 4 parameter lines **plus one trailing `F` line**
  (the "Common Eg" flag). **Gaussian / Lorentz** = 4 lines `Amp, iAmp, Br, En`
  (iAmp is hidden in the GUI). Drude(RT) = `Resistivity (Ohm·cm)`, `Scat. Time (fs)`.
* `Gen-Osc Grade Parms` holds exactly one 4-line group per oscillator
  **parameter**; a mismatch gives "no model loaded".
* The `lo`/`hi` fields on each parameter line are enforced fit limits — use
  them to keep an oscillator from escaping to a degenerate solution
  (a Gaussian sliding to 0.3 eV, outside the data, produced ±458 error bars).
* Library materials can be **decoded out of a .mod**: `Mat Table` arrays are
  base64 of gzip of big-endian float32, wavelength in Å (`scripts/lib/ce_mat.py`).
  `SI_JAW` = Herzinger Si; `NTVE_JAW` has n = 1.734 at 633 nm, not 1.46.
* Add a layer by duplicating the whole `start_LayerN … end_LayerN` block,
  renumbering it and its `'Thickness # N'`, and bumping the layer count on the
  first `Model Parms` line.

Every generated file is **round-trip verified** (`check_mod.py`,
`thin_ag_260819/ag_roundtrip.py`): the `.mod` is parsed back, the stack
rebuilt from what the file says, and its MSE against the data must match the
fit.

## 3. Sputtered ITO / IZO, second pass (`scripts/ito_izo_v2/`)

Oscillator forward model in CompleteEASE's own conventions
(`scripts/lib/ce_osc.py`, analytic ε₁ for Drude, Tauc-Lorentz and Gaussian,
validated against brute-force Kramers-Kronig to ≤ 6e-4). Analytic ε₁ matters:
a truncated KK integral dumps missing high-energy weight into Einf, so the Einf
written to a `.mod` would not mean what CompleteEASE's Einf means.

Process (`ce_v11.py` → `ce_v12.py` → `ce_v13.py`):

1. **Valid window = 340–1080 nm.** Below 340 nm an exhaustive free-(n,k) grid
   — the best any homogeneous layer can do — still leaves residual 0.06–0.09
   versus a 0.003 floor, unchanged for roughness 0–15 nm (`ce_grid.py`).
   Above 1080 nm Si is transparent (1107 nm) and wafer backside reflection
   enters. Fitting the full 192–1688 nm range is what drove the runaway
   negative-Gaussian fits.
2. Geometry scan (d, roughness, angle offset) scored by the model-free chain
   residual (`ce_chain2.py` uses a coarse (n,k) grid initialiser so no
   wavelength is trapped at k = 0).
3. Gen-Osc fit at the scanned geometry; oscillator count chosen by MSE.
4. **Identifiability** (`ce_cov.py`, `ce_lock2.py`): Jacobian covariance and
   correlation matrix; parameters with > 50 % relative CI are locked in the
   `.mod`. Also confine the sub-gap Gaussian to the measured window
   (`ce_gscan.py`).

Angle offset: scan jointly with oxide thickness and roughness (`ce_dth.py`);
if neither can absorb it the offset is real. Locking d and the offset is what
prevents the fit sliding onto the n·d alias branch.

## 4. Thin Ag on HATCN vs MoOx seeds (`scripts/thin_ag_260819/`)

Stack air / roughness (Bruggeman Ag+void) / Ag / seed / NTVE 2 nm / Si,
5 angles, 260–1080 nm.

1. `ag_seed.py` — bare-seed samples: seed dispersion (Einf + Tauc-Lorentz) and
   thickness, scanned against the assumed native-oxide thickness (the two are
   optically interchangeable; only their sum is measured).
2. `ag_v3.py` — pin every Ag thickness at its deposited value and scan the
   seed thickness; the minimum gives the underlayer the *Ag-covered* samples
   require (≈ +1 nm over the bare-seed value). With that seed, free Ag
   thicknesses land within ±10 % of QCM (mean 1.01 ± 0.05).
3. `ag_final.py` — Ag = Einf + Drude + **3 Gaussians** (residuals otherwise pile
   up at 260–350 nm where Ag's interband structure is). Cutting the UV instead
   lets d drift +2.5 nm away from nominal — the UV carries thickness
   information, keep it.
4. `ag_polish.py` — restart every sample from every other sample's solution
   (local-minimum check).
5. Effective Ag density = (ħω_p / ħω_p,ref)² from the Drude plasma energy,
   ref = the thickest film. Reliable to ±0.02 for continuous films; qualitative
   for island films. Cross-checked against Bruggeman in n,k space (`ag_ema.py`).
6. **EMA as the layer model was tested and rejected** (`ag_emafit.py`): holding
   bulk-Ag n,k fixed and fitting only (d, f) gives MSE 9–43 (Bruggeman) or
   65–209 (Maxwell-Garnett) versus 4–9 for a free Gen-Osc layer. Thin Ag is not
   diluted bulk Ag — its damping is larger — and EMA lacks that degree of freedom.
7. Angle-by-angle inversion (`ag_anis.py`) and an (n,k) residual map at 633 nm
   quantify how tightly n is pinned (0.07–0.11 for Ag ≥ 7 nm on HATCN;
   n = 0.28 is 12× worse).

Quartz-monitor mass is reproduced to 3 %: (d_SE × density)/nominal = 0.98 ± 0.03
regardless of the seed assumption, although d alone trades 1:1 with the seed.

## 5. Cross-check against absolute T/R on glass (`scripts/tr_crosscheck/`)

Coherent front stack + incoherent 1 mm substrate (`tr_check2.py`; bare-glass
check T = 91.76 %, R = 8.24 % vs measured 91.66 / ~8.5). With SE n,k fixed,
thickness inverted from measured T+R agrees with SE to 0.3–0.5 nm, and k agrees
to 1–4 % (both at bulk Ag). **n disagrees for Ag ≥ 7 nm** (SE 0.07–0.10 vs T/R
0.23–0.28 at 550 nm): SE is specular, a 2-observable T/R inversion counts any
non-specular light as absorption and cannot detect it. Use T/R values for
conservative device-loss estimates, SE values as the material constants; an
off-specular scan or SE on the same glass pieces resolves it.

## 6. Device optics (`scripts/device_optics/`)

Angle-resolved single-pass electrode loss for s and p, Lambertian-weighted.
`dev_sub.py` is the correct stack (organic 1.8 / seed / Ag / capping / **matched
n = 1.8 substrate** with an MLA above): the 50°, s-polarised 15 % peak seen with
air above the capping is a TE guided mode of the capping and disappears with
the matched substrate; the profile is then flat at ≈ 2 % per pass. Capping
index (n 2.1 → 2.7) still cuts the angle-averaged loss ~40 %. A thick Ag mirror
loses 1.9 % per bounce — no better than one pass through the 8 nm electrode.

## 7. Layout

```
scripts/
  ellipsometry_fit.py        shared library (TMM FIXED 2026-08-20), EMA, KK, loaders
  lib/       tmm_fix.py      corrected TMM + the three regression tests
             ce_osc.py       CompleteEASE-convention oscillators with analytic eps1
             ce_mat.py       decode SI_JAW / NTVE_JAW tables out of a .mod
  completeease_mod/          make_mod_v13.py (1 Gen-Osc layer), make_mod_ag.py (4 layers), check_mod.py
             archive/        earlier generator versions (superseded)
  ito_izo_v2/                sputtered ITO/IZO pipeline, section 3
  thin_ag_260819/            Ag on HATCN / MoOx, section 4
  tr_crosscheck/             build_tra.py (raw Cary exports -> consolidated CSV,
                             documents the back-surface R correction),
                             tr_check.py / tr_check2.py - section 5
  device_optics/             angle-resolved electrode loss, section 6
  archive/                   ellipsometry_fit_prebugfix.py (do not use for new work)
```

Requires `numpy`, `scipy`, `matplotlib`, `openpyxl`.
