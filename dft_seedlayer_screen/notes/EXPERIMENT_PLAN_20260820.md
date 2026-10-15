# Measurement plan, 2026-08-20 — Cary UMA G6875A

Samples already deposited (16):

| # | stack |
|---|---|
| 1 | glass / HATCN 5 |
| 2 | glass / MoOx 5 |
| 3–9 | glass / HATCN 5 / Ag 4, 5, 6, 7, 8, 10, 12 |
| 10–16 | glass / MoOx 5 / Ag 4, 5, 6, 7, 8, 10, 12 |

Sheet resistance is already in `data/TR_20260819/rs.csv`.

## What is actually unresolved

The Rs data alone leaves the absorption undetermined over a factor of two,
because it cannot say how much of bulk Ag's eps2 grows with the size effect:

| assumption | A on glass, HATCN/Ag5 | A in the device |
|---|---|---|
| only the Drude part scales (36 %) | 4.4 % | 1.5 % |
| all of eps2 scales (100 %) | 9.2 % | 3.3 % |

Measuring k removes the question entirely — that is the point of tomorrow.

A second, larger risk is not modellable at all: if the film is not a uniform
slab, localized plasmons add absorption no Drude model predicts, and the
Bruggeman effective medium gets the sign wrong near percolation.  Only the
spectral shape and the scatter scan can see it.

## Instrument setup — do once, then do not touch

| | setting |
|---|---|
| range | 350–850 nm |
| mode | **%T and %R — never Abs** (Abs = log10(1/T), a different quantity) |
| SBW | 2 nm |
| scan rate | <= 300 nm/min |
| averaging | 2–3 scans |
| geometry, T | sample 0 deg, detector 180 deg |
| geometry, R | sample 6 deg, detector **12 deg** (detector = 2 x sample) |

1. Run the UMA alignment routine.
2. Baseline: **no sample**, 0 deg / 180 deg → 100 % line.
3. Dark: block the beam → 0 % line.
4. From here on, do not move the accessory or re-align.

The UMA gives absolute %R with no reference mirror — the detector arm rotates
to catch the reflected beam and there are no fold mirrors in the path, so the
direct beam is a valid 100 % reference at every detector angle.

Illuminate from the **film side** for R, and keep it the same for every sample.
T is direction-independent; R is not.

## Step 1 — one sample decides the day (30 min)

Measure **glass / HATCN 5 / Ag 5** only.  Absolute T, then absolute R without
remounting.  Compute A = 1 - T - R at 550 nm.

| A measured | meaning | action |
|---|---|---|
| < 0 anywhere | a baseline is not absolute | stop, fix before anything else |
| 4–5 % | optimistic end; the old 10–13 % was a baseline artefact | continue, story holds |
| 5.5–7.5 % | middle; most likely | continue |
| ~9 % | pessimistic end; defect absorption scales too | continue, revise the device target |
| > 11 % | outside the model — islands or plasmons | check the spectral shape before continuing |

Also measure **bare glass** once: A should be < 0.5 %.  If not, the baseline is
wrong.

## Step 2 — the full series (2–3 h)

All 16 samples, T then R, no remount between them.  Suggested order so that a
drift is visible: HATCN 4 → 12, then MoOx 4 → 12, then the two bare seeds,
then repeat HATCN/Ag5 at the end.  **The repeat is the drift check** — if it
differs from the morning value by more than 0.5 %p, something moved.

## Step 3 — scatter, only the thin ones (30 min)

Sample fixed at 6 deg, sweep the **detector** from 8 to 40 deg, on:
- HATCN 5 / Ag 4
- MoOx 5 / Ag 4

Specular light sits at 12 deg.  Anything outside it is scattered light that a
fixed-angle measurement would have counted as absorption.  Roughness alone
predicts < 0.2 % at rms 1–2 nm, so a large off-specular signal means islands,
not roughness.

## Step 4 — analysis

Drop the CSVs into `data/TR_20260819/` as

```
HATCN5_T.csv  HATCN5_R.csv          wavelength_nm, value
HATCN5_Ag4_T.csv ... MoOx5_Ag12_R.csv
```

then

```
python scripts/83_analyze_TR_series.py
```

It reports A vs thickness for both seeds, one specularity per seed from a
global fit over the whole series, an independent specularity from the Rs
series, and it flags any wavelength where A < 0.

## What tomorrow cannot answer — put these in the next deposition

| sample | what it settles |
|---|---|
| glass / HATCN 5 / Ag 5 / CPL, and / Ag 8 / CPL | validates the glass → device transfer, and pins the CPL's own n, k, d, which is currently the single largest uncertainty (+0.73 / -0.39 %p) |
| glass / HATCN 5 / Ag 2, Ag 3 | locates the percolation threshold; both seeds already conduct at 4 nm, so d_c is currently an extrapolation |
| glass / Ag 6, Ag 8 (no seed) | isolates how much of the seed's benefit is the seed rather than the substrate |

Measure the capped samples from the **CPL side**, where A is ~5 % rather than
~1.8 % and the check has real discriminating power.

## Files to retrieve

- `nk_JH_total.mat` — the measured ITO/IZO n,k.  The Ag-vs-TCO comparison
  currently rests on my estimated k and the ranges overlap, so this decides it.
- the CPL material's identity, thickness and optical constants.
- XRR or a calibrated QCM tooling factor: k carries a ~5 % systematic unless
  the thickness is known to +/-0.3 nm.
