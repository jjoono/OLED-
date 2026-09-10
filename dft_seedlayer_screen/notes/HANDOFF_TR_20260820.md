# Absolute T/R campaign, 2026-08-20 — handoff for the ellipsometry work

Ultrathin Ag top electrodes for top-emitting OLEDs, seeded with HATCN or MoOx.
Everything below comes from one Cary 6000i + UMA session plus the four-point
probe series taken the day before.

**Data files**

| file | contents |
|---|---|
| `data/TR_20260820/ALL_SAMPLES_TRA.csv` | 350–850 nm at 2 nm, 15 samples: T, R as measured, R corrected, A |
| `data/TR_20260820/raw/*.csv` | the 28 untouched Cary exports |
| `data/nk/Ag{5,7,8}nm_on_HATCN5_measured.csv` | n, k inverted from T/R, 400–800 nm at 1 nm |
| `data/nk/*.csv` | McPeak/Palik Ag, and the lab's HATCN, MoO3, ITO, IZO, organics, from `nk_JH_total.mat` |
| `runs/TR_20260820_measured_nk.json` | per-sample inversion output |

Scripts `88`, `90`, `91`, `92` regenerate all of it.

---

## 1. Samples

Soda-lime glass / seed 5 nm / Ag. Seed thickness fixed at 5 nm throughout.

| id | stack | Rs (Ω/sq) |
|---|---|---|
| 1-2 | HATCN 5 / Ag 4 | 52.2 |
| 1-3 | HATCN 5 / Ag 5 | 23.3 |
| 1-4 | HATCN 5 / Ag 6 | 18.9 |
| 2-1 | HATCN 5 / Ag 7 | 12.7 |
| 2-2 | HATCN 5 / Ag 8 | 9.1 |
| 2-3 | HATCN 5 / Ag 10 | 7.5 |
| 2-4 | HATCN 5 / Ag 12 | 5.3 |
| 1-9 | MoOx 5, bare | — |
| 1-10 | MoOx 5 / Ag 4 | 138.0 |
| 1-11 | MoOx 5 / Ag 5 | 49.1 |
| 1-12 | MoOx 5 / Ag 6 | 32.0 |
| 2-9 | MoOx 5 / Ag 7 | 14.9 |
| 2-10 | MoOx 5 / Ag 8 | 10.8 |
| 2-11 | MoOx 5 / Ag 10 | 9.6 |
| 2-12 | MoOx 5 / Ag 12 | 6.6 |

Two exports are missing: **1-9 has no T** and **2-12 has no R**. Everything
else is complete. Ag thicknesses are nominal QCM values — see §6.

## 2. Measurement conditions

| | |
|---|---|
| instrument | Agilent Cary 6000i with the G6875A UMA |
| absolute T | sample 0°, detector 180° |
| absolute R | sample 6°, detector 12° (detector angle = 2 × sample angle) |
| baseline | direct beam, no sample; dark taken with a Si wafer at the sample position |
| SBW | 8 nm |
| data interval | 2 nm |
| averaging time | 1.0 s → 120 nm/min |
| source changeover | 350 nm; detector and grating changeover 800 nm |
| valid range | **380–800 nm** |

The UMA reports absolute reflectance with no reference mirror: the detector arm
rotates to catch the reflected beam and no fold mirror sits between the baseline
path and the sample path, so the direct beam is a valid 100 % reference at any
detector angle.

**Below 355 nm the data is unusable.** R breaks its trend at 350 nm and both T
and R jump discontinuously at 348 nm — the deuterium/tungsten source changeover.

Noise after settling: **σ(A) = 0.15 %p** over 500–700 nm. Raising SBW from 2 to
8 nm was what bought this; the 700–800 nm region was the worst before, because
the PMT's quantum efficiency collapses there.

## 3. The one correction that matters — R is 0.6 %p low

The accessory collects only part of the substrate's back-surface reflection. At
6° that beam is displaced about 0.14 mm from the front-surface one and its edge
misses the aperture.

Measured on bare glass:

| λ | R theory | R measured | deficit | back-surface beam collected |
|---|---|---|---|---|
| 450 | 8.68 % | 8.16 % | −0.52 | 87.4 % |
| 550 | 8.53 % | 7.87 % | −0.65 | 84.0 % |
| 650 | 8.44 % | 7.63 % | −0.81 | 79.9 % |
| 750 | 8.38 % | 7.78 % | −0.60 | 85.0 % |

**Collection factor 84.3 ± 2.4 %, flat across the visible** — which is itself
evidence the interpretation is right, since a wavelength-dependent cause would
not be flat.

Validation: applying 84.3 % drives the bare-glass absorptance over 450–700 nm to
**+0.01 ± 0.15 %p**, i.e. to zero, which is correct for soda-lime in its
transparent window. It simultaneously leaves the physically expected Fe³⁺ and
Fe²⁺ features intact (§4). Getting three independent things right at once is why
this number is trusted.

`Rcorr` in the consolidated CSV already has this applied. The correction is
smaller for the Ag samples because the back-surface beam must cross the film
twice: 0.44 %p at Ag 5 nm, 0.20 %p at Ag 12 nm.

**T needs no correction.** At normal incidence the internally reflected
components stay collinear, so all of them are collected. Bare glass measured
T = 91.66 % against a lossless-slab prediction of 91.47 %.

## 4. Substrate

Soda-lime, 1 mm. After the R correction the real absorptance is:

| λ | A | what it is |
|---|---|---|
| 350 nm | 7.69 % | Fe³⁺ UV edge |
| 360 | 3.88 % | Fe³⁺ |
| 380 | 1.20 % | Fe³⁺ tail |
| 450–700 | **+0.01 ± 0.15 %** | transparent window |
| 800 | 0.89 % | Fe²⁺ NIR tail (band centred ~1050 nm) |

The internal transmittance used in the inversions is in
`scripts/91_nk_library_series.py` as `GLASS_TI`.

## 5. Results

### 5.1 Absorptance, A = 1 − T − R_corr (%)

| sample | 450 | 500 | 550 | 600 | 650 | 700 |
|---|---|---|---|---|---|---|
| HATCN / Ag 4 | 13.35 | 13.77 | 14.59 | 16.15 | 17.86 | 19.51 |
| HATCN / Ag 5 | 11.81 | 11.31 | 11.17 | 11.73 | 12.34 | 13.00 |
| HATCN / Ag 6 | 11.39 | 10.19 | 9.72 | 9.79 | 10.12 | 10.53 |
| HATCN / Ag 7 | 10.72 | 9.47 | 8.09 | 7.67 | 7.52 | 7.46 |
| HATCN / Ag 8 | 10.31 | 8.56 | 7.32 | 7.05 | 6.92 | 6.98 |
| HATCN / Ag 10 | 10.28 | 8.26 | 7.04 | 6.70 | 6.60 | 6.63 |
| HATCN / Ag 12 | 8.75 | 6.94 | 6.27 | 6.17 | 6.10 | 6.27 |
| MoOx / Ag 4 | 19.62 | 24.92 | 27.63 | 28.26 | 27.33 | 26.10 |
| MoOx / Ag 5 | 20.35 | 26.01 | 27.93 | 27.47 | 26.29 | 25.06 |
| MoOx / Ag 6 | 19.39 | 25.77 | 27.43 | 26.81 | 25.87 | 25.04 |
| MoOx / Ag 7 | 18.35 | 22.69 | 22.79 | 21.74 | 20.00 | 18.48 |
| MoOx / Ag 8 | 16.45 | 19.25 | 19.54 | 20.46 | 21.52 | 22.00 |
| MoOx / Ag 10 | 18.61 | 21.53 | 21.60 | 20.97 | 20.49 | 19.59 |

### 5.2 Extracted optical constants, silver on HATCN

Point-by-point inversion of (T, R_corr) with d fixed at the QCM value and the
seed taken as the lab's lossless `l_HATCN` (n = 1.849 at 550, k ≡ 0). Marching
from long wavelength keeps the solver on the metallic branch. Residuals are zero
to printed precision at every wavelength.

| λ | Ag 5 nm | Ag 7 nm | Ag 8 nm | McPeak bulk |
|---|---|---|---|---|
| 450 | 0.585 + 2.696i | 0.392 + 2.650i | 0.333 + 2.675i | 0.041 + 2.676i |
| 550 | 0.528 + 3.528i | 0.282 + 3.479i | 0.228 + 3.574i | 0.044 + 3.610i |
| 650 | 0.592 + 4.331i | 0.267 + 4.339i | 0.227 + 4.457i | 0.051 + 4.460i |
| 750 | 0.662 + 5.116i | 0.261 + 5.177i | 0.235 + 5.306i | — |

**k matches bulk silver to a few percent at every thickness. The entire
difference is in n**, which runs 5–12× bulk. Since ε₂ = 2nk, n carries all of
the excess absorption. Physically this is bulk-density silver with the Drude
damping raised: a Kramers-Kronig-consistent Drude+Lorentz fit to the Ag 5 nm
spectrum returns ħω_p = 9.37 eV against bulk 9.0–9.2, with ε₁ and ε₂ residuals
of 2.2 % and 4.5 %.

### 5.3 ε₁ — the seed comparison

| d_Ag | HATCN ε₁ | MoOx ε₁ |
|---|---|---|
| 4 | −11.57 | −4.72 |
| 5 | −12.17 | −5.15 |
| 6 | −11.23 | −6.06 |
| 7 | −12.02 | −8.10 |
| 8 | **−12.72** | −10.55 |
| 10 | −11.53 | −8.94 |
| 12 | −12.29 | — |

Bulk silver is −13.03 (McPeak, 550 nm).

**On HATCN, ε₁ sits within 8 % of bulk at every thickness from 4 to 12 nm** — the
film is already a proper metal at 4 nm and thickness does not change its
character. **On MoOx, ε₁ climbs monotonically** and only reaches the metallic
range above 7 nm.

### 5.4 Closure thickness — two independent measurements agree

A closed film obeys ρ = ρ₀ + C/d. Fitting that line to the thicknesses where
closure is certain (d ≥ 7 nm) and reading where the thin points fall back onto
it:

| | line | closes at |
|---|---|---|
| HATCN | ρ = 3.57 + 34.9/d µΩ·cm | **5 nm** (4 nm sits +70 % above) |
| MoOx | ρ = 5.75 + 30.1/d µΩ·cm | **7 nm** (6 nm sits +78 % above) |

The optics say the same thing without any transport model: MoOx absorptance is
flat at 27.6–27.9 % for 4, 5 and 6 nm, then **steps down to 22.8 % at 7 nm**.
HATCN shows no such step — it decreases monotonically from 4 nm.

Fitting Rs ∝ (d − d_c)^−t instead returns d_c = 3.3 nm for HATCN against 2.8 nm
for MoOx, i.e. the wrong ordering. That fit extrapolates outside the data, since
both seeds already conduct at 4 nm, and should not be used.

### 5.5 Device numbers

Organic n = 1.8 / Ag / capping n = 2.1 / air, capping thickness re-optimised per
entry, 550 nm:

| d_Ag | A_device | best cap |
|---|---|---|
| 5 | 3.69 % | 65 nm |
| 7 | 2.76 % | 65 nm |
| **8** | **2.60 %** | 65 nm |
| 10 | 2.64 % | 64 nm |

For scale, McPeak silver at 8 nm in the same stack gives 0.52 %. The measured
electrode is **4.7× above that floor**, and since ε₂/ε₂_bulk = ρ/ρ_bulk holds
here to within 20 %, that factor is exactly the sheet resistance ratio.

## 6. What is still open — and where ellipsometry helps most

### 6.1 Thickness is assumed, not measured

Every n and k above is conditional on the QCM thickness. The inversion has two
observables (T, R) and cannot also solve for d. A ±0.3 nm error moves k by about
5 %, and k·d is **not** invariant here — over d = 4.5 to 5.5 nm, k drifts −20 %
while k·d drifts +19 %, because 4πkd/λ = 0.41 at 5 nm is too large for the
optically-thin limit.

**This is the single most valuable thing ellipsometry can add.** XRR would also
do it. Either fixes the largest systematic in the whole dataset.

### 6.2 Scattering is not separated from absorption

T and R capture only the specular channels; anything scattered out of them is
counted as absorption. With two observables and three unknowns (damping,
scattering, thickness) the problem is underdetermined, and this is not academic:

| model | HATCN / Ag 5 | MoOx / Ag 5 | fits (T,R) |
|---|---|---|---|
| free n, k, no scattering | ε₂ 8.6× bulk | ε₂ 29.4× bulk | exactly |
| bulk ε₁ + damping + scattering | damping 5.9×, scatter 3.4 % | damping 21.4×, scatter 7.3 % | exactly |
| what Rs independently says | ρ 7.3× bulk | ρ 15.4× bulk | — |

Both models reproduce the data exactly. The second agrees better with the
independent sheet resistance for both seeds, which is why the films are read as
continuous-but-scattering rather than as islands.

An earlier reading of MoOx as ε₁ > 0 (islands) was a solver-branch artefact and
has been retracted; with a metallic starting guess the same data gives ε₁ < 0.

**Resolving this needs a direct scattering measurement**: sample fixed at 6°,
detector swept 8°→40°, single wavelength. Specular sits at 12°; anything off it
is scattered light. Roughness alone predicts under 0.2 % at 1–2 nm rms, so a
large off-specular signal would mean islands rather than roughness.

### 6.3 Smaller open items

- Ag 10 nm is anomalous in **both** measurements — ρ rises against 8 nm (7.50 vs
  7.28 µΩ·cm) and so does absorptance. Two independent measurements agreeing on
  the anomaly points at the sample, not the instrument. Worth re-depositing.
- 1-9 (bare MoOx) has no T file and 2-12 (MoOx/Ag 12) has no R file.
- No bare HATCN sample was measured, so the seed's own contribution is taken
  from the library rather than from this session.
- The `HATCN` and `HATCN_NIR` entries in `nk_JH_total.mat` both have k rising
  toward the red, which a wide-gap organic cannot do; **use `l_HATCN`**, which
  has k ≡ 0 over 430–830 nm.

## 7. Where the measurement is fragile

Worth knowing before trusting any repeat run:

- The 100 % baseline must be the **empty beam**. Putting a substrate there makes
  everything relative and silently breaks A = 1 − T − R.
- The dark must be taken **at the sample position** with an opaque piece the same
  size as the samples, so light bypassing the sample is subtracted. A Si wafer
  works below 900 nm and is transparent above 1000 nm.
- Any change to SBW, wavelength range, interval, scan rate or geometry
  invalidates the baseline.
- Do not use Abs mode. Cary's Abs is log₁₀(1/T), a different quantity.
- The session ended when the lamp failed, so the remaining samples (Ag 4, 6 and
  the MoOx series repeats) were never re-measured with the final settings.
