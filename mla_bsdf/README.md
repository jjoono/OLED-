# Microlens-array BSDF by geometric Monte-Carlo ray tracing

Open-source, license-free reproduction of a **LightTools**-style BSDF calculation
for a hemispherical microlens array (MLA), for validation against LightTools
results.

## Device modelled

| quantity | value |
|---|---|
| lens shape | full hemisphere |
| lens radius `R` | 20 µm |
| arrangement | hexagonal close-packed, pitch = 2R = 40 µm (lenses touch) |
| lens index | 1.5 |
| substrate | glass, **same index 1.5**, semi-infinite below |
| ambient | air (1.0), semi-infinite above |
| wavelength | single (index is non-dispersive, so λ does not enter geometric tracing) |
| illumination | plane wave inside the glass, from below, incident angle θᵢ = 0…90° |

Because the lens and the substrate share the same refractive index, the **only
optical interfaces** are:

1. the **spherical caps** (glass ⇄ air), and
2. the **flat inter-lens gaps** (glass ⇄ air) — the ~9.3 % of the area not
   covered by the close-packed circles.

The flat gaps produce the specular `θr = θi` diagonal seen in the reflection
BSDF; the caps produce the spread transmission/reflection lobes.

### Convex vs concave (`concave` flag)

The same array can be built two ways (sphere centres stay on `z = 0`):

* **convex** (`concave=False`, default): glass hemispheres bulge **up** into the
  air; the real spherical surface is the **upper** hemisphere (`z ≥ 0`), inside
  the sphere is glass.
* **concave** (`concave=True`): hemispherical **air cavities** are carved **down**
  into the glass; the real surface is the **lower** hemisphere (`z ≤ 0`), inside
  the sphere is the air pit.

Illumination (plane wave from below, in the glass) and the flat gaps are
identical for both. Result: below the 41.8° critical angle the **concave** pits
transmit more (≈65 %→53 %) than the convex lenses (≈45 %→48 %); above ~50° the
convex lenses win, and the two transmittance curves cross near the critical
angle (see `compare_total_RT.png`).

## Simplifications (equivalent to the full 25 mm array + 1 mm source)

The physical setup (25×25 mm array, 1 mm-radius surface source just below it) is
reduced — **without loss of generality** — to:

* **One hexagonal unit cell** with periodic surroundings. The array is periodic
  and enormous compared with the 40 µm pitch, and the finite source merely
  illuminates many cells uniformly, so edge effects are negligible.
* **One plane wave per incident angle.** Ray start points are sampled uniformly
  over the unit cell, reproducing uniform illumination and the exact
  90.7 % / 9.3 % lens / gap area split. Sweeping θᵢ gives the BSDF directly.
* **Azimuthal averaging**: each ray gets a random incidence azimuth φ, and only
  the polar output angle is binned (output azimuth integrated) → the θᵢ-vs-θ_out
  maps that match the LightTools plots.

## Method

Monte-Carlo non-sequential ray tracing (same principle as LightTools):

* At every interface, one uniform deviate chooses **reflection** (probability =
  unpolarised Fresnel R, including **total internal reflection** where R = 1) or
  **refraction**. Ray count is constant; multiple bounces and re-entry into
  neighbouring lenses are handled naturally.
* A ray leaving **upward into the air** half-space is **transmitted** (θt from
  +z); a ray returning **downward into the glass** half-space is **reflected**
  (θr from −z).
* No absorption anywhere ⇒ **T + R = 1** for every incident angle (verified).

Output normalisation matches LightTools "Power ratio (%)": each map cell is the
fraction of total incident power landing in that 1°-wide output-angle bin
(azimuth integrated), in percent.

## Validation (in `raytrace.py` limits)

* **Index-matched** (`n_glass = n_air`): T = 1 at all angles → interface
  bookkeeping is correct.
* **Flat limit** (`R → 0`): integrated T(θᵢ) reproduces the analytic unpolarised
  Fresnel transmittance and the 41.8° TIR cut-off to 3–4 digits.
* **Energy**: T + R = 1 to within the (≈3×10⁻⁵) fraction of rays exceeding
  `MAX_BOUNCE`.

## Run

```bash
pip install numpy matplotlib numba
python run_bsdf.py            # convex  -> convex_*   (~10 s, 90 angles x 3e5 rays)
python run_bsdf.py --concave  # concave -> concave_*
python compare.py            # export all CSVs + convex-vs-concave overlay
```

Each `run_bsdf.py` produces `PREFIX_result.npz`, `PREFIX_bsdf_maps.png`,
`PREFIX_total_RT.png` (PREFIX = `convex_` / `concave_`).

## Files

| file | purpose |
|---|---|
| `raytrace.py` | numba Monte-Carlo tracer (`concave` flag) + `sweep()` driver |
| `run_bsdf.py` | full incident-angle sweep, saves data + plots |
| `compare.py`  | export CSVs for both cases + convex-vs-concave R/T overlay |

## Comparing with your LightTools export

`PREFIX_result.npz` holds `bsdf_T`, `bsdf_R` (shape `[n_theta_i, 90]`, percent),
`theta_i`, `theta_out_centers`, and the integrated `T_total`/`R_total`; the same
data is written as `PREFIX_bsdf_T.csv`, `PREFIX_bsdf_R.csv`, `PREFIX_totals.csv`.
Load your LightTools BSDF onto the same 1° grid to overlay / difference the two.
