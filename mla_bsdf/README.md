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
python run_bsdf.py          # ~10 s for 90 angles × 3×10^5 rays (numba parallel)
```

Produces `mla_bsdf_result.npz`, `bsdf_maps.png`, `total_RT.png`.

## Files

| file | purpose |
|---|---|
| `raytrace.py` | numba Monte-Carlo tracer + `sweep()` driver |
| `run_bsdf.py` | full incident-angle sweep, saves data + plots |

## Comparing with your LightTools export

`mla_bsdf_result.npz` holds `bsdf_T`, `bsdf_R` (shape `[n_theta_i, 90]`, percent),
`theta_i`, `theta_out_centers`, and the integrated `T_total`/`R_total`. Load your
LightTools BSDF onto the same 1° grid to overlay / difference the two.
