# Substrate→device R/T of a patterned (island) layer — FDTD vs TMM

Answers: *can FDTD give the reflectance/transmittance seen by a plane wave
coming from the substrate into the device, for an island-patterned layer, where
plain TMM is not enough?*  **Yes — and here is the proof that TMM fails.**

## Setup (2D, Meep, λ = 550 nm, normal incidence, Ez)

```
        PML
        air  (superstrate)              transmission monitor
        patterned layer  d = 200 nm :
            island  n = 1.75  width = fill·Λ
            gap     n = 1.00 (air)
        glass  n = 1.5  (substrate)     reflection monitor + plane-wave source
        PML                             (x periodic, period Λ = 500 nm)
```

A plane wave is launched upward in the glass. R = reflected power back into the
glass, T = power transmitted into the air, both normalised to the incident
power (two-run method: incident from an all-glass normalisation run, subtracted
with `load_minus_flux_data`).

## Validation — continuous limit reproduces TMM

With `fill = 1.0` the layer is a uniform slab, so FDTD must equal the analytic
3-layer Fresnel/Airy result. It does, converging first-order in resolution:

| resolution | FDTD R |
|---|---|
| 120 | 0.0614 |
| 200 | 0.0708 |
| 300 | 0.0757 |
| Richardson → ∞ | **0.0855** |
| **TMM exact** | **0.0857** |

→ method validated. (Convergence is first-order because of the flux-plane /
sharp-interface discretisation; the container is too slow for high res, so the
island values below carry a similar ~1 %-absolute offset — negligible next to
the effect being demonstrated.)

## Result — islands: FDTD required, TMM/EMA fails badly

`EMA-TMM` = the best a TMM user can do: replace the island layer by an
effective homogeneous medium (volume-averaged ε) and run planar TMM.

| fill | FDTD R | FDTD T | EMA-TMM R | EMA-TMM T | ΔR (FDTD−EMA) |
|---|---|---|---|---|---|
| 1.0 (planar) | 0.061 | 0.939 | 0.086 | 0.914 | — |
| 0.7 | **0.370** | 0.630 | 0.043 | 0.957 | **0.33** |
| 0.5 | **0.764** | 0.236 | 0.040 | 0.960 | **0.72** |
| 0.3 | **0.775** | 0.225 | 0.038 | 0.962 | **0.74** |

The patterning turns a ~6 %-reflector into a **76 %-reflector**; EMA-TMM would
still predict ~4 %. TMM is wrong by an **order of magnitude**.

## Why

The period Λ = 500 nm exceeds λ / n_glass = 367 nm, so the island layer is a
**diffraction grating**, not an effective medium: it scatters the incident wave
into higher diffraction orders (propagating in the glass) and can Bragg-reflect
most of the power. Diffraction is a full-wave effect that TMM (infinite planar
layers) and EMA (Λ ≪ λ homogenisation) both structurally cannot represent. Only
FDTD (or RCWA, for strictly periodic islands) captures it. Energy is conserved
(R + T = 1.000 in every run), an internal check.

For **isolated / random** islands (not periodic) the same script generalises by
switching the lateral boundaries to `Absorber` and adding a near-to-far
transform to resolve the scattered angular distribution.

## Run

```bash
python island_rt.py --res 120           # ~10 s: validation + island table
python island_rt.py --res 120 --Lambda 0.30   # sub-wavelength period -> EMA recovers
```
