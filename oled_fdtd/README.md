# Minimal OLED FDTD with Meep (open source)

Single-wavelength (550 nm), 2D FDTD power budget of dipole emission in the
simplest credible bottom-emitting OLED stack, using
[Meep](https://meep.readthedocs.io/) (MIT, free, conda-forge `pymeep`).

## Stack (units: µm)

```
   PML
   glass     n=1.50   semi-infinite (1.2 µm slab ending in PML)
   ITO       n=1.90   100 nm   (lossless — simplest)
   organic   n=1.75   100 nm   <- point dipole at layer centre
   Al        Rakić Drude–Lorentz (meep.materials.Al)
   PML
```

Lateral size 6 µm, PML 0.5 µm all around, resolution 100 px/µm (10 nm).

## What is computed

For each dipole orientation (Ez = out-of-plane/TE, Ex = in-plane horizontal,
Ey = vertical/TM) and their isotropic ⅓ average:

| output | meaning |
|---|---|
| `frac_glass` | power crossing a plane just above the ITO into the glass — substrate-coupled light (upper bound for what can be extracted to air) |
| `frac_wg`    | power flowing laterally through vertical planes spanning the ITO+organic — waveguided modes |
| `frac_metal` | remainder — absorption in the Al cathode, incl. SPP (dominant for the vertical dipole) |

Total emitted power = Poynting flux through a small box enclosing the dipole.
Narrowband Gaussian pulse at f₀ = 1/0.55 µm⁻¹; run until fields decay (1e-6).

## Run

```bash
conda create -n mp -c conda-forge pymeep
conda run -n mp python oled_meep.py            # ~a minute on a laptop, 2D
conda run -n mp python oled_meep.py --res 150  # convergence check
```

Writes `oled_meep_result.json`.

## Physics expectations (sanity checks)

* Vertical dipole (Ey, TM): strong SPP coupling → large `frac_metal`.
* Out-of-plane dipole (Ez, TE): no SPP (TM-only) → most power split between
  glass and TE waveguide mode.
* Typical planar OLED literature: ~20-30 % substrate-coupled, ~30-50 %
  waveguide + SPP, rest absorbed — our 2D single-λ numbers should be in this
  ballpark but not identical (2D vs 3D dipole radiation patterns differ).

## Notes / simplifications

* 2D (not 3D): dipole radiation patterns and mode densities differ from 3D;
  fractions are qualitative. 3D is a straight extension (same script, add z).
* ITO lossless, organic lossless, single wavelength, single dipole position
  (EML centre). Each is easy to refine later (Drude ITO, dipole-position
  average, wavelength sweep).
* `frac_metal` is computed as a remainder, so it also absorbs the (small)
  numerical flux-bookkeeping error.
