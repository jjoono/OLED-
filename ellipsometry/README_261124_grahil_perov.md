# Session 2026-11-24 — GraHIL and the perovskite-NC film of the green PeLED

`scripts/grahil_perov_261124/` fits `GraHIL_perov.xlsx`: two Si pieces measured
at 65 / 70 / 75°, 192.4–1688.3 nm, with Ψ/Δ, N/C/S, intensity and
depolarization.

- **layer1** — GraHIL (PEDOT:PSS : PFI) on Si
- **layer2** — the same GraHIL with the perovskite nanocrystal film on top

The device they belong to is the green PeLED of Kim *et al.*, Nature 611, 688
(2022): glass / FTO 350 / BufHIL 75 / perovskite 270 / ZADN 45 / LiF 1.2 /
Al 100 nm. Its Methods, SEM cross-section (SI Fig. 5), absorbance-derived k
(SI Fig. 6a), EL spectrum (Extended Data Fig. 7c) and angular emission
(Extended Data Fig. 7g) are used below; the two curves and the angular points
were digitised.

```bash
export GRAHIL_XLSX=/path/to/GraHIL_perov.xlsx
export ELLIPS_DATA=/path/to/data      # the lab .mod with SI_JAW / NTVE_JAW, as in the 261117 session
export ELLIPS_OUT=/path/to/output     # also holds ed7g_digitized.json, ed7c_digitized.json (DATA.md)
cd scripts/grahil_perov_261124
python gp_load.py            # inventory; layer1's repeated scans matched to angles and averaged
python gp_grahil.py          # layer1: GraHIL thickness and n,k
python gp_explore.py         # layer2: which structure - nine alternatives on one footing
python gp_perov.py           # layer2: full-range fit and its variants   (python gp_perov.py "soft boundary" reruns one)
python gp_device_check.py    # angular emission of the device with these n,k
python gp_report.py          # -> $ELLIPS_OUT/grahil_perov/  n,k tables, structure.json, figure
```

Si and its native oxide (fixed, 1.5 nm) are the CompleteEASE SI_JAW / NTVE_JAW
tables; the transfer matrix is the corrected `lib/tmm_fix`. Fits stop at
1000 nm: beyond it the wafer is transparent, and layer2's exported
depolarization goes to −50 % and lower, which no sample can do.

## Result in one table

| | this work | assumed in the paper's optical simulation |
|---|---|---|
| GraHIL thickness (on Si) | 108.3 nm | 75 nm |
| GraHIL n (550 nm) | 1.381, k ≈ 0.006 | 1.5 |
| perovskite, total | 284 nm | 270 nm |
| perovskite, depth profile | 122 nm dense on 162 nm of effective n ≈ 1.5–1.65 that scatters | one homogeneous layer |
| dense perovskite n (540 nm) | 2.34 (2.34–2.37 over the fit variants) | 2.3, constant |
| dense perovskite k (540 nm) | 0.076 (0.065–0.089) | 0.035 (absorbance) |
| device angular emission, rms vs measured | 0.019 | 0.12–0.20 |

## 1. GraHIL (layer1)

One homogeneous isotropic layer with a Kramers-Kronig-consistent dispersion
(Gaussians on a fixed energy grid + UV pole + ε∞) fits to **MSE 1.37**:

| variant | MSE | thickness (nm) | what the extra freedom does |
|---|---|---|---|
| isotropic | 1.37 | 108.34 ± 0.06 | — |
| graded through the depth | 1.35 | 108.35 | grading −0.6 %: none |
| uniaxial | 1.26 | 108.65 | out-of-plane (ε−1) × 0.992 |
| graded + uniaxial | 1.23 | 108.81 | 1.6 % and 0.986 |

n(550 nm) = 1.381, k rising from 0.002 (300 nm) through 0.006 (550 nm) to
0.007–0.008 (633–800 nm), the polaron tail of PEDOT. The self-organised PFI gradient is not visible to the
optics: PFI and PSS are too close in index. The thickness moves by ±1.2 nm when
the native oxide is varied over 0.5–3 nm. It is 108 nm on this piece, not the
75 nm the paper quotes for the device.

## 2. Perovskite film (layer2): the structure

The nominal stack — one perovskite layer, with or without surface roughness, on
the 108 nm GraHIL — does not fit at any thickness or index. `gp_explore.py`
compares the alternatives in the transparent range (570–1000 nm, Sellmeier
perovskite, 26–80 random starts each):

| model | MSE | rms R ×10³ | best parameters |
|---|---|---|---|
| one layer on GraHIL 108 nm | 428 | 109 | — |
| one layer, GraHIL thickness free | 94 | 42 | GraHIL → 1 nm, n(633) 2.75 |
| one layer + height spread (coherent average) | 70 | 61 | σ 44 nm, n(633) 2.64 |
| one layer + scattering loss | 45 | 37 | n(633) 2.45 |
| dense / low-index lower part, lossless | 148 | 156 | — |
| dense / lower part, absorption band in the GraHIL | 11.3 | 36 | GraHIL k(600) ≈ 0.3 |
| **dense / lower part that scatters** | **10.0** | **25** | 122 / 164 nm, 51 % perovskite, k_sc(633) 0.26, p 3.8 |
| the same, host index of the lower part free | 10.1 | 24 | 22 % perovskite in n = 1.35 |
| the same, dense part allowed to scatter | 10.0 | 25 | k_top → 0 |

So two things are needed together, and neither alone is close:

1. **The perovskite is not uniform through its depth.** The upper ~120 nm is
   dense (n(633) ≈ 2.1); the lower ~160 nm has an effective index of only
   ≈ 1.5. Half perovskite and half void, or a quarter perovskite in an n ≈ 1.35
   host, describe it equally well — the optics see the effective index, not
   what is between the grains.
2. **The lower part removes light from the specular beam**, with
   k_sc ∝ λ^−3.6…−3.8. That exponent is Rayleigh-like: scattering by grains and
   gaps smaller than the wavelength, not absorption. The dense part does not
   scatter (k_top → 0).

Putting the loss into the GraHIL instead fits the ellipsometric quantities
nearly as well (11.3) but the reflectance worse (36 against 25), and it needs
k = 0.3–0.7 at 500–600 nm in a film whose k is 0.006 on layer1. A dedoped
PEDOT minority phase cannot supply that; the alternative is rejected on
physical grounds, not on MSE alone.

This matches the cross-section SEM of the device (SI Fig. 5): a featureless
upper band on blocky ~100–200 nm grains with dark gaps, and the paper's grain
statistics (in situ particles 123 ± 34 nm, 10 ± 2 nm nanocrystals inside
~100 nm macro-grains).

R is the unpolarised reflectance obtained from the two sheets' intensities,
R₂ = (I₂/I₁)·R₁(model of layer1); the lamp cancels because both pieces were
measured in the same run (layer1's own angle ratios agree with its model to
0.26 % rms). The model reproduces its shape, not its amplitude (rms 0.024, up
to 0.06 locally); dropping R from the fit leaves the structure and the n,k
unchanged.

## 3. Perovskite film: full-range fit

`gp_perov.py`, 205–1000 nm, three angles:

```
air
EMA roughness          0.6 nm     dense / void 50:50
dense perovskite     121.7 nm     KK-consistent basis (fine grid across the band edge)
lower layer          161.9 nm     Bruggeman(dense, void, 0.51) + i·0.258·(633/λ)^3.60
GraHIL               108.3 nm     n,k of layer1, fixed
native oxide 1.5 nm / Si
```

**MSE 9.38** (3.6 / 5.0 / 7.6 / 11.4 for 205–300 / 300–450 / 450–570 /
570–1000 nm). Each variant below changes one assumption and refits everything:

| variant | MSE | dense / lower (nm) | perovskite fraction below (top → bottom) | note |
|---|---|---|---|---|
| base | 9.38 | 122 / 162 | 0.51 | |
| GraHIL thickness free | 8.29 | 120 / 183 | 0.52 | GraHIL 82 nm |
| surface height spread | 8.63 | 121 / 166 | 0.50 | σ 8.8 nm |
| host of the lower part free | 8.38 | 119 / 168 | 0.23 | host n 1.34 |
| dense part may scatter | 8.45 | 122 / 163 | 0.51 | k_top 0.006 |
| no reflectance in the fit | 9.44 | 122 / 162 | 0.51 | |
| 65° + 70° only | 9.13 | 122 / 159 | 0.52 | |
| 70° + 75° only | 8.40 | 123 / 157 | 0.52 | |
| from 230 nm | 9.38 | 122 / 161 | 0.51 | |
| lower part graded | 8.88 | 122 / 163 | 0.53 → 0.49 | no gradient |
| soft dense/lower boundary | 8.44 | 123 / 163 | 0.50 | 21 nm transition |
| graded + GraHIL free | 8.32 | 120 / 185 | 0.52 → 0.53 | GraHIL 79 nm |
| graded + height spread | 7.28 | 115 / 190 | 0.15 → 0.66 | σ 22 nm |
| all of the above | 5.50 | 133 / 142 | 0.15 → 0.88 | 104 nm transition, GraHIL 118 nm |

What is stable across them, and what is not:

- **Stable:** a dense upper part of 115–133 nm; a lower part whose mean index is
  ≈ 1.5; the scattering law (k_sc(633) 0.22–0.27, exponent 3.5–4.0; only
  "graded + height spread" trades half of it for the spread); the dense n,k
  from 480 nm to 1000 nm (table below).
- **Not determined:** where the lower part ends and the GraHIL begins. Their
  indices are 1.5 and 1.38, so 25 nm can move from one to the other for 1 MSE
  unit (GraHIL 79–118 nm when freed). The fit therefore does not show the
  GraHIL was changed by the perovskite deposition — and cannot exclude a
  change of that size either.
- **Surface roughness** comes out at 0.5–6 nm of effective-medium layer, which
  is small for a nanocrystal film; the UV data bound it below 15 nm. A height
  spread of σ ≈ 9–22 nm with a lateral scale too large to act as an effective
  medium is allowed and slightly preferred.
- The most flexible variant (5.50) turns the two-zone step into a continuous
  profile: 81 nm dense, the perovskite fraction then falling to 0.15 over
  104 nm and rising again to 0.88 at the GraHIL over the last 90 nm, i.e. the
  emptiest zone sits just under the dense cap. "Graded + height spread" does
  the same (0.15 → 0.66); the lower part graded on its own shows no gradient.
  The data do not choose between the step and such a profile. All of them
  have a dense top with low-index material beneath it.

### Dense perovskite n,k

| λ (nm) | n | range over variants | k | range over variants |
|---|---|---|---|---|
| 450 | 2.333 | 2.281–2.346 | 0.324 | 0.226–0.341 |
| 500 | 2.283 | 2.250–2.289 | 0.208 | 0.191–0.232 |
| 520 | 2.303 | 2.295–2.322 | 0.179 | 0.165–0.202 |
| 530 | 2.320 | 2.313–2.348 | 0.156 | 0.145–0.178 |
| 540 | 2.342 | 2.335–2.368 | 0.076 | 0.065–0.089 |
| 550 | 2.284 | 2.272–2.313 | 0.016 | 0.005–0.037 |
| 600 | 2.144 | 2.140–2.175 | 0 | |
| 633 | 2.106 | 2.103–2.141 | 0 | |
| 800 | 2.022 | 2.021–2.066 | 0 | |
| 1000 | 1.985 | 1.985–2.034 | 0 | |

The band edge and the excitonic bump in n at 540 nm are resolved without any
anchor: the absorbance-derived k of SI Fig. 6a was used only as a starting
point for 510–565 nm. Averaged over the film (dense part + 51 % of the lower
part = 0.72 of the volume) the fitted k is 0.129 at 520 nm, against 0.125 from
the paper's absorbance. Below 450 nm the dense k depends on how the few
nanometres of surface roughness are described (0.36–0.54 at 400 nm).

The lower layer as an effective medium: 1.645 + 0.49i (540 nm), 1.535 + 0.26i
(633 nm), 1.495 + 0.11i (800 nm), the imaginary part being scattering. That
law is determined where light reaches the layer (≳ 480 nm); the delivered
table holds it constant below 450 nm instead of extrapolating λ^−3.6 into
the UV.

## 4. The device's angular emission

`gp_device.py` computes the far field by reciprocity (|E|² in the emitting
layer of a plane wave incident from the glass; isotropic dipoles spread over
the perovskite volume; incoherent glass with its multiple reflections).
`gp_device_check.py` scales the ellipsometric structure to the device's 270 nm,
puts the lower part on the HIL side, uses the GraHIL n,k at 75 nm, and solves
the internal emission spectrum so that the normal-direction output equals the
measured EL spectrum. Intensity relative to normal, rms over 0–70°:

| | 20° | 40° | 60° | 70° | rms |
|---|---|---|---|---|---|
| measured (ED Fig. 7g) | 0.918 | 0.725 | 0.435 | 0.273 | |
| Lambertian | 0.940 | 0.766 | 0.500 | 0.342 | 0.043 |
| fully scattered light (isotropic in the glass) | 0.939 | 0.761 | 0.474 | 0.295 | 0.028 |
| paper optics: n = 2.3, k from absorbance, one layer, HIL n = 1.5 | 1.025 | 0.927 | 0.562 | 0.355 | 0.130 |
| ellipsometric dense n,k as one 270 nm layer | 1.011 | 0.907 | 0.565 | 0.354 | 0.121 |
| **ellipsometric structure (base fit)** | 0.912 | 0.726 | 0.470 | 0.297 | **0.019** |

(perovskite thickness averaged over a 15 nm lateral spread; without it the
one-layer descriptions are worse, 0.17–0.20, and the ellipsometric structure
is 0.018.)

A single high-index layer, with either set of n,k, predicts a lobe at 20–40°
that the device does not have. The two-part structure removes it. Sensitivity
of that agreement, one change at a time from the base fit:

| change | rms |
|---|---|
| lower part made as dense as the top (no low-index region) | 0.066 |
| scattering loss off / ×0.5 / ×1.5 | 0.076 / 0.028 / 0.025 |
| dense n ×0.95 / ×1.05 | 0.035 / 0.014 |
| perovskite fraction of the lower part 0.35 / 0.70 | 0.009 / 0.020 |
| dense part 90 nm / 150 nm of 284 | 0.066 / 0.020 |
| emission only in the dense part / only in the lower part | 0.079 / 0.072 |
| BufHIL 108 nm instead of 75 nm | 0.035 |
| ZADN n 1.85 instead of 2.0 | 0.023 |
| fit variants of §3 | 0.011–0.036, 0.098 for "graded + height spread" |

How far this goes: the angular data exclude the one-layer picture and are
consistent with the ellipsometric structure, but they do not refine it. A
featureless emitter (fully scattered light) is already at 0.028, and the
digitised points carry about 0.01, so differences between 0.01 and 0.03 are
not significant. The device model also treats scattered light as removed;
if it leaves through the glass it adds a component of the "fully scattered"
shape, which is close to the model's own curve and would not change the
comparison.

## 5. Open points

- **The lower layer is predicted to be hazy.** Its scattering optical depth is
  0.8 at 633 nm and 1.7 at 540 nm at normal incidence, i.e. a film on glass
  should transmit only ~45 % of 633 nm light specularly, the rest diffusely.
  The paper's absorbance (SI Fig. 6a) shows no sub-gap tail, so either that
  spectrum was baseline-corrected, or this Si piece differs from the paper's
  films, or the loss is not ordinary scattering. A total-versus-specular
  transmittance of a film on glass at 600–800 nm settles it.
- The reflectance amplitude is not reproduced (§2).
- The angular data are the paper's figure, digitised, for a device whose HIL
  and perovskite thicknesses are taken from its Methods. Angle-resolved
  spectra of the actual device (spectrum against angle, not only intensity)
  would constrain n,k far more than the integrated intensity does.
- `perovfit*.py`, `k_anchored.py`, `final_window.py` and `Perov_FINAL_nk.csv`
  (2026-08-16) predate the transfer-matrix fix and assumed one homogeneous
  270 nm layer with k tied to SI Fig. 6a. This session supersedes them.

## Files written to `$ELLIPS_OUT/grahil_perov/`

| file | content |
|---|---|
| `Perov_dense_nk.csv` | dense perovskite n,k, 205–1000 nm, 1 nm; min/max over the fit variants; lower-layer effective n,k |
| `GraHIL_nk.csv` | GraHIL n,k |
| `lib_Perov_dense.csv`, `lib_Perov_lower_with_scattering.csv`, `lib_Perov_lower_no_scattering.csv`, `lib_GraHIL.csv` | the same media as `wavelength_nm,n,k` |
| `structure.json` | thicknesses, every variant, the tables above |
| `grahil_perov_summary.png` | fit at three angles, n, k, device angular emission |

For a device simulation use two perovskite layers, not one: `lib_Perov_dense`
on the ETL side (0.43 of the perovskite thickness) and a `lib_Perov_lower`
table on the HIL side (0.57). "with scattering" reproduces the specular
(coherent) field; "no scattering" is the same effective index without the
loss, for codes that handle scattered light separately.
