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
python gp_genosc.py          # both dispersions refitted as compact Gen-Osc materials
python gp_make_mod.py        # -> $ELLIPS_OUT/grahil_perov_mod/  CompleteEASE .mod files, read back and re-scored
python gp_rough.py           # the layer2 fit with the surface roughness imposed (§9)
python gp_device_model.py    # dipole model of the device: where every photon goes, per description (§8)
python gp_device_fit.py      # best angular pattern of each description within the micrograph's thicknesses
python gp_stitch.py          # luminescence-stitched band edge (§11); DASGUPTA_CSV = published set, optional
python gp_report2.py         # -> grahil_perov_device.png
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
| device angular emission, rms vs measured, nominal thicknesses (§4; not decisive, §8) | 0.019 | 0.12–0.20 |
| EQE ceiling of the dipole model (η_rad = 1, substrate light lost, §8) | 24–26 % | 20–24 % |

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

The device cross-section (SI Fig. 5) has the same orientation but not the same
proportions: a low-relief band under the ETL, 45 nm thick at the median
(24–97 nm), over blocky 80–150 nm grains that reach the HIL, with no visible
boundary at 43 % of the depth. The paper's grain statistics (in situ particles
123 ± 34 nm, 10 ± 2 nm nanocrystals inside ~100 nm macro-grains) give the
scatterers their size. §10 goes through what the micrographs do and do not
support.

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
  is small for a nanocrystal film. A height spread of σ ≈ 9–22 nm with a
  lateral scale too large to act as an effective medium is allowed and
  slightly preferred. §9 imposes both kinds and compares with the micrographs:
  the small value is required by the data.
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

At these nominal thicknesses a single high-index layer, with either set of
n,k, predicts a lobe at 20–40° that the device does not have, and the two-part
structure removes it. (§8 repeats the comparison with the thicknesses the
device micrograph allows: there the single layer fits as well. The numbers of
this section are for the nominal stack only.) Sensitivity
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

## 6. CompleteEASE `.mod` files

The fits above describe each material with a fixed basis of 32 / 44 Gaussians.
That is the right tool for finding the n,k without assuming a line shape, but
not something to edit in a Gen-Osc dialog. `gp_genosc.py` therefore refits
both, against the data, as ordinary Gen-Osc materials, and `gp_make_mod.py`
writes them into the lab's own model file frame:

| file | stack (substrate up) | MSE, basis fit | MSE, Gen-Osc fit | MSE from the file |
|---|---|---|---|---|
| `GraHIL_108.4nm_on_NTVE15A_Si_GenOsc.mod` | SI_JAW / NTVE_JAW 15 Å / **GraHIL** 108.36 nm | 1.37 | 1.74 | 1.74 |
| `PerovNC_on_GraHIL_on_NTVE15A_Si_GenOsc.mod` | SI_JAW / NTVE_JAW 15 Å / **GraHIL** 108.3 nm (fit off) / **Perov_lower** 163.8 nm / **Perov_dense** 121.4 nm / Roughness 1.36 nm | 9.38 | 9.89 | 9.89 |

"MSE from the file": each `.mod` is parsed back, its tables decoded and its
oscillators re-evaluated, the stack rebuilt from that alone and scored against
the data (200–1000 nm for layer1, 205–1000 nm for layer2). The files were not
opened in CompleteEASE here.

**GraHIL** — Einf 1.244, UV pole 89.28 at 12 eV (fixed energy), IR pole 0.0223,
four Gaussians:

| # | Amp | Br (eV) | En (eV) |
|---|---|---|---|
| 1 | 0.0440* | 3.123 | 0.600* |
| 2 | 0.0145 | 1.612 | 4.883 |
| 3 | 0.0622 | 0.448 | 5.479 |
| 4 | 0.2188 | 0.889 | 6.500 |

With Einf fixed at 1, as in the tooling files, this film cannot be fitted (MSE
12–14): its dispersion is too flat for a pole alone. The IR pole carries the
free-carrier dispersion of PEDOT; a Gaussian placed there instead slides below
the data range. n,k agree with the basis result to 0.0014 / 0.0023 over
250–1000 nm.

**Perov_dense** — Einf 2.067, UV pole 1.48* at 6.34 eV, seven Gaussians:

| # | Amp | Br (eV) | En (eV) | |
|---|---|---|---|---|
| 1 | 0.344 | 0.080 | 2.333 | exciton / band edge, 531 nm |
| 2 | 0.549* | 0.164 | 2.416 | |
| 3 | 0.515* | 0.255* | 2.560 | |
| 4 | 0.982* | 0.417 | 2.825 | |
| 5 | 1.626 | 0.658 | 3.228 | |
| 6 | 3.605 | 1.236 | 4.100 | |
| 7 | 2.137 | 1.939 | 6.007 | |

(* written with the fit flag off: 90 % interval wider than half the value, or
on a bound.) The structure refitted with it is 121.4 / 163.8 nm, 50.7 %
perovskite below, k_sc(633) 0.252, exponent 3.52 — the base fit within its
variant spread — and n,k follow the basis result to 0.012 / 0.010 over
450–1000 nm (540 nm: 2.336 + 0.080i). Left to themselves the broad UV
Gaussians put an absorption tail under the gap (k = 0.007 at 633 nm, with a
slightly lower MSE, 9.41); ε₂ above 580 nm is penalised in this fit, which
brings it to k < 0.0012 there.

**Perov_lower** is a tabulated layer, not an EMA layer: CompleteEASE has no
effective-medium layer with a scattering loss. The table is
Bruggeman(Perov_dense, void) plus k_sc·(633/λ)^p, held constant below 450 nm,
computed from the same fit on a 0.005 eV grid (188–1771 nm; an extrapolation
outside 205–1000 nm). It does not follow if Perov_dense is refitted in
CompleteEASE; thicknesses and roughness can be.

To reproduce the MSE in CompleteEASE, limit the fit range to 205–1000 nm
(200–1000 nm for layer1): beyond 1000 nm the wafer is transparent and
layer2's exported depolarization is unusable.

Format notes, beyond those of `README_260820_session.md` §2 — all read off
files CompleteEASE saved:

* a material with its own name uses that name in the section tags
  (`start_<name> Fit Parms`, `… Misc Parms`, `… Grade Parms`,
  `… Permanent Poles`), and on the layer line twice
  (`'Layer' '<name>' 'GENOSC' '' '<name>'`);
* a tabulated film is the `NTVE_JAW` block with another name: `Mat Table`
  count, then `Wvl` (Å, ascending), `e1`, `e2`, each base64 of gzip of
  big-endian float32. The gzip header is ten bytes `1f 8b 08` + zeros and the
  stream is deflate level 6: re-encoding the template's own tables this way
  returns its bytes exactly, which `gp_make_mod.py` asserts before writing;
* `IR Pole Amp.` is A in −A/E²;
* numbers are written as CompleteEASE prints them (`1.0E-8`, not `1e-08`).

## 7. Using these n,k in a device model

What the description of the perovskite layer does to the amount of light a
model predicts (`gp_device_check.py`; far-field power into air per unit dipole
strength — no Purcell normalisation — with one internal spectrum for all rows):

| perovskite described as | normal direction | total into air |
|---|---|---|
| one layer, paper optics (n 2.3, k from absorbance, HIL n 1.5) | 1.52 | 1.70 |
| one layer, ellipsometric dense n,k | 1.29 | 1.37 |
| **two layers, lower part without the scattering loss** | 1.00 | 1.00 |
| two layers, scattering loss × 0.5 | 0.49 | 0.46 |
| two layers, scattered light counted as lost | 0.26 | 0.24 |

A factor of seven between the first and the last row. The structure matters
more than the exact n,k, and how the scattered light is treated matters most.

1. **Two perovskite layers, not one.** `lib_Perov_dense` on the ETL side,
   0.43 of the thickness; a `lib_Perov_lower` table on the HIL side, 0.57.
   Use the tables over the whole EL band, not one index: n goes 2.28 → 2.34 →
   2.28 → 2.14 between 500 and 600 nm.
2. **Scattered light is not absorbed light.** At 540 nm 93 % of the lower
   part's extinction is scattering, and one pass through it leaves 17 % of
   the light in the coherent beam. `lib_Perov_lower_with_scattering` gives the
   coherent field — right for the shape of the specular interference
   pattern, wrong as a loss: it removes three quarters of the output, which a
   device with the reported 28.9 % EQE cannot have lost. Use
   `…_no_scattering` (only the perovskite's own absorption) for efficiency
   and energy balance, as the upper bound; the coherent-only number is the
   lower bound. The angular pattern is then the coherent part plus a diffuse
   part of the "fully scattered" shape (§4), its share being the one number
   the ellipsometry cannot supply.
3. **Re-absorption.** The EL peak sits on the absorption edge: in the dense
   part k = 0.179 / 0.076 / 0.003 at 520 / 540 / 560 nm, so one pass through
   its 116 nm absorbs 39 / 18 / 1 %. The internal emission spectrum is
   therefore bluer than the measured EL (peak 536–538 nm against 540 nm
   here) and must be the model's source, not the EL; and an EQE model needs
   photon recycling. The paper's k at 540 nm is half this value.
4. **Emitters** spread over both parts, weighted by the perovskite fraction
   (1 and 0.51); the angular data reject emission from either part alone
   (§4). The depth profile of recombination is the free parameter to fit.
5. **GraHIL** n = 1.38, k = 0.006, not 1.5; and its thickness in the device
   needs its own measurement — 108 nm on this Si piece, 75 nm in the paper's
   Methods, 50 ± 11 nm in the paper's device micrograph — and the angular
   pattern is sensitive to the difference (§4, §8).
6. Carry the spread, not one curve: the `n_min … k_max` columns, a lower-part
   perovskite fraction of 0.35–0.70 and a dense share of 0.40–0.47 cover the
   fit variants.

§8 carries this out with a dipole model in which every photon is followed.

Two measurements would replace most of these assumptions: total against
specular transmittance of the film on glass (is the scattering real, and how
strong), and the reflectance spectrum of the finished device from the glass
side, which tests the stack's optics with no assumption about the emitters.

## 8. Device model: every photon followed

`gp_dipole.py` is a dipole-emission solver for planar stacks (plane-wave
expansion of the dipole field, exact multiple reflection at the dipole plane,
power flow through every interface). `gp_device_model.py` uses it to follow
each emission event to its end: into air, trapped in the glass, re-absorbed
by the perovskite (and re-emitted: photon recycling), absorbed elsewhere, or
scattered in the lower zone — in which case it restarts as an isotropic,
incoherent source inside that zone at the same wavelength, and is followed
again. EQE = η_rad·air / (1 − η_rad·A_perovskite).

Checks of the solver: a dipole in a homogeneous medium (exact); in front of a
perfect mirror (the analytic decay-rate curves, to 10⁻⁴); the far field
against the independent reciprocity calculation of `gp_device.py` (four
digits, absorbing emitter layer included). And it reproduces the paper's own
simulation when given the paper's optics:

| one layer, n = 2.3, k from absorbance | this code | paper (SI Fig. 6b, Methods) |
|---|---|---|
| direct outcoupling | 8.5 % | 9 % |
| re-absorbed by the perovskite | 72.5 % | ~69 % |
| parasitic | 19.1 % | ~22 % |
| EQE at η_rad = 1, substrate light recycled | 30.8 % | 29.2 % |
| escape probability of the bare film → η_rad from the 46 % PLQE | 0.104 → 0.89 | ~0.11 → 0.88 |

### Descriptions compared

Nominal thicknesses (FTO 350 / BufHIL 75 / perovskite 270 / ZADN 45), 15 nm
thickness spread, emission by perovskite volume, substrate light counted as
lost (the EQE was measured with a luminance meter):

| perovskite described as | angular rms | air | glass | re-absorbed | parasitic | EQE at η_rad = 1 | η_rad needed for 25.2 % |
|---|---|---|---|---|---|---|---|
| one layer, the paper's optics | 0.112 | 8.0 | 8.2 | 66.2 | 17.6 | 23.7 % | 1.02 |
| one layer, ellipsometric dense n,k | 0.104 | 6.1 | 6.0 | 73.3 | 14.6 | 22.7 % | 1.03 |
| two layers, no scattering | 0.080 | 7.8 | 13.3 | 56.4 | 22.6 | 17.8 % | 1.15 |
| two layers, scattered light counted as lost | 0.018 | 1.8 | 1.8 | 29.2 | 8.6 (+58.5 scattered) | 2.6 % | 2.7 |
| **two layers, scattered light followed** | 0.064 | 9.7 | 10.9 | 62.1 | 17.3 | **25.6 %** | 0.99 |
| the same, substrate light recycled | 0.066 | 11.1 | — | 68.8 | 20.1 | 35.6 % | 0.89 |

(shares in % of emitted photons, averaged over the internal spectrum; angular
rms over 0–70° against Extended Data Fig. 7g, where a Lambertian source
scores 0.043.)

What this says:

- **Scattered light counted as lost** gives the best angular pattern and an
  EQE of 3 %. It is not a model of the device; its pattern is that of the
  coherent part alone.
- **The paper's flat optics**, at the nominal thicknesses, misses its own
  angular measurement by the widest margin: a lobe at 20–40° (0.90 at 40°
  against the measured 0.725) that a lateral thickness spread reduces but does
  not remove (rms 0.19 / 0.11 / 0.08 for 0 / 15 / 25 nm). Other thicknesses
  change that — see below.
- **Two layers without scattering** lowers the EQE ceiling to 18 %: the
  low-index zone sends light into the glass (13 %), and nothing brings it
  back.
- **Two layers with the scattered light followed** conserves energy, has no
  lobe, and reaches the mean measured EQE — at η_rad ≈ 1 and unit charge
  balance, as the paper's own model also requires. The champion 28.9 % is out
  of reach for every description here unless substrate light is recovered
  (η_rad 0.93 would then do).

### Thicknesses within what the micrograph allows

The Methods give nominal thicknesses. `gp_device_fit.py` searches the ranges
the device cross-section shows (FTO 342–396, BufHIL 34–75, perovskite 231–273,
ZADN 35–51 nm, with the coat under the FTO; dense share 0.15–0.53; x_d), 140
random geometries per family, the best re-evaluated at full resolution:

| description | best angular rms | at FTO / HIL / perovskite / ZADN (nm), dense share, x_d | EQE at η_rad = 1: substrate light lost / recycled | over all geometries: 10th percentile / median / worst |
|---|---|---|---|---|
| one layer (dense n,k) | 0.010 | 357 / 34 / 258 / 47 | 20.4 % / 23.9 % | 0.021 / 0.044 / 0.094 |
| two layers, no scattering | 0.016 | 383 / 44 / 269 / 39, 0.52, 1.0 | 20.4 % / 23.6 % | 0.038 / 0.128 / 0.278 |
| two layers, scattered light followed | 0.011 | 393 / 40 / 266 / 47, 0.52, 1.0 | 21.1 % / 28.8 % | 0.025 / 0.046 / 0.058 |
| **the same, emission by perovskite volume** | **0.015** | 366 / 54 / 272 / 42, 0.53, 0.69 | **23.8 % / 28.8 %** | — |

- **The angular distribution does not choose between the descriptions.** With
  the thicknesses free inside the micrograph's ±10 %, every one of them
  reaches the noise of the digitised points (~0.01) — the single layer too.
  The lobe of the one-layer optics, and the statement of §4 that the angular
  data exclude that picture, belong to the nominal thicknesses only (a thinner
  BufHIL, 34–45 nm, removes the lobe).
- What does differ is how much depends on getting the thicknesses right.
  With the scattering followed the pattern hardly moves over the whole range
  (worst 0.058); without it the same structure ranges from 0.016 to 0.28.
- And the EQE ceiling: 20 % for one layer, 24 % for two layers with the
  scattering followed and emission spread by volume (29 % if substrate light
  is recovered), against a measured mean of 25.2 % and a best of 28.9 %.

The description to use is therefore the last row: two perovskite layers, the
scattered light followed, emission by perovskite volume, photon recycling,
incoherent glass, and the device's own thicknesses (FTO ≈ 366 + coat, BufHIL
≈ 54, perovskite ≈ 272, ZADN ≈ 42 nm). It is chosen because the ellipsometry
says the film is not one layer and because it is the least sensitive to what
is not known — not because the angular data prefer it.

### The photoluminescence yield

The paper turns the film's external PLQE of 46 % into an internal efficiency
of 88 % with its flat-film optics. The same step with each description
(bare film on a transparent substrate, substrate-guided light re-absorbed):

| film described as | escape probability | implied η_rad | EQE the device model then predicts |
|---|---|---|---|
| one layer, the paper's optics | 0.104 | 0.89 | 17.4 % |
| two layers, no scattering | 0.135 | 0.86 | 13.1 % |
| two layers, scattered light followed | 0.287 | 0.75 | 13.5 % |

(0.20 / 0.32 / 0.38 escape and η_rad 0.81 / 0.73 / 0.69 if substrate-guided
light is counted as escaping through the edges inside the sphere.) No
description connects the 46 % PLQE to a 25–29 % EQE: the device at its
efficiency peak must run at a higher radiative efficiency than the film under
62.5 mW cm⁻² of 325 nm light, by any of these optics. A scattering film makes
that gap larger, not smaller, because it lets more of the photoluminescence
out. This is the one observation that speaks against scattering as strong
as on the Si piece; a haze measurement decides it (§5).

## 9. Surface roughness: is the small fitted value an artefact?

The fit gives 0.6 nm of effective-medium roughness, while the device
cross-section (SI Fig. 5) shows a visibly undulating perovskite surface. Two
different things are being compared. `gp_rough.py` repeats the fit of §3 with
the roughness imposed, everything else refitted:

| imposed | MSE | dense / lower (nm) | n(633) of the dense part |
|---|---|---|---|
| effective-medium layer fitted: 0.6 nm | 9.38 | 122 / 162 | 2.106 |
| effective-medium layer 3 nm | 10.4 | 121 / 154 | 2.125 |
| 6 nm | 14.4 | 120 / 144 | 2.150 |
| 10 nm | 23.6 | 120 / 130 | 2.189 |
| 15 nm | 38.3 | 120 / 116 | 2.232 |
| 20 nm | 55.0 | 119 / 106 | 2.286 |
| height spread σ = 5 nm (long lateral scale) | 9.02 | 121 / 163 | 2.106 |
| σ = 10 nm | 8.62 | 121 / 167 | 2.107 |
| σ = 13.5 nm | 9.35 | 119 / 172 | 2.111 |
| σ = 20 nm | 14.7 | 116 / 181 | 2.122 |
| 6 nm layer and σ = 13.5 nm | 13.2 | 118 / 162 | 2.139 |

- **An effective-medium layer** stands for roughness much finer than the
  wavelength. The data reject more than ~3–5 nm of it: the MSE grows
  quadratically, and the fit can only follow by raising the index of the dense
  part above that of bulk lead bromide perovskite (2.107 at 633 nm in the
  published data of §11). The small value is required by the data, not an
  over-fit.
- **A height spread with a long lateral scale** (amplitudes averaged
  coherently) is what a film undulating over hundreds of nanometres is. To
  first order it leaves r_p/r_s unchanged, so ellipsometry barely sees it;
  it shows in the fringe contrast and in the reflected intensity. 5–13.5 nm
  fits as well as or better than none, with the same n,k.

The micrographs say the same. Along 3 µm of the device cross-section the
perovskite top surface has rms 13.5 nm (peak-to-valley ~55 nm), but split by
lateral period it is 1.7 nm below 100 nm, 3.2 nm at 100–300 nm and 11.4 nm at
0.3–1 µm; part of it is inherited from the FTO (rms 14 nm), which a film on
Si does not have. The top-view SEM of the core/shell film (SI Fig. 1c) shows
close-packed, rounded cobbles of 40–130 nm with dark grooves and no pinholes.
A top view gives the lateral size, not the height: full hemispheres would
mean rms 8–11 nm, caps a quarter as high as wide 5–6 nm, an eighth 2–3 nm,
and groove contrast in an SEM does not tell which. The ellipsometry does, for
this Si piece: an equivalent layer of at most ~5 nm, i.e. shallow caps (rms
≲ 3 nm). The paper has no AFM; one scan of the film would settle it.

So: effective-medium roughness 1–5 nm, plus a long-scale height spread of
5–13 nm on Si (8–15 nm in the device). The "surface spread" fits of §3 are the
preferred description; the device model averages over a thickness spread of
15 nm.

One caution on how the roughness is determined. Below 450 nm an imposed
roughness layer costs nothing (the fitted n,k of the dense part absorb it:
n(300 nm) goes from 1.88 to 2.09 at 6 nm); the penalty comes from 450–1000 nm,
where the same, Kramers-Kronig-smooth dispersion must also give the fringes.
The roughness is therefore fixed by the consistency of n between the opaque
and the transparent range, not by the UV data alone.

## 10. Assumptions checked against the paper, its SI and its micrographs

| assumption | evidence | verdict | used |
|---|---|---|---|
| perovskite 270 nm | Methods: 270; cross-section: 254 ± 14 (231–273) | supported | 254–270 |
| BufHIL 75 nm | Methods: 75; cross-section: 50 ± 11 (34–67); our Si piece: 108 | the device is thinner than nominal, and not 108 | 35–75 scanned |
| FTO 350 nm, flat | cross-section: 367 ± 15, rms 14 nm, plus a ~47 nm coat under it | thickness supported; coat and texture not in any model so far | 367 + SnO₂/SiO₂ coat in the "micrograph" geometry |
| ZADN 45 nm, n = 2.0 | cross-section 42 ± 5; the index is the paper's assumption | thickness supported, index unverified | 42–45; 1.85 tested |
| dense zone on the ETL side | device image: a low-relief band under the ZADN, 45 nm median (24–97), blocky grains below to the HIL | orientation supported; the 43 % share is from the Si piece, the image suggests less | share 0.15–0.53 scanned |
| lower zone = 51 % perovskite + **void** | cross-section dark fraction 10–30 % (fracture surface, low trust); through-film gaps ≤ 3 % in plan-view TEM; the film is 10 nm nanocrystals in an organic (BPA) matrix inside ~100 nm macro-grains | 50 % void is not supported; the same effective index with an organic host and 10–30 % void is | the "organic matrix" fit (23 % perovskite in a host of n = 1.34) alongside the base fit |
| lower zone scatters | the paper has no haze, transmittance or full absorbance spectrum; 80–150 nm grains, rms 13 nm interfaces | plausible, unverified | followed, never counted as absorption; ×0.5 tested |
| dense zone is dense perovskite | n in the transparent range equals published Cs₀.₁₇FA₀.₈₃PbBr₃ to 0.02 (§11) | supported | — |
| isotropic dipoles, emission spread in depth | the paper's own simulation: both orientations, 20 positions | shared assumption | isotropic; x_d free |
| internal spectrum by reverse calculation | the paper does the same (peak 535 nm) | supported; ours peaks at 535–537.5 nm | — |
| internal radiative efficiency | external PLQE 46 % measured; the 88 % is the paper's flat-film model | model-dependent (§8) | free, ≤ 1 |
| EQE counts forward light only | luminance meter + angular distribution; the glass edges glow in the photographs | supported | substrate light counted as lost |
| target EQE | maximum 28.9 %, mean 25.2 ± 1.6 % (40 devices), ~24 % in the second laboratory | — | both quoted |

Not verifiable from the paper: BPA concentration, film haze, any AFM, the
glass thickness, and whether the Si piece measured here followed the same
recipe as the device films.

On the void fraction specifically. The ellipsometry fixes the effective index
of the lower zone (1.65 at 540 nm, 1.53 at 633 nm), not its composition:
51 % perovskite + 49 % void and 23 % perovskite in a host of n = 1.34 fit
equally (§3). The micrographs rule out the first reading — there is no 50 %
of empty space — and the paper describes nanocrystals bound by excess BPA, so
the second, with part of the host being void (10–30 %), is the physical one.
For the device model the difference is the perovskite content of the lower
zone, i.e. how much of the emission and of the re-absorption happens there
(x_d by volume 0.60 against 0.76); the angular pattern and the EQE barely
move (§8).

## 11. Band edge and the method of Dasgupta *et al.* (arXiv:2601.11793)

"Deriving a comprehensive dataset of optical constants for metal halide
perovskites" (Dasgupta, Hu, Seo, Yuan, Boeije, Johnston, Stranks, Snaith,
January 2026) makes three points that apply here.

1. *Fit one layer at a time, with an intensity measurement in the fit.* They
   characterise the substrate and each underlayer separately and add
   transmittance to the ellipsometric data. Here: GraHIL from its own piece,
   then the perovskite on it; the reflected-intensity ratio of the two pieces
   plays the role of their transmittance.
2. *Ellipsometry leaves an unphysical extinction below the gap, caused by
   roughness, and an effective-medium roughness layer does not remove it.*
   That is the scattering term of §2–3 seen from the other side: here it is
   kept out of the perovskite's own k by giving it its own layer and
   wavelength law, and §7–8 treat it as redistribution, never as absorption.
3. *The band-edge shape must come from a measurement that sees absorption
   only; scale it to the ellipsometric k just above the edge, extend the tail
   log-linearly, recalculate n by Kramers-Kronig.* They use PDS. This film has
   none, but its luminescence is tied to its absorption by detailed balance,
   α(E) ∝ φ_int(E)·e^{E/kT}/(n²E²), with φ_int the internal spectrum of the
   device model. `gp_stitch.py` does the stitching with that:

| λ (nm) | k, ellipsometry | k, stitched | k, published set shifted to this edge |
|---|---|---|---|
| 530 | 0.156 | 0.156 | 0.176 |
| 540 | 0.076 | 0.060 | 0.059 |
| 545 | 0.036 | 0.018 | — |
| 550 | 0.016 | 0.0048 | 0.010 |
| 560 | 0.0027 | 0.00033 | 0.0022 |
| 580 | 0 | 2×10⁻⁶ | — |
| 600 | 0 | 2×10⁻⁸ | 3×10⁻⁶ |

   Stitch point 537.5 nm (k = 0.10); Urbach energy 15.4 meV (theirs, for
   Cs₀.₁₇FA₀.₈₃PbBr₃: 21 meV); n changes by at most 0.019 (at 540 nm:
   2.342 → 2.357). The ellipsometric fit alone has an edge that is too soft
   in its last decade — its shape there is that of the basis functions, not
   of the data — which is exactly their warning. Above the stitch point the
   ellipsometric k is kept: there the luminescence-derived α depends on the
   self-absorption correction of the device model and is not independent.
   `Perov_dense_nk_stitched.csv` is the table to use where re-absorption of
   the red side of the emission matters.

Their public dataset (github.com/akashdasgupta/Perovskite-Dielectric-Constants-Repository)
has Cs₀.₁₇FA₀.₈₃PbBr₃, close to this film's
(FA₀.₇MA₀.₁GA₀.₂)₀.₈₇Cs₀.₁₃PbBr₃:

| λ (nm) | n, dense part (this work) | n, published | k, this work | k, published | k, published shifted |
|---|---|---|---|---|---|
| 400 | 2.296 | 2.175 | 0.525 | 0.423 | 0.368 |
| 450 | 2.333 | 2.144 | 0.324 | 0.332 | 0.313 |
| 500 | 2.283 | 2.161 | 0.208 | 0.299 | 0.307 |
| 520 | 2.303 | 2.186 | 0.179 | 0.303 | 0.287 |
| 633 | 2.106 | 2.107 | 0 | 5×10⁻⁷ | — |
| 800 | 2.022 | 2.006 | 0 | — | — |
| 1000 | 1.985 | 1.968 | 0 | — | — |

- In the transparent range the dense part has the refractive index of bulk
  lead bromide perovskite, to 0.02. It is dense perovskite, and an
  effective-medium roughness layer that pushes it higher (§9) is unphysical.
- Their band edge lies 15–23 nm to the red (half-height of k at 554.5 nm
  against 539 nm): this film's composition and 10 nm nanocrystals put it at
  2.30 eV.
- Just above the edge this film absorbs less (k = 0.18–0.21 at 500–520 nm
  against 0.30) and towards the UV more. Either that is the nanocrystal
  spectrum, which rises steadily where the bulk film has a plateau, or the
  ellipsometry cannot fully separate absorption in the dense part from
  extinction in the lower part over 480–535 nm. The film-average k of the
  paper's absorbance (0.125–0.132 at 520–526 nm) is met either way: 0.129
  with the fitted values (§3), and equally by a bulk-like 0.30 in the dense
  part with almost none below. With η_rad = 1 the EQE of §8 does not depend
  on which (k × 1.5: 25.1 % against 24.8 %).

## Files written to `$ELLIPS_OUT/grahil_perov/`

| file | content |
|---|---|
| `Perov_dense_nk.csv` | dense perovskite n,k, 205–1000 nm, 1 nm; min/max over the fit variants; lower-layer effective n,k |
| `GraHIL_nk.csv` | GraHIL n,k |
| `lib_Perov_dense.csv`, `lib_Perov_lower_with_scattering.csv`, `lib_Perov_lower_no_scattering.csv`, `lib_GraHIL.csv` | the same media as `wavelength_nm,n,k` |
| `structure.json` | thicknesses, every variant, the tables above |
| `grahil_perov_summary.png` | fit at three angles, n, k, device angular emission |
| `Perov_dense_nk_stitched.csv`, `lib_Perov_dense_stitched.csv`, `stitch.json` | dense perovskite with the luminescence-stitched band edge (§11) |
| `grahil_perov_device.png` | angular emission and photon budget per description, roughness scan, band edge |
| `../grahil_perov_mod/*.mod`, `mod_check.json` | the two CompleteEASE models of §6 and their read-back check |

For a device simulation use two perovskite layers, not one: `lib_Perov_dense`
on the ETL side (0.43 of the perovskite thickness) and a `lib_Perov_lower`
table on the HIL side (0.57). "with scattering" reproduces the specular
(coherent) field; "no scattering" is the same effective index without the
loss, for codes that handle scattered light separately.
