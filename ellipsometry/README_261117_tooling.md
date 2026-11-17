# Session 2026-11-17 — tooling set: NPB, TAPC, Bphen on Si and glass

`scripts/tooling_261117/` fits `tooling.xlsx`: 16 pieces, AnE Ψ/Δ at 45–65°,
245.8–1688.1 nm, no depolarization export. `*_21–24` are Si wafer, `*_31–34`
glass. Material and thickness were not given; a guessed map was supplied
afterwards and is checked here rather than assumed.

```bash
export ELLIPS_DATA=/path/to/data      # tooling.xlsx
export ELLIPS_OUT=/path/to/output
cd scripts/tooling_261117
python tl_load.py         # inventory
python tl_cauchy_si.py    # Si pieces: exhaustive (d, A) grid, then refine - n.d alias check
python tl_pointnk.py      # model-free n,k at that thickness -> absorption fingerprints
python tl_osc.py          # KK-consistent Gaussian dispersion per material
python tl_npb_sep.py      # NPB fitted per piece, to see the shared-dispersion compromise
python tl_si_id.py        # every Si piece against every material, d only
python tl_thin_diag.py; python tl_thin_diag2.py    # why 2_23 / 2_24 fit nothing
python tl_glass.py        # glass pieces, front surface only
python tl_glass_sub.py    # backside reflection test
python tl_glass_global.py # one glass + one backside fraction for all glass pieces
python tl_tail.py         # is the k ~ 0.01 visible tail of TAPC / Bphen required?
python tl_report.py; python tl_notail_csv.py       # -> $ELLIPS_OUT/tooling_nk/
```

Reuses the corrected TMM and reference constants of `joint_nk_260910/`
(Green 2008 Si, Malitson native oxide fixed at 2 nm) and the CompleteEASE
Gaussian of `lib/ce_osc.py`.

## 1. Identification is from the data

On the thick Si pieces the model-free k shows one absorption band per material,
and the bands do not overlap:

| piece | main band | refractive index n(550) | → |
|---|---|---|---|
| 1_21, 1_22 | 346–353 nm, weaker 280–295 nm | 1.78 | NPB |
| 1_24 | 308 nm | 1.66 | TAPC |
| 2_21, 2_22 | 277 nm, shoulder ~315 nm | 1.73 | Bphen |

The NPB dispersion lands on the lab library `data/nk/NPB.csv` across its whole
400–700 nm range (n(450) 1.882 vs 1.886, n(550) 1.782 vs 1.784). Every
thickness search is an exhaustive grid first; the runner-up thickness branch is
MSE > 100 for all five thick pieces, so no n·d alias is possible. Scoring every
Si piece against every dispersion, the right material wins by a factor 15–90 in
MSE on the thick pieces and 3× on 1_23 (NPB 10 nm). **The guessed map holds for
every piece that carries a continuous film.**

## 2. Thickness

| pieces | material | nominal | Si | glass |
|---|---|---|---|---|
| 1_21 / 1_31 | NPB | 30 | 27.3 | 29.3 |
| 1_22 / 1_32 | NPB | 50 | 47.5 | 49.3 |
| 1_23 / 1_33 | NPB | 10 | 9.9 | 9.4 |
| 1_24 / 1_34 | TAPC | 50 | 48.6 | 49.2 |
| 2_21, 2_22 / 2_31, 2_32 | Bphen | 50 | 47.0, 46.0 | 45.5, 46.4 |
| 2_23, 2_24 / 2_33, 2_34 | Bphen? | 10 | — | (0.7, 3.8) |

Si: ±1 nm for a ±1 nm native-oxide assumption, which is the dominant
uncertainty. Glass: the spread over the three depolarization conventions below
is up to 1.8 nm. Si and glass pieces agree to ±2 nm once the glass model is
right (§4).

## 3. The Bphen 10 nm pieces are not a continuous film

2_23 and 2_24 fit **no** homogeneous layer — any of the three dispersions, free
oxide, a Bruggeman void fraction, an angle offset: MSE 34–42 against 4–10 for
everything else, residual in every band including 500–1000 nm. A free
absorbing Cauchy does best (MSE ~21) only by going to n = 1.2 with a flat
k ≈ 0.08 — a low-density, light-losing layer, which is what a dewetted or
crystallized island film looks like to a specular model. The glass twins say
the same from the other side: 2_33 reads as bare glass (0.7 nm equivalent),
2_34 as 3.8 nm, both with a Bphen-like 277 nm feature.

NPB at the same nominal 10 nm (1_23) is a perfectly ordinary film (MSE 4). Bphen
(T_g ≈ 62 °C) is known to crystallize where NPB (T_g ≈ 95 °C) does not, so the
likely story is morphology, not deposition — but that is an inference; an AFM or
optical-microscope look at 2_23 would settle it. Either way, no tooling factor
can be taken from these four pieces.

## 4. Glass needs its back surface

Front-surface-only fits leave MSE 13–20 and put every glass piece ~10 % thinner
than its Si twin, with the fitted glass index scattering 1.51–1.57 from piece to
piece. Adding the incoherent back-surface beams (a fraction f reaching the
detector, summed as Mueller matrices in N,C,S) brings the near-bare piece 2_33
from MSE 15 to 9 and the transparent-region residual from 13 to 4.

Fitting all eight glass pieces together — one glass Cauchy, one f — gives
n(550) = 1.517, f = 0.52, and Si/glass agreement within ±2 nm. Without a
depolarization export it is not known how the instrument turned the
depolarized N,C,S into Ψ/Δ; three readings were tried and the data prefer
normalising the N,C,S vector (MSE 10.5, vs 10.9 raw and 13.1 for Ψ from N
alone). The glass thicknesses carry that spread.

## 5. Dispersions, and what to trust in them

* **NPB** is fitted to the 30 nm and 50 nm pieces together. Fitted alone, the
  30 nm piece has a ~15 % weaker 350 nm band; the shared fit is the compromise
  (MSE 7.5 / 10.6). Both per-piece versions are delivered.
* **TAPC / Bphen** show k ≈ 0.01 in 400–550 nm. Forcing that tail to zero costs
  only 9 % / 17 % in MSE and moves thickness ≤ 0.5 nm, and k ≈ 0.01 is the
  detection limit for these films — so the tail is not established. Both
  versions are delivered (`*_no_visible_tail.csv`); use that one for device
  absorption unless there is independent evidence of sub-gap loss.
* Isotropic throughout. NPB in particular can be weakly uniaxial; nothing in
  45–65° data on 10–50 nm films constrains that.

Outputs (`$ELLIPS_OUT/tooling_nk/`): `thickness_summary.csv`,
`nk_{NPB,TAPC,Bphen}.csv` and variants on a 1 nm grid 250–1000 nm,
`n_glass_substrate.csv`, `fit_parameters.json`, two figures.
