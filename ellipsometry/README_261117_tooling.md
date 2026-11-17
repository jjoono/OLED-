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

# Gen-Osc .mod files (§6) - ELLIPS_DATA also holds the two lab .mod files below
python tl_genosc.py       # constrained Gen-Osc in this pipeline's own frame (method check)
python tl_ce_genosc.py    # the same in CompleteEASE's frame: SI_JAW + NTVE_JAW 33.68 A
python tl_ce_thickness.py # every continuous Si piece, at the lab oxide and at 15 A
python tl_ce_glass.py     # glass pieces with those dispersions; glass as Einf 1 + UV pole
python tl_make_mod.py     # Si .mod files, each read back and re-scored
python tl_make_mod_glass.py   # glass .mod files - unverified format, see §6
python tl_mod_report.py   # parameter table, n,k CSVs, figure -> $ELLIPS_OUT/tooling_mod/
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

## 6. Gen-Osc model files for CompleteEASE

`tl_make_mod.py` writes one `.mod` per continuous-film piece. Structure is the
lab's own CompleteEASE-saved Gen-Osc file (`ADNC_ITO_GenOSC_v2_claude.mod`,
SI_JAW / NTVE_JAW / Gen-Osc, model v2); only the Gen-Osc layer, thickness,
oxide, roughness and angle offset are rewritten. Fits are redone in that frame
— SI_JAW and NTVE_JAW decoded out of the file, oxide fixed at 33.68 Å from
`Nativeoxide_180sO2plasma_Si_wafer.mod` — so Generate alone reproduces the data.
Every file is read back and re-scored; all twelve reproduce their fit MSE.

**What "realistic for an organic" was taken to mean**, and enforced:

| | rule | why |
|---|---|---|
| Einf | fixed at 1 | the UV pole carries the transitions above the data; a free Einf trades with it and means nothing |
| UV pole | 6–12 eV, amplitude ≥ 0 | σ→σ*/π→π* above the measured 5.06 eV |
| Gaussians | amplitude ≥ 0, FWHM 0.10–1.20 eV, centre above the measured absorption onset | vibronic bands; no negative oscillators cancelling each other |
| sub-gap | ε₂ penalised below the onset during the fit | a broad band centred above the gap can still leak a visible tail; these are wide-gap films |
| count | fewest Gaussians within 10 % of the best MSE | NPB 5, TAPC 4, Bphen 4 |
| locking | fit flag OFF if the 90 % CI > 50 % or on a bound | only TAPC G4 En (on its 5.6 eV bound) |

Results, CompleteEASE frame:

| material | UV pole | bands (nm) | n(550) | k(450) | MSE |
|---|---|---|---|---|---|
| NPB | 76 eV² @ 8.27 eV | 384, 359, 337, 281, 223 | 1.772 | 5e-4 | 9.3, 9.7 |
| TAPC | 84 eV² @ 8.21 eV | 320, 311, 285, 221 | 1.649 | 0 | 3.7 |
| Bphen | 95 eV² @ 8.54 eV | 315, 281, 276, 228 | 1.733 | 0 | 3.8, 4.1 |

The lab oxide beats a 15 Å bare-wafer oxide on every Si piece (e.g. 1_23 MSE
3.50 vs 4.39, 1_24 3.66 vs 4.65), so the tooling wafers look plasma-treated like
the reference wafer; at 15 Å every film would be ~1.9 nm thicker. Roughness was
tested per material: NPB 0, TAPC 0.3 nm (no gain), Bphen ~2 nm (11 % better MSE)
— written only for Bphen, which is also the material that crystallizes. The
lab's BsBCZ models put roughness at 9–10 nm, pinned on their 100 Å bound;
nothing here asks for that.

**Glass `.mod` files are unverified.** No glass-substrate model was supplied, so
the substrate (a Gen-Osc layer with no oscillators: Einf 1 + UV pole, |Δn| ≤ 0.003
from the fitted glass Cauchy) and the back-surface switch (the flag before
`'# Back Reflections'`, with `'% 1st Reflection'` = 46 and free) are inferred.
They round-trip in this pipeline's model; whether CompleteEASE reads them the
same way has to be checked by loading one.
