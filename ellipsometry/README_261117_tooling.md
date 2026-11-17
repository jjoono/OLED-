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
python tl_thin_diag.py; python tl_thin_diag2.py    # why 2_23 / 2_24 fit nothing (first pass)
python tl_bphen10_signal.py; python tl_bphen10.py; python tl_bphen10_sub.py; python tl_bphen10_sub2.py
python tl_bphen10_back.py; python tl_bphen10_plot.py   # Bphen 10 nm, in depth (section 3)
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
python tl_mod_report.py   # parameter table, n,k CSVs, figure -> $ELLIPS_OUT/tooling_mod/

# The lab's tooling template (§7) - ELLIPS_DATA also holds BPhen_S12_tooling.SEsnap
python tl_s5_extract.py   # the model inside the snapshot -> $ELLIPS_OUT/S5_template.mod
python tl_s5.py           # its Bphen and no5 glass against this set
python tl_ce_glass2.py    # glass pieces with a backside fraction per piece
python tl_ce_thickness_ntve0.py   # Si thickness in the template's NTVE 0 A convention
python tl_make_mod_glass2.py      # glass .mod files on the template's structure
python tl_make_mod_ntve0.py       # Si .mod files at NTVE 0 A, for like-with-like tooling
```
(`tl_make_mod_glass.py`, the first, inferred-format glass writer, is superseded.)

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
| 2_23, 2_24 / 2_33, 2_34 | Bphen | 10 | none detectable | none, 3.4 |

These are this pipeline's own frame (2 nm Malitson oxide). The numbers in the
`.mod` files - CompleteEASE frame, measured 33.68 Å oxide, and the lab
template's NTVE 0 Å convention - are tabled in §7.

Si: ±1 nm for a ±1 nm native-oxide assumption, which is the dominant
uncertainty. Glass: the spread over the three depolarization conventions below
is up to 1.8 nm. Si and glass pieces agree to ±2 nm once the glass model is
right (§4).

## 3. The Bphen 10 nm pieces: no 10 nm film on any of the four

Re-examined in depth (`tl_bphen10*.py`) because the data matter. Every picture
below was fitted to all four pieces, in CompleteEASE's frame for the Si pieces.

**The measured signal is not a Bphen signal.** Subtracting the bare-substrate
model from 2_23 leaves a Ψ change of up to −6.4° with its extrema at Si's E1
(365 nm) and E2 (~290 nm) critical points (`tl_bphen10_signal.py`). Any Bphen
film raises Ψ there; NPB at the same nominal 10 nm (1_23) follows its 10 nm
prediction closely.

| 2_23 / 2_24 on Si | MSE | what it says |
|---|---|---|
| Bphen film, ± roughness, ± void fraction | 41–42 | roughness → 0, void → 0, d ≈ 3.6 nm |
| uniaxial Bphen (order parameter) | 40–41 | S pinned at −0.5, no help |
| islands next to bare substrate (incoherent area mix) | 30 | fixes 500–1000 nm only |
| oxide thickness free | 41–43 | no help |
| **rough / damaged Si surface, no Bphen** | **13–14** | ~8 nm Si + 34 % void |
| the same + a Bphen overlayer | 13–14 | overlayer → 0 |
| rough Si with the valleys filled by Bphen | 18–19 | worse than air-filled |
| islands + Si damage | 17 | 3.7 % coverage → ~2 nm Bphen-equivalent |

The Bphen 50 nm pair from the same run (2_21, 2_22) does **not** want that
substrate layer (it fits to 0.04–0.06 nm), so 2_23/2_24 sit on a different
surface - another wafer, or the other side of one.

| glass | best | front film | back-surface film only |
|---|---|---|---|
| 2_33 | bare glass | 0.5 nm | 0 (MSE 11) |
| 2_34 | Bphen 3.4 nm on the measured face (MSE 7.1; 3–6 nm across all models) | 3.4 nm | 0 (MSE 18) |

The glass test rules out the tidy explanation (deposit on the far face, the
measured face bare): the back-film fits are worse on both pieces.

**Reading:** these four pieces did not receive a 10 nm Bphen film. Best
estimates - 2_23, 2_24: none detectable, ≤ 2-4 nm Bphen-equivalent depending on
the substrate model; 2_33: none; 2_34: ~3.4 nm. The 50 nm Bphen pieces came in at
0.86-0.92 of nominal, so this is not a tooling factor; it points at the
deposition or the sample handling (shutter, mask, plume position, labels). What
would settle it: which face of 2_23/2_24 was measured and whether it is the
polished one; the depolarization export of these four from the .SE files (rough
Si and island scattering both depolarize); an optical-microscope or AFM look for
crystallites; a bare piece of the 2_23/2_24 wafer.

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

One f for all pieces is too stiff: the lab template carries `% 1st Reflection`
per sample, and freeing it per piece (`tl_ce_glass2.py`, glass still shared)
gives f = 0.24–0.77 on the 1_3x pieces and 0.68–1.0 on the 2_3x pieces, with
2_31 going from MSE 8.2 to 3.5 and 2_34 from 7.1 to 4.7. The collected fraction
depends on how each piece sat on the stage (backside contact, size against the
beam), so a per-piece value is physical.

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

**Glass `.mod` files**: the first set inferred the glass-substrate format and
was marked unverified. It is replaced by files built on the lab's own
CompleteEASE-saved glass-capable model - see §7.

## 7. The lab's tooling template (`BPhen_S12_tooling.SEsnap`)

A `.SEsnap` is a zip of the data file and `_FitLog`; the fit log carries the
whole model ('Bphen_S5') as CompleteEASE saved it. `tl_s5_extract.py` takes it
out unchanged (one level of tab nesting removed) as `S5_template.mod`. No
depolarization export came with it.

**What the template is.** Substrate: a Bruggeman EMA of SI_JAW and a Gen-Osc
glass ('no5_glass_sub_genosc': Einf 1.42, one Tauc-Lorentz with Eg 3.27 eV, UV
and IR poles), 1 mm thick, so one model serves Si (`EMA % (Mat 2)` = 0) and
glass (100). Then NTVE_JAW at 0 Å, then Bphen as Gen-Osc (Einf 0.66, Gaussians
at 315 / 282 / 225 nm, UV pole 186 eV² at 11 eV). Back reflections switched
on, `% 1st Reflection` per sample. Two samples fitted together, 54.8 nm and
13.8 nm of Bphen, both at EMA % 0, i.e. both on Si.

**Bphen.** The lab dispersion and this set's agree to 0.002 in n across the
visible (n 1.778 / 1.731 / 1.712 vs 1.780 / 1.733 / 1.713 at 450 / 550 / 633 nm;
k(350) 0.020 vs 0.013). That is an independent check of §5/§6. On 2_21 / 2_22
the lab Bphen gives MSE 7.8 vs 6.4 / 6.9 with this set's, and thicknesses ~0.4 nm apart.

**Bphen 10 nm.** The template model reads 2_23 / 2_24 the same way §3 does: MSE
41–43 at any thickness, its best d 3.5–6.8 nm depending on the oxide, and a
roughness layer goes to 0. The same model was used for the 13.8 nm film in the
snapshot, so a real ~10 nm Bphen film is within what it reads; §3 stands.

**Glass.** The lab's no5 glass has n lower by 0.014 at 550 nm and a weak UV
absorption. With the lab Bphen it fits 2_31 / 2_32 / 2_33 nearly as well as
this set's glass (MSE 4.5 / 4.4 / 9.0 vs 3.8 / 4.5 / 8.8, thickness within
1 nm), but not 2_34 (8.3 vs 4.7). This set's shared glass is kept. Dropping the
no5 absorption (`no5 k=0`) makes 2_31 / 2_32 worse and leaves 2_33 / 2_34 about
the same, so a weak UV absorption in this glass is possible but not
established; this set's glass is left lossless.

> Correction: the first run of this comparison reported MSE > 370 for no5.
> The backside model handed the absorbing glass to `tmm()` as the incidence
> medium. The Snell invariant is then complex, and `tmm()`'s Im(q) ≥ 0 branch
> rule flips q in air, so even k = 1e-9 gave a back-surface |r|² of ~80.
> `tl_glass.ncs_glass` now uses Re(n) at the interfaces and Beer-Lambert
> through the 1 mm bulk. The shared `lib/tmm_fix.py` is unchanged. Every other
> glass result here used a lossless Cauchy (k exactly 0), so none of them change.

**Glass `.mod` files** (`tl_make_mod_glass2.py`, `$ELLIPS_OUT/tooling_mod/glass/`)
are `S5_template.mod` with values changed:

| field | template | written |
|---|---|---|
| `EMA % (Mat 2)` | 0, fit on | 100 (glass), fit off |
| glass material | no5: Einf 1.42, TL Amp 0.0137, Eg 3.27, UV pole 99.3 @ 11, IR pole 0.016 | renamed `tooling2611_glass`: Einf 1, UV pole 138.5 eV² @ 10.55 eV (≤ 0.0023 from the fitted Cauchy), IR pole 0; TL kept but parked: Amp 0.001 (its lower bound), Eg 6.0, Eo 7.0 |
| NTVE_JAW | 0 Å, fit on | 0 Å, fit off |
| film | Bphen Gen-Osc | this set's NPB / TAPC / Bphen Gen-Osc (§6), d ± 10 nm free |
| `# Back Reflections` | 0 | 5, fit off |
| `% 1st Reflection` | 100 | 100·f of the piece, fit on |
| MultiSamp | 2 samples | single-sample block, copied from the ADNC file |
| angle offset, roughness | 0 | 0, fit off, bounds narrowed |

The TL is parked rather than left at Amp 0.001, Eg 3.27: even that gives k ~ 1e-6
below 379 nm, which over 1 mm takes ~16 % off the 300 nm backside beam, and
these fits assumed a lossless glass. With Eg above the data (≤ 5.04 eV) the
glass has k = 0 everywhere measured. Each file is read back and re-scored; all
eight reproduce their fit MSE within 0.15:

| piece | material | d (nm) | % 1st Reflection | MSE |
|---|---|---|---|---|
| 1_31 | NPB | 29.18 | 41.4 | 13.98 |
| 1_32 | NPB | 48.22 | 23.6 | 18.80 |
| 1_33 | NPB | 10.35 | 76.7 | 6.22 |
| 1_34 | TAPC | 48.94 | 27.0 | 8.97 |
| 2_31 | Bphen | 46.53 | 71.1 | 3.52 |
| 2_32 | Bphen | 47.25 | 67.8 | 4.06 |
| 2_33 | Bphen | 0.74 (bare) | 100 | 8.27 |
| 2_34 | Bphen | 4.41 | 98.9 | 4.77 |

Not verifiable here: the template's own samples were on Si, where the
backside does nothing, so how CompleteEASE counts `# Back Reflections` and
applies `% 1st Reflection` has not been seen on glass. Load 2_31 and press
Generate. If the MSE is far from 3.5, those two fields are the first suspects.

**Oxide convention, and what the glass twins say.** The template fits Si with
NTVE 0 Å; §6 uses the measured 33.68 Å. Both are given
(`tl_ce_thickness_ntve0.py`; `.mod` files in
`tooling_mod/Si_NTVE0_template_convention/`, all six re-scored OK):

| Si / glass | material | Si, NTVE 33.68 Å: d (MSE) | Si, NTVE 0 Å: d (MSE) | glass twin: d (MSE) |
|---|---|---|---|---|
| 1_21 / 1_31 | NPB | 25.94 (9.7) | 29.28 (14.3) | 29.18 (14.0) |
| 1_22 / 1_32 | NPB | 46.44 (9.3) | 49.86 (12.5) | 48.22 (18.8) |
| 1_23 / 1_33 | NPB | 8.33 (3.5) | 11.55 (5.9) | 10.35 (6.2) |
| 1_24 / 1_34 | TAPC | 47.16 (3.7) | 50.58 (6.1) | 48.94 (9.0) |
| 2_21 / 2_31 | Bphen | 43.94 + 1.9 rough (3.8) | 47.76 + 0.8 rough (6.4) | 46.53 (3.5) |
| 2_22 / 2_32 | Bphen | 42.87 + 2.1 rough (4.1) | 46.65 + 1.1 rough (6.7) | 47.25 (4.1) |

The two pulls disagree. MSE prefers 33.68 Å on every Si piece (by 35–67 %).
The glass twins, which have no oxide question, land on average 0.9 nm below the
NTVE 0 Å values and 2.6 nm above the 33.68 Å ones. If each Si / glass pair saw
the same flux, that points to a thinner oxide (the 2 nm of §2 sits between).
The Si data themselves never ask for 0 Å: fitting d and the oxide together
(the NPB UV bands make that possible on a thin film, §8) returns 3.41 nm on
1_23 - the lab value - and 4.2-4.6 nm on 1_21 / 1_22, where the
shared-dispersion compromise of §5 leaks in. For
the 10 nm NPB the choice moves the tooling factor by 39 %, so the oxide on these
wafers is the open question. A bare piece of the same wafer lot, measured with
these settings, would settle it. Until then, compare tooling factors with earlier
template runs using the NTVE 0 Å column, and quote the glass twins as the
oxide-free check.

## 8. Src10_10nm: a nominal 10 nm organic film on Si, material not given

`tl_src10.py` (`SRC10_XLSX=.../Src10_10nm.xlsx`; same export layout, one
sheet). Same CompleteEASE frame and footing as the tooling set; the three NPB
tooling pieces are refitted alongside as the like-for-like reference.

**It is NPB, at NPB band strength.** Its Ψ at 65° tracks 1_23 (8.3 nm NPB)
point for point from 260 to 900 nm. Against every picture tried:

| model | d (nm) | MSE |
|---|---|---|
| NPB film (oxide 33.68 Å) | 9.31 | 4.00 |
| NPB film, band amplitudes × s free | 9.41 (s = 0.91) | 3.91 |
| NPB film, oxide free | 7.98 (oxide 4.76 nm) | 3.47 |
| NPB film + Si/void EMA layer | 7.39 (+ 2.75 nm at 81 % void: a low-index interlayer) | 2.99 |
| NPB film, NTVE 0 Å (template convention) | 12.54 | 7.47 |
| TAPC / Bphen film, best variant | 5.1 / 5.4 | 8.9 / 10.1 |
| transparent Cauchy film | 8.6 | 11.2 |
| no film: bare Si / Si damage only | - | 13.2 / 10.9 |

Only three dispersions are on hand, so "NPB" means "absorbs where and as
strongly as NPB" - α-NPD-type hole transporters would read the same.

**Thickness: 9.3 nm** in the convention of the tooling set (§6 files, where
1_23 is 8.3 nm). Statistical 90 % intervals are ±0.01-0.06 nm; the spread is
the interface. Unlike the tooling NPB pieces, Src10 asks for a little more
than the lab oxide: freed, the oxide goes to 4.8 nm (on 1_23 the same fit
returns 3.41 nm, the lab value), and a Si/void EMA layer, which comes out
81 % void (oxide-like, not damaged Si), does the same job.
Either way the film reads 7.4-8.0 nm. Read 9.3 nm as the comparable number
and 7.4-9.4 nm as the range. A different wafer batch or clean would explain a
thicker oxide; a bare piece of the same wafer would settle it.

`tl_make_mod_src10.py` writes the NPB Gen-Osc `.mod` at 33.68 Å (9.3 nm, MSE
4.00) and at NTVE 0 Å (12.5 nm, MSE 7.47) to `tooling_mod/Src10/`; both read
back and re-score to the same MSE.
