# Session 2026-09-10 — n,k re-extracted from the SE and T/R data together

`scripts/joint_nk_260910/` re-does the thin-Ag optical constants from the two
2026-08 measurement sets at once:

* **SE** — `summary.xlsx`, 16 Si pieces, 5 angles 45–65°, 675 points
  245.8–1688.1 nm, Ψ/Δ + `% Depolarization`.
* **T/R** — the tracked 260820 campaign,
  `dft_seedlayer_screen/data/TR_20260820/raw/`, 2 nm grid, read through
  `tr_crosscheck/build_tra.py` so the reflectance correction has one owner.

Only the SE workbook needs `ELLIPS_DATA`; the T/R data is in the repository.

```bash
export ELLIPS_DATA=/path/to/data       # summary.xlsx
export ELLIPS_OUT=/path/to/output
cd scripts/joint_nk_260910
python jnk_data.py      # inventory + the sample map, run this first
python jnk_glass.py     # substrate from its own T/R
python jnk_seed.py      # seed dispersions from the Ag-free pieces
python jnk_ag.py        # seed-thickness scan, then the Ag fits
python jnk_tr.py        # thickness and model-free n,k from the glass side
python jnk_joint.py     # one dielectric function for both measurements
python jnk_report.py    # n,k CSVs, summary table, figures
```

---

## 1. Why this is a new pipeline and not a re-run of `thin_ag_260819/`

The 260819 scripts fit in CompleteEASE's own frame: `SI_JAW` and `NTVE_JAW` are
decoded out of a `.mod` file (`lib/ce_mat.py`) and `ce_fit._mat` serves them to
every fit. No `.mod` is available here, so the substrate is carried by published
tabulated data instead — Green 2008 for Si, Malitson for the native oxide
(`jnk_ref.py`, both CC0 and embedded in the source). Everything else — the
corrected TMM, the CompleteEASE-convention oscillators of `lib/ce_osc.py`, the
N/C/S residual, the Bruggeman roughness layer — is shared with the earlier work.

The frames turn out to agree. Fitting the Ag-free HATCN piece here gives

| | 260819 (CompleteEASE frame) | this pass (Green + Malitson) |
|---|---|---|
| d_seed | 6.42 nm | 6.68 nm |
| Einf | 2.65 | 2.69 |
| TL Br | 0.54 | 0.544 |
| TL Eo | 3.95 | 3.965 |
| TL Eg | 2.54 | 2.65 |

so the substrate choice is not what limits this analysis. `NTVE_JAW` is not
SiO₂ (n = 1.734 vs 1.457 at 633 nm), so only the *sum* of oxide and seed is
comparable between the two frames, not the oxide thickness itself.

## 2. Sample map — the two sets use different piece numbers

The glass piece and the Si piece of the same deposition are **four apart**
(glass `1-2` ↔ Si `1-6`). The glass column below is the campaign handoff's own
mapping, carried by `tr_crosscheck/build_tra.py`; the Si column follows from
`thin_ag_260819/ag_load.py`. The T/R series agrees with it independently —
absorption at 550 nm falls monotonically 14.6 → 6.3 % across HATCN Ag 4 → 12 nm
and the Ag-free piece (`1-9`) reflects like bare glass.

| seed | Ag (nm) | 0 | 4 | 5 | 6 | 7 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|---|---|---|
| HATCN | Si (SE) | 1-5 | 1-6 | 1-7 | 1-8 | 2-5 | 2-6 | 2-7 | 2-8 |
| HATCN | glass (T/R) | *(1-1, not measured)* | 1-2 | 1-3 | 1-4 | 2-1 | 2-2 | 2-3 | 2-4 |
| MoOx | Si (SE) | 1-13 | 1-14 | 1-15 | 1-16 | 2-13 | 2-14 | 2-15 | 2-16 |
| MoOx | glass (T/R) | 1-9 *(R only)* | 1-10 | 1-11 | 1-12 | 2-9 | 2-10 | 2-11 | 2-12 *(T only)* |

The workbook is not uniform — sheet `1-5` has an extra leading column in the
depolarization block and `1-8` an extra header row — so `jnk_data` locates the
header row and every block by name. (`thin_ag_260819/ag_load.py` hard-codes
columns 51–55, right only for `1-5`; every consumer there discards the
depolarization array, so its fits are unaffected and only its own diagnostic
table is wrong.)

## 3. Reflectance is always the corrected one — and the correction now has a direct check

The UMA collects only f ≈ 0.84 of the substrate's back-surface beam at 6°, so
the raw R export is low by (1−f)·R_g·(T/100)² — 0.59 %p on bare glass, 0.2 %p at
12 nm of Ag. `build_tra.py` owns that constant and this pipeline imports it
rather than repeating it; the corrected R computed here matches the delivered
`ALL_SAMPLES_TRA.csv` to **0.012 %p** at worst.

f had been recovered by least squares against the delivered file, because the
bare-substrate exports were not part of the delivered set. They are now
(`raw/glass.csv`, `raw/glassR.csv`), and they check f directly — a transparent
slab must show zero absorptance:

| bare glass, 450–700 nm | 1 − T − R |
|---|---|
| raw R | **+0.581 %p** |
| corrected R | **−0.007 %p** (campaign noise 0.15 %p) |

So the recovered f = 0.837 is right, and the substrate is lossless over the
working window. Fitting the raw R instead makes the glass look 0.5 % absorbing
and biases every intensity-derived thickness.

With that settled the substrate index follows from the bare scan, a two-term
Cauchy fitted to both channels over 430–780 nm:

```
n(λ) = 1.52354 + 0.002989/λ²   (λ in µm)   →   n(550) = 1.5334
```

reproducing the bare-glass T to 0.157 %p and R to 0.085 %p. That is the floor on
every absolute-intensity claim below.

Usable T/R window is **420–780 nm**: soda-lime absorbs below ~400 nm, and the
instrument changes source/detector near 800 nm, where T drops ~0.7 %p and
1−T−R jumps by 0.5 %p.

## 4. Seeds

Ag-free pieces, air / seed (Einf + Tauc-Lorentz) / oxide / Si, 260–1080 nm.
The seed **must** be allowed to absorb: with a transparent Cauchy the residual
below 400 nm is 4–8× the visible one and looks exactly like a bad substrate,
and with the Tauc-Lorentz it is flat across the whole range (MSE 1.5–2.5).

| seed | d (nm) at 2 nm oxide | n(450) | n(550) | n(633) |
|---|---|---|---|---|
| HATCN | 6.68 | 1.943 | 1.877 | 1.855 |
| MoOx | 7.20 | 1.948 | 1.904 | 1.887 |

HATCN n(550) = 1.877 sits between the 260819 SE fit (1.867) and the vendor
library value (1.849). Seed and oxide thickness are interchangeable at this
thinness — over d_ox 0→4 nm the seed moves 8.2→5.4 nm at essentially constant
MSE — so a seed thickness is only meaningful together with the oxide it assumes.

## 5. Ag from the SE data

Stack air / roughness (Bruggeman 50 % Ag + void) / Ag / seed / oxide 2 nm / Si,
Ag = Einf + Drude(RT) + 3 Gaussians, 260–1080 nm.

With Ag pinned at the deposited thickness and one common seed thickness scanned
per seed, the **HATCN pieces want +1.0 nm of seed over the bare-seed fit**
(6.68 → 7.68 nm) — the same +1 nm the 260819 pass found. MoOx wants none
(7.20 nm), but its MSE is 3–4× worse at every offset, which is the first sign
that a single homogeneous layer is the wrong description there.

Freeing everything at that seed thickness:

| seed | Ag nom (nm) | 4 | 5 | 6 | 7 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|---|---|
| HATCN | d fitted | 4.02 | 5.37 | 6.18 | 6.84 | 7.87 | 8.97 | 10.71 |
| HATCN | n(633) | 0.839 | 0.519 | 0.164 | 0.112 | 0.111 | 0.103 | 0.180 |
| HATCN | ħω_p (eV) | 8.66 | 8.68 | 8.81 | 9.28 | 9.42 | 9.49 | 9.69 |
| MoOx | d fitted | 5.90 | 7.44 | 7.92 | 8.14 | 9.18 | 10.34 | 12.92 |
| MoOx | n(633) | 2.401 | 1.611 | 1.278 | 0.850 | 0.554 | 0.522 | 0.222 |
| MoOx | ħω_p (eV) | 6.35 | 6.51 | 6.39 | 7.48 | 7.86 | 8.01 | 7.85 |

Taking the densest film (HATCN 12 nm) as the reference, the effective density
(ħω_p/ħω_p,ref)² is **0.80 → 1.00 across the HATCN series but only 0.43 → 0.68
across MoOx**, and MoOx thicknesses run 8–50 % over nominal while HATCN sits
within ±10 %. Both say the same thing: Ag on HATCN is continuous and near-bulk
by 6–7 nm, Ag on MoOx is still porous at 12 nm. MSE follows — 1.8–2.8 for
HATCN, 4.3–7.3 for MoOx, where a single homogeneous layer is the wrong model.

## 6. One dielectric function for both measurements

`jnk_joint.py` fits one Einf + Drude + 3 Gaussians to the Si-piece Ψ/Δ **and**
the glass-piece absolute T and R at once, giving each piece its own thickness
and roughness. Every sample is fitted three ways and each block's residual is
reported against all three. In the SE-only row the glass piece is forced to the
Si piece's geometry — i.e. what the ellipsometry predicts for the other
substrate if the two pieces really got the same film.

| Ag on HATCN | SE MSE | mean \|ΔT\| | mean \|ΔR\| |
|---|---|---|---|
| SE-only solution, same geometry on both | 1.8–2.8 | 0.9–6.0 %p | 0.3–1.1 %p |
| joint | 2.3–3.2 | **0.16–0.24 %p** | **0.10–0.18 %p** |
| T/R-only | 54–196 | 0.05–0.98 %p | 0.04–1.79 %p |
| *Ag on MoOx, joint* | *5.2–11.1* | *0.4–1.2 %p* | *0.4–0.9 %p* |

Three things come out of this.

**For HATCN the two measurements are consistent.** One dispersion reproduces
the ellipsometry within ~0.5 MSE of the SE-only fit *and* the absolute T/R to
0.16–0.24 %p — which is the bare-substrate calibration floor of §3
(0.157 %p in T, 0.085 %p in R). The joint n(550) is also monotone across the
series (0.48 → 0.10 for Ag 4 → 12 nm) where the SE-only values wobble
(0.53, 0.47, 0.14, 0.09, 0.09, 0.08, 0.14).

**The SE-vs-T/R gap in n is real but smaller than it looked, and it is a
compromise, not a contradiction.** Against the repository's own T/R inversion
in `data/nk/Ag<d>nm_on_HATCN5_measured.csv`, at 550 nm:

| Ag (nm) | n: SE / joint / T/R inversion | k: SE / joint / T/R inversion |
|---|---|---|
| 5 | 0.472 / 0.272 / 0.528 | 3.475 / 3.243 / 3.528 |
| 7 | 0.088 / 0.184 / 0.282 | 3.588 / 3.360 / 3.479 |
| 8 | 0.087 / 0.153 / 0.228 | 3.675 / 3.450 / 3.574 |

k has always agreed to a few percent. n did not, and the joint solution lands
between the two — close enough to each that both data sets are reproduced at
their own noise level. The 260820 notes read the gap as light scattered out of
the specular beam, which a 2-observable inversion must book as absorption; the
joint fit shows the geometry difference between the two pieces (0.2–1.5 nm)
accounts for most of it, without needing that.

**T/R alone does not determine Ag n,k.** The T/R-only column fits its own data
as well as the joint one and still returns k(550) of 1.5 at Ag 7–8 nm against
3.4–3.6 from either measurement done properly, with SE MSE two orders of
magnitude out. Two observables per wavelength over 420–780 nm do not constrain
a 16-parameter model. That column is a demonstration, not a result — and it is
the reason the campaign's own inversion needs the ellipsometry beside it.

**MoOx is not one homogeneous layer.** The joint fit costs a factor ~2 in SE MSE
and still leaves 0.4–1.2 %p in T/R, and the two pieces come out with very
different thicknesses (Ag 4 nm: 7.0 nm on Si vs 3.0 nm on glass). For island
films that is expected — morphology depends on what is under the seed — but it
means the MoOx numbers here are effective-medium values for one particular
piece, not material constants. Quote them with the sample.

For the seed itself the agreement with the lab library is good: HATCN
n(550) = 1.877 here against 1.849 in `data/nk/l_HATCN.csv`, both with k ≈ 0,
with this fit's dispersion ~2 % flatter toward the red.

## 7. Deliverables

`jnk_report.py` writes to `$ELLIPS_OUT/jnk_nk/`:

* `seed_HATCN_nk.csv`, `seed_MoOx_nk.csv` — 300–1000 nm, 5 nm grid
* `Ag_<seed>_<sheet>_<nominal>nm_nk.csv` — n,k from all three weightings side
  by side, so the spread between them is visible in the delivered file. **Use
  the `joint` columns**; `SE` and `TR` are there to show the spread.
* `nk_by_seed.png` — n and k against bulk Ag (McPeak), both seeds
* `n550_se_vs_tr.png` — the three weightings against nominal thickness
* `tr_joint_fit.png` — measured vs joint-model T and R
* `seed_nk.png`

## 8. What is still open

* `1-1` (HATCN, no Ag) was never measured on glass, `1-9` has no T and `2-12`
  no R, so those three have no full T/R pair. The consolidated
  `ALL_SAMPLES_TRA.csv` blanks `MoOx5_bare` and `MoOx5_Ag12` entirely; this
  pipeline reads the raw exports, so `2-12` keeps its transmittance.
* Above 1080 nm the SE data is unused (Si turns transparent at 1107 nm and the
  wafer backside enters); above 780 nm the T/R data is unused (instrument
  changeover). Nothing here constrains the NIR beyond those limits.
* The MoOx series wants an inhomogeneous model — a graded layer or an explicit
  island description — not a thicker single layer.
* If a CompleteEASE `.mod` with `SI_JAW`/`NTVE_JAW` turns up, `lib/ce_mat.py`
  decodes it and the whole pipeline can be re-run in that frame for a direct
  comparison with the 260819 numbers; §1 suggests the difference will be small.

## 9. How far do these constants carry into a TMM / CPS device model?

`jnk_devsens.py` propagates the remaining ambiguity into the quantities a
device calculation asks for, at a fixed thickness so only the dispersion
varies. Read the **SE-vs-joint** gap as the bracket: both fit the ellipsometry
to MSE 2–3, only the joint one also reproduces the absolute T/R. (The T/R-only
column is there to show that intensity alone does not determine n,k.)

| Ag on HATCN, 550 nm | SE-only | joint | gap |
|---|---|---|---|
| ε₁ = n²−k² (8 nm) | −13.50 | −11.88 | 14 % |
| ε₂ = 2nk (8 nm) | 0.639 | 1.056 | 39 % |
| electrode A at 0° in a device stack (8 nm) | 2.19 % | 3.63 % | 40 % |
| same at 60° internal | 2.50 % | 4.21 % | 41 % |
| quenching FOM Im[(ε_m−ε_d)/(ε_m+ε_d)] (8 nm) | 0.039 | 0.090 | 57 % |
| short-range SPP length L (8 nm) | 152 nm | 69 nm | 2.2× |

Across the series ε₁ holds to 2–14 % but ε₂ moves 25–86 %, and the device-level
losses move with ε₂. The asymmetry is structural: for a metal with n ≪ k,
ε₁ ≈ −k² is carried by the well-determined half of the dispersion and
ε₂ = 2nk is linear in n, the half that only the absolute intensity pins down.
**Where a resonance sits is safe; how lossy it is inherits the whole n
uncertainty.**

What is genuinely validated: the joint fit reproduces T and R at normal
incidence on glass to 0.16–0.24 %p, so A = 1−T−R over 420–780 nm is good to
~0.3 %p absolute — a few percent relative on a film absorbing 6–15 %. Inside
that box a TMM calculation is on solid ground.

What is extrapolation, in order of how much it should worry you:

1. **The evanescent region is not measured at all.** SE at 45–65° in air
   reaches u = sin θ ≤ 0.91 and the T/R is at normal incidence; the CPS
   integral is dominated by u > 1. Every SPP and lossy-surface-wave number
   above comes from the oscillator model continued past its data.
2. **The layer was fitted isotropic.** A columnar or island film is uniaxial,
   and p-polarised large-u response depends on ε_z, which neither data set
   constrains. `ellipsometry_fit._tmm_uniaxial` and `thin_ag_260819/ag_anis.py`
   are the starting point if this matters.
3. **The adjacent media differ.** These films were measured with air above and
   Si or glass below; in a device the Ag sits between the seed and an organic
   capping. For the near-bulk HATCN films (≥ 7 nm) that is a material constant
   and travels; for the 4–6 nm films and for all of MoOx it is an
   effective-medium value that does not.
4. **Roughness is a coherent EMA layer, not scattering.** The fits put it at
   0–1.3 nm, implausibly smooth for a 4–6 nm Ag film, so some real roughness is
   being absorbed into n,k. A device model will book that as absorption.
5. **MoOx should not be used as material constants at all** (§6).

So: for electrode transmission, reflection and absorbed fraction in the visible
at moderate internal angles, on HATCN, at 7 nm and above — yes, and the joint
fit is markedly better for this than an SE-only fit, because an absolute
intensity measurement is exactly the constraint such a calculation needs. For
a CPS mode split, quote the SPP and quenching channels with a factor ~2, or
measure the evanescent region first.

The measurement that would close it is one that reaches u > 1 — Kretschmann /
ATR or prism-coupled ellipsometry on the same films. Angle- and
polarisation-resolved absolute T/R would be cheaper and would at least extend
the absolute anchor from u = 0 to u ≈ 0.9; a haze or integrating-sphere
measurement would separate the scattering that item 4 currently hides in k.
