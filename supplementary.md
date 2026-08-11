# Supplementary Information

Practical Saturation of Freeform Microlens Arrays on Extended OLED Emitters
and Design Routes beyond Lens Shape

All simulations share one source model, one substrate, and one textured patch
size. Values that differ between runs are listed explicitly below; anything not
listed is identical across all campaigns.

---

## Table S1 | Common platform

| Quantity | Value |
|---|---|
| Source | disc, radius $r_{\mathrm{OLED}} = 1$ mm, CPS microcavity dipole distribution $I_{\mathrm{sub}}(\theta,\lambda)$ |
| Substrate | thickness $d_{\mathrm{sub}} = 1.295$ mm, $n = 1.51$ |
| Textured patch | 25 × 25 mm (all campaigns) |
| Lenslet | ~10 μm radius, hexagonal placement (X 0.0866 mm, Y 0.1 mm) |
| Emission window | 453–753 nm |
| Emitter | $\eta_{\mathrm{rad}} = 0.98$, horizontal dipole ratio 0.865 |
| Stack | Al / ETL / EML / HTL / ITO / glass, ETL and HTL thicknesses free in 10–150 nm |
| Shape parameterization | 7 spline control points, endpoints fixed at (0,1) and (1,0); five free points give 10 shape variables with $x_2 \ldots x_6$ constrained monotonic |
| Design vector | 13 variables = 10 shape + $d_{\mathrm{ETL}}$ + $d_{\mathrm{HTL}}$ + stretch$_Z$ |

Angular bands are polar bins of the far field over the full azimuth:
0–20°, 20–40°, 40–60°, 60–80°. Band selectivity is
$S_j = \mathrm{EQE}_{\mathrm{band},j}/\mathrm{EQE}_{\mathrm{total}}$.
The substrate-mode-free Lambertian partition,
$S_j^{\mathrm{Lam}} = \sin^2\theta_{\mathrm{hi}} - \sin^2\theta_{\mathrm{lo}}$,
gives 0.117 / 0.296 / 0.337 / 0.220.

---

## Table S2 | Per-campaign settings

| Campaign | Script | Free variables | Search budget | Search fidelity | Final fidelity | Log marker |
|---|---|---|---|---|---|---|
| Convex, weighted sweep | `pareto_front_freeform.m` | 13 | 150 random + 120/weight + 15 polish, $w \in \{0,\,0.25,\,0.5,\,0.75,\,1\}$ | 10,000 rays, 31 λ | 50,000 rays, 151 λ, ×3 | phase 1–3 |
| Convex, per-band | `opt_4band_freeform.m` | 13 | 60 + 15 polish per arm, 5 arms | 10,000 rays, 31 λ | 50,000 rays, 151 λ, ×3 | phase 4 |
| Hemispherical reference | `opt_hemisphere_arms.m` | 3 (cavity + height; shape fixed on a quarter circle) | 30 + 10 polish per arm, 5 arms | 10,000 rays, 31 λ | 50,000 rays, 151 λ, ×3 | phase 8 |
| Inverted (concave) | `opt_4band_inverted.m` | 13 | 60 + 15 polish per arm, 4 arms | 10,000 rays, 31 λ | 50,000 rays, 151 λ, ×3 | phase 6 |
| Randomly assembled | `stress_random_mla.m` | 6 assembly statistics | 50 random + 60 + 15 polish | 5,000 rays, 16 λ | 10,000 rays, 151 λ, ×3 seeds | phase 7 |
| Convergence check | `convergence_check.m` | — (re-evaluation only) | 20 designs re-evaluated | — | 200,000 rays ×3, and 450–750 nm broadband | — |

The randomly assembled family is the only campaign with reduced sampling. Its
supercell carries 6 × 6 lenslets on a 141 × 141 height grid, so one evaluation
traces far more geometry than a single-lens texture; Table S3 records what that
reduction costs in accuracy. The inverted family uses a model file identical to
the convex one except for the texture relief setting
(`LibraryElementUnitCell.Bumps`, "Yes" → "No"), so the two differ in
construction and in nothing else.

---

## Table S3 | Cost calibration for the randomly assembled family

One fixed geometry (seed 7777, mid-range assembly statistics) re-evaluated under
each setting. The acceptance criterion is band selectivity, since the reported
quantities are ratios; total EQE may shift by a few tenths of a percent without
consequence.

| Setting | lenslets | grid | rays | λ | time (s) | speed-up | EQE | Δ selectivity (pp) |
|---|---|---|---|---|---|---|---|---|
| reference | 8 × 8 | 201 | 10,000 | 31 | 789 | 1.00 | 0.4975 | — |
| half rays | 8 × 8 | 201 | 5,000 | 31 | 544 | 1.45 | 0.4992 | 0.38 |
| coarser λ | 8 × 8 | 201 | 10,000 | 16 | 453 | 1.74 | 0.5012 | 0.66 |
| coarser grid | 8 × 8 | 141 | 10,000 | 31 | 725 | 1.09 | 0.5016 | 0.28 |
| fewer lenslets | 6 × 6 | 201 | 10,000 | 31 | 559 | 1.41 | 0.4948 | 0.54 |
| **adopted** | **6 × 6** | **141** | **5,000** | **16** | **191** | **4.13** | **0.4950** | **0.39** |

The adopted combination is 4.13× faster while shifting selectivity by 0.39
percentage points — the same size as the Monte-Carlo spread measured by
repeating one design (0.5% relative on total EQE). Two individual reductions
exceed 0.5 pp on their own; that they do not accumulate confirms the deviations
are noise rather than bias.

---

## Table S4 | Convergence of the selectivity–efficiency drift

The correlation $R(\mathrm{EQE}_{\mathrm{total}}, S_j)$ was re-measured at
twenty-fold ray count and over the full emission band to establish that it is
neither Monte-Carlo noise nor a narrowband artifact.

| Band | baseline | 200,000 rays, ×3 repeats | broadband 450–750 nm |
|---|---|---|---|
| 0–20° | +0.59 | +0.61 | +0.64 |
| 20–40° | +0.67 | +0.67 | +0.72 |
| 40–60° | +0.07 | +0.04 | +0.05 |
| 60–80° | −0.70 | −0.71 | −0.76 |

---

## Table S5 | Family comparison

Natural composition is the median selectivity of the twenty designs with the
highest total EQE in each campaign.

| Family | best total EQE | natural composition (0–20 / 20–40 / 40–60 / 60–80) |
|---|---|---|
| Hemispherical reference | 0.5468 | — |
| Convex freeform | 0.5481 | 0.096 / 0.278 / 0.360 / 0.237 |
| Randomly assembled | 0.5216 ± 0.0011 | 0.095 / 0.281 / 0.360 / 0.233 |
| Inverted (concave) | 0.5165 | 0.113 / 0.303 / 0.340 / 0.215 |
| Lambertian partition | — | 0.117 / 0.296 / 0.337 / 0.220 |

The uncertainty quoted for the randomly assembled family is the spread over
three independent disorder realisations of the same assembly statistics
(coefficient of variation 0.2%).

Relative gain of the freeform optimum over the equally optimized hemisphere,
$G_j = \max[\mathrm{EQE}_j \mid \mathrm{freeform}] / \max[\mathrm{EQE}_j \mid \mathrm{hemisphere}]$:

| Objective | freeform | hemisphere | $G_j$ |
|---|---|---|---|
| 0–20° | 0.06318 | 0.06682 | 0.945 |
| 20–40° | 0.16032 | 0.16591 | 0.966 |
| 40–60° | 0.19802 | 0.18503 | 1.070 |
| 60–80° | 0.13863 | 0.13589 | 1.020 |
| total EQE | 0.5481 | 0.54679 | 1.002 |

---

## Reproducing the figures

| Figure | Script | Inputs |
|---|---|---|
| Fig. 2a | `plot_fig2a_selectivity.m` | `opt_4band_result.mat` |
| Fig. 3 | `pareto_front_freeform.m` (final block) | its own log |
| Fig. 5 | `plot_fig5_families.m` | `opt_4band_result.mat`, `opt_4band_inverted_result.mat`, `stress_random_result.mat` |
| Table S3 | `calibrate_random_cost.m` | — |
| Table S4 | `convergence_check.m` | `pareto_front_result.mat` |

Result archives: `pareto_front_result.mat`, `opt_4band_result.mat`,
`opt_hemisphere_result.mat`, `opt_4band_inverted_result.mat`,
`stress_random_result.mat`, `calibrate_random_cost.mat`,
`convergence_check_result.mat`.
