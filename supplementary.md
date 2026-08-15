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
| Convex, total EQE only | `freeform_EQEtotal.m` | 13 | 140 per start, 3 independent starts | 10,000 rays, 31 λ | 50,000 rays, 151 λ, ×3 | — |
| Hemisphere-seeded control | `warmstart_from_hemisphere.m` | 13 | 40 + 15 + 15 polish per arm, 4 arms | 10,000 rays, 31 λ | 50,000 rays, 151 λ, ×3 | phase 9 |
| Hemispherical reference | `opt_hemisphere_arms.m` | 3 (cavity + height; shape fixed on a quarter circle) | 30 + 10 polish per arm, 5 arms | 10,000 rays, 31 λ | 50,000 rays, 151 λ, ×3 | phase 8 |
| Inverted (concave) | `opt_4band_inverted.m` | 13 | 60 + 15 polish per arm, 4 arms (identical protocol to the convex benchmark) | 10,000 rays, 31 λ | 50,000 rays, 151 λ, ×3 | phase 6 |
| Randomly assembled | `stress_random_mla.m` | 6 assembly statistics | 50 random + 60 + 15 polish | 5,000 rays, 16 λ | 10,000 rays, 151 λ, ×3 seeds | phase 7 |
| Convergence check | `convergence_check.m` | — (re-evaluation only) | 20 designs re-evaluated | — | 200,000 rays ×3, and 450–750 nm broadband | — |
| Patch-size check | `check_patch_convergence.m` | — (re-evaluation only) | 1 design × 3 patch sizes × 3 repeats | — | 50,000 rays, 151 λ | phase 11 |
| 20–40° confirmation | `reeval_confirm_2040.m` | — (re-evaluation only) | 2 designs × 5 repeats | — | 50,000 rays, 151 λ | phase 10 |

The randomly assembled family is the only campaign with reduced sampling. Its
supercell carries 6 × 6 lenslets on a 141 × 141 height grid, so one evaluation
traces far more geometry than a single-lens texture; Table S3 records what that
reduction costs in accuracy. The inverted family uses a model file identical to
the convex one except for the texture relief setting
(`LibraryElementUnitCell.Bumps`, "Yes" → "No"), so the two differ in
construction and in nothing else.

---

## Table S3 | Cost calibration for the randomly assembled family (Fig. S3)

One fixed geometry (seed 7777, mid-range assembly statistics) re-evaluated under
each setting. The acceptance criterion is band selectivity, since the reported
quantities are ratios; total EQE may shift by a few tenths of a percent without
consequence.

Δ selectivity is the largest absolute deviation from the reference across the
four bands, in percentage points.

| Setting | lenslets | grid | rays | λ step | λ | time (s) | speed-up | EQE | Δ selectivity (pp) |
|---|---|---|---|---|---|---|---|---|---|
| reference | 8 × 8 | 201 | 10,000 | 10 nm | 31 | 777.5 | 1.00 | 0.4975 | — |
| half rays | 8 × 8 | 201 | 5,000 | 10 nm | 31 | 529.3 | 1.47 | 0.4992 | 0.38 |
| coarser λ | 8 × 8 | 201 | 10,000 | 20 nm | 16 | 446.5 | 1.74 | 0.5012 | 0.66 |
| coarser grid | 8 × 8 | 141 | 10,000 | 10 nm | 31 | 704.6 | 1.10 | 0.5016 | 0.28 |
| fewer lenslets | 6 × 6 | 201 | 10,000 | 10 nm | 31 | 543.3 | 1.43 | 0.4947 | 0.54 |
| **adopted** | **6 × 6** | **141** | **5,000** | **20 nm** | **16** | **181.0** | **4.30** | **0.4950** | **0.39** |

The adopted combination is 4.3× faster while shifting selectivity by 0.39
percentage points — the same size as the Monte-Carlo spread measured by
repeating one design (0.5% relative on total EQE). Two individual reductions
exceed 0.5 pp on their own; that they do not accumulate confirms the deviations
are noise rather than bias.

Because the seed is fixed (7777), the accuracy columns are reproducible exactly;
the wall-clock times are not, and vary by a few percent with machine load.

---

## Table S4 | Convergence of the selectivity–efficiency drift (Fig. S4)

Twenty designs stratified across the efficiency range were re-evaluated at
twenty-fold ray count and over the full emission band, to establish that the
drift is neither Monte-Carlo noise nor a narrowband artifact. Because the subset
is stratified rather than representative, its correlations differ in magnitude
from the full-population values quoted in Section 2.3
($+0.60 / +0.56 / -0.12 / -0.57$ over all 606 usable of the 691 logged evaluations); the sign pattern and
band ordering agree.

| Band | baseline | 200,000 rays, ×3 repeats | broadband 450–750 nm |
|---|---|---|---|
| 0–20° | +0.59 | +0.61 | +0.64 |
| 20–40° | +0.67 | +0.67 | +0.72 |
| 40–60° | +0.07 | +0.04 | +0.05 |
| 60–80° | −0.70 | −0.71 | −0.76 |

---

## Table S5 | Family comparison

Two statistics are given for each family, both computed from that family's own
evaluation log: the **top-20 median**, i.e. the median selectivity of the twenty
designs with the highest total EQE, and the **population mean** over every
usable evaluation of the campaign. Reporting both matters, because the gap
between them within one family is comparable to the gap between families under
either one — a single-statistic table would create apparent family differences
that the data do not support.

| Family | best total EQE | top-20 median (0–20 / 20–40 / 40–60 / 60–80) | population mean |
|---|---|---|---|
| Hemispherical reference | 0.5468 | — | — |
| Convex freeform (per-band campaign) | 0.5481 | 0.094 / 0.278 / 0.361 / 0.236 | 0.097 / 0.278 / 0.350 / 0.242 |
| Convex freeform (weighted sweep) † | 0.5556 | 0.102 / 0.276 / 0.354 / 0.239 | 0.095 / 0.274 / 0.354 / 0.244 |
| Convex freeform (total EQE only) | 0.5539 | — | — |
| Randomly assembled | 0.5216 ± 0.0011 | 0.111 / 0.299 / 0.347 / 0.213 | 0.102 / 0.289 / 0.355 / 0.224 |
| Inverted (concave) | 0.5165 | 0.113 / 0.303 / 0.340 / 0.215 | 0.100 / 0.288 / 0.357 / 0.231 |
| Lambertian partition | — | 0.117 / 0.296 / 0.337 / 0.220 | — |

† Large-patch normalization: this campaign's absolute totals correspond to the
infinite-patch limit (Table S7); its best design re-measures at 0.5428 on the
fixed 25 mm patch. Its compositions, being ratios, are comparable with the
other rows.

Every entry lies within 0.094–0.113 (0–20°), 0.274–0.303 (20–40°),
0.340–0.361 (40–60°) and 0.213–0.244 (60–80°). For the randomly assembled
family the mean over its unbiased random-sampling phase alone is
0.095 / 0.281 / 0.360 / 0.233, and the three high-precision winning realisations
give 0.110 / 0.298 / 0.349 / 0.213; both fall inside the same window.

The uncertainty quoted for the randomly assembled family is the spread over
three independent disorder realisations of the same assembly statistics
(coefficient of variation 0.2%).

Relative gain of the freeform optimum over the equally optimized hemisphere,
$G_j = \max[\mathrm{EQE}_j \mid \mathrm{freeform}] / \max[\mathrm{EQE}_j \mid \mathrm{hemisphere}]$:

| Objective | freeform | hemisphere | $G_j$ |
|---|---|---|---|
| 0–20° | 0.06699 | 0.06682 | 1.003 |
| 20–40° | 0.16628 | 0.16591 | 1.002 |
| 40–60° | 0.19802 | 0.18503 | 1.070 |
| 60–80° | 0.14075 | 0.13589 | 1.036 |
| total EQE | 0.5539 | 0.54679 | 1.013 |

Every entry in this table is from a fixed-patch (25 × 25 mm) campaign, so the
two columns share one normalization. The 0–20°, 20–40° and 60–80° freeform entries
come from the hemisphere-seeded control (Table S6) — in the 60–80° arm the
control surpassed the per-band campaign's own optimum (0.14075 against
0.13863); the 40–60° entry is from the per-band campaign, which the control did
not surpass there (0.19355 against 0.19802); and the total is from the
dedicated total-EQE campaign (0.5495 / 0.5530 / 0.5539 over three independent
starts), which the control's same-session value of 0.5523 confirms to within
0.3%.

The weighted-sweep campaign's 0.5556 does not appear in this table because its
normalization is the large-patch limit (Table S7), not the fixed patch; its
coarse-fidelity log maximum of 0.5591 at 10,000 rays is likewise not used
anywhere.

---

## Table S6 | Hemisphere-seeded control (Fig. S2)

Each arm's search restarted at that arm's hemispherical optimum, with the
hemisphere point in the seed set, eight perturbations within 8% of each
variable's range, and a pattern search launched directly from the hemisphere
point. The hemisphere baseline is re-measured in the same session at final
fidelity; the archived value is a cross-check only.

| Arm | hemisphere (archived) | hemisphere (re-measured) | dev. | warm start | gain | $t$ | winning branch |
|---|---|---|---|---|---|---|---|
| 0–20° | 0.06682 | 0.06677 ± 0.00017 | 0.07% | 0.06699 ± 0.00013 | +0.32% | 1.7 | polish from hemisphere |
| 20–40° | 0.16591 | 0.16595 ± 0.00017 | 0.02% | 0.16628 ± 0.00020 | +0.20% ‡ | 2.1 ‡ | surrogateopt |
| 40–60° | 0.18503 | 0.18486 ± 0.00008 | 0.09% | 0.19355 ± 0.00010 | **+4.70%** | **114.6** | polish from surrogate |
| 60–80° | 0.13589 | 0.13581 ± 0.00014 | 0.06% | 0.14075 ± 0.00025 | **+3.64%** | **29.6** | surrogateopt |
| total EQE | 0.54679 | 0.54679 ± 0.00013 | 0.00% | 0.55231 ± 0.00004 | **+1.01%** | **68.7** | polish from surrogate |

Uncertainties are the standard deviation over three high-precision repeats. The
threshold is the pooled one-sided 95% $t$ value at four degrees of freedom,
2.13.

‡ The 20–40° arm was the one borderline case, and was settled by a dedicated
confirmation run: both stored designs re-measured with five fresh repeats each,
search excluded. This gives hemisphere 0.16583 ± 0.00015 against warm start
0.16631 ± 0.00021, a residue of +0.29% at $t = 4.1$ (five-repeat threshold
1.86) — statistically real, three tenths of a percent in size
(`reeval_confirm_2040_result.mat`).

Two properties of this control matter for how it should be read. It is
one-sided: the hemisphere is among the screened candidates, so the result cannot
fall meaningfully below it. And it is biased toward reporting a gain, since the
winner is the maximum over roughly seventy noisy search evaluations per arm. The
near-zero outcomes are therefore conservative, and the separation in scale
between the two groups—residues of +0.29% to +0.32% against margins of +1.01%
to +4.70% at $t = 30$ to $115$—is the evidence that the residues reflect an
essentially absent margin rather than a weak search. In the 60–80° arm the
control exceeded the original per-band campaign's own optimum (0.14075 against
0.13863), so restarting at the hemisphere is also simply a better search
protocol in this design space.

The first four arms took 7.6 h (108–118 min per arm) with 1 of 192 evaluations
rejected on geometry; the 60–80° arm and the confirmation run completed in a
second overnight session.

---

## Table S7 | Patch-size dependence (Fig. S1)

One fixed high-efficiency design (the weighted-sweep optimum at $w = 0.75$),
re-evaluated at final fidelity (50,000 rays, 151 λ) with only the textured
patch size varied. The 15/25/35 mm rows are 3 repeats
(`patch_convergence_result.mat`); the 100 mm row is 2 repeats from a follow-up
run (`patch_convergence_100.mat`).

| Patch (mm) | total EQE | s.d. | vs 25 mm | selectivity (0–20 / 20–40 / 40–60 / 60–80) |
|---|---|---|---|---|
| 15 × 15 | 0.51681 | 0.00016 | −4.79% | 0.108 / 0.283 / 0.346 / 0.233 |
| 25 × 25 | 0.54282 | 0.00008 | — | 0.106 / 0.280 / 0.346 / 0.238 |
| 35 × 35 | 0.55128 | 0.00003 | +1.56% | 0.106 / 0.279 / 0.346 / 0.239 |
| 100 × 100 | 0.5636 | 0.00003 | +3.83% | 0.106 / 0.277 / 0.347 / 0.240 |

Three observations follow.

1. The initial rise saturates with a decay length of about 9 mm — four
   critical-angle round trips of 2.32 mm each — but the far tail does not:
   between 35 and 100 mm the total still gains 2.2%, about three times what a
   single-exponential tail fitted to the first three points would allow.
   (An earlier draft of this table extrapolated an infinite-patch limit of
   0.5554 from the first three points alone; the 100 mm measurement falsifies
   that extrapolation, and it is withdrawn.) Multiply recycled light migrates
   laterally over many tens of round trips before escaping, so the total EQE
   of a finite film increases slowly with its extent, and the 25 mm value is a
   lower bound sitting 3.7% below the 100 mm value.
2. The angular composition is patch-converged: while the total gains 3.8%
   from 25 to 100 mm, no band selectivity changes by more than 0.24 pp over
   the same fourfold enlargement (0.11 pp between 25 and 35 mm). All ratios,
   correlations and class comparisons, which are evaluated at equal patch,
   are unaffected by the patch choice.
3. The weighted-sweep campaign's archived 0.5556 falls between the 35 mm and
   100 mm values of the same design (0.5513 and 0.5636), while the design
   re-measures at 0.5428 on the fixed 25 mm patch. That campaign's absolute
   normalization therefore corresponds to a larger effective patch than the
   unified 25 mm one, which is why its totals are excluded from the
   fixed-patch $G_j$ comparison of Table S5.

---

## Reproducing the figures

All five manuscript figures are produced by `make_figures.py`; the
supplementary figures S1–S4 and the raw-data workbooks by
`make_supp_figures_and_data.py`, run from the repository root in that order:

```
python3 make_figures.py
python3 make_supp_figures_and_data.py
```

The first script reads only the archived result files, prints every number
quoted in the captions to `figure_numbers.txt`, and writes each figure as both
PNG (300 dpi) and PDF. The second exports, for **every** figure (main and
supplementary), an Excel workbook named after the figure — e.g.
`fig2_achievable_region.xlsx` alongside `fig2_achievable_region.png` — with one
sheet per panel containing exactly the plotted arrays, so any figure can be
re-plotted or restyled from its workbook without touching the `.mat` archives.
Because the export script imports the plotting script, the workbooks cannot
drift from the figures.

| Figure | Output | Inputs |
|---|---|---|
| Fig. 1 | `fig1_platform.png` / `.pdf` | `opt_4band_result_25by25.mat`, `opt_hemisphere_result.mat` |
| Fig. 2 | `fig2_achievable_region.png` / `.pdf` | `pareto_front_result.mat`, `opt_4band_result_25by25.mat`, `opt_hemisphere_result.mat`, `warmstart_hemisphere_result.mat`, `freeform_EQEtotal_result.mat` |
| Fig. 3 | `fig3_selectivity_map.png` / `.pdf` | `pareto_front_result.mat` |
| Fig. 4 | `fig4_recycling_routes.png` / `.pdf` | `angular_recycling_result.npz`, `angular_recycling_bandwidth.npz` |
| Fig. 5 | `fig5_families.png` / `.pdf` | `opt_4band_result_25by25.mat`, `opt_4band_inverted_result.mat`, `stress_random_result.mat` |
| Fig. S1 | `figS1_patch_dependence.png` / `.pdf` | `patch_convergence_result.mat`, `patch_convergence_100.mat` |
| Fig. S2 | `figS2_warmstart_control.png` / `.pdf` | `warmstart_hemisphere_result.mat`, `reeval_confirm_2040_result.mat` |
| Fig. S3 | `figS3_cost_calibration.png` / `.pdf` | `calibrate_random_cost.mat` |
| Fig. S4 | `figS4_convergence.png` / `.pdf` | `convergence_check_result.mat` |
| Table S3 | — | `calibrate_random_cost.m` |
| Table S4 | — | `convergence_check.m`, `pareto_front_result.mat` |

The MATLAB plotting scripts `plot_fig2a_selectivity.m` and
`plot_fig5_families.m` remain in the repository as the in-session views used
during the campaigns; the published figures come from `make_figures.py`.

Result archives: `pareto_front_result.mat`, `opt_4band_result_25by25.mat`,
`freeform_EQEtotal_result.mat`, `opt_hemisphere_result.mat`,
`opt_4band_inverted_result.mat`,
`stress_random_result.mat`, `calibrate_random_cost.mat`,
`warmstart_hemisphere_result.mat`, `reeval_confirm_2040_result.mat`,
`patch_convergence_result.mat`, `patch_convergence_100.mat`,
`convergence_check_result.mat`,
`angular_recycling_result.npz`,
`angular_recycling_bandwidth.npz`.
