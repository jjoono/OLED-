# Supplement 1

An Angular Mixing Number for Passive Light-Extraction Structures:
Why Microlens Shape Saturates on Extended OLED Emitters

This supplement collects the material that supports the main text but is not needed to follow its argument. Notes S1–S9 document the freeform-versus-hemisphere benchmark (platform and campaign settings, statistical controls, normalization checks, convergence and cost calibrations, the inverted and randomly assembled families, the separately constrained asymmetric study, and the pass-by-pass analysis of angle retention). Notes S10 and S11 document the angular mixing number and the library of structures and return matrices. Main-text references are cited as "main ref. n".

---

## Note S1 | Common platform and campaign settings

All simulations share one source model, one substrate, and one textured patch size. Table S1 lists the common platform. Table S2 lists every setting that differs between campaigns. Anything not listed is identical across campaigns.

**Table S1 | Common platform.**

| Quantity | Value |
|---|---|
| Source | disc, radius 1 mm, CPS microcavity dipole distribution $I_{\mathrm{sub}}(\theta,\lambda)$ |
| Substrate | thickness $d_{\mathrm{sub}} = 1.295$ mm, $n = 1.51$ |
| Textured patch | 25 × 25 mm (all campaigns; see Note S3) |
| Lenslet | ~10 μm radius, hexagonal placement (X 0.0866 mm, Y 0.1 mm) |
| Lens material | N-BK7, $n = 1.517$ at 589 nm, index-matched to the substrate, identical for every lens class and family |
| Emission window | 453–753 nm |
| Emitter | $\eta_{\mathrm{rad}} = 0.98$, horizontal dipole ratio 0.865 |
| Stack | Al / ETL / EML / HTL / ITO / glass, ETL and HTL thicknesses free in 10–150 nm |
| Shape parameterization | 7 spline control points, endpoints fixed at (0,1) and (1,0); five free points give 10 shape variables with $x_2 \ldots x_6$ constrained monotonic |
| Design vector | 13 variables = 10 shape + $d_{\mathrm{ETL}}$ + $d_{\mathrm{HTL}}$ + stretch$_Z$ |

Angular bands are polar bins of the far field over the full azimuth: 0–20°, 20–40°, 40–60°, and 60–80°. The band selectivity is $S_j = \mathrm{EQE}_{j}/\mathrm{EQE}_{\mathrm{total}}$. The Lambertian partition, $S_j^{\mathrm{Lam}} = \sin^2\theta_{\mathrm{hi}} - \sin^2\theta_{\mathrm{lo}}$, gives 0.117 / 0.296 / 0.337 / 0.220.

**Table S2 | Per-campaign settings.**

| Campaign | Script | Free variables | Search budget | Search fidelity | Final fidelity |
|---|---|---|---|---|---|
| Convex, weighted sweep | `pareto_front_freeform.m` | 13 | 150 random + 120/weight + 15 polish, $w \in \{0,\,0.25,\,0.5,\,0.75,\,1\}$ | 10,000 rays, 31 λ | 50,000 rays, 151 λ, ×3 |
| Convex, per-band | `opt_4band_freeform.m` | 13 | 60 + 15 polish per arm, 5 arms | 10,000 rays, 31 λ | 50,000 rays, 151 λ, ×3 |
| Convex, total EQE only | `freeform_EQEtotal.m` | 13 | 140 per start, 3 independent starts | 10,000 rays, 31 λ | 50,000 rays, 151 λ, ×3 |
| Restart control | `warmstart_from_hemisphere.m` | 13 | 40 + 15 + 15 polish per arm, 5 arms | 10,000 rays, 31 λ | 50,000 rays, 151 λ, ×3 |
| Hemispherical reference | `opt_hemisphere_arms.m` | 3 (cavity + height) | 30 + 10 polish per arm, 5 arms | 10,000 rays, 31 λ | 50,000 rays, 151 λ, ×3 |
| Inverted (concave) | `opt_4band_inverted.m` | 13 | 60 + 15 polish per arm, 4 arms | 10,000 rays, 31 λ | 50,000 rays, 151 λ, ×3 |
| Randomly assembled | `stress_random_mla.m` | 6 assembly statistics | 50 random + 60 + 15 polish | 5,000 rays, 16 λ | 10,000 rays, 151 λ, ×3 seeds |
| Convergence check | `convergence_check.m` | re-evaluation only | 20 designs | — | 200,000 rays ×3; 450–750 nm broadband |
| Patch-size check | `check_patch_convergence.m` | re-evaluation only | 1 design × 4 patch sizes | — | 50,000 rays, 151 λ |
| 20–40° confirmation | `reeval_confirm_2040.m` | re-evaluation only | 2 designs × 5 repeats | — | 50,000 rays, 151 λ |

The inverted family uses a model file identical to the convex one except for the texture relief setting (`LibraryElementUnitCell.Bumps`, "Yes" → "No"). The two families therefore differ in construction and nothing else. The randomly assembled family is the only campaign with reduced sampling (Note S6).

**Usable evaluations and trace failures.** An evaluation counts as usable when the trace returns a positive total EQE; no outliers are removed. The weighted-sum campaign logged 691 evaluations (150 random, 526 optimizer, 15 polish), of which 606 are usable (114 random, 492 optimizer and polish). The convex per-band campaign logged 238 evaluations with 25 trace failures (10.5%); the hemispherical campaign logged 224 with none. Candidates rejected by the geometric feasibility rules (main-text Table 1) are not traced; they are returned to `surrogateopt` as infeasible through its inequality-constraint output, so the optimizer learns the feasible region rather than spending evaluations on it.

---

**Fig. S1 | Platform and benchmark design.** (a) OLED–substrate–MLA architecture: CPS microcavity dipole source (emitting radius 1 mm), glass substrate ($d_{\mathrm{sub}} = 1.295$ mm), hexagonal lenslet array (~10 μm lens radius, 25 × 25 mm). (b) The two lens classes under identical constraints, in units of the lens radius: the hemispherical reference (3 free variables) and the axisymmetric freeform (13 variables). Open circles are the five free spline control points; the freeform shown is the 40–60° optimum. (c) The four polar bands and the band selectivity $S_j$, taken over the full azimuth. (d) Optimization workflow: CPS source, LightTools ray tracing through the MATLAB COM interface, surrogate-based global search, pattern-search polishing, and high-ray re-evaluation. *(`fig1_platform.png`)*

---

## Note S2 | Hemispherical benchmark and restart control

**Table S3** lists the entries behind $G_j$ in main-text Fig. 5b. Every entry comes from a fixed-patch (25 × 25 mm) campaign, so the two columns share one normalization.

**Table S3 | Relative gain of the freeform optimum over the equally optimized hemisphere.**

| Objective | freeform | hemisphere | $G_j$ | source of freeform entry |
|---|---|---|---|---|
| 0–20° | 0.06699 | 0.06682 | 1.003 | restart control |
| 20–40° | 0.16628 | 0.16591 | 1.002 | restart control |
| 40–60° | 0.19802 | 0.18503 | 1.070 | per-band campaign |
| 60–80° | 0.14075 | 0.13589 | 1.036 | restart control |
| total EQE | 0.5539 | 0.54679 | 1.013 | dedicated total-EQE campaign |

Each freeform entry is the better of the matched freeform campaigns. In the 40–60° arm the restart control did not surpass the per-band optimum (0.19355 against 0.19802). In the 60–80° arm it did (0.14075 against 0.13863). The total is from the dedicated campaign (0.5495 / 0.5530 / 0.5539 over three starts), and the restart control's same-session value of 0.5523 confirms it to within 0.3%.

**Test statistic.** For each objective the two designs are re-evaluated $n = 3$ times. The statistic is $t = \Delta/(s_p\sqrt{2/n})$, where $\Delta$ is the absolute difference of the means and $s_p$ the pooled standard deviation; the one-sided 95% threshold at four degrees of freedom is 2.13. The five-repeat re-measurement of the 20–40° arm uses eight degrees of freedom (threshold 1.86; 2.90 at 1%). With a Holm correction over the five objectives, the step-wise one-sided thresholds are 3.75, 3.50, 3.19, 2.78, and 2.13; the 40–60°, 60–80°, and total-EQE gains ($t$ = 115, 30, 69) remain significant. The large $t$ values reflect the small re-evaluation spread, about $10^{-4}$ in absolute EQE.

**Why the restart control was needed.** In the original per-band campaign, which was randomly seeded, the 0–20° and 20–40° gains came out below unity, at 0.946 and 0.966. Because the hemisphere is a point of the 13-variable feasible set, $G_j \ge 1$ holds for an exhaustive search. A value below unity is therefore a statement about the search, not about optics. It was also the strongest objection to the main result: a search that fails to recover a solution it contains cannot be trusted when it reports that nothing better exists. The restart control removes this ambiguity.

**Protocol.** Each arm's search restarts at that arm's hemispherical optimum: the fixed quarter-circle control points with that arm's optimized cavity thicknesses and lens height. The hemisphere point and eight perturbations within 8% of each variable's range form the seed set. The search runs three ways: `surrogateopt` from the seed set, `patternsearch` from the surrogate winner, and `patternsearch` launched directly from the hemisphere point. The last is the most direct test of whether the hemisphere is a local optimum in thirteen dimensions. The winner is re-evaluated at final fidelity. The hemisphere baseline is re-measured in the same session at final fidelity. The archived value serves only as a cross-check. An arm counts as improved when the difference exceeds the pooled one-sided 95% $t$ value at $2(N_{\mathrm{rep}}-1) = 4$ degrees of freedom, 2.13.

**Table S4 | Restart control.**

| Arm | hemisphere (archived) | hemisphere (re-measured) | dev. | restarted | gain | $t$ | winning branch |
|---|---|---|---|---|---|---|---|
| 0–20° | 0.06682 | 0.06677 ± 0.00017 | 0.07% | 0.06699 ± 0.00013 | +0.32% | 1.7 | polish from hemisphere |
| 20–40° | 0.16591 | 0.16595 ± 0.00017 | 0.02% | 0.16628 ± 0.00020 | +0.20% ‡ | 2.1 ‡ | surrogateopt |
| 40–60° | 0.18503 | 0.18486 ± 0.00008 | 0.09% | 0.19355 ± 0.00010 | **+4.70%** | **114.6** | polish from surrogate |
| 60–80° | 0.13589 | 0.13581 ± 0.00014 | 0.06% | 0.14075 ± 0.00025 | **+3.64%** | **29.6** | surrogateopt |
| total EQE | 0.54679 | 0.54679 ± 0.00013 | 0.00% | 0.55231 ± 0.00004 | **+1.01%** | **68.7** | polish from surrogate |

Uncertainties are the standard deviation over three high-precision repeats.

‡ The 20–40° arm was the only borderline case. It was settled by re-measuring both stored designs with five fresh repeats each, with the search excluded. This gives hemisphere 0.16583 ± 0.00015 against restarted 0.16631 ± 0.00021, a residue of +0.29% at $t = 4.1$ (five-repeat threshold 1.86). The residue is statistically real but three tenths of a percent in size.

**Reading the control.** The control is one-sided: the hemisphere is among the screened candidates, so the result cannot fall meaningfully below it. It is also biased toward reporting a gain, because the winner is the maximum over roughly seventy noisy search evaluations per arm. The near-zero outcomes are therefore conservative. The two groups differ in scale: residues of +0.29% to +0.32%, against margins of +1.01% to +4.70% at $t = 30$ to 115. This separation shows that the residues reflect an essentially absent margin rather than a weak search. In the 60–80° arm the restart also beat the original campaign's own optimum, so starting at the hemisphere is simply a better search protocol in this design space. The two original sub-unity values were a search artifact, and the corrected values are 1.003 and 1.002.

**Fig. S2 | Restart control.** (a) Hemisphere baseline (re-measured in the same session) and restarted optimum for every arm, with standard deviations over three repeats and the $t$ value of each difference. (b) Five-repeat confirmation of the 20–40° arm: individual repeats and means for both designs ($t = 4.1$, +0.29%). *(`figS1_restart_control.png`)*

---

## Note S3 | Patch-size dependence and normalization

One fixed high-efficiency design, the weighted-sweep optimum at $w = 0.75$, was re-evaluated at final fidelity (50,000 rays, 151 λ) with only the textured patch size varied. The 15/25/35 mm rows are 3 repeats (`patch_convergence_result.mat`). The 100 mm row is 2 repeats from a follow-up run (`patch_convergence_100.mat`).

**Table S5 | Patch-size dependence.**

| Patch (mm) | total EQE | s.d. | vs 25 mm | selectivity (0–20 / 20–40 / 40–60 / 60–80) |
|---|---|---|---|---|
| 15 × 15 | 0.51681 | 0.00016 | −4.79% | 0.108 / 0.283 / 0.346 / 0.233 |
| 25 × 25 | 0.54282 | 0.00008 | — | 0.106 / 0.280 / 0.346 / 0.238 |
| 35 × 35 | 0.55128 | 0.00003 | +1.56% | 0.106 / 0.279 / 0.346 / 0.239 |
| 100 × 100 | 0.5636 | 0.00003 | +3.83% | 0.106 / 0.277 / 0.347 / 0.240 |

Three observations follow.

1. The initial rise saturates with a decay length of about 9 mm, four critical-angle round trips of 2.32 mm each. The far tail does not saturate. Between 35 and 100 mm the total still gains 2.2%, about three times what a single exponential fitted to the first three points would allow. Light that makes many round trips migrates laterally over many tens of round trips before escaping. The total EQE of a finite film therefore increases slowly with its extent, and the 25 mm value is a lower bound, 3.7% below the 100 mm value.
2. The angular composition is converged with respect to patch size. While the total gains 3.8% from 25 to 100 mm, no band selectivity changes by more than 0.24 pp (0.11 pp between 25 and 35 mm). All ratios, correlations, and class comparisons are evaluated at equal patch and are unaffected.
3. The weighted-sweep campaign's archived best total, 0.5556, falls between the 35 mm and 100 mm values of the same design, while the design re-measures at 0.5428 on the 25 mm patch. That campaign's absolute normalization therefore corresponds to a larger effective patch. Its absolute totals are excluded from the $G_j$ comparison (Table S3). Its compositions and correlations, being ratios, are used. Its coarse-fidelity log maximum (0.5591 at 10,000 rays) is not used anywhere.

**Fig. S3 | Patch-size dependence.** (a) Total EQE of one fixed design versus textured patch size, with the hemisphere baseline at 25 mm for scale. (b) Band selectivity of the same design versus patch size: the composition does not move. *(`figS2_patch_dependence.png`)*

---

## Note S4 | Angular-composition statistics

**Natural composition.** The natural composition is the median selectivity of the twenty designs with the highest total EQE in a campaign log. For the convex per-band log, whose top designs were located by four different single-band objectives, it is 0.094 / 0.278 / 0.361 / 0.236, with 10–90% spreads of 0.092–0.096, 0.277–0.281, 0.358–0.364, and 0.233–0.238. For the weighted-sum log it is 0.102 / 0.276 / 0.354 / 0.239.

**Sensitivity to the choice of twenty.** For the weighted-sum log:

**Table S6a | Natural composition for different top-k choices (weighted-sum log).**

| Statistic | 0–20° | 20–40° | 40–60° | 60–80° | bootstrap 95% half-width |
|---|---|---|---|---|---|
| top-10 median | 0.100 | 0.277 | 0.353 | 0.239 | 0.003–0.005 |
| top-20 median | 0.102 | 0.276 | 0.354 | 0.239 | 0.002–0.003 |
| top-50 median | 0.100 | 0.277 | 0.356 | 0.237 | 0.001–0.002 |
| top-20 mean | 0.101 | 0.276 | 0.353 | 0.239 | — |

The choices differ by at most 0.003 in any band.

**Dedicated versus natural.** The band-dedicated optima of the per-band campaign reach selectivities of 0.119 / 0.300 / 0.364 / 0.283. These are +27%, +8%, +1%, and +20% above the natural medians. The same designs lose total EQE: 3%, 2%, 1%, and 11% below that campaign's best total of 0.548. The net band-power gain is the product of the selectivity gain and the total-EQE ratio: 1.22, 1.05, 1.00, and 1.07 (Fig. S4). A stricter benchmark compares each dedicated optimum with the best value of the same band reached incidentally while another band was being optimized. Against that benchmark the dedicated optima win by only 6.3%, 0.6%, 4.5%, and 5.8%.

**Fig. S4 | Gain decomposition of the band-dedicated optima.** Selectivity gain $S_{\mathrm{win}}/S_{\mathrm{nat}}$, total-EQE ratio $E_{\mathrm{win}}/E_{\max}$, and their product, the net band-power gain (labelled). *(`figS3_gain_decomposition.png`)*

**Selectivity across the design space.** Fig. S5 shows the four bands separately, over all 606 usable evaluations of the weighted-sum campaign. The log holds 691 evaluations. The remaining 85 traces failed to return a positive total EQE. The correlations between total EQE and $S_j$ are +0.60, +0.55, −0.12, and −0.57. With the bands normalized to 0–80° instead of 0–90° they become +0.60, +0.55, −0.23, and −0.59. The 80–90° remainder is 3.1% for the median design (2.1–8.2% over all designs) and falls as efficiency rises (correlation −0.63).

**Collapse statistics.** Least-squares fits of $\mathrm{EQE}_j$ against $\mathrm{EQE}_{\mathrm{total}}$ over the 606 usable designs (95% intervals):

**Table S6b | Band-EQE collapse.**

| Band | slope | intercept | $R^2$ (all 606) | $R^2$ (114 random only) |
|---|---|---|---|---|
| 0–20° | 0.125 ± 0.004 | −0.014 ± 0.002 | 0.848 | 0.828 |
| 20–40° | 0.318 ± 0.007 | −0.020 ± 0.004 | 0.931 | 0.925 |
| 40–60° | 0.348 ± 0.005 | 0.003 ± 0.003 | 0.968 | 0.967 |
| 60–80° | 0.189 ± 0.008 | 0.026 ± 0.004 | 0.769 | 0.724 |

The random designs alone reproduce the collapse, so it is not produced by correlated points along optimizer trajectories. Because each band is a part of the total, the substantive content is the narrow spread of $S_{40\text{–}60}$, 0.336–0.375 (5–95%) over all designs.

**Band ordering.** The tilt toward lower polar angles does not reorder the bands among efficient designs:

- In all 520 designs with total EQE ≥ 0.40, the 40–60° band carries the largest share.
- Among the 421 designs with total EQE ≥ 0.50, 416 follow the order 40–60° > 20–40° > 60–80° > 0–20°.
- The 60–80° band is the largest in 8 of the 606 designs, all with total EQE ≤ 0.354.
- The 60–80° share exceeds the 20–40° share in 85 of the 606 designs (14%). Of these, 50 have total EQE < 0.40 and only 5 have total EQE ≥ 0.50; these 5 are the exceptions to the full ordering above.

The same count in the other logs follows. The 40–60° band is the largest in 100% of the convex per-band designs, 99.4% of the concave designs, and 100% of the randomly assembled designs. The 60–80° share exceeds the 20–40° share in 19.2%, 28.8%, and 4.6% of them, respectively.

**Fig. S5 | Band selectivity across the design space.** $S_j$ versus total EQE over all 606 usable evaluations of the weighted-sum campaign, separated into random feasible designs (grey) and optimizer-visited designs (colored), one panel per band. Black lines are least-squares fits; red dashed lines mark the Lambertian partition. Each vertical axis spans only a few percentage points. *(`figS4_selectivity_map.png`)*

**Family comparison.** Two statistics are given for each family, both computed from the family's own log: the top-twenty median and the population mean over every usable evaluation. Both are needed. The gap between them within one family is comparable to the gap between families under either statistic. A table with only one statistic would suggest family differences that the data do not support.

**Table S6c | Family comparison.**

| Family | best total EQE | top-20 median (0–20 / 20–40 / 40–60 / 60–80) | population mean |
|---|---|---|---|
| Hemispherical reference | 0.5468 | — | — |
| Convex freeform (per-band campaign) | 0.5481 | 0.094 / 0.278 / 0.361 / 0.236 | 0.097 / 0.278 / 0.350 / 0.242 |
| Convex freeform (weighted sweep) † | 0.5556 | 0.102 / 0.276 / 0.354 / 0.239 | 0.095 / 0.274 / 0.354 / 0.244 |
| Convex freeform (total EQE only) | 0.5539 | — | — |
| Randomly assembled | 0.5216 ± 0.0011 | 0.111 / 0.299 / 0.347 / 0.213 | 0.102 / 0.289 / 0.355 / 0.224 |
| Inverted (concave) | 0.5165 | 0.113 / 0.303 / 0.340 / 0.215 | 0.100 / 0.288 / 0.357 / 0.231 |
| Lambertian partition | — | 0.117 / 0.296 / 0.337 / 0.220 | — |

† Large-patch normalization (Note S3); this design re-measures at 0.5428 on the 25 mm patch. Its compositions are ratios and are comparable with the other rows.

Every entry lies within 0.094–0.113 (0–20°), 0.274–0.303 (20–40°), 0.340–0.361 (40–60°), and 0.213–0.244 (60–80°). For the randomly assembled family, the mean over the unbiased random-sampling phase alone is 0.095 / 0.281 / 0.360 / 0.233. Its three high-precision winning realizations give 0.110 / 0.298 / 0.349 / 0.213. Both lie inside the same window. The uncertainty quoted for this family is the spread over three independent disorder realizations of the same assembly statistics (coefficient of variation 0.2%). In the concave family, band-dedicated optimization gives net band-power gains of 0.98, 0.99, 1.05, and 1.19. The largest, at 60–80°, raises selectivity by 35% while losing 12% of total EQE.

---

**Fig. S6 | Total efficiency fixes the angular composition.** (a) Total EQE versus 40–60° band EQE for random feasible designs (grey) and optimizer-visited designs (green) of the weighted-sum campaign. All points collapse onto one near-linear locus (black line, $R^2 = 0.968$, $n = 606$). (b) Selectivity of each band-dedicated optimum (red) against the natural composition of the twenty most efficient designs (grey, 10–90% range) and the Lambertian partition (dashed), labelled with the net band-power gain. (c) Selectivity versus total EQE for all four bands, with least-squares fits and the correlation $R$ in the legend. Among designs with total EQE ≥ 0.40 the 40–60° band is always the largest; the 20–40° and 60–80° fits cross only in the inefficient region. *(`fig3_efficiency_composition.png`)*

---

## Note S5 | Convergence of the selectivity–efficiency drift

Twenty designs stratified across the efficiency range were re-evaluated under three conditions of increasing strictness: the baseline ray count; twenty-fold ray count (200,000 rays) with three repeats; and broadband evaluation over 450–750 nm. Because the subset is stratified rather than representative, its correlations differ in magnitude from the full-population values (+0.60 / +0.55 / −0.12 / −0.57). The outer-band signs agree, and the 40–60° value is near zero in both.

**Table S7 | Convergence of the drift correlations.**

| Band | baseline | 200,000 rays, ×3 repeats | broadband 450–750 nm |
|---|---|---|---|
| 0–20° | +0.59 | +0.61 | +0.64 |
| 20–40° | +0.67 | +0.67 | +0.72 |
| 40–60° | +0.07 | +0.04 | +0.05 |
| 60–80° | −0.70 | −0.71 | −0.76 |

The drift is therefore neither Monte-Carlo noise nor a narrowband artifact.

**Fig. S7 | Convergence of the drift.** Correlation between total EQE and $S_j$ for the stratified 20-design subset under the three evaluation conditions of Table S7. *(`figS5_convergence.png`)*

---

## Note S6 | Cost calibration for the randomly assembled family

**Supercell construction.** Lenslets sit on a jittered hexagonal lattice. Each lenslet radius is drawn uniformly within $r_{\mathrm{mean}}(1 \pm r_{\mathrm{jitter}})$ and each aspect ratio within $\mathrm{AR}_{\mathrm{mean}}(1 \pm \mathrm{AR}_{\mathrm{jitter}})$. Positions are displaced uniformly within a disc whose radius is a fraction (`posJitter`) of the largest displacement that keeps neighbours from overlapping, $(p - 2r_{\max})/2$. Each profile is a convex blend, with weight `profileMix`, of the hemispherical control points and an independently drawn valid set of freeform control points ($x$ sorted in 0–1, $y$ in 0–1.5, the same feasibility rules as Table 1). Lenslets are kept inside the supercell so that tiling produces no seams. A realization is fixed by its random seed.

One fixed geometry (seed 7777, mid-range assembly statistics) was re-evaluated under each setting. Because the reported quantities are ratios, the acceptance criterion is band selectivity. Δ selectivity is the largest absolute deviation from the reference across the four bands, in percentage points.

**Table S8 | Cost calibration.**

| Setting | lenslets | grid | rays | λ step | λ | time (s) | speed-up | EQE | Δ selectivity (pp) |
|---|---|---|---|---|---|---|---|---|---|
| reference | 8 × 8 | 201 | 10,000 | 10 nm | 31 | 777.5 | 1.00 | 0.4975 | — |
| half rays | 8 × 8 | 201 | 5,000 | 10 nm | 31 | 529.3 | 1.47 | 0.4992 | 0.38 |
| coarser λ | 8 × 8 | 201 | 10,000 | 20 nm | 16 | 446.5 | 1.74 | 0.5012 | 0.66 |
| coarser grid | 8 × 8 | 141 | 10,000 | 10 nm | 31 | 704.6 | 1.10 | 0.5016 | 0.28 |
| fewer lenslets | 6 × 6 | 201 | 10,000 | 10 nm | 31 | 543.3 | 1.43 | 0.4947 | 0.54 |
| **adopted** | **6 × 6** | **141** | **5,000** | **20 nm** | **16** | **181.0** | **4.30** | **0.4950** | **0.39** |

The adopted combination is 4.3× faster and shifts selectivity by 0.39 pp. This equals the Monte-Carlo spread measured by repeating one design (0.5% relative on total EQE). Two individual reductions exceed 0.5 pp on their own, but their effects do not accumulate. This indicates that the deviations are noise rather than bias. The accuracy columns are exactly reproducible because the seed is fixed. The wall-clock times vary by a few percent with machine load. The supercell carries 6 × 6 lenslets on a 141 × 141 height grid. The disorder correlation length is one lenslet pitch, a sixth of the supercell.

**Fig. S8 | Supercell cost calibration.** Evaluation time (bars, with speed-up) and largest selectivity deviation from the reference (red markers) for each sampling reduction and for the adopted combination. *(`figS6_cost_calibration.png`)*

---

## Note S7 | Inverted and randomly assembled MLA families

We next ask whether the result depends on the array construction. We repeat the analysis for the inverted (concave) and the randomly assembled families (Fig. S9), using the same statistics, each computed from the family's own evaluation log. The concave family uses exactly the variable bounds and feasibility rules of main-text Table 1. Three signatures of saturation recur in all three families:

- **Collapse.** The (total EQE, 40–60° EQE) points again collapse onto a line. $R^2$ is 0.94 for the convex per-band campaign, 0.94 for the concave family, and 0.84 for the random family (Fig. S9, left column). The logs differ in sampling density and, for the random family, in fidelity, so we read these values only as evidence that the collapse recurs, not as a ranking of the families.
- **Narrow composition.** The concave family's best total EQE is 0.517, 6% below the convex 0.548, as expected for a recessed rather than protruding surface. Its natural composition is 0.113 / 0.303 / 0.340 / 0.215. The random family reaches 0.522 ± 0.001 over three disorder realizations, a coefficient of variation of 0.2%. Only the assembly statistics matter, not the particular realization. Its natural composition is 0.111 / 0.299 / 0.347 / 0.213. Band-dedicated optimization of the concave family gives a largest net band-power gain of 1.19, in the 60–80° band, close to the convex family's 1.22.
- **Same direction of tilt.** Selectivity again rises with efficiency at 0–40° and falls at 60–80° (Fig. S9, right column). The 40–60° correlation is weak in the two periodic families and negative in the random one. In every family the 40–60° band remains the largest in at least 99% of sampled designs.

The families also bound how far the composition can be moved at all. Every high-efficiency composition we obtained lies within 0.09–0.11 (0–20°), 0.28–0.30 (20–40°), 0.34–0.36 (40–60°), and 0.21–0.24 (60–80°). This spans three constructions, four single-band objectives, and two averaging conventions (top-twenty median and population mean). The change from convex to concave lenslets moves the composition to the edge of this window. The shift is at most 0.02 in any band, about the same as switching the averaging convention within one family (Note S4, Table S6c). Positional disorder and per-lenslet shape disorder leave the composition inside the window. The composition is therefore set by the platform—extended source, thick substrate, coextensive external film—not by any design variable within it.

The randomly assembled MLA should not be confused with a volumetric scattering layer. Each lenslet remains a refractive surface with deterministic local normals and a single index step. A scattering layer instead redirects light by multiple volumetric scattering. Its behaviour depends on the particle–matrix index contrast and loading, which are not parameters of this study. Our conclusions apply to refractive external films and make no comparison with scattering layers.

**Fig. S9 | Generality across external-film MLA families.** Rows: convex freeform (reference, per-band campaign), inverted (concave), and randomly assembled (pseudo-random supercell). Every statistic is computed from each family's own evaluation log. (a1–c1) Total EQE versus 40–60° band EQE with a least-squares line; the near-linear collapse recurs in all three families. (a2–c2) Natural composition (top-twenty median, colored) against the convex reference (grey) and the Lambertian partition (dashed outline). The open circle is the population mean of the same family; the spread between statistics is comparable to the spread between families. (a3–c3) Correlation $R$ between total EQE and $S_j$, per band, against the convex reference. *(`fig4_families.png`)*

---

## Note S8 | Asymmetric three-dimensional freeform exploration

To check whether axial symmetry itself limits the outcome, three-dimensional asymmetric freeform lenslets were optimized directly against a restricted $(\theta,\phi)$ window, with the objective

$$
\mathrm{EQE}_{\mathrm{win}}=\int_{\theta_1}^{\theta_2}\!\!\int_{\phi_1}^{\phi_2} I_{\mathrm{air}}(\theta,\phi)\,\sin\theta\,d\phi\,d\theta .
$$

This was a separate, differently constrained study. It used a different OLED stack (Ag rather than ITO electrode), an anisotropic emitter cell, and 52 shape variables rather than 13. It is therefore not part of the same-constraints benchmark.

The asymmetric shapes could shift the centroid and fine structure of the far field. They did not reliably achieve an absolute window power appreciably above the hemispherical MLA. We report this comparison qualitatively and give no ratio. The stack, emitter cell, and parameter count all differ from the controlled benchmark, so a numerical gain would invite a side-by-side reading that the scope does not support. What carries over is the sign of the effect, not its magnitude. Periodic arrays can steer light under other conditions: asymmetric prisms, diffractive elements, metasurfaces, or sufficient aperture expansion can produce directional emission (see, e.g., main ref. 10). Within coextensive refractive MLAs on a uniform extended source, such redistribution did not become a useful output channel beyond the hemispherical reference.

---

## Note S9 | Angle retention over many passes and the effect of a randomizing rear

**Purpose.** This note quantifies why reshaping the lens buys little, by following the light pass by pass, and tests whether breaking the angle retention of the returned light would restore leverage to the lens shape.

**Method.** The bidirectional scattering distribution functions (BSDFs) of the optimized lenses of this work were computed with a Monte-Carlo height-field ray tracer for a hexagonally close-packed array, index-matched to the substrate ($n = 1.51$). The tracer reproduces the LightTools BSDF of a hemispherical MLA: escape probability 0.4275 against 0.4274, and angle-resolved transmittance within 0.01. The profiles were rebuilt from the optimized control points with a PCHIP spline; LightTools uses its own spline, so the freeform BSDFs approximate the traced geometry. The 0–20° freeform optimum is excluded because it ends in a vertical sidewall at the rim, which a height-field tracer cannot represent. The BSDFs enter the angle-resolved multi-pass model of main ref. 11: the light in the substrate is a vector over polar angle, the film transmits $B_T$ and returns $B_R$, and the device reflects with $R(\theta)$ from a CPS/transfer-matrix model of the same stack at 550 nm (Al or Ag cathode, ITO 150 nm, horizontal dipole ratio 0.865). Two device reflections are compared: angle-preserving (the actual device) and a Lambertian reset, in which the returned light is redistributed as $\cos\theta\sin\theta$ at every pass. To separate lens from cavity, all lenses were also evaluated with one common cavity (ETL 61.3 nm, HTL 97.5 nm, the hemispherical total-EQE optimum).

**Depletion.** With the angle-preserving device, the escape probability of the hemispherical total-EQE design falls pass by pass: 0.43, 0.39, 0.35, 0.31, 0.27, 0.23 (Al) and 0.41, 0.37, 0.34, 0.30, 0.27, 0.24 (Ag). Light at angles the lens transmits well leaves first, and the remainder sits at angles it transmits poorly. With the Lambertian reset the escape probability stays at 0.42.

**Table S9 | Extraction efficiency at a common cavity (550 nm).**

| Lens | Al, angle-preserving | Al, reset | Ag, angle-preserving | Ag, reset |
|---|---|---|---|---|
| hemisphere (total-EQE height) | 0.808 | 0.834 | 0.926 | 0.944 |
| hemispheres, all five heights (range) | 0.800–0.809 | 0.832–0.834 | 0.921–0.927 | 0.943–0.944 |
| freeform 20–40°, 40–60°, 60–80°, total (range) | 0.789–0.807 | 0.827–0.835 | 0.915–0.931 | 0.942–0.945 |

The reset raises the extraction by 2.7–5.5% (Al) and 1.2–3.0% (Ag), and the spread between lenses falls from 2.5% to 1% (Al) and from 1.7% to 0.3% (Ag). The output composition converges as well: across the nine lenses the range of $S_{40\text{–}60}$ falls from 0.325–0.354 to 0.335–0.346. The gain disappears if the randomizing element adds round-trip loss of more than about 3.1 percentage points (median over lenses; 2.3–4.8) with Al, or about 1.4 percentage points (0.9–2.2) with Ag.

**Lens versus cavity.** With each design's own cavity, the 60–80° hemisphere design reaches $S_{60\text{–}80}$ = 0.303, against 0.247 for the total-EQE design. With the common cavity the same two lenses give 0.249 and 0.247, and the 60–80° freeform lens gives 0.252. In this model, therefore, almost all of the band steering of the dedicated designs comes from the cavity, and the lens shape adds about 0.005. A fixed-cavity LightTools control is in progress to confirm this at full fidelity.

**Consequences.** Breaking the angle retention of the returned light does not restore leverage to the lens shape; it removes most of what remains, because every lens then sees the same Lambertian input on every pass.

**Fig. S10 | Angle-preserving versus Lambertian-reset device reflection.** (a) Escape probability per pass for the hemispherical total-EQE design with the Al cathode. (b) Extraction efficiency of nine optimized lenses at a common cavity, with Al and Ag cathodes. (c) Range of the output composition over the nine lenses (Al), with the Lambertian partition for reference. *(`figS8_angle_reset.png`)*

---

## Note S10 | Angular mixing number: definition, robustness, and validation

**Definition.** The return matrix $B_R(\theta_i\to\theta_r)$ of a structure gives the fraction of the power incident in polar bin $i$ that re-enters the substrate in bin $r$. Normalizing every column by the returned fraction $1-p(\theta_i)$ gives the column-stochastic kernel $M$ of the conditional Markov chain of the returned light. When transmission or absorption removes light, the surviving light relaxes to a quasi-stationary distribution; normalizing the columns separates this escape from the redistribution itself. For a reciprocal structure $w_i B_R(i\to r) = w_r B_R(r\to i)$ with the étendue weight $w_k = \sin^2\theta_{k+1}-\sin^2\theta_k$, so $D^{1/2} M D^{-1/2}$ is symmetric for a suitable diagonal $D$, and the spectrum of $M$ is real. We define $N_{\mathrm{mix}} = 1/(1-|\lambda_2|)$, with $\lambda_2$ the eigenvalue of second-largest modulus (the absolute spectral gap of main ref. 20). For every structure in this work the second eigenvalue was real and positive and also the second largest in modulus. $N_{\mathrm{mix}}$ is the relaxation time, not the e-folding pass number $-1/\ln|\lambda_2|$; for $\lambda_2$ = 0.5 the two are 2 and 1.44.

**Grazing-angle truncation.** Bins near 90° carry almost no étendue but can return light at their own angle; a few such states can then dominate $\lambda_2$ without affecting extraction. We therefore evaluate $\lambda_2$ on $\theta < 85°$, which excludes 0.76% of the étendue. Table S10 shows the effect.

**Table S10 | Sensitivity of $|\lambda_2|$ to the angular range.**

| Structure | 0–90° | 0–85° | 0–80° | 0–70° | weight of second eigenvector above 80° |
|---|---|---|---|---|---|
| hemispherical MLA (total-EQE height) | 0.876 | 0.859 | 0.846 | 0.840 | 0.23 |
| freeform MLA (40–60° optimum) | 0.941 | 0.867 | 0.689 | 0.532 | 0.49 |
| random lenslet shapes, s = 0.2 | 0.679 | 0.628 | 0.514 | 0.419 | 0.30 |
| random lenslet heights, spread 0.2 | 0.935 | 0.859 | 0.826 | 0.787 | 0.46 |

For the hemisphere the value is insensitive to the range. For the freeform 40–60° lens and for height-randomized hemispheres, the full-range value is set by grazing states, and truncation at 85° brings it to the hemisphere level. Truncation at 85° leaves the mean error of the reduced model unchanged (0.57 pp) and lowers the largest error for freeform lenses from 4.1 to 2.2 pp; truncation at 80° raises the mean error to 0.80 pp.

**Validation.** The reduced model of main-text Sec. 2.3 was compared with the full series for every structure of the library (Note S11), with the Al, Ag, and low-loss devices at 500, 550, and 600 nm. Table S11 summarizes the errors by class. [Table S11 and the spectrally averaged comparison to be completed with the final library.]

---

## Note S11 | Library of structures and return matrices

**Height-field tracer.** Microlens arrays, random lens arrays, and wavy surfaces are traced with a Monte-Carlo height-field ray tracer: the surface $z = h(x, y)$ separates the index-matched medium ($n = 1.51$) from air, rays are marched with a Lipschitz-safe step and refined by bisection, and each interface event chooses reflection or refraction by the unpolarized Fresnel coefficients. Rays start at the base plane with a uniform position over the periodic cell and a uniform azimuth; 20,000 rays per 1° incidence bin are used (10,000 for wavy surfaces). The transmitted rays are histogrammed by their polar angle in air ($B_T$), the returned rays by their polar angle in the substrate ($B_R$). For a single hemispherical lenslet the tracer reproduces the LightTools return matrix of main ref. 11: escape probability 0.4275 against 0.4274, angle-resolved transmittance within 0.008, and less than 0.1% of rays undetermined.

**Random lens arrays.** A 6 × 6 supercell of lenslets is tiled periodically. Every lenslet may carry its own spline profile, its own height, and a jittered centre and radius. Shape spread $s$ perturbs the five control points of the quarter circle by Gaussian noise of standard deviation $0.25s$ (radial) and $0.4s$ (height); $s$ = 1 approaches independently drawn profiles. Where lenslets overlap, the higher surface is taken; gaps are flat. Profiles with near-vertical walls, which a height field cannot represent, are rejected.

**Wavy surfaces.** Periodic, Gaussian-correlated random height fields $h = h_0 + \sum_m a_m\cos(\mathbf{g}_m\cdot\mathbf{x}+\phi_m)$ on a 12 × 12 box with correlation length one, scaled to a prescribed rms slope. Ray optics is scale-free, so only the slope statistics matter; the physical correlation length must exceed a few wavelengths.

**Scattering layers.** Non-absorbing, index-matched slabs with a Henyey–Greenstein phase function, Monte Carlo with 20,000 rays per 0.5° incidence bin, combined with the flat exit face by the adding method.

**Table S12 | Randomness sweep.** [To be completed when the sweep finishes: $p$, $|\lambda_2|$, $N_{\mathrm{mix}}$, $\eta_{\mathrm{sta}}$ (Al, Ag, low loss), $S_{40\text{–}60}$ for every random structure.]

---

## Reproducing the figures and data

Main-text Figs. 1–4 are produced by `selectivity/oe_figs.py`, which builds the structure library (Note S11), evaluates the full and reduced models, and writes each figure to `figures_oe/` as PNG (300 dpi), PDF, and an Excel workbook with one sheet per panel holding exactly the plotted arrays. The random-structure return matrices are generated by `multipass_reset/run_random_sweep.py` (tracer `hftrace_gen.py`). Main-text Fig. 5 and Figs. S1–S10 come from the freeform campaign archives through `make_figures.py` and `make_supp_figures_and_data.py`.

```
cd multipass_reset && python3 run_random_sweep.py && cd ..
cd selectivity && python3 extract_inputs.py && python3 run_devices.py && python3 oe_figs.py && cd ..
python3 make_figures.py
python3 make_supp_figures_and_data.py
```

**Table S13 | Figure files and their inputs.**

| Figure | Output | Inputs |
|---|---|---|
| Fig. 1 | `figures_oe/fig1_mixing_concept` | `selectivity/data/devices.npz`, `lt_hemisphere_bsdf.mat`, `multipass_reset/random_sweep/` |
| Fig. 2 | `figures_oe/fig2_reduced_model` | structure library (Note S11), `devices.npz` |
| Fig. 3 | `figures_oe/fig3_structure_map` | structure library, `devices.npz` |
| Fig. 4 | `figures_oe/fig4_bound_window` | structure library, `devices.npz`, `opt_4band_result_25by25.mat` |
| Fig. 5 | `figures_oe/fig5_hemisphere_benchmark` | `opt_4band_result_25by25.mat`, `opt_hemisphere_result.mat`, `warmstart_hemisphere_result.mat`, `freeform_EQEtotal_result.mat` |
| Fig. S1 | `fig1_platform` | `opt_4band_result_25by25.mat`, `opt_hemisphere_result.mat` |
| Fig. S2 | `figS1_restart_control` | `warmstart_hemisphere_result.mat`, `reeval_confirm_2040_result.mat` |
| Fig. S3 | `figS2_patch_dependence` | `patch_convergence_result.mat`, `patch_convergence_100.mat` |
| Fig. S4 | `figS3_gain_decomposition` | `pareto_front_result.mat`, `opt_4band_result_25by25.mat` |
| Fig. S5 | `figS4_selectivity_map` | `pareto_front_result.mat` |
| Fig. S6 | `fig3_efficiency_composition` | `pareto_front_result.mat`, `opt_4band_result_25by25.mat` |
| Fig. S7 | `figS5_convergence` | `convergence_check_result.mat` |
| Fig. S8 | `figS6_cost_calibration` | `calibrate_random_cost.mat` |
| Fig. S9 | `fig4_families` | `opt_4band_result_25by25.mat`, `opt_4band_inverted_result.mat`, `stress_random_result.mat` |
| Fig. S10 | `figS8_angle_reset` | `multipass_reset/` (`hftrace.py`, `multipass.py`, `bsdfs.npz`, `fig_reset.py`) |

Each output exists as `.png`, `.pdf`, and `.xlsx`.

Result archives: `pareto_front_result.mat`, `opt_4band_result_25by25.mat`, `freeform_EQEtotal_result.mat`, `opt_hemisphere_result.mat`, `opt_4band_inverted_result.mat`, `stress_random_result.mat`, `calibrate_random_cost.mat`, `warmstart_hemisphere_result.mat`, `reeval_confirm_2040_result.mat`, `patch_convergence_result.mat`, `patch_convergence_100.mat`, `convergence_check_result.mat`, `multipass_reset/random_sweep/*.npz`.
