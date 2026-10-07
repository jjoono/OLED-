# Practical Saturation of Freeform Microlens Arrays on Extended OLED Emitters and Design Routes beyond Lens Shape

**Authors** — [TBD]

---

## Abstract

Freeform microlens arrays (MLAs) are widely expected to control both the light-extraction efficiency and the angular emission of organic light-emitting diodes (OLEDs). We test this expectation for MLAs tiled over an extended emitter. We couple an experimentally validated dipole-microcavity source model to three-dimensional ray tracing and compare an axisymmetric freeform MLA with an equally optimized hemispherical MLA under identical constraints. Total external quantum efficiency (EQE) and four polar bands are optimized as separate objectives. The freeform optimum exceeds the hemisphere by at most 7% in any band and by 1.3% in total EQE. A control that restarts every search at the hemispherical optimum shows that these small margins come from the design space, not from a weak search. All sampled designs collapse onto one line relating band power to total EQE. Among efficient designs the angular composition stays within a narrow window. This holds for convex, concave, and randomly assembled arrays alike. We trace the saturation to three causes: the output/source area ratio of a coextensive array is fixed, the thick substrate mixes light laterally, and a passive film cannot raise radiance. We then map the design routes that remain productive.

---

## 1. Introduction

A planar OLED loses much of its generated light. Waveguided modes in the organic layers, losses at the metal electrode, and total internal reflection (TIR) at the substrate/air interface all trap photons [1,2]. External MLAs extract the substrate modes without touching the electrical stack, which makes them attractive for large-area and flexible OLEDs [3–5]. Hemispherical and near-hemispherical MLAs already give high efficiencies [4,5]. High total extraction has also been reached without lenses, by combining an external scattering layer with horizontally oriented emitters [12].

More complex lens shapes are still expected to do better. Freeform illumination optics and inverse design are powerful tools for producing a prescribed intensity distribution. They have recently been applied to micro-LED and OLED packaging [6,7]. Tuning lenslet curvature, height, asymmetry, and the microcavity together should, in principle, give finer control over both total extraction and emission into chosen directions.

That expectation comes from optics built around point-like sources, and it may not carry over to an OLED. Freeform design methods were developed for sources of near-zero étendue. A thin-film OLED is an extended source, and its light reaches the array only after travelling, and partly recycling, inside a millimetre-thick substrate. Whether shape freedom pays off in this regime is therefore an open question, and the answer has a practical cost. Freeform molding raises fabrication and metrology effort. If it buys little over a hemisphere, the next design step should be a change of source, aperture, or recycling path rather than a more complex surface.

We therefore ask: **how much is lens-shape freedom worth for a tiled refractive MLA on an extended OLED?** We answer with a controlled numerical benchmark. We use one source model that has been validated against fabricated devices [15,16] and one set of manufacturable constraints. Under these we compare a hemispherical MLA and an axisymmetric freeform MLA, each optimized for total EQE and for four polar bands. A restart control separates the limits of the search from the limits of the physics. Two further external-film families, inverted (concave) and randomly assembled arrays, test whether the outcome depends on how the array is built. The result is a reproducible numerical observation for this platform, not an impossibility theorem: a *practical saturation* of the explored design space. We explain why it occurs and map which levers remain productive beyond it.

---

## 2. Methods

### 2.1 Source model and experimental validation

Dipole emission in the OLED microcavity is computed with the Chance–Prock–Silbey (CPS) formalism. The result is the substrate-side angular–spectral intensity $I_{\mathrm{sub}}(\theta,\lambda)$, which serves as the source for three-dimensional ray tracing in LightTools, driven from MATLAB through the COM interface (Fig. 1). The stack is Al / ETL / EML / HTL / ITO / glass. The emitting region is a disc of radius 1 mm, and the glass substrate is $d_{\mathrm{sub}} = 1.295$ mm thick with $n = 1.51$. All CPS inputs are listed in Supplement 1, Note S1 (Table S1).

The same pipeline has been validated against fabricated devices with almost the same geometry. In ref. [16], OLEDs with a 1 mm emitting-aperture radius on a 1.4 mm substrate with $n_{\mathrm{sub}} = 1.51$ were made in three configurations: bare, MLA-attached, and with a structured outcoupler. Measured angular profiles and current efficiencies agreed with the simulations in all three cases, with measured maximum EQEs of 35.6%, 35.4%, and 48.0% (Fig. 3 and Table 2 of ref. [16]). In ref. [15], the pipeline reproduced the emission pattern and near-Lambertian angular electroluminescence of an inverted-MLA device. That device reached 58.0% EQE, against 30.5% for its planar control. The validation therefore covers the source model, the ray-optical treatment of microtextured exit faces, and both attached and inverted lenses, on the same emitter platform and substrate geometry used here.

### 2.2 Lens classes and constraints

The exit face carries a hexagonal array of lenslets of about 10 μm radius over a 25 × 25 mm patch. The emitting disc is therefore about 100 lenslet radii across, and the array extends well beyond it. The lens material is N-BK7 ($n = 1.517$ at 589 nm), index-matched to the substrate. The film/air boundary alone sets the escape condition. The lens index is held fixed throughout. It is a property of the platform here, not a design variable.

Two lens classes are compared (Fig. 1b):

- **Hemispherical reference.** The profile is a fixed quarter circle. The lens height and the ETL and HTL thicknesses are free.
- **Axisymmetric freeform.** The profile is a spline with endpoints fixed at (0,1) and (1,0) and five free control points, whose radial coordinates are constrained monotonic. Together with the two cavity thicknesses and the lens height, this gives 13 design variables.

Both classes share the same material, pitch, fill factor, height range, and maximum draft angle. The hemisphere is a member of the freeform set. The freeform constraints are never more permissive than the hemisphere's, so any difference measures shape freedom alone. Cavity thicknesses are optimized together with the lens in every campaign, so each class is compared at its own best cavity. Candidates with self-intersection, negative thickness, or unmanufacturable draft angles are rejected before evaluation.

### 2.3 Figures of merit

The far field is divided into four polar bands over the full azimuth: 0–20°, 20–40°, 40–60°, and 60–80°. The band EQE, $\mathrm{EQE}_j$, is the fraction of generated photons emitted into band $j$. The band selectivity is

$$
S_j=\frac{\mathrm{EQE}_{j}}{\mathrm{EQE}_{\mathrm{total}}} .
$$

As an analytic reference we use the Lambertian partition $S_j^{\mathrm{Lam}} = \sin^2\theta_{\mathrm{hi}} - \sin^2\theta_{\mathrm{lo}}$, which gives 0.117, 0.296, 0.337, and 0.220. The value of shape freedom for objective $j$ is the relative gain

$$
G_j=\frac{\max\left[\mathrm{EQE}_j\mid\mathrm{freeform}\right]}{\max\left[\mathrm{EQE}_j\mid\mathrm{hemisphere}\right]} .
$$

Each maximum is the mean of three high-precision re-evaluations of the best design. Because the hemisphere lies inside the freeform set, $G_j \ge 1$ holds for an exhaustive search.

### 2.4 Optimization and the restart control

Each objective is optimized by surrogate-based global search (MATLAB `surrogateopt`), followed by `patternsearch` polishing from the surrogate winner. The campaigns are:

- one optimization per objective, for total EQE and each of the four bands;
- a dedicated total-EQE campaign with three independent starts;
- a weighted-sum sweep, $J_w = w\,\hat{\eta}_{\mathrm{ext}} + (1-w)\,\hat{P}_{40\text{–}60}$ with $w \in \{0, 0.25, 0.5, 0.75, 1\}$, together with 150 random feasible designs that populate the design space without optimizer bias.

Searches use 10,000 rays and 31 wavelengths. Every reported optimum is re-evaluated three times with 50,000 rays and 151 wavelengths. Budgets for every campaign are listed in Supplement 1, Table S2.

A search that misses part of its own design space could report a false absence of gain. To exclude this, we ran a restart control. For each objective, the search was restarted from that objective's hemispherical optimum, using the identical objective, constraints, and fidelities. The hemisphere point and eight perturbations within 8% of each variable's range were placed in the seed set. A pattern search was also launched directly from the hemisphere point. The hemisphere baseline was re-measured in the same session. An objective counts as improved when the gain exceeds the pooled one-sided 95% threshold, $t = 2.13$ at four degrees of freedom. The control is biased toward reporting a gain, because the winner is the maximum over many noisy evaluations. A null result under it is therefore conservative (Supplement 1, Note S2).

The sampling uncertainty of the selectivity–efficiency correlations was checked on 20 designs spread across the efficiency range. These were re-evaluated with 200,000 rays (three repeats) and over a 450–750 nm broadband spectrum (Note S5). The dependence of absolute EQE on patch size was measured separately (Note S3).

### 2.5 Additional external-film families

Two further manufacturable families test whether the outcome depends on array construction:

- **Inverted (concave) MLA.** The concave counterpart of the same profile class, as realized by dimples in an ultrathin substrate [15]. It uses the same per-band protocol as the convex benchmark: four single-band optimizations of 60 surrogate evaluations plus 15 polish steps each.
- **Randomly assembled MLA.** A 6 × 6-lenslet pseudo-random supercell, tiled over the patch. Lenslet positions are jittered off the hexagonal lattice, and each lenslet carries an independently drawn profile from the same class. Six assembly statistics replace the lens shape as design variables: fill factor, radius jitter, position jitter, mean aspect ratio, aspect jitter, and a hemisphere-to-random profile blend. The campaign samples 50 random realizations and then optimizes these statistics. Because the disorder correlation length is one pitch, a sixth of the supercell, the supercell is statistically equivalent to a fully random array.

A supercell trace is expensive, so the random family uses reduced sampling (5,000 rays and 16 wavelengths during search). A calibration shows that this shifts band selectivity by at most 0.39 percentage points (pp), comparable to the Monte-Carlo spread (Note S6).

All three families are external outcoupling films: they act on light that has already entered the substrate. Internal arrays beneath the electrode, which address waveguided rather than substrate modes [5], change the emitting stack itself. They lie outside this comparison.

**Fig. 1 | Platform and benchmark design.** (a) OLED–substrate–MLA architecture: CPS microcavity dipole source (emitting radius 1 mm), glass substrate ($d_{\mathrm{sub}} = 1.295$ mm), hexagonal lenslet array (~10 μm lens radius, 25 × 25 mm). (b) The two lens classes under identical constraints, drawn in units of the lens radius: the hemispherical reference (3 free variables) and the axisymmetric freeform (13 variables). Open circles are the five free spline control points; the freeform shown is the 40–60° optimum. (c) The four polar bands and the band selectivity $S_j$, taken over the full azimuth. (d) Optimization workflow: CPS source, LightTools ray tracing through the MATLAB COM interface, surrogate-based global search, pattern-search polishing, and high-ray re-evaluation. *(`fig1_platform.png`)*

---

## 3. Results

### 3.1 Freeform versus an equally optimized hemisphere

Figure 2 compares the two lens classes. The band-optimized freeform profiles differ visibly from the hemisphere and from one another (Fig. 2a). These shape differences buy little performance. The relative gains are $G_j$ = 1.003, 1.002, 1.070, and 1.036 for the four bands, and $G$ = 1.013 for total EQE (Fig. 2b). In no objective does the freeform exceed the hemisphere by more than 7%, and in none does it fall below.

The total-EQE ceiling of the freeform class was found twice, by two different procedures. The dedicated campaign reaches 0.5539, with a 0.8% spread across three independent starts. The restart control reaches 0.5523 in the same session as the hemisphere baseline of 0.5468. The two estimates agree within 0.3%.

The restart control also shows how much of the remaining margin is real (Fig. 2c). Gains above the hemisphere are significant in three objectives: +4.70% at 40–60° ($t = 115$), +3.64% at 60–80° ($t = 30$), and +1.01% for total EQE ($t = 69$). In the two lower bands the residual gains are +0.32% at 0–20° ($t = 1.7$) and +0.20% at 20–40° ($t = 2.1$), both below the threshold. A five-repeat re-measurement settles the 20–40° case as a real but tiny residue of +0.29% ($t = 4.1$). The same procedure that resolves a +0.29% residue does not miss gains ten times larger. The small values at 0–40° therefore show that shape freedom is worth little there, not that the search is weak. The original, randomly seeded per-band campaign under-recovered the hemisphere in these two bands. The restart corrects this, and the restart values are the ones reported (Note S2).

The hemisphere is therefore not merely a convenient baseline. It is a practical near-optimum of this structural class. This does not mean the optimum is exactly hemispherical. The best profiles differ from one another, and the low-efficiency part of the design space is diverse. Near the top of the performance range, however, this diversity does not convert into meaningful extra efficiency or band power.

**Fig. 2 | Freeform versus an equally optimized hemisphere.** (a) Band-optimized freeform profiles (colored) against the hemispherical reference at its own optimized height (black), in units of the lens radius. (b) Best EQE reached by each class at each objective, labelled with the relative gain $G_j$. All entries come from campaigns at the same 25 × 25 mm patch. (c) Restart control: gain over the hemisphere when each search is restarted from the hemispherical optimum, labelled with $t$ against three high-precision repeats. Green marks gains above the one-sided 95% threshold ($t = 2.13$). A five-repeat re-measurement resolves the 20–40° residue as +0.29% ($t = 4.1$; Note S2). *(`fig2_hemisphere_benchmark.png`)*

### 3.2 Total efficiency fixes the angular composition

If shape freedom offered an independent angular degree of freedom, optimizing for one band would trade against total EQE and trace out a Pareto front. The data show no such front. Across all 606 usable evaluations of the weighted-sum campaign, the 40–60° band EQE is almost proportional to total EQE ($R^2 = 0.968$, Fig. 3a). This holds for random and optimized designs alike, across total EQE from 0.12 to 0.56. Every weight $w$ returns a design in the same high-efficiency cluster. For this design space, maximizing the band and maximizing total extraction are the same goal.

The angular composition of efficient designs is correspondingly narrow. We take the median selectivity of the twenty designs with the highest total EQE as the *natural composition*: what a design acquires when only efficiency is pushed. It is 0.094 / 0.278 / 0.361 / 0.236 for the four bands. The 10–90% spread is a few thousandths in each band. These twenty designs were located by four different single-band objectives, so the agreement is not an artifact of a shared objective. Compared with the Lambertian partition, the platform depletes the 0–20° band and enriches the 40–60° band.

Band-dedicated optimization does move the composition, and we measure by how much (Fig. 3b). Selectivity rises by +27%, +8%, +1%, and +20% relative to the natural composition, but total EQE falls by 3%, 2%, 1%, and 11%. The two effects largely cancel. Referred to the natural composition at the best total EQE, the net band-power gains are 1.22, 1.05, 1.00, and 1.07 (decomposition in Note S4, Fig. S3). Two points follow. First, the 40–60° band, the natural target of a directional film, gains nothing: the most efficient design already maximizes it. Second, the largest gain is at 0–20°, the band the platform under-serves. Shape freedom restores a depleted direction better than it creates a new one.

The composition is not strictly constant, and we report the deviation. Over the full sampled population, the correlation $R$ between total EQE and $S_j$ is +0.60, +0.55, −0.12, and −0.57 for the four bands (Fig. 3c; per-band panels in Fig. S4). As efficiency rises, the distribution tilts toward lower polar angles. The tilt is robust to sampling: a stratified subset of 20 designs, re-evaluated at 20 times the ray count, gives +0.61, +0.67, +0.04, and −0.71, with the same outer-band signs and a near-zero 40–60° value (Note S5). It is also not a patch artifact, because band selectivity changes by at most 0.24 pp when the patch is enlarged from 25 to 100 mm (Note S3).

The tilt does not reorder the bands among efficient designs. In all 520 designs with total EQE ≥ 0.40, the 40–60° band carries the largest share. Among the 421 designs with EQE ≥ 0.50, 416 follow the order 40–60° > 20–40° > 60–80° > 0–20°. Reordering appears only in the inefficient part of the space. In 8 of the 606 designs, all with EQE ≤ 0.354, the 60–80° band is largest. This is why the 20–40° and 60–80° fits cross at low efficiency in Fig. 3c.

Passive optics can trade position for angle, so étendue conservation alone does not fix the angular distribution, and we claim no absolute invariance. The data support a narrower, practical statement. **Within the explored manufacturable class, angular redistribution is limited to a weak tilt of a few percentage points in selectivity, and its net value for band power is at most about 20%.**

**Fig. 3 | Total efficiency fixes the angular composition.** (a) Total EQE versus 40–60° band EQE for random feasible designs (grey) and optimizer-visited designs (green) of the weighted-sum campaign. All points collapse onto one near-linear locus (black line, $R^2 = 0.968$, $n = 606$). (b) Selectivity of each band-dedicated optimum (red) against the natural composition of the twenty most efficient designs (grey, 10–90% range) and the Lambertian partition (dashed), labelled with the net band-power gain. (c) Selectivity versus total EQE for all four bands, with least-squares fits and the correlation $R$ in the legend. Among designs with total EQE ≥ 0.40 the 40–60° band is always the largest; the 20–40° and 60–80° fits cross only in the inefficient region. *(`fig3_efficiency_composition.png`)*

### 3.3 Generality across external-film families

We next ask whether the result depends on the array construction. We repeat the analysis for the inverted (concave) and the randomly assembled families (Fig. 4), using the same statistics, each computed from the family's own evaluation log. Three signatures of saturation recur in all three families:

- **Collapse.** The (total EQE, 40–60° EQE) points again collapse onto a line. $R^2$ is 0.94 for the convex per-band campaign, 0.94 for the concave family, and 0.84 for the random family (Fig. 4, left column).
- **Narrow composition.** The concave family's best total EQE is 0.517, 6% below the convex 0.548, as expected for a recessed rather than protruding surface. Its natural composition is 0.113 / 0.303 / 0.340 / 0.215. The random family reaches 0.522 ± 0.001 over three disorder realizations, a coefficient of variation of 0.2%. Only the assembly statistics matter, not the particular realization. Its natural composition is 0.111 / 0.299 / 0.347 / 0.213. Band-dedicated optimization of the concave family gives a largest net band-power gain of 1.19, in the 60–80° band, close to the convex family's 1.22.
- **Same direction of tilt.** Selectivity again rises with efficiency at 0–40° and falls at 60–80° (Fig. 4, right column). The 40–60° correlation is weak in the two periodic families and negative in the random one. In every family the 40–60° band remains the largest in at least 99% of sampled designs.

The families also bound how far the composition can be moved at all. Every high-efficiency composition we obtained lies within 0.09–0.11 (0–20°), 0.28–0.30 (20–40°), 0.34–0.36 (40–60°), and 0.21–0.24 (60–80°). This spans three constructions, four single-band objectives, and two averaging conventions (top-twenty median and population mean). The change from convex to concave lenslets moves the composition to the edge of this window. The shift is at most 0.02 in any band, about the same as switching the averaging convention within one family (Note S4, Table S6). Positional disorder and per-lenslet shape disorder leave the composition inside the window. The composition is therefore set by the platform—extended source, thick substrate, coextensive external film—not by any design variable within it.

The randomly assembled MLA should not be confused with a volumetric scattering layer. Each lenslet remains a refractive surface with deterministic local normals and a single index step. A scattering layer instead redirects light by multiple volumetric scattering. Its behaviour depends on the particle–matrix index contrast and loading, which are not parameters of this study. Our conclusions apply to refractive external films and make no comparison with scattering layers.

**Fig. 4 | Generality across external-film MLA families.** Rows: convex freeform (reference, per-band campaign), inverted (concave), and randomly assembled (pseudo-random supercell). Every statistic is computed from each family's own evaluation log. (a1–c1) Total EQE versus 40–60° band EQE with a least-squares line; the near-linear collapse recurs in all three families. (a2–c2) Natural composition (top-twenty median, colored) against the convex reference (grey) and the Lambertian partition (dashed outline). The open circle is the population mean of the same family; the spread between statistics is comparable to the spread between families. (a3–c3) Correlation $R$ between total EQE and $S_j$, per band, against the convex reference. *(`fig4_families.png`)*

---

## 4. Discussion

### 4.1 Why the design space saturates

Three considerations explain the saturation.

**No area leverage.** Angular compression requires an output aperture larger than the source. A small source under a large lens can narrow its emission cone because the lens supplies extra output area. In a tiled MLA on a large-area OLED, each lenslet's aperture and the source area it serves grow together. The output/source area ratio is pinned near unity, whatever the lenslet size, so enlarging the lenslets does not restore the point-source collimation gain. The same relation explains why an index-matched hemispherical macroextractor works well on a small emitter, and why it cannot be tiled over a large area.

**Lateral mixing in the substrate.** A ray at the critical angle moves sideways by $d_{\mathrm{sub}}\tan\theta_c = 1.16$ mm in one pass through the substrate, and by 2.32 mm per recycling round trip. That is about 100 and 200 lenslet radii. Each lenslet therefore receives light launched from a substrate region two orders of magnitude wider than itself. The light reaching the lens surface has already been spatially averaged before the lens shape can act on it. The patch-size series measures this mixing directly (Note S3). One fixed design gives total EQE 0.517, 0.543, 0.551, and 0.564 at patch sizes of 15, 25, 35, and 100 mm. Total EQE keeps rising slowly up to 100 mm, so light travels laterally over many round trips before escaping. The 25 mm values are therefore lower bounds. The angular composition changes by at most 0.24 pp over the same range, so every ratio, correlation, and class comparison is unaffected. Ref. [16] shows the same loss channel experimentally: an MLA film only 2 mm in radius on a 1 mm emitter gave no EQE gain (35.4% against 35.6% bare), because light escaped past the film edge before it could be extracted.

**The radiance limit of a passive film.** For a given source radiance and exit area, the power a passive external layer can deliver into a chosen solid angle is bounded [10]. We use this bound as a benchmark, not as an impossibility proof. Once a hemisphere operates near the bound, a more complex profile can only reshuffle small amounts of power; it cannot create radiance.

### 4.2 Design routes beyond lens shape

Saturation of the lens shape does not end light extraction. It changes which quantity is worth varying. Table 1 maps each target figure of merit to the lever that remains productive once the hemispherical benchmark has confirmed saturation.

**Table 1 | Design-route map after refractive-MLA saturation.**

| Target figure of merit | Productive lever | Evidence | Practical constraint |
|---|---|---|---|
| Total EQE from substrate modes | Hemispherical MLA at its optimum; freeform shape adds ≤ 1.3% | Sec. 3.1, Fig. 2 | Molding tolerance only |
| Total EQE beyond the external film | Source/cavity engineering: dipole orientation, cavity design, internal extraction [11,12] | Cavity co-optimized in every campaign | Modifies the emitting stack |
| Power in a chosen polar band | Follows total EQE; band-dedicated shape gives a net gain ≤ 1.22 | Sec. 3.2, Fig. 3 | Selectivity is paid for in total EQE |
| Angular compression (collimation) | Aperture expansion: output/source area ratio > 1 (macroextractor, pixel-level optics) | Composition window fixed in coextensive films, Fig. 4 | Loses planar tiling |
| Selectivity beyond the single-pass limit | Angle-selective light recycling: filter plus reflective electrode [8,13,14] | Ideal filter 62.1% vs 29.1% planar at 10% loss; DBR 48.2% → 32.5% over 100 nm (Note S7) | Round-trip loss and source bandwidth |

**Source and cavity engineering** changes the substrate-side source distribution itself, through cavity thickness, dipole orientation, or resonant structures. It is the established lever for OLED angular emission [11,12] and the first to examine once the MLA has saturated.

**Aperture expansion** is needed when collimation is the main requirement. Everything we varied within the coextensive geometry moved the natural composition by at most 0.02 in any band: shape freedom, convex versus concave lenslets, and positional and shape disorder. Angular compression therefore requires changing the output/source area ratio, which means leaving the coextensive geometry. Small emitters and pixel-level optics allow this, at the cost of planar scalability.

**Angle-selective light recycling** can provide the selectivity that a non-selective refractive film cannot. Here a filter reflects out-of-band light back toward the reflective electrode for another attempt. In a single pass, non-selective layers deliver 33.7% of the escaping light into the 40–60° band. A Markov recycling model with an ideal angular filter and 10% round-trip loss raises the delivery into that band to 62.1% of generated power, against 29.1% for a planar non-selective reference. A realizable eight-pair dielectric multilayer gives 48.2% for a monochromatic source but only 32.5% over a 100 nm bandwidth (Note S7, Fig. S7). Round-trip loss and source bandwidth are therefore the practical limits of this route. This is a reference calculation that shows why selective recycling is the natural next lever, not a device proposal.

### 4.3 Scope and limitations

The saturation statement covers three manufacturable external-film families at substrate-matched lens index. It does not cover internal arrays, scattering layers, or high-index films. A separate, differently constrained exploration with three-dimensional asymmetric lenslets is consistent with the result. That exploration used an Ag electrode, an anisotropic emitter cell, and 52 shape variables. Its asymmetric shapes moved the far-field centroid but did not reliably raise the absolute power in a target $(\theta,\phi)$ window above the hemisphere. Directional emission is possible with other elements, such as resonant structures, metasurfaces, and asymmetric prisms [8,9]; the result here concerns refractive lenslets only. Because its stack and constraints differ, we report it only qualitatively (Note S8). Absolute EQE values at the 25 mm patch are lower bounds (Section 4.1), but ratios and comparisons are not affected. Finally, the simulated surfaces are ideal: there is no form error, misalignment, or lens-to-lens variation. Fabrication imperfections can only reduce whatever advantage a freeform profile holds over the hemisphere, so the near-null result is conservative.

---

## 5. Conclusion

We tested the practical value of lens-shape freedom for manufacturable, tiled refractive MLAs on an extended OLED. Against an equally optimized hemisphere, the freeform optimum gains at most 7% in any polar band and 1.3% in total EQE. A restart control confirms that this is a property of the design space and not of the search. All designs collapse onto a single line relating band power to total EQE. Among efficient designs the angular composition is fixed to within a few percentage points, whatever objective located them. Optimizing for one band raises its selectivity by up to 27%, but costs up to 11% of total EQE, leaving a net band-power gain of at most 1.22.

Convex, concave, and randomly assembled films all keep the composition within one narrow window. That window is therefore a property of the platform: an extended source, a thick substrate, and a coextensive external film. The finding is not a claim that the hemisphere is optimal in every optical system. It provides a decision criterion. Once the hemispherical benchmark confirms saturation, further improvement should come from the source and cavity, from aperture expansion, or from angle-selective light recycling, not from more lens-shape freedom. The same reasoning applies to other extended thin-film emitters, such as perovskite and quantum-dot LEDs.

---

## Funding

[TBD — 과제/기관 정보 입력 필요]

## Acknowledgments

[TBD]

## Disclosures

The authors declare no conflicts of interest.

## Data availability

The result archives underlying all figures and tables are available from the corresponding author upon reasonable request. These comprise the evaluation logs and optima of every optimization campaign and the recycling-model outputs, listed in Supplement 1. The archives come with the scripts that generated them. `make_figures.py` and `make_supp_figures_and_data.py` reproduce every figure, together with an Excel workbook of its plotted data. Data for the separately constrained asymmetric exploration are also available on request. Running the CPS source model and LightTools project files requires a LightTools license.

---

## References

1. Brütting, W.; Frischeisen, J.; Schmidt, T. D.; Scholz, B. J.; Mayr, C. **Device Efficiency of Organic Light-Emitting Diodes: Progress by Improved Light Outcoupling.** *Phys. Status Solidi A* **2013**, *210*, 44–65.

2. Yablonovitch, E. **Statistical Ray Optics.** *J. Opt. Soc. Am.* **1982**, *72*, 899–907.

3. Möller, S.; Forrest, S. R. **Improved Light Out-Coupling in Organic Light Emitting Diodes Employing Ordered Microlens Arrays.** *J. Appl. Phys.* **2002**, *91*, 3324–3327.

4. Wrzesniewski, E.; et al. **Enhancing Light Extraction in Top-Emitting Organic Light-Emitting Devices Using Molded Transparent Polymer Microlens Arrays.** *Small* **2012**, *8*, 2647–2651. https://doi.org/10.1002/smll.201102662.

5. Qu, Y.; Kim, J.; Coburn, C.; Forrest, S. R. **Efficient, Nonintrusive Outcoupling in Organic Light Emitting Devices Using Embedded Microlens Arrays.** *ACS Photonics* **2018**, *5*, 2453–2458. https://doi.org/10.1021/acsphotonics.8b00255.

6. Kim, S.; Shin, J. M.; Lee, J.; Park, C.; Lee, S.; Park, J.; Seo, D.; Park, S.; Park, C. Y.; Jang, M. S. **Inverse Design of Organic Light-Emitting Diode Structure Based on Deep Neural Networks.** *Nanophotonics* **2021**, *10*, 4533–4541. https://doi.org/10.1515/nanoph-2021-0434.

7. Ni, Y.; Feng, D.; Ma, D. **Design of Freeform Microlens Arrays with Prescribed Luminance Distributions for MicroLED Optical Packaging.** *Appl. Opt.* **2025**, *64*, 7875–7884. https://opg.optica.org/ao/abstract.cfm?uri=ao-64-27-7875.

8. Buhl, M.; et al. **Resonance-Based Directional Light Emission from Organic Light-Emitting Diodes.** *Adv. Photonics Res.* **2023**, *4*, 2200143. https://doi.org/10.1002/adpr.202200143.

9. Abdelkhalik, M. S.; Garcia-Santiago, X.; van Raaij, T.-J.; López, T.; Berghuis, A. M.; de Jong, L. M. A.; Gómez Rivas, J. **Enhanced and Directional Electroluminescence from MicroLEDs Using Metallic or Dielectric Metasurfaces.** *Commun. Eng.* **2025**, *4*, 63. https://doi.org/10.1038/s44172-025-00401-w.

10. Winston, R.; Jiang, L.; Ricketts, M. **Nonimaging Optics: A Tutorial.** *Adv. Opt. Photon.* **2018**, *10*, 484–511.

11. Xiang, C.; Koo, W.; So, F.; Sasabe, H.; Kido, J. **A Systematic Study on Efficiency Enhancements in Phosphorescent Green, Red and Blue Microcavity Organic Light-Emitting Devices.** *Light: Sci. Appl.* **2013**, *2*, e74. https://doi.org/10.1038/lsa.2013.30.

12. Song, J.; et al. **Lensfree OLEDs with over 50% External Quantum Efficiency via External Scattering and Horizontally Oriented Emitters.** *Nat. Commun.* **2018**, *9*, 3207. https://doi.org/10.1038/s41467-018-05671-x.

13. Liao, P.-H.; Lee, W.-K.; Lee, C.-C.; Huang, C.-W.; Wen, S.-W.; Chen, Y.-T.; Chen, C.-C.; Lin, W.-Y.; Kwak, B. L.; Visser, R. J.; Wu, C.-C. **Using Angle-Selective Optical Film to Enhance the Light Extraction of a Thin-Film Encapsulated 3D Reflective Pixel for OLED Displays.** *Opt. Express* **2022**, *30*, 46435–46449. https://doi.org/10.1364/OE.477797.

14. Kim, H.-J.; et al. **High Efficient OLED Displays Prepared with the Air-Gapped Bridges on Quantum Dot Patterns for Optical Recycling.** *Sci. Rep.* **2017**, *7*, 43063. https://doi.org/10.1038/srep43063.

15. Kim, J.; Kim, E.; Park, J.; Song, J.; Kim, S.; Moon, H.; Yoo, S. **Toward Near-Foldable Surface Light Sources with Ultimate Efficiency: Ultrathin Substrates Embedded with Micron-Scale Inverted Lens Arrays.** *ACS Photonics* **2023**, *10*, 1775–1782. https://doi.org/10.1021/acsphotonics.3c00017.

16. Kim, M.; Kim, J.; Yoo, S. **Near-Planar Light Outcoupling Structures with Finite Lateral Dimensions for Ultra-Efficient and Optical Crosstalk-Free OLED Displays.** *Nat. Commun.* **2025**, *16*, 11606. https://doi.org/10.1038/s41467-025-66538-6.
