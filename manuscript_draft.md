# Practical Saturation of Freeform Microlens Arrays on Extended OLED Emitters and Design Routes beyond Lens Shape

**Authors** — [TBD]

---

## Abstract

Freeform microlens arrays (MLAs) and numerical inverse design are widely expected to give simultaneous control over the light-extraction efficiency and the angular emission profile of organic light-emitting diodes (OLEDs). This expectation, however, largely originates in illumination optics built around point-like sources far smaller than the lens aperture, and it does not automatically transfer to an MLA tiled coextensively over a uniform, extended emitter. Here we combine a dipole microcavity source model based on the classical power spectrum (CPS) formalism with three-dimensional ray tracing to compare hemispherical reference lenses and axisymmetric freeform lenses under identical geometric, material, and process constraints, supplemented by a separately constrained asymmetric, off-axis freeform exploration. Total external quantum efficiency (EQE) and the 40–60° polar band were optimized independently; the remaining polar bands (0–20°, 20–40°, 60–80°) are covered by a weighted-sum sweep between the two objectives and by a three-band optimization [3-band optimization results: TO BE INSERTED]. Across 150 randomly sampled feasible designs and all optimized designs, spanning total EQE from 0.11 to 0.56, the (total EQE, band EQE) points collapse onto a single near-linear frontier: efficiency and band power rise together, and no efficiency-versus-directionality trade-off is observed. At the frontier, the 40–60° selectivity is approximately 0.355, only marginally above the substrate-mode-free Lambertian value of 0.337, and no explored design achieved meaningful angular steering.

We explain this practical saturation by three factors. First, in a coextensively tiled array the lens aperture and the effective source patch grow in the same proportion, so the area leverage that underlies point-source collimation does not arise. Second, lateral propagation in the substrate mixes light between neighboring lenslets and erases the spatial information a local surface shape could exploit. Third, a passive external layer operating on a given source radiance and exit area faces a practical ceiling on the power it can deliver into a chosen angular channel. These results do not render freeform MLAs meaningless; they provide a quantitative criterion for deciding when further investment in shape complexity on extended OLEDs should stop. Finally, we organize the three routes that remain productive after MLA saturation—source/cavity engineering, aperture expansion, and angular-selective recycling—in the same physical language, and present a design-route map for choosing the next lever. In an idealized recycling model with 10% round-trip loss, an angular filter delivers 62.1% of the generated power into the 40–60° band, whereas the planar non-selective reference stays at 29.1%; a realizable eight-pair dielectric multilayer reaches 48.2% for a monochromatic source but degrades to 32.5% over a 100 nm bandwidth, identifying source bandwidth and loss as the practical constraints of this route.

---

## 1. Introduction

External light extraction remains a central limit on the efficiency and high-luminance lifetime of OLEDs. Waveguided modes in the organic layers, loss channels associated with the metal electrode, and total internal reflection at the substrate/air interface together trap a large fraction of the photons generated in a planar OLED [1,2]. External MLAs can extract the substrate modes without perturbing the electrical structure of the device, and have therefore been an attractive solution for large-area and flexible OLEDs in particular [3–5]. Hemispherical or near-hemispherical MLAs already reach high efficiencies, and very high EQE has been reported for OLEDs with embedded hemispherical MLAs [4,5].

Nevertheless, the expectation attached to more complex shapes persists. Freeform illumination optics and inverse design are powerful tools for producing prescribed intensity distributions, and have recently been applied to the optical packaging of micro-LEDs and OLEDs [6,7]. Intuitively, tuning the asymmetry, curvature distribution, height, pitch, or the microcavity condition of the lenses simultaneously should offer finer control over total extraction and over emission into specific directions than a hemispherical MLA does.

Two observations sharpen this intuition into a testable question. On one hand, freeform design methods were largely developed for sources of near-zero étendue, whereas a thin-film emitter is an extended source whose light reaches the array only after propagation, and partial recycling, inside the substrate. On the other hand, high total extraction demonstrably does not require lenses at all: Song et al. reported lens-free OLEDs exceeding 50% EQE using an external random scattering layer combined with horizontally oriented emitters [13]. If total efficiency is reachable without shaped surfaces, the question that remains—and the one this paper addresses—is whether freeform lens shape buys angular control *on top of* efficiency on an extended emitter.

This question matters in practice. Freeform design and precision molding raise fabrication cost and metrology burden; if the achievable gain barely exceeds a hemispherical reference, the next design step should not be a more complex surface but a change of source, aperture, or recycling path. Our aim is therefore not "can a better freeform lens be found?" but: **how much is shape freedom actually worth for a tiled refractive MLA on an extended OLED?**

We compare, under one OLED source model and one set of manufacturable geometric constraints, (i) a hemispherical MLA and (ii) an axisymmetric freeform MLA, and we additionally examine (iii) an asymmetric three-dimensional freeform exploration performed as a separate, differently constrained study (Section 2.4). We use total EQE and polar-band EQE as independent objectives, employ surrogate-based multi-start optimization with high-precision re-evaluation to separate optimizer dependence from physics, and connect the outcome to a radiance/étendue picture and to lateral mixing in the substrate. We emphasize the epistemic status of the result throughout: it is not an impossibility theorem but a reproducible numerical observation—**a practical saturation of the explored freeform design space**—together with a quantitative account of why it occurs and of which optical levers remain productive beyond it.

---

## 2. Results and Discussion

### 2.1 A fairly aligned benchmark: isolating the effect of shape complexity itself

Figure 1 shows the optical system under comparison. Dipole emission in the OLED microcavity is computed with the CPS formalism to obtain the substrate-side angular–spectral distribution $I_{\mathrm{sub}}(\theta,\lambda)$, which serves as the source for three-dimensional ray tracing. The emitting region is a disc of radius $r_{\mathrm{OLED}} = 1$ mm on a substrate of thickness $d_{\mathrm{sub}} = 1.295$ mm; the exit face carries a two-dimensional hexagonal array of lenslets of approximately 10 μm radius extending over 15 × 15 mm, so that laterally spreading substrate light remains within the array. The source is thus extended by two orders of magnitude relative to an individual lenslet—the regime of interest for tiled outcoupling films. All lenses satisfy the same substrate index, lens material, pitch, fill factor, maximum height, and maximum draft-angle constraints. Cases with fixed cavity thickness and with joint lens/cavity optimization are reported separately.

This alignment is essential. If the freeform class were allowed a wider variable range, a greater lens height, a larger aperture, or a more favorable cavity, shape effects would be confounded with other effects. Here the hemispherical MLA is itself optimized over the same outer radius and the same manufacturable height range, and the freeform class is subject to constraints that are identical—never more permissive. Subsequent performance differences can therefore be read, as far as possible, as the pure effect of shape freedom.

**Fig. 1 | Platform and benchmark design.** (a) OLED–substrate–MLA architecture: CPS microcavity dipole source ($r_{\mathrm{OLED}} = 1$ mm), glass substrate ($d_{\mathrm{sub}} = 1.295$ mm), hexagonal lenslet array (~10 μm lens radius, 15 × 15 mm extent). (b) Lens classes compared under identical constraints: hemispherical reference and axisymmetric freeform (13 spline variables). (c) Definition of the polar bands and of the source/aperture geometry. (d) Optimization workflow: surrogate-based global search, pattern-search polishing, and high-ray re-evaluation.

### 2.2 Changing the objective does not carry designs far beyond the hemispherical reference

Figure 2 summarizes the best attained performance when total EQE and individual polar bands are maximized. For each objective $j$ we define the relative gain

$$
G_j=\frac{\max\left[\mathrm{EQE}_j\mid\mathrm{freeform}\right]}
{\max\left[\mathrm{EQE}_j\mid\mathrm{hemisphere}\right]},
$$

where each maximum is defined after independent initializations, multiple starts, and high-ray re-evaluation, as the mean of the re-evaluated repeats. Across all explored polar objectives, $G_j$ remained within **[FINAL VALUE: range of $G_j$ with uncertainty]**. Some freeform designs showed local gains in a particular angular bin, but once the accompanying reduction of total EQE, the displacement of power into other bins, and the Monte-Carlo uncertainty are accounted for, these gains did not form an independent performance axis that would displace the hemispherical MLA.

This does not mean that the optimum "must be exactly hemispherical." The best shapes differed from one another in slope and curvature distribution, and the low-efficiency region of the design space contains great shape diversity. The important point is that near the top of the performance range this diversity did not convert into meaningful additional efficiency or band power. The hemisphere therefore operates not as a mere convenience baseline but as a **practical near-optimum** of this structural class.

**Fig. 2 | Achievable region and Pareto collapse.** (a) Best hemispherical and freeform values for each objective, with re-evaluation error bars, and the relative gains $G_j$. (b) Total EQE versus 40–60° band EQE for 150 random feasible designs (grey) and all weighted-sum optima $w \in \{0, 0.25, 0.5, 0.75, 1\}$ (colored): all points collapse onto a near-linear frontier spanning total EQE 0.11–0.56, with no observed trade-off between efficiency and band power. (c) Best freeform profiles compared with the hemispherical reference.

### 2.3 Pareto analysis: polar shaping is not an independent design axis

To probe total extraction and the 40–60° band power jointly, we optimized the weighted objective

$$
J_w=w\,\hat{\eta}_{\mathrm{ext}}+(1-w)\,\hat{P}_{40\text{–}60}
$$

for $w = 0, 0.25, 0.5, 0.75, 1$ using surrogate-based global optimization with pattern-search polishing, and separately collected $N = 150$ random feasible freeform designs to populate the interior of the achievable region without optimizer sampling bias.

The outcome is a collapse rather than a front. All (total EQE, band EQE) points—random and optimized alike—fall on a single near-linear locus spanning total EQE from 0.11 to 0.56 (Fig. 2b). Every weight returns a design in the same high-efficiency cluster: optimizing for the band and optimizing for total extraction are not competing objectives within this design space.

The angular composition is best expressed through the band selectivity $S_j = \mathrm{EQE}_{\mathrm{band},j}/\mathrm{EQE}_{\mathrm{total}}$, compared against the substrate-mode-free Lambertian reference $S_j^{\mathrm{Lam}} = \sin^2\theta_{\mathrm{hi}} - \sin^2\theta_{\mathrm{lo}}$: 0.117 (0–20°), 0.296 (20–40°), 0.337 (40–60°), and 0.220 (60–80°). At the best designs the measured 40–60° selectivity is approximately 0.355, only slightly above the Lambertian value.

The selectivity is, however, not strictly constant, and we report the deviation as a result rather than suppressing it. Across all designs, the Pearson correlation $R$ between $\mathrm{EQE}_{\mathrm{total}}$ and $S_j$ is approximately $+0.6$ for 0–20°, $+0.7$ for 20–40°, $+0.05$ for 40–60°, and $-0.7$ for 60–80° (Fig. 3). As total extraction improves, the angular distribution drifts systematically toward lower polar angles; yet the ordering of the bands never inverts—the 40–60° band carries the largest share in every explored design. This drift survives a dedicated convergence check with a twenty-fold ray count (200,000 rays), three independent repeats, and broadband 450–750 nm evaluation, which reproduce the same correlations (Methods 4.3; Supplementary Table S2 and Fig. S1). It is therefore a real, if small, physical trend, not Monte-Carlo noise or a narrowband artifact.

A conceptual distinction matters here. Passive optics can redistribute power between position and angle, so étendue conservation by itself does not fix the angular distribution, and we make no claim of absolute invariance. What the data support is a more limited and more practical statement: **within the explored manufacturable freeform class, the magnitude of angular redistribution is small—a weak, monotonic tilt of a few percentage points in selectivity—so the polar-shaping freedom available to a designer is effectively saturated.** This is presented not as a universal theorem but as a reproducible benchmark with quantified uncertainty. It is also consistent with the lens-free result of Song et al. [13]: if a random scattering layer already reaches the efficiency frontier without shaped surfaces, shape can hardly be the carrier of a large independent angular degree of freedom on an extended emitter.

**Fig. 3 | Selectivity map across the design space.** Four panels showing the band selectivity $S_j$ versus total EQE for all 150 random and all optimized designs, for the bands 0–20°, 20–40°, 40–60°, and 60–80°. Horizontal lines mark the substrate-mode-free Lambertian references (0.117, 0.296, 0.337, 0.220). Annotated correlation coefficients ($R \approx +0.6$, $+0.7$, $+0.05$, $-0.7$) quantify the systematic drift toward lower polar angles with rising efficiency; the band ordering does not invert in a single explored design. Convergence-check values at twenty-fold ray count, three repeats, and 450–750 nm broadband evaluation are overlaid (Supplementary Fig. S1).

### 2.4 Asymmetric and off-axis freeform lenses as a stress test

To check whether symmetry itself limits the outcome, we optimized three-dimensional asymmetric freeform lenslets directly against a restricted $(\theta,\phi)$ window, with the objective

$$
\mathrm{EQE}_{\mathrm{win}}=\int_{\theta_1}^{\theta_2}\!\!\int_{\phi_1}^{\phi_2}
I_{\mathrm{air}}(\theta,\phi)\,\sin\theta\,d\phi\,d\theta .
$$

This test is important: if an asymmetric freeform lens delivered substantially higher window power or contrast at comparable total EQE, the saturation of the preceding sections would merely reflect the limits of an axisymmetric parameterization.

We state the scope of this test explicitly. The asymmetric exploration was a **separate, differently constrained study**: it used a different OLED stack (Ag rather than ITO electrode), an anisotropic emitter cell, and 52 shape variables rather than the 13 of the controlled benchmark. It is therefore *not* part of the same-constraints comparison of Section 2.1 and is used here only as supporting evidence that an appreciably richer parameterization, under its own favorable conditions, likewise produced no stable angular steering.

In our calculations the asymmetric shapes could shift the centroid and fine structure of the far field, but did not stably achieve an **absolute window power** appreciably above the hemispherical MLA **[FINAL VALUE: window-power comparison with uncertainty]**. This is not a claim that periodic arrays are incapable of steering under all circumstances: asymmetric prisms, diffractive elements, metasurfaces, or sufficient aperture expansion can produce directional light distributions [8,9]. Within the scope of this work—coextensive refractive MLAs on a uniform extended OLED source—such redistribution simply did not convert into a useful output channel beyond the hemispherical reference.

### 2.5 Physical origin of the saturation: area, mixing, and the radiance envelope

Three mutually complementary considerations account for the observed saturation.

First, angular compression generally requires expansion of the output aperture. A small source under a large macro-lens can reduce the emitted solid angle because it uses an output area larger than the source area. In a tiled MLA covering a large-area OLED, by contrast, each lenslet's aperture and the source patch it serves grow in the same proportion: what matters is the ratio of output to source area, not the absolute lenslet size, and for a coextensive array that ratio is pinned near unity. Enlarging the lenslets therefore does not recover the point-source collimation gain. The same relation explains, from the other side, why an index-matched hemispherical macroextractor on a small emitter works so well—and why that configuration does not tile over a large area.

Second, substrate thickness and source size set the degree of lateral mixing. Light that has propagated through the substrate does not originate solely beneath the lenslet it strikes; it is blended with light from neighboring cells, each lenslet drawing from an effective collection radius far exceeding its own. Sweeping source radius, substrate thickness, and pitch shows the attainable steering falling as the angular extent of the source grows **[FINAL VALUE: sweep quantification]** (Supplementary Fig. S2). The input phase space presented to the freeform surface is thus already averaged before shape can act on it.

Third, for a given source radiance and exit area there is a radiance/étendue envelope on the power a passive external layer can deliver into a chosen solid angle [10]. We use this envelope not as a global impossibility statement but as an interpretive benchmark for the numerical frontier: once a flat interface or a hemisphere operates close to the envelope, a more complex profile can only re-shuffle small amounts of power rather than create new radiance.

### 2.6 Design routes after MLA saturation

That shape optimization yields diminishing returns does not mean light extraction is finished; it means the quantity to change is no longer the lens profile. Figure 4 organizes the three subsequent routes.

**Source/cavity engineering.** Microcavity thickness, dipole orientation, reflective electrodes, and resonant structures alter the substrate-side source distribution itself. This is already the established lever for OLED angular-emission control [12,13]—Song et al.'s combination of horizontal emitter orientation with a non-selective external layer [13] is precisely a source-side success—and we position it as the first lever to examine once the MLA has saturated.

**Aperture expansion.** Where a small emitter or a pixel-level architecture genuinely permits an output/source area ratio well above unity, a macroextractor or non-coplanar optical element enables angular compression. This route lacks the planar scalability of a tiled MLA but is appropriate when collimation is the dominant requirement.

**Angular-selective recycling.** An angular filter that selectively reflects light outside the target band, returning it for further attempts, can create the selectivity that a non-selective refractive MLA does not provide. In our single-pass calculations the 40–60° band selectivity of non-selective layers saturates at 33.7% (analytic) and 33.8% (numeric), independent of the asymmetry parameter of an adjacent scattering layer (Methods 4.4)—recycling with angular selection is what breaks this wall. In the idealized filter model with 10% round-trip loss, delivery into the 40–60° band reaches 62.1% of the generated power, whereas the planar non-selective (class-A) reference stays at 29.1%. A transfer-matrix calculation of a realizable eight-pair dielectric multilayer yields 48.2% for a monochromatic source, degrading to 32.5% over a 100 nm bandwidth, showing that source bandwidth and round-trip loss are the practical constraints of this route. Photon recycling as such has a firm quantitative footing—in photovoltaics through the reciprocity relation between quantum efficiency and electroluminescence [11]—and angle-selective films and recycling structures have been demonstrated in OLED contexts [8,14,15]. Our calculation is not a new device proposal; it is a reference computation that isolates *why* selective recycling is the natural next lever after refractive-MLA saturation.

**Fig. 4 | The recycling route and the design-route map.** (a) Markov recycling model: power delivered into the 40–60° band versus round-trip loss $a$ for the ideal angular filter (62.1% at $a = 0.1$) and the planar non-selective reference (29.1%). (b) Transfer-matrix result for an eight-pair dielectric multilayer: band delivery versus source bandwidth, from 48.2% (monochromatic) to 32.5% (100 nm). (c) Single-pass selectivity wall: 33.7% analytic / 33.8% numeric, independent of scattering asymmetry $g$. (d) Design-route map: recommended lever (source/cavity engineering, aperture expansion, angular-selective recycling) as a function of the target figure of merit, entered once the hemispherical benchmark confirms saturation.

---

## 3. Conclusion

We systematically tested the practical value of freeform shape freedom for manufacturable tiled refractive MLAs on extended OLEDs. Total EQE and the 40–60° polar band were optimized independently, a weighted-sum sweep and 150 random designs mapped the achievable region, and a separately constrained asymmetric off-axis exploration served as a stress test; the remaining bands are covered by a three-band optimization [3-band optimization results: TO BE INSERTED]. Freeform MLAs showed no stable performance gain substantially beyond a hemispherical MLA under identical constraints: all designs collapsed onto a near-linear efficiency–band-power frontier with selectivity close to the Lambertian partition, modulated only by a small, systematic, convergence-verified drift toward lower polar angles at higher efficiency. This is not a universal claim that the hemisphere is optimal in all optical systems. It is a quantitative finding that when coextensive tiling, an extended source, substrate-mediated lateral mixing, and a fixed radiance/area budget hold simultaneously, the additional freedom purchasable by shape complexity is small.

The study accordingly does not conclude that freeform MLAs should be abandoned; it supplies a criterion that accelerates the design decision. Once saturation is confirmed against the hemispherical benchmark, the next improvement should be sought not in more lens-shape degrees of freedom but in source/cavity engineering, aperture expansion, or angular-selective recycling. This design-route map applies beyond OLEDs, offering a practical starting point for PeLED, QLED, and micro-LED optical packaging built on extended thin-film sources.

---

## 4. Methods

### 4.1 Trans-scale source model

Dipole emission in the OLED stack was computed with the CPS formalism. The wavelength- and substrate-angle-resolved intensity $I_{\mathrm{sub}}(\theta,\lambda)$ served as the spectral/angular source of the macroscopic ray tracing. The emitting disc has radius 1 mm; the substrate thickness is 1.295 mm. CPS inputs—material dispersion, dipole orientation, internal radiative efficiency, emission spectrum, layer thicknesses, and wavelength sampling—are collected in Table S1. The framework has previously been compared against measured angular profiles and EQE in our group's earlier MLA–OLED work [5].

### 4.2 Lens classes and manufacturing constraints

The hemispherical reference and the axisymmetric freeform class were compared under one set of constraints: identical lens material, pitch, fill factor, nominal height, boundary continuity, and maximum draft angle. Lenslets of approximately 10 μm radius are tiled hexagonally over a 15 × 15 mm area on the substrate exit face. The axisymmetric freeform surface is parameterized by 13 design variables (spline control points, with the radial control points $x_2 \ldots x_6$ constrained monotonic) and implemented as a LightTools native freeform entity. Candidates exhibiting geometric self-intersection, negative thickness, or unmanufacturable draft angles were rejected before optimization. The asymmetric off-axis study of Section 2.4 used a different stack (Ag electrode), an anisotropic emitter cell, and 52 variables, and is reported as a separate exploration, not as part of the controlled benchmark.

### 4.3 Optimization and validation

Objectives were total EQE, polar-band EQE, or $(\theta,\phi)$-window EQE. Each objective was optimized with MATLAB surrogate-based global optimization (`surrogateopt`) followed by `patternsearch` polishing, with multiple independent starts; low-ray search results were re-evaluated with independent high-ray repeats. The achievable region was populated with $N = 150$ random feasible designs in addition to the weighted-sum optima at $w \in \{0, 0.25, 0.5, 0.75, 1\}$. Optimization budgets, ray numbers, independent-run counts, feasible-sample counts, and standard deviations are reported in Table S2. Freeform superiority was judged from high-precision means against the hemispherical reference, never from a single best run.

**Convergence check for the selectivity drift.** Because the correlations between total EQE and band selectivity are small in magnitude for some bands, we verified that they are not sampling artifacts. The full design set was re-evaluated under three progressively stricter conditions: (i) the baseline ray count, (ii) a twenty-fold ray count (200,000 rays) with three independent repeats, and (iii) broadband evaluation over 450–750 nm. The Pearson correlations were stable across all three: $R(0\text{–}20^\circ) = 0.59 \rightarrow 0.61 \rightarrow 0.64$; $R(20\text{–}40^\circ) = 0.67 \rightarrow 0.67 \rightarrow 0.72$; $R(40\text{–}60^\circ) = 0.07 \rightarrow 0.04 \rightarrow 0.05$; $R(60\text{–}80^\circ) = -0.70 \rightarrow -0.71 \rightarrow -0.76$ (Supplementary Table S2 and Fig. S1). The systematic drift of Section 2.3 is therefore robust against Monte-Carlo noise and narrowband bias.

### 4.4 Angular-selective recycling model

The external layer is represented by its angular transmission/reflection $T(\theta,\lambda)$ and the reflective electrode by a round-trip loss $a$, combined in a Markov recycling model implemented in Python. The ideal angular filter is defined by $T = 1$ inside the target band and $T = 0$ outside. In the single-pass limit the 40–60° band selectivity of non-selective layers is 33.7% analytically and 33.8% numerically, and is unchanged over the scattering asymmetry range $g \in [-0.5, +0.5]$ of a two-stream scattering layer. Realizable filters were computed by transfer-matrix calculation of an alternating high/low-index dielectric multilayer (eight pairs). This model is not a proposal of a combined MLA–DBR device; it is a reference calculation separating the saturation of non-selective refractive MLAs from the behavior of selective recycling.

---

## References

1. Brütting, W.; Frischeisen, J.; Schmidt, T. D.; Scholz, B. J.; Mayr, C. **Device Efficiency of Organic Light-Emitting Diodes: Progress by Improved Light Outcoupling.** *Phys. Status Solidi A* **2013**, *210*, 44–65.

2. Yablonovitch, E. **Statistical Ray Optics.** *J. Opt. Soc. Am.* **1982**, *72*, 899–907.

3. Möller, S.; Forrest, S. R. **Improved Light Out-Coupling in Organic Light Emitting Diodes Employing Ordered Microlens Arrays.** *J. Appl. Phys.* **2002**, *91*, 3324–3327.

4. Wrzesniewski, E.; et al. **Enhancing Light Extraction in Top-Emitting Organic Light-Emitting Devices Using Molded Transparent Polymer Microlens Arrays.** *Small* **2012**, *8*, 2647–2651. https://doi.org/10.1002/smll.201102662.

5. Qu, Y.; et al. **Efficient, Nonintrusive Outcoupling in Organic Light Emitting Devices Using Embedded Microlens Arrays.** *ACS Photonics* **2018**. https://doi.org/10.1021/acsphotonics.8b00255.

6. Kim, Y.; et al. **Inverse Design of Organic Light-Emitting Diode Structure Based on Deep Neural Networks.** *Nanophotonics* **2021**, *10*.

7. **Design of Freeform Microlens Arrays with Prescribed Luminance Distributions for MicroLED Optical Packaging.** *Appl. Opt.* **2025**, *64*, 7875.

8. Buhl, M.; et al. **Resonance-Based Directional Light Emission from Organic Light-Emitting Diodes.** *Adv. Photonics Res.* **2023**. https://doi.org/10.1002/adpr.202200143.

9. **Enhanced and Directional Electroluminescence from MicroLEDs Using Metallic or Dielectric Metasurfaces.** *Commun. Eng.* **2025**. https://doi.org/10.1038/s44172-025-00401-w.

10. Winston, R.; Jiang, L.; Ricketts, M. **Nonimaging Optics: A Tutorial.** *Adv. Opt. Photon.* **2018**, *10*, 484–511.

11. Rau, U. **Reciprocity Relation between Photovoltaic Quantum Efficiency and Electroluminescent Emission of Solar Cells.** *Phys. Rev. B* **2007**, *76*, 085303.

12. Xiang, C.; Koo, W.; So, F.; Sasabe, H.; Kido, J. **A Systematic Study on Efficiency Enhancements in Phosphorescent Green, Red and Blue Microcavity Organic Light-Emitting Devices.** *Light: Sci. Appl.* **2013**, *2*, e74. https://doi.org/10.1038/lsa.2013.30.

13. Song, J.; et al. **Lensfree OLEDs with over 50% External Quantum Efficiency via External Scattering and Horizontally Oriented Emitters.** *Nat. Commun.* **2018**, *9*, 3207. https://doi.org/10.1038/s41467-018-05671-x.

14. **Using Angle-Selective Optical Film to Enhance the Light Extraction of a Thin-Film Encapsulated 3D Reflective Pixel for OLED Displays.** **2022**. https://pubmed.ncbi.nlm.nih.gov/36558597/.

15. Kim, H.-J.; et al. **High Efficient OLED Displays Prepared with the Air-Gapped Bridges on Quantum Dot Patterns for Optical Recycling.** *Sci. Rep.* **2017**, *7*, 43063. https://doi.org/10.1038/srep43063.
