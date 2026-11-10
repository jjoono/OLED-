# Reciprocity Constrains Angular Mixing in Passive Light-Extraction Structures: Two Spectral Metrics for OLED Outcoupling

**Authors** — [TBD]

---

## Abstract

Microlens arrays, scattering layers, and textured films are widely used to extract substrate modes from organic light-emitting diodes (OLEDs), but there is no common metric for comparing them. We treat multi-pass extraction as a Markov chain on the polar angle in the substrate. Reciprocity makes the return matrix of a passive structure obey detailed balance with respect to the étendue, which yields two spectral metrics: an angular mixing number $N_{\mathrm{mix}} = 1/(1-\lambda_2)$, the number of passes needed to forget the angle, and a trapping ratio $\xi = p_\infty/p$, where $p_\infty$ is the escape probability of long-circulating light and is bounded by the escape probability $p$ for uniform radiance. Reciprocity also fixes the redistribution in any reduced model. With it, $p(\theta)$ and $N_{\mathrm{mix}}$ reproduce the full angle-resolved series to within 0.30 percentage points (pp) on average and 3.0 pp at most for 342 combinations of 38 planar, scattering, microlens, random-lens, and wavy structures, three devices, and three wavelengths; the memoryless formula deviates by up to 60 pp. Random lens arrays and wavy surfaces combine the escape probability of microlenses with the mixing of scattering layers, yet gain less than 1 pp, because reciprocity prevents mixing from opening angles that a structure does not transmit. This explains why freeform lens optimization gains at most 7% in any polar band over an equally optimized hemisphere.

---

## 1. Introduction

A planar OLED loses much of its generated light. Waveguided modes, losses at the metal electrode, and total internal reflection (TIR) at the substrate/air interface all trap photons [1,2]. External structures on the substrate extract the substrate modes without touching the electrical stack: microlens arrays (MLAs) [3–5], scattering layers [6,7], and textured films [8]. Ray-tracing and experimental studies have mapped how their extraction depends on lens profile, aspect ratio, and fill factor [3,9–11], and their shapes are increasingly optimized by inverse design and machine learning [12,13]. Periodic, disordered, and quasi-periodic textures have also been compared, mainly for internal, wave-optical structures [14]. With such structures, OLED external quantum efficiencies (EQEs) have generally plateaued near 55–60%, a ceiling recently traced to the loss that light returned by the structure suffers on each round trip through the device [15].

What is missing is a common way to compare the structures themselves. Two descriptions are in use. The first is the statistical ray-optics picture, in which the light in the substrate is assumed to be fully randomized in angle on every pass [16,17]. The extraction efficiency then depends on one escape probability $p$ and one round-trip loss $A'$, $\eta = p/[p+(1-p)A']$ [15]. The second propagates the angle-resolved power distribution from pass to pass with redistribution matrices, as in the OPTOS formalism for solar cells [18], in models of textured LEDs [19], and in the trans-scale model of ref. [15]. The matrix models are accurate but numerical, and the degree of randomization they contain is not extracted. Where randomization has been quantified, it has been done empirically: by a randomization fraction [20], by the number of passes after which the internal distribution looks Lambertian [21], or not at all, as noted for partial randomization in thin-film solar cells [22].

The distinction matters for OLEDs. A hemispherical MLA returns unextracted light largely at its original angle [15]. Light at angles it transmits well leaves in the first passes, and later passes see mostly light at angles it transmits poorly. The memoryless formula then overestimates extraction, and the error depends on the structure. It also matters for the expected value of shape optimization. If a structure acts mainly through how it redistributes returned light, a freeform profile may change the redistribution without changing how much light escapes.

Here we introduce an angular mixing number, $N_{\mathrm{mix}} = 1/(1-\lambda_2)$, defined by the second eigenvalue of a structure's return matrix. It is the per-pass counterpart of the transport factor $1/(1-g)$ of radiative transfer [23] and of the relaxation time of a reversible Markov chain [24]. We show that reciprocity fixes the form of the mixing in any reduced description, that $p(\theta)$ and $N_{\mathrm{mix}}$ then suffice to predict multi-pass extraction for planar, scattering, microlens, random-lens, and wavy structures, and that a second spectral quantity, the escape probability of the quasi-stationary distribution [25,26], identifies light trapped at dead angles. Both spectral quantities are standard objects of Markov-chain theory and of the theory of open, leaking dynamical systems, where the escape rate and the conditionally invariant measure describe the long-time survivors of chaotic billiards and dielectric microcavities [27–30]. What is new here is their identification for passive extraction structures, the constraints that reciprocity places on them, and their use as reporting metrics for OLED outcoupling. We place the structure classes on a common map and connect it to the radiance-ordered extraction bound [31–33]. As a case study, we use a controlled freeform-versus-hemisphere benchmark on an experimentally validated OLED platform [34,35] to show why lens-shape freedom saturates.

---

## 2. Theory

### 2.1 Multi-pass extraction as a Markov chain on angle

Light delivered to the substrate is described by a vector $\mathbf{v}$ of power in polar-angle bins $\theta_k$. On each encounter with the structure a fraction $p(\theta_k)$ escapes into air, and the remainder returns into the substrate with the redistribution matrix $B_R(\theta_i \to \theta_r)$. The device reflects the returned light with its angle-resolved reflectance $R(\theta)$, without changing its angle. The power reaching the structure after $m$ round trips is $\mathbf{v}^{(m)} = (\mathrm{diag}\,R\; B_R)^m\,\mathbf{v}^{(0)}$, and the substrate-to-air extraction efficiency is

$$
\eta_{\mathrm{sta}} = \mathbf{p}^{\mathsf T}\left(I - \mathrm{diag}\,R\; B_R\right)^{-1}\mathbf{v}^{(0)} ,
$$

the matrix form of the trans-scale model [15]. Normalizing each column of $B_R$ gives a column-stochastic matrix $M$, the transition matrix of a Markov chain on the polar angle of the returned light.

### 2.2 Reciprocity, detailed balance, and the angular mixing number

For a passive, reciprocal structure, the power scattered from bin $i$ into bin $j$ and from $j$ into $i$ obey $w_i\,B_R(i\to j) = w_j\,B_R(j\to i)$, where $w_k = \sin^2\theta_{k+1}-\sin^2\theta_k$ is the étendue of bin $k$; this follows from Helmholtz reciprocity of the bidirectional scattering distribution [36] and has been noted as a symmetry of the redistribution matrices of textured sheets [18]. This is detailed balance with respect to the étendue weight. $M$ is therefore similar to a symmetric matrix, its eigenvalues are real, and for a lossless structure the uniform-radiance (Lambertian) distribution is stationary. The largest eigenvalue is one. The second, $\lambda_2$, sets how fast any other initial distribution relaxes towards the stationary one: the deviation shrinks by $\lambda_2$ per pass. We define the angular mixing number

$$
N_{\mathrm{mix}} = \frac{1}{1-|\lambda_2|},
$$

the number of passes the structure needs to forget the angle of the incident light. A planar interface has $\lambda_2 = 1$ and $N_{\mathrm{mix}} = \infty$; a Lambertian scatterer has $\lambda_2 = 0$ and $N_{\mathrm{mix}} = 1$. $N_{\mathrm{mix}}$ is the relaxation time of the chain [24] and the per-pass analogue of the $1/(1-g)$ scattering events that a photon needs to lose its direction in radiative transfer [23]. Because states at near-grazing angles carry almost no étendue but can return light at the same angle, we evaluate $\lambda_2$ on $\theta < 85°$, which excludes 0.76% of the étendue (Note S10).

### 2.3 Reduced model and the reciprocity constraint

The reduced model keeps the angle-resolved escape probability $p(\theta)$ and replaces $B_R$ by a single mixing parameter $m$: a fraction $1-m$ of the returned light keeps its angle, and a fraction $m$ is redistributed with a fixed distribution $q$,

$$
B_R(\theta_i\to\theta_r) = \left[1-p(\theta_i)\right]\left[(1-m)\,\delta_{ir} + m\,q(\theta_r)\right].
$$

Reciprocity fixes $q$. Detailed balance, $w_i B_R(i\to r) = w_r B_R(r\to i)$, holds for the mixing term only if $q(\theta) \propto w(\theta)\,[1-p(\theta)]$: light is redistributed into each angle in proportion to how much light that angle itself returns. The common choice $q \propto w$, a Lambertian redistribution, violates reciprocity whenever $p$ depends on angle, because it sends returned light into the escape cone of a structure that returns almost nothing from there. With the reciprocal $q$, the normalized kernel has eigenvalues 1 and $1-m$, so we set $m = 1-\lambda_2$ and the reduced model relaxes at the same rate as the full matrix. With $m = 1$, $q \propto w$, and a scalar $p$, the model reduces to the memoryless formula of ref. [15]; with $m = 0$ every ray returns at its own angle.

### 2.4 Late-pass escape probability and trapping ratio

$\lambda_2$ describes how fast the angle of the returned light is forgotten, not where the light ends up. The second quantity we use is the largest eigenvalue $\rho$ of the sub-stochastic return operator $B_R$ itself. Light that has circulated many times without escaping approaches the corresponding eigenvector, the quasi-stationary distribution of the chain [25,26], and escapes with probability

$$
p_\infty = 1-\rho(B_R)
$$

per encounter. In the language of open dynamical systems, $p_\infty$ is the escape rate per encounter and the eigenvector is the conditionally invariant measure [27]. The same construction describes the long-time survival of rays in dielectric microcavities [28], in systems with partial absorption [29], and in systems with partial, Fresnel-type escape [30], to which the device reflectance $R(\theta)$ adds a further, angle-dependent loss. Because $B_R$ is self-adjoint in the étendue-weighted inner product, the Rayleigh quotient with the uniform-radiance vector gives $\rho \ge 1-p$, and the column sums give $\rho \le \max_\theta\,[1-p(\theta)]$. Hence

$$
\min_\theta p(\theta) \le p_\infty \le p .
$$

We call $\xi = p_\infty/p$ the trapping ratio. It equals one when long-lived light escapes as easily as Lambertian light, and zero when some angles never escape and never mix (planar interface). Like $\lambda_2$, $p_\infty$ is evaluated on $\theta < 85°$, and $p$ in $\xi$ is taken on the same range.

### 2.5 Band-resolved radiance-ordered bound

The power that a passive structure can send into a set of output directions in one encounter is bounded by filling the corresponding output étendue with the brightest internal radiance first [31–33]. For a uniform internal radiance the bound reduces to the étendue ratio, $1/n^2$ for the full hemisphere [16] and $(\sin^2\theta_{\mathrm{hi}}-\sin^2\theta_{\mathrm{lo}})/n^2$ for a polar band in air. An OLED microcavity, however, launches an anisotropic radiance into the substrate. We resolve the bound into the four polar bands used below, with the actual first-pass radiance of the device, and use it to rate how much of the available phase space each structure uses.

---

## 3. Methods

### 3.1 Devices

Dipole emission and the device reflectance are computed with the Chance–Prock–Silbey model [37] extended to birefringent layers and with transfer matrices, using measured optical constants [15,34]. Three bottom-emitting stacks are used. The first is the platform of the freeform benchmark: Al (100 nm) / B3PyMPM / TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd (25 nm) / TCTA (10 nm) / TAPC / ITO (150 nm) / glass ($n = 1.51$), with a horizontal dipole ratio of 0.865. Its round-trip loss for Lambertian light is $A' = 0.17$ at 530 nm. The second replaces Al by Ag ($A' = 0.049$). The third is a low-loss stack, Ag / B3PyMPM (200 nm) / EML / TAPC (180 nm) / ITO (50 nm), on glass ($A' = 0.027$) [15]. The angular distribution of the light delivered to the substrate and the device reflectance are computed per wavelength.

### 3.2 Structures and their return matrices

All structures are index-matched to the substrate, so the structure/air boundary sets the escape condition. Their transmission and return matrices are computed on 1° (microlenses) or 0.5° (others) polar bins.

- **Planar interface.** Fresnel reflection, returned at the incident angle.
- **Scattering layers.** Index-matched, non-absorbing slabs with a Henyey–Greenstein phase function [38] (asymmetry $g$ = 0–0.9) and transport optical thickness $\tau^* = S(1-g)$ = 0.25–4, combined with the flat exit face by the adding method; Monte Carlo with 20,000 rays per incidence bin.
- **Hexagonal MLAs.** Close-packed axisymmetric lenslets of the spline family of the benchmark (Table 1), traced with a height-field Monte-Carlo tracer that reproduces the LightTools return matrix of the hemispherical MLA (escape probability 0.4275 against 0.4274; Note S11). The set includes the hemispherical and freeform optima of the benchmark.
- **Random MLAs.** A 6 × 6 supercell in which every lenslet carries its own profile (shape spread $s$ about the quarter circle), its own height (aspect spread), or a jittered position and radius; overlapping lenslets take the higher surface.
- **Wavy surfaces.** Periodic Gaussian-correlated random height fields with rms slope 0.2–1.5. Ray optics is scale-free, so only the slope distribution matters. Such surfaces resemble the spontaneously formed buckles that have been used to extract light from OLEDs [39].

The random structures (22 in total) are listed with their parameters and results in Table S12. All return matrices are traced with 20,000 rays per incidence bin (10,000 for wavy surfaces); the fraction of rays lost to the iteration limit is below $6\times10^{-4}$ for every structure. Monte Carlo uncertainty is small: re-tracing with other ray seeds changes $N_{\mathrm{mix}}$ by at most 0.01 and $\eta_{\mathrm{sta}}$ by at most 0.06 pp, and drawing other geometric realizations with the same statistics changes $N_{\mathrm{mix}}$ by at most 0.24, the trapping ratio by at most 0.021, and $\eta_{\mathrm{sta}}$ by at most 0.23 pp (Table S13). The traced return matrices satisfy detailed balance within Monte Carlo noise: the relative asymmetry of $w_i B_R(i\to j)$ is 5–13% on 1° bins and falls to 0.2–0.7% on 15° bins (Note S10).

### 3.3 Freeform benchmark

The case study of Section 4.4 compares an axisymmetric freeform MLA with an equally optimized hemispherical MLA on the Al platform. LightTools ray tracing is driven from MATLAB, with the CPS source on a 1 mm-radius disc and a 1.295 mm glass substrate. The pipeline has been validated against fabricated devices with nearly the same geometry, including attached and inverted MLAs [34,35]. Both lens classes share the material, pitch, height range, and cavity bounds (Table 1); the hemisphere is a point of the freeform design space. Total EQE and the EQE in four polar bands (0–20°, 20–40°, 40–60°, 60–80°) are optimized as separate objectives by surrogate-based global search followed by pattern search. Every optimum is re-evaluated three times with 50,000 rays and 151 wavelengths. A restart control repeats every search from the hemispherical optimum and tests the gain with a pooled $t$ statistic, with Holm correction [40] for the five objectives (Note S2).

**Table 1 | Design variables, bounds, and feasibility rules of the benchmark.**

| Quantity | Freeform | Hemispherical reference |
|---|---|---|
| Radial control points $x_2$–$x_6$ (lens radii) | 0–1, sorted ascending | fixed on the quarter circle |
| Height control points $y_2$–$y_6$ (lens radii) | 0–1.5 | fixed on the quarter circle |
| End points | (0, 1) and (1, 0), fixed | same |
| Height scale (stretch$_Z$) | 0.1–3 | 0.1–3 |
| ETL and HTL thickness | 10–150 nm each | 10–150 nm each |
| Feasibility rules | no self-intersection or collinear triples; control-point spacing 0.05–1.0; turning angle between adjacent segments ≤ 120° | always satisfied |
| Free variables | 13 | 3 |

**Fig. 1 | Multi-pass extraction as a Markov chain on angle.** (a) Light delivered to the substrate meets the passive structure, which transmits a fraction $p(\theta)$ and returns the rest with the redistribution $B_R(\theta_i\to\theta_r)$; the device reflects the returned light with $R(\theta)$ and a round-trip loss $A'$. (b–d) Column-normalized return matrices $M$ (log scale) of a planar interface, the hemispherical MLA, and a scattering layer ($g$ = 0.8, $\tau^*$ = 1), with $\lambda_2$ and $N_{\mathrm{mix}}$. *(`fig1_mixing_concept.png`)*

---

## 4. Results

### 4.1 Two parameters predict multi-pass extraction

Figure 2 compares four models with the full angle-resolved series for 342 combinations: three devices, three wavelengths (500, 550, and 600 nm), and 38 structures (planar, six scattering layers, nine hemispherical and freeform MLAs, 17 random MLAs, and five wavy surfaces). The memoryless formula fails for structures that retain the angle (Fig. 2a). It overestimates the planar interface by up to 60 pp, because it lets returned light reach the escape cone that a flat face never lets it reach, and the hemispherical MLA by 2–3 pp; over all cases its mean error is 3.9 pp. Assuming that every ray returns at its own angle fails in the opposite direction for scattering layers (Fig. 2b).

With $m = 1-\lambda_2$ and Lambertian redistribution ($q \propto w$), the reduced model is accurate for lens arrays but overestimates the shallow wavy surfaces by up to 10.5 pp (Fig. 2c). Their $\lambda_2$ is low, yet light trapped beyond the critical angle is returned mostly to other trapped angles: on the slope-0.2 surface, only 3–8% of the light returned from 70–80° lands at angles with $p(\theta) > 0.3$, against 50% for Lambertian redistribution. The reciprocal redistribution $q \propto w(1-p)$ of Section 2.3 builds this in, since angles that return much light also receive much of the mixed light. It reduces the mean error to 0.30 pp and the maximum to 3.0 pp over all 342 cases, with no additional parameter (Fig. 2d, Table S11). The mean error is 0.29 pp for scattering layers, 0.41 pp for hemispherical and freeform MLAs, 0.16 pp for random MLAs, 0.19 pp for wavy surfaces of rms slope ≥ 0.7, and 1.4 pp for the two shallowest wavy surfaces, which remain the hardest case.

A structure is therefore described, for the purpose of multi-pass extraction, by its angle-resolved escape probability and one number for how fast it forgets the angle, provided that the mixing respects reciprocity.

**Fig. 2 | Validation of the reduced model.** Extraction efficiency $\eta_{\mathrm{sta}}$ predicted by (a) the memoryless formula $\eta = p/[p+(1-p)A']$, (b) the angle-preserving limit ($m = 0$), and the reduced model with $m = 1-\lambda_2$ and (c) Lambertian redistribution $q \propto w$ or (d) reciprocal redistribution $q \propto w(1-p)$, against the full angle-resolved series. Symbols: planar (squares), scattering layers (circles), hemispherical (diamonds), freeform (triangles), and random (inverted triangles) MLAs, and wavy surfaces (crosses); Al, Ag, and low-loss devices at 500, 550, and 600 nm; 342 cases. The outliers in (c) are the wavy surfaces of rms slope 0.2 and 0.4. *(`fig2_reduced_model.png`)*

### 4.2 The structure map

Figure 3a places all structures on the $(p, N_{\mathrm{mix}})$ plane, with contours of $\eta_{\mathrm{sta}}$ for the Al device. The classes separate:

- **Planar interface:** $p$ = 0.40 and $N_{\mathrm{mix}} = \infty$. Light outside the escape cone never reaches it, and $\eta_{\mathrm{sta}}$ = 0.38.
- **Scattering layers:** $N_{\mathrm{mix}}$ = 1.1–2.3, but $p$ = 0.19–0.39, falling with thickness because more light is back-scattered before it reaches the exit face. $\eta_{\mathrm{sta}}$ = 0.58–0.75.
- **Hemispherical MLAs:** the highest escape probability, $p$ = 0.42, but $N_{\mathrm{mix}}$ = 6–8. $\eta_{\mathrm{sta}}$ = 0.78–0.79.
- **Freeform MLAs:** $p$ = 0.42–0.43, like the hemisphere, with $N_{\mathrm{mix}}$ = 3–7.5.

Randomness lowers $N_{\mathrm{mix}}$ without lowering $p$ (Fig. 3b, Table S12). Drawing every lenslet profile independently reduces $N_{\mathrm{mix}}$ from 7.0 to 1.9 ± 0.1 (spread over geometric realizations) while $p$ stays at 0.42–0.43. Jittering the lenslet positions by up to 30% of the pitch reduces it to 3.5 ± 0.1, and adding a 20% radius spread to 3.2. Varying only the lenslet heights is the least effective (7.1 → 4.1 ± 0.2 at ±60%). A wavy surface of rms slope 0.4–1.5 has $\lambda_2$ = 0.12–0.18 ($N_{\mathrm{mix}}$ = 1.1–1.2), as low as a scattering layer, and its escape probability rises with slope from 0.41 to 0.43, reaching that of the hemispherical MLA at slope 1.0. At slope 1.5 the wavy surface gives the highest $\eta_{\mathrm{sta}}$ of all structures on the Al device, 0.801, against 0.792 for the hemispherical MLA.

The trapping ratio adds what $N_{\mathrm{mix}}$ cannot see (Fig. 3c). For lens arrays the two are tied: $\xi$ falls from 0.80 to 0.39 as $N_{\mathrm{mix}}$ rises from 2 to 8, because light that keeps its angle for many passes ends up at the angles the lenses transmit worst. The hemispherical MLA has $\xi$ = 0.47: light that survives many passes escapes at less than half the rate of fresh Lambertian light. Randomizing the lenslets raises $\xi$ to 0.7–0.8, wavy surfaces of rms slope ≥ 0.7 reach 0.81–0.88, and thick scattering layers 0.93–0.99. Two structures fall well below the lens-array trend at low $N_{\mathrm{mix}}$: the wavy surface of slope 0.2 ($N_{\mathrm{mix}}$ = 2.2, $\xi$ = 0.35) and the thin isotropic scattering layer ($N_{\mathrm{mix}}$ = 2.2, $\xi$ = 0.54). Both mix quickly but locally, so trapped light is recycled among trapped angles. These are the structures with the largest residual errors in Fig. 2d. $\xi$ thus flags the cases in which a single mixing rate is least reliable.

The contours of Fig. 3a show what lower $N_{\mathrm{mix}}$ and higher $\xi$ are worth, and the answer is: little, once $p$ is high. The contours run almost vertically, so $\eta_{\mathrm{sta}}$ is set mainly by $p$. At fixed $p(\theta)$, the span available to mixing lies between full angle retention ($m = 0$) and complete reciprocal mixing ($m = 1$). For the hemispherical MLA on the Al device this span is 0.780–0.806, and the real structure sits at 0.792. Randomizing the lenslet profiles moves $N_{\mathrm{mix}}$ from 7.0 to 1.9 but changes $\eta_{\mathrm{sta}}$ by less than 1 pp (0.794 → 0.785), because the same randomness also reshapes $p(\theta)$. The planar interface makes the reciprocity constraint explicit. Even complete mixing at its own $p(\theta)$ raises $\eta_{\mathrm{sta}}$ only from 0.38 to 0.49, because a face that transmits nearly all light inside the escape cone also returns almost none from there, and reciprocity sends the mixed light back in the same proportion. Mixing cannot open angles at which a structure does not transmit; only a change of $p(\theta)$ can (Section 5).

**Fig. 3 | The structure map.** (a) Escape probability for Lambertian incidence and angular mixing number of planar, scattering, hemispherical, freeform, random-lens, and wavy structures, with contours of $\eta_{\mathrm{sta}}$ from the reciprocal reduced model for the Al device at 550 nm, computed with the angular shape of $p(\theta)$ of the hemispherical MLA scaled to each $p$. (b) $N_{\mathrm{mix}}$ against the randomness parameter: lenslet shape spread $s$, position jitter (fraction of the pitch), lenslet height spread, and rms slope of a wavy surface; the blue dotted line marks the hemispherical MLA and the green dotted line a scattering layer ($g$ = 0.8, $\tau^*$ = 1). (c) Trapping ratio $\xi = p_\infty/p$ against $N_{\mathrm{mix}}$; the planar interface ($N_{\mathrm{mix}} = \infty$, $\xi$ = 0) is drawn at the right edge. *(`fig3_structure_map.png`)*

### 4.3 How much of the radiance bound passive structures use

On the first pass, the microcavity delivers light to the substrate with a radiance that is far from uniform. For the Al platform at 550 nm, filling each air band with the brightest internal radiance allows 0.085, 0.205, 0.230, and 0.156 of the substrate power into the 0–20°, 20–40°, 40–60°, and 60–80° bands, and 0.540 in total (Fig. 4a). The uniform-radiance bound gives 0.051, 0.130, 0.148, and 0.096, and 0.439 in total. Every passive structure we examined stays well below the band-resolved bound. The MLAs reach 76–80% of it in total and 55–70% per band, a wavy surface of rms slope 1.0 77%, the planar interface 65%, and the scattering layers 35–68%. The lenses transmit almost exactly as much as a structure facing uniform radiance would (0.43 against 0.44), so the anisotropy that the cavity provides is not used.

The output composition, finally, is nearly independent of the structure type (Fig. 4b). Over planar, scattering, hemispherical, freeform, random, and wavy structures, the share of the output in the 40–60° band stays between 0.33 and 0.36, close to the Lambertian partition of 0.337. Within the class of angle-insensitive passive structures, the angular composition of multi-pass output is therefore set by the platform rather than by the structure.

**Fig. 4 | Radiance bound and composition window.** (a) First-pass power into each polar band for the Al platform at 550 nm: uniform-radiance bound, radiance-ordered bound for the actual cavity emission, hemispherical MLA, wavy surface (rms slope 1.0), scattering layer ($g$ = 0.8, $\tau^*$ = 1), and planar interface. (b) 40–60° selectivity of the multi-pass output against $\eta_{\mathrm{sta}}$ for all structures; the band marks 0.33–0.36 and the dotted line the Lambertian partition. *(`fig4_bound_window.png`)*

### 4.4 Case study: freeform versus an equally optimized hemisphere

The map explains the result of a controlled shape-optimization benchmark. On the Al platform, band-optimized freeform profiles differ visibly from the hemisphere and from one another (Fig. 5a), but the gains are small: $G_j$ = 1.003, 1.002, 1.070, and 1.036 for the four bands, and $G$ = 1.013 for total EQE (Fig. 5b, Table 2). The restart control shows that these margins come from the design space and not from a weak search (Fig. 5c). Gains are significant at 40–60° (+4.70%, $t$ = 115), 60–80° (+3.64%, $t$ = 30), and for total EQE (+1.01%, $t$ = 69), and remain significant under Holm correction. At 0–40° the residues are +0.32% and +0.20%; a five-repeat re-measurement resolves the 20–40° residue as +0.29% ($t$ = 4.1). The same procedure that resolves a +0.29% residue does not miss gains ten times larger.

**Table 2 | Freeform versus hemisphere for each objective.**

| Objective | Hemisphere (archived) | Best freeform | Source of freeform entry | $G_j$ | Restart control: hemisphere → restarted (mean ± SD) | Gain, $t$ |
|---|---|---|---|---|---|---|
| 0–20° | 0.06682 | 0.06699 | restart control | 1.003 | 0.06677 ± 0.00017 → 0.06699 ± 0.00013 | +0.32%, 1.7 |
| 20–40° | 0.16591 | 0.16628 | restart control | 1.002 | 0.16595 ± 0.00017 → 0.16628 ± 0.00020 | +0.20%, 2.1 (five repeats: +0.29%, 4.1) |
| 40–60° | 0.18503 | 0.19802 | per-band campaign | 1.070 | 0.18486 ± 0.00008 → 0.19355 ± 0.00010 | +4.70%, 115 |
| 60–80° | 0.13589 | 0.14075 | restart control | 1.036 | 0.13581 ± 0.00014 → 0.14075 ± 0.00025 | +3.64%, 30 |
| Total EQE | 0.54679 | 0.5539 | dedicated total-EQE campaign | 1.013 | 0.54679 ± 0.00013 → 0.55231 ± 0.00004 | +1.01%, 69 |

In the language of Section 4.2, every optimized freeform lens has the same escape probability as the hemisphere ($p$ = 0.42–0.43) and differs from it only in $N_{\mathrm{mix}}$ and in the shape of $p(\theta)$. Shape freedom thus moves a structure vertically on the map, where the $\eta_{\mathrm{sta}}$ contours are flat, and cannot move it to the right, where they are steep. The band selectivities of all designs follow the same picture: across 606 sampled designs the 40–60° band EQE is proportional to the total EQE ($R^2$ = 0.968), and band-dedicated optimization raises the band power by at most a factor of 1.22 net of the loss in total EQE (Note S4). The same collapse holds for inverted and randomly assembled arrays (Note S7).

**Fig. 5 | Freeform versus an equally optimized hemisphere.** (a) Band-optimized freeform profiles against the hemispherical reference at its optimized height. (b) Best EQE of each class at each objective, labelled with $G_j$. (c) Restart control: gain over the hemisphere when each search restarts from the hemispherical optimum, labelled with $t$; green marks gains above the one-sided 95% threshold. *(`fig5_hemisphere_benchmark.png`)*

---

## 5. Discussion

**Why lens shape saturates.** The escape probability of a lens array for uniform radiance is close to the étendue bound $1/n^2$ already for a hemisphere, so profile changes cannot raise $p$ appreciably. What they can change is how the returned light is redistributed, that is, $N_{\mathrm{mix}}$, and the shape of $p(\theta)$. Along that direction the efficiency contours are flat. The anisotropy of the cavity emission would allow more, but a refractive surface that is coextensive with an extended source cannot map the bright part of the internal radiance into a chosen band; it behaves as if the radiance were uniform.

**Why mixing matters little.** Reciprocity ties mixing to escape. A structure that sends trapped light into the escape cone must, by detailed balance, return light from the cone to trapped angles in the same étendue-weighted proportion, which lowers its escape probability for fresh light. This is the trade-off seen in scattering layers: they reach $\xi \approx 1$ but lose $p$ with thickness. For a structure with a given $p(\theta)$, the span between full retention and complete reciprocal mixing is 1.8–4.6 pp on the Al device, 1.7–7.0 pp on the Ag device, and 0.1–1.7 pp on the low-loss device for the lens arrays we examined, and real structures already capture part of it. The random lens arrays show that the remainder is hard to collect: the randomness that lowers $N_{\mathrm{mix}}$ and raises $\xi$ also flattens $p(\theta)$, and the effects nearly cancel. Lowering $N_{\mathrm{mix}}$ further therefore buys at most a few percentage points, and only in lossy devices. In low-loss devices the choice of outcoupling structure, and its shape in particular, matters even less, in line with ref. [15].

**Design implications.** Three numbers are worth reporting for any new outcoupling structure: its escape probability for uniform radiance $p$, its angular mixing number $N_{\mathrm{mix}}$, and its trapping ratio $\xi$. The first sets the efficiency, the second the accuracy of any memoryless description, and the third whether light that survives many passes is still extractable. A structure that combines a high $p$ with a low $N_{\mathrm{mix}}$ and a high $\xi$, such as a random MLA or a wavy surface of rms slope one or more, avoids both the low escape probability of thick scattering layers and the angle retention of periodic lenses, without introducing absorption. Its advantage in efficiency over a periodic MLA is small (≤ 1 pp); whether it brings other practical benefits, such as weaker periodicity-related angular features, is outside the present analysis. Optimizers, including machine-learning ones [12,13], can then be benchmarked against the $(p, N_{\mathrm{mix}})$ map and the band-resolved bound rather than against each other.

**Scope.** The analysis covers passive, angle-insensitive external structures index-matched to the substrate, in the ray-optical regime. Elements that select the angle of transmission, such as angular filters, are also passive and also obey the radiance bound per encounter, but in combination with recycling they can change the output composition; they are outside the present scope. Internal structures that act on waveguided modes, absorbing scatterers, and wavelength-converting layers are not covered either. The reduced model assumes that the device reflectance preserves the angle and that $p(\theta)$ is azimuthally symmetric. Its largest residual errors (up to 3 pp) occur for shallow surfaces with local mixing, which the trapping ratio identifies. The three wavelengths (500, 550, and 600 nm) bracket the green emission band; the conclusions hold at each of them (Table S11).

---

## 6. Conclusion

We introduced an angular mixing number, $N_{\mathrm{mix}} = 1/(1-\lambda_2)$, and a trapping ratio, $\xi = p_\infty/p$ with $p_\infty = 1-\rho(B_R)$, from the spectrum of the return matrix of a passive light-extraction structure. Reciprocity makes the spectrum real, bounds the late-pass escape probability by the escape probability for uniform radiance, and fixes how a reduced model must redistribute the returned light. With this reciprocal redistribution, the angle-resolved escape probability and $N_{\mathrm{mix}}$ predict multi-pass extraction to within 0.30 pp on average and 3.0 pp at most across 38 planar, scattering, microlens, random-lens, and wavy structures, where the memoryless formula deviates by up to 60 pp. On the resulting map, periodic microlenses combine high escape probability with long angle retention and a low trapping ratio, scattering layers the opposite, and random lens arrays and wavy surfaces combine the escape probability of the former with the mixing of the latter. Because reciprocity prevents mixing from opening angles that a structure does not transmit, the efficiency is set mainly by the escape probability, and this combination is worth less than 1 pp. Every passive structure uses only part of the band-resolved radiance bound, and angle-insensitive structures share one output composition. Freeform optimization of microlens shape mainly changes $N_{\mathrm{mix}}$, which explains why it gains at most 7% in any band over an equally optimized hemisphere on an extended OLED emitter.

---

## Funding

[TBD]

## Acknowledgments

[TBD]

## Disclosures

The authors declare no conflicts of interest.

## Data availability

Return matrices, scripts, and the data underlying all figures are available in Supplement 1 and from the authors upon reasonable request.

## Supplemental document

See Supplement 1 for supporting content.

---

## References

1. Brütting, W.; Frischeisen, J.; Schmidt, T. D.; Scholz, B. J.; Mayr, C. **Device Efficiency of Organic Light-Emitting Diodes: Progress by Improved Light Outcoupling.** *Phys. Status Solidi A* **2013**, *210*, 44–65.

2. Greiner, H. **Light Extraction from Organic Light Emitting Diode Substrates: Simulation and Experiment.** *Jpn. J. Appl. Phys.* **2007**, *46*, 4125–4137. https://doi.org/10.1143/JJAP.46.4125.

3. Möller, S.; Forrest, S. R. **Improved Light Out-Coupling in Organic Light Emitting Diodes Employing Ordered Microlens Arrays.** *J. Appl. Phys.* **2002**, *91*, 3324–3327.

4. Wrzesniewski, E.; et al. **Enhancing Light Extraction in Top-Emitting Organic Light-Emitting Devices Using Molded Transparent Polymer Microlens Arrays.** *Small* **2012**, *8*, 2647–2651. https://doi.org/10.1002/smll.201102662.

5. Qu, Y.; Kim, J.; Coburn, C.; Forrest, S. R. **Efficient, Nonintrusive Outcoupling in Organic Light Emitting Devices Using Embedded Microlens Arrays.** *ACS Photonics* **2018**, *5*, 2453–2458. https://doi.org/10.1021/acsphotonics.8b00255.

6. Song, J.; et al. **Lensfree OLEDs with over 50% External Quantum Efficiency via External Scattering and Horizontally Oriented Emitters.** *Nat. Commun.* **2018**, *9*, 3207. https://doi.org/10.1038/s41467-018-05671-x.

7. Levell, J. W.; Harkema, S.; Pendyala, R. K.; Rensing, P.; Senes, A.; Bollen, D.; MacKerron, D.; Wilson, J. S. **Device Reflectivity as a Simple Rule for Predicting the Suitability of Scattering Foils for Improved OLED Light Extraction.** *Proc. SPIE* **2013**, *8829*, 88291L. https://doi.org/10.1117/12.2023021.

8. Kovačič, M. **Design of External Microtextures for Efficient Light Outcoupling in OLEDs with Different Preferential Orientation of Emission Dipoles.** *Opto-Electron. Rev.* **2022**, *30*, e141542. https://doi.org/10.24425/opelre.2022.141542.

9. Sun, Y.; Forrest, S. R. **Organic Light Emitting Devices with Enhanced Outcoupling via Microlenses Fabricated by Imprint Lithography.** *J. Appl. Phys.* **2006**, *100*, 073106. https://doi.org/10.1063/1.2356904.

10. Peng, H.; Ho, Y. L.; Yu, X.-J.; Wong, M.; Kwok, H.-S. **Coupling Efficiency Enhancement in Organic Light-Emitting Devices Using Microlens Array—Theory and Experiment.** *J. Display Technol.* **2005**, *1*, 278–282. https://doi.org/10.1109/JDT.2005.858944.

11. Hwang, D. H.; Kwon, O. T.; Lee, W. J.; Hong, J. W.; Kim, T. W. **Outcoupling Efficiency of Organic Light-Emitting Diodes Depending on the Fill Factor and Size of the Microlens Array.** *Phys. Status Solidi A* **2014**, *211*, 1773–1777. https://doi.org/10.1002/pssa.201330613.

12. Kim, S.; Shin, J. M.; Lee, J.; Park, C.; Lee, S.; Park, J.; Seo, D.; Park, S.; Park, C. Y.; Jang, M. S. **Inverse Design of Organic Light-Emitting Diode Structure Based on Deep Neural Networks.** *Nanophotonics* **2021**, *10*, 4533–4541. https://doi.org/10.1515/nanoph-2021-0434.

13. Zhao, Z.-R.; et al. **Deep Q-Network-Guided TiO2/NOA61 Microstructures for Light Outcoupling and Angular Emission Redistribution in Blue Top-Emitting OLEDs.** *Adv. Funct. Mater.* **2026**, e78939. https://doi.org/10.1002/adfm.78939.

14. Zhang, Y.; Biswas, R.; Shinar, R.; Shinar, J. **Simulation of Enhanced Light Extraction from Periodic, Disordered, and Quasi-Periodic OLED Structures.** *J. Opt. Soc. Am. B* **2021**, *38*, C144. https://doi.org/10.1364/JOSAB.430593.

15. Kim, J.; Lee, S.; Yoo, S.; Gather, M. C. **Minimizing Round-Trip Loss for Near-Unity Quantum Efficiency in Organic Light-Emitting Diodes.** *arXiv* **2026**, 2610.12065.

16. Yablonovitch, E. **Statistical Ray Optics.** *J. Opt. Soc. Am.* **1982**, *72*, 899–907.

17. Schnitzer, I.; Yablonovitch, E.; Caneau, C.; Gmitter, T. J.; Scherer, A. **30% External Quantum Efficiency from Surface Textured, Thin-Film Light-Emitting Diodes.** *Appl. Phys. Lett.* **1993**, *63*, 2174–2176. https://doi.org/10.1063/1.110575.

18. Eisenlohr, J.; Tucher, N.; Höhn, O.; Hauser, H.; Peters, M.; Kiefel, P.; Goldschmidt, J. C.; Bläsi, B. **Matrix Formalism for Light Propagation and Absorption in Thick Textured Optical Sheets.** *Opt. Express* **2015**, *23*, A502–A518. https://doi.org/10.1364/OE.23.00A502.

19. David, A. **Surface-Roughened Light-Emitting Diodes: An Accurate Model.** *J. Display Technol.* **2013**, *9*, 301–316. https://doi.org/10.1109/JDT.2013.2240373.

20. Deckman, H. W.; Roxlo, C. B.; Yablonovitch, E. **Maximum Statistical Increase of Optical Absorption in Textured Semiconductor Films.** *Opt. Lett.* **1983**, *8*, 491–493. https://doi.org/10.1364/OL.8.000491.

21. Manzoor, S.; Filipič, M.; Onno, A.; Topič, M.; Holman, Z. C. **Visualizing Light Trapping within Textured Silicon Solar Cells.** *J. Appl. Phys.* **2020**, *127*, 063104. https://doi.org/10.1063/1.5131173.

22. Battaglia, C.; Boccard, M.; Haug, F.-J.; Ballif, C. **Light Trapping in Solar Cells: When Does a Lambertian Scatterer Scatter Lambertianly?** *J. Appl. Phys.* **2012**, *112*, 094504. https://doi.org/10.1063/1.4761988.

23. van Rossum, M. C. W.; Nieuwenhuizen, T. M. **Multiple Scattering of Classical Waves: Microscopy, Mesoscopy, and Diffusion.** *Rev. Mod. Phys.* **1999**, *71*, 313–371. https://doi.org/10.1103/RevModPhys.71.313.

24. Levin, D. A.; Peres, Y. **Markov Chains and Mixing Times**, 2nd ed.; American Mathematical Society: Providence, RI, **2017**. https://doi.org/10.1090/mbk/107.

25. Darroch, J. N.; Seneta, E. **On Quasi-Stationary Distributions in Absorbing Discrete-Time Finite Markov Chains.** *J. Appl. Probab.* **1965**, *2*, 88–100. https://doi.org/10.2307/3211876.

26. Collet, P.; Martínez, S.; San Martín, J. **Quasi-Stationary Distributions: Markov Chains, Diffusions and Dynamical Systems**; Springer: Berlin, Heidelberg, **2013**. https://doi.org/10.1007/978-3-642-33131-2.

27. Altmann, E. G.; Portela, J. S. E.; Tél, T. **Leaking Chaotic Systems.** *Rev. Mod. Phys.* **2013**, *85*, 869–918. https://doi.org/10.1103/RevModPhys.85.869.

28. Ryu, J.-W.; Lee, S.-Y.; Kim, C.-M.; Park, Y.-J. **Survival Probability Time Distribution in Dielectric Cavities.** *Phys. Rev. E* **2006**, *73*, 036207. https://doi.org/10.1103/PhysRevE.73.036207.

29. Altmann, E. G.; Portela, J. S. E.; Tél, T. **Chaotic Systems with Absorption.** *Phys. Rev. Lett.* **2013**, *111*, 144101. https://doi.org/10.1103/PhysRevLett.111.144101.

30. Clauß, K.; Altmann, E. G.; Bäcker, A.; Ketzmerick, R. **Structure of Resonance Eigenfunctions for Chaotic Systems with Partial Escape.** *Phys. Rev. E* **2019**, *100*, 052205. https://doi.org/10.1103/PhysRevE.100.052205.

31. Bortz, J. C.; Shatz, N. E.; Ries, H.; Winston, R. **Consequences of Étendue and Skewness Conservation for Nonimaging Devices with Inhomogeneous Sources and Targets.** *Proc. SPIE* **1997**, *3139*, 59–75. https://doi.org/10.1117/12.290222.

32. Muschaweck, J.; Wiesmann, C. **Extracting Light out of LEDs.** *Adv. Opt. Technol.* **2013**, *2*, 291–300. https://doi.org/10.1515/aot-2013-0027.

33. Lenef, A.; Piquette, A.; Kelso, J. **Thermodynamics of Light Extraction from Luminescent Materials.** *ECS J. Solid State Sci. Technol.* **2018**, *7*, R3211–R3226. https://doi.org/10.1149/2.0271801jss.

34. Kim, M.; Kim, J.; Yoo, S. **Near-Planar Light Outcoupling Structures with Finite Lateral Dimensions for Ultra-Efficient and Optical Crosstalk-Free OLED Displays.** *Nat. Commun.* **2025**, *16*, 11606. https://doi.org/10.1038/s41467-025-66538-6.

35. Kim, J.; Kim, E.; Park, J.; Song, J.; Kim, S.; Moon, H.; Yoo, S. **Toward Near-Foldable Surface Light Sources with Ultimate Efficiency: Ultrathin Substrates Embedded with Micron-Scale Inverted Lens Arrays.** *ACS Photonics* **2023**, *10*, 1775–1782. https://doi.org/10.1021/acsphotonics.3c00017.

36. Clarke, F. J. J.; Parry, D. J. **Helmholtz Reciprocity: Its Validity and Application to Reflectometry.** *Light. Res. Technol.* **1985**, *17*, 1–11. https://doi.org/10.1177/14771535850170010301.

37. Chance, R. R.; Prock, A.; Silbey, R. **Molecular Fluorescence and Energy Transfer Near Interfaces.** *Adv. Chem. Phys.* **1978**, *37*, 1–65. https://doi.org/10.1002/9780470142561.ch1.

38. Henyey, L. G.; Greenstein, J. L. **Diffuse Radiation in the Galaxy.** *Astrophys. J.* **1941**, *93*, 70–83. https://doi.org/10.1086/144246.

39. Koo, W. H.; Jeong, S. M.; Araoka, F.; Ishikawa, K.; Nishimura, S.; Toyooka, T.; Takezoe, H. **Light Extraction from Organic Light-Emitting Diodes Enhanced by Spontaneously Formed Buckles.** *Nat. Photonics* **2010**, *4*, 222–226. https://doi.org/10.1038/nphoton.2010.7.

40. Holm, S. **A Simple Sequentially Rejective Multiple Test Procedure.** *Scand. J. Stat.* **1979**, *6*, 65–70.
