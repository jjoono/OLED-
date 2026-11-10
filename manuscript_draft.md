# Reciprocity Constrains Angular Mixing in Passive Light-Extraction Structures: Two Spectral Metrics for OLED Outcoupling

**Authors** — [TBD]

---

## Abstract

Microlens arrays, scattering layers, and textured films are widely used to extract substrate modes from organic light-emitting diodes (OLEDs), but there is no common metric for comparing them. We treat multi-pass extraction as a Markov chain on the polar angle in the substrate. Building on the reciprocity of passive scattering, which makes such transport operators self-adjoint in the étendue inner product, we characterize each structure by two spectral metrics: an angular mixing number $N_{\mathrm{mix}} = 1/(1-|\lambda_2|)$, the number of passes needed to forget the angle, and a late-pass escape ratio $\xi$, the escape probability $p_\infty$ of long-circulating light relative to that of uniform radiance, by which it is bounded. A reduced model whose redistribution respects reciprocity, a constraint known from reflectance modelling but not imposed in multi-pass extraction models, then predicts extraction from two inputs: $p(\theta)$ and $N_{\mathrm{mix}}$ reproduce the full angle-resolved series to within 0.30 percentage points (pp) on average and 3.0 pp at most for 342 combinations of 38 planar, scattering, microlens, random-lens, and wavy structures, three devices, and three wavelengths; the memoryless formula deviates by up to 60 pp. Random lens arrays and wavy surfaces combine the escape probability of microlenses with the mixing of scattering layers, yet gain less than 1 pp, because reciprocity prevents mixing from opening angles that a structure does not transmit. Consistently, freeform lens optimization gains at most 7% in any polar band over an equally optimized hemisphere.

---

## 1. Introduction

A planar OLED loses much of its generated light. Waveguided modes, losses at the metal electrode, and total internal reflection (TIR) at the substrate/air interface all trap photons [1,2]. External structures on the substrate extract the substrate modes without touching the electrical stack: microlens arrays (MLAs) [3–5], scattering layers [6,7], and textured films [8]. Ray-tracing and experimental studies have mapped how their extraction depends on lens profile, aspect ratio, and fill factor [3,9–11], and their shapes are increasingly optimized by inverse design and machine learning [12,13]. Periodic, disordered, and quasi-periodic textures have also been compared, mainly for internal, wave-optical structures [14]. With external and combined structures, the highest reported EQEs of single-unit OLEDs lie near 55–62% [6,15–17], a ceiling recently traced to the loss that light returned by the structure suffers on each round trip through the device [18]. Recent reviews summarize internal and external extraction strategies [19,20].

What is missing is a compact way to compare the structures themselves. Two descriptions are in use. The first is the statistical ray-optics picture, in which the light in the substrate is assumed to be fully randomized in angle on every pass [21,22]; closed-form multi-pass solutions exist for this ideal Lambertian limit, including light-emitting devices [23]. The extraction efficiency then depends on one escape probability $p$ and one round-trip loss $A'$, $\eta = p/[p+(1-p)A']$ [18]. The second propagates the angle-resolved power distribution from pass to pass with redistribution matrices, as in angle-discretized radiative-transfer models of OLED substrates [24], the OPTOS formalism for solar cells [25], models of textured LEDs [26], net-radiation models of OLEDs with scattering layers [27], pass-resolved analyses of external textures [28], and the trans-scale model of ref. [18]. Greiner characterized structures by tabulated multi-pass escape probabilities versus emission angle and stack reflectance [2]. These descriptions are accurate but numerical, and they do not separate escape from redistribution or extract the degree of randomization. Where randomization has been quantified, it has been done empirically: by a randomization fraction [29], by the number of passes after which the internal distribution looks Lambertian [30], or by a single-event Lambertianity factor, with the approach to a Lambertian distribution over successive passes left open [31].

The distinction matters for OLEDs. A hemispherical MLA returns unextracted light largely at its original angle [18]. Light at angles it transmits well leaves in the first passes, and later passes see mostly light at angles it transmits poorly. The memoryless formula then overestimates extraction, and the error depends on the structure. It also matters for the expected value of shape optimization. If a structure acts mainly through how it redistributes returned light, a freeform profile may change the redistribution without changing how much light escapes.

Here we introduce an angular mixing number, $N_{\mathrm{mix}} = 1/(1-|\lambda_2|)$, defined by the second eigenvalue of a structure's return matrix. It is the per-pass counterpart of the transport factor $1/(1-g)$ of radiative transfer [32] and of the relaxation time of a reversible Markov chain [33]. We impose reciprocity on the mixing of a reduced description and show that $p(\theta)$ and $N_{\mathrm{mix}}$ then suffice to predict multi-pass extraction for planar, scattering, microlens, random-lens, and wavy structures, and that a second spectral quantity, the escape probability of the quasi-stationary distribution [34,35], identifies light trapped at dead angles. Both spectral quantities, and the reciprocity properties that make them real and bounded, are standard objects of Markov-chain theory, of open dynamical systems, where the escape rate and the conditionally invariant measure describe the long-time survivors of chaotic billiards and dielectric microcavities [36–39], and of light-transport theory [40,41]; the reciprocal coupling of a diffuse lobe to angle-dependent reflection is likewise standard in reflectance modelling [42,43]. What is new here is their use as structure metrics for passive light extraction, the reciprocal reduced model of multi-pass extraction built on them, and the resulting map of OLED outcoupling structures, which we connect to the radiance-ordered extraction bound [44–46]. As a case study, we use a controlled freeform-versus-hemisphere benchmark on an experimentally validated OLED platform [15,16] to show that lens-shape freedom saturates as the map predicts.

---

## 2. Theory

### 2.1 Multi-pass extraction as a Markov chain on angle

The substrate angle is divided into 90 polar bins with edges $\theta_k = k°$ ($k = 0, \dots, 90$). The state is a vector $\mathbf{v}$ whose element $v_k$ is the power travelling towards the structure in bin $k$, integrated over the bin and over the full azimuth. The étendue (projected solid angle) of bin $k$ is proportional to $w_k = \sin^2\theta_{k+1}-\sin^2\theta_k$, so uniform radiance corresponds to $v_k \propto w_k$. On each encounter with the structure a fraction $p_k = p(\theta_k)$ of the power in bin $k$ escapes into air. The remainder returns into the substrate, and $B_R(i\to r)$ is the fraction of the power incident in bin $i$ that returns in bin $r$. Hence $\sum_r B_R(i\to r) = 1-p_i$, and $B_R$ factorizes as

$$
B_R = M\,\mathrm{diag}(1-\mathbf{p}),
$$

where $M(i\to r)$ is column-stochastic: it is the distribution of the returned angle conditioned on return. $M$ is therefore the transition matrix of a Markov chain on the polar angle of the light that is returned; escape is not an explicit absorbing state of $M$ but is carried by the survival factor $1-p_i$ in $B_R$. The escape probability for uniform radiance is $p = \sum_k w_k p_k / \sum_k w_k$.

The device reflects the returned light with its angle-resolved reflectance $R_k = R(\theta_k)$, without changing its angle. The power reaching the structure after $m$ round trips is $\mathbf{v}^{(m)} = (\mathrm{diag}\,\mathbf{R}\; B_R)^m\,\mathbf{v}^{(0)}$, where $\mathbf{v}^{(0)}$ is the power that the OLED launches into the substrate, and the substrate-to-air extraction efficiency is

$$
\eta_{\mathrm{sta}} = \mathbf{p}^{\mathsf T}\left(I - \mathrm{diag}\,\mathbf{R}\; B_R\right)^{-1}\mathbf{v}^{(0)} ,
$$

the matrix form of the trans-scale model [18]. All losses enter through $R(\theta)$. The scalar round-trip loss $A' = 1-\sum_k w_k R_k/\sum_k w_k$, the uniform-radiance average of $1-R(\theta)$, is a derived quantity used only in the memoryless comparison formula $\eta = p/[p+(1-p)A']$ [18]; it is not an additional loss.

The transport model is scalar. $R(\theta)$ is the polarization average $(|r_s|^2+|r_p|^2)/2$ of the stack reflectance, $\mathbf{v}^{(0)}$ is the sum of the s- and p-polarized power launched into the substrate, and the ray tracers apply the unpolarized Fresnel coefficients $(R_s+R_p)/2$ at every interface encounter without tracking polarization. Textures that are not axisymmetric in a single realization (hexagonal arrays, random supercells, wavy height fields) are reduced to a polar-angle kernel by averaging: rays are launched at the centre angle of each incidence bin, uniformly over the unit cell or supercell and uniformly in azimuth, and their exit directions are binned by polar angle only. $B_R$ is therefore the position- and azimuth-averaged return operator. $\eta_{\mathrm{sta}}$ is computed by a direct linear solve on the full 0–90° grid. Over all cases the spectral radius of $\mathrm{diag}\,\mathbf{R}\,B_R$ stays below 0.9991, the solve residual below $10^{-16}$, and a truncated pass-by-pass sum agrees with the solve to $10^{-14}$.

### 2.2 Reciprocity, detailed balance, and the angular mixing number

For a passive, reciprocal structure, the power scattered from bin $i$ into bin $j$ and from $j$ into $i$ obey $w_i\,B_R(i\to j) = w_j\,B_R(j\to i)$, where $w_k = \sin^2\theta_{k+1}-\sin^2\theta_k$ is the étendue of bin $k$; this follows from Helmholtz reciprocity of the bidirectional scattering distribution [47] and has been noted as a symmetry of the redistribution matrices of textured sheets [25]. It is the discrete form of the self-adjointness of reciprocal scattering operators in the projected-solid-angle (étendue) inner product [40], whose spectra are therefore real [41]; eigenmode expansions of interreflection operators are classical in radiometry [48]. This is detailed balance with respect to the étendue weight. $M$ is therefore similar to a symmetric matrix, its eigenvalues are real, and its stationary distribution is $\propto w(1-p)$, the distribution of returned light under uniform illumination. The largest eigenvalue is one. The eigenvalue of second-largest modulus, $\lambda_2$, sets how fast any other initial distribution relaxes towards the stationary one: the deviation shrinks by $|\lambda_2|$ per pass. We define the angular mixing number

$$
N_{\mathrm{mix}} = \frac{1}{1-|\lambda_2|},
$$

the number of passes the structure needs to forget the angle of the incident light. Using the second-largest eigenvalue modulus covers oscillatory relaxation; for every structure examined here $\lambda_2$ is real and positive (0.09–0.88 for $\theta < 85°$), so $|\lambda_2| = \lambda_2$. A planar interface has $\lambda_2 = 1$ and $N_{\mathrm{mix}} = \infty$; a Lambertian scatterer has $\lambda_2 = 0$ and $N_{\mathrm{mix}} = 1$. $N_{\mathrm{mix}}$ is the relaxation time of the chain [33] and the per-pass analogue of the $1/(1-g)$ scattering events that a photon needs to lose its direction in radiative transfer [32]. Because states at near-grazing angles carry almost no étendue but can return light at the same angle, we evaluate the spectral metrics on $\theta < 85°$, which excludes 0.76% of the étendue and less than 0.1% of $\mathbf{v}^{(0)}$. Concretely, $B_R$ is restricted to the bins below 85°; $M$ is formed by renormalizing the columns of this sub-block, and $\rho$ in Section 2.4 is taken from the sub-block itself. The efficiencies $\eta_{\mathrm{sta}}$, in the full series and in the reduced model, are always computed on the full 0–90° grid. The absolute values of the metrics depend on this cutoff, but their ranking across structures does not: between 80° and 90° the rank correlation with the 85° values is 0.96–0.99 for $N_{\mathrm{mix}}$ and 0.90–0.98 for $\xi$ (Note S12). The traced matrices obey detailed balance only within Monte Carlo noise (Section 3.2). Enforcing it exactly, by symmetrizing $w_i B_R(i\to j)$ while keeping each column sum, changes $N_{\mathrm{mix}}$ by at most 2.7%, $\xi$ by at most 0.005, and $\eta_{\mathrm{sta}}$ by at most 0.04 pp (Note S12).

### 2.3 Reduced model and the reciprocity constraint

The reduced model keeps the angle-resolved escape probability $p(\theta)$ and replaces $B_R$ by a single mixing parameter $m$: a fraction $1-m$ of the returned light keeps its angle, and a fraction $m$ is redistributed with a fixed distribution $q$,

$$
B_R(\theta_i\to\theta_r) = \left[1-p(\theta_i)\right]\left[(1-m)\,\delta_{ir} + m\,q(\theta_r)\right].
$$

As in reciprocal, energy-conserving reflectance models that couple a diffuse lobe to angle-dependent specular reflection [42,43], reciprocity fixes $q$. Detailed balance, $w_i B_R(i\to r) = w_r B_R(r\to i)$, holds for the mixing term only if $q(\theta) \propto w(\theta)\,[1-p(\theta)]$: light is redistributed into each angle in proportion to how much light that angle itself returns. To our knowledge this constraint has not been imposed in multi-pass light-extraction or light-trapping models. Reduced models usually redistribute returned light isotropically, as in the photon-gas picture [21–23] and in geometrical OLED models with isotropic scattering [49]. Inside a separate scattering layer this is reciprocal, but applied to the combined return of a structure whose escape depends on angle, the choice $q \propto w$ violates reciprocity, because it sends returned light into the escape cone of a structure that returns almost nothing from there. With the reciprocal $q$, the normalized kernel has eigenvalues 1 and $1-m$, so we set $m = 1-\lambda_2$ and the reduced model relaxes at the same rate as the full matrix. With $m = 1$, $q \propto w$, and a scalar $p$, the model reduces to the memoryless formula of ref. [18]; with $m = 0$ every ray returns at its own angle.

### 2.4 Late-pass escape probability and escape ratio

$\lambda_2$ describes how fast the angle of the returned light is forgotten, not where the light ends up. The second quantity we use is the largest eigenvalue $\rho$ of the sub-stochastic return operator $B_R$ itself. Light that has circulated many times without escaping approaches the corresponding eigenvector, the quasi-stationary distribution of the chain [34,35], and escapes with probability

$$
p_\infty = 1-\rho(B_R)
$$

per encounter. In the language of open dynamical systems, $p_\infty$ is the escape rate per encounter and the eigenvector is the conditionally invariant measure [36]. The same construction describes the long-time survival of rays in dielectric microcavities [37], in systems with partial absorption [38], and in systems with partial, Fresnel-type escape [39], to which the device reflectance $R(\theta)$ adds a further, angle-dependent loss; that late-order escape is set by the spectral radius of the transport recursion, with a limiting angular distribution, is also known in radiative transfer [50,51]. Because $B_R$ is self-adjoint in the étendue-weighted inner product, the Rayleigh quotient with the uniform-radiance vector gives $\rho \ge 1-\tilde p$, and the column sums give $\rho \le \max_\theta\,[1-p(\theta)]$, the angular counterpart of the row-sum bounds on the principal eigenvalue of interreflection kernels [52]. Hence

$$
\min_\theta p(\theta) \le p_\infty \le \tilde p , \qquad \tilde p = 1-\frac{\sum_{i} w_i \sum_{r} B_R(i\to r)}{\sum_i w_i},
$$

with both sums over the bins below the cutoff. $\tilde p$ is the probability that uniform-radiance light leaves the truncated chain on one encounter, by escape or by return beyond the cutoff; it exceeds $p$ on the same range by less than 0.008 for every structure. Numerically, all 37 non-planar structures satisfy both bounds, with $p_\infty \le \tilde p - 0.0014$ and $p_\infty \ge \min_\theta p(\theta)$.

We call $\xi = p_\infty/\tilde p$ the late-pass escape ratio. It equals one when long-lived light escapes as easily as Lambertian light, and zero when some angles never escape and never mix (planar interface). For a Lambertian integrating sphere the redistribution has rank one [48], so that $N_{\mathrm{mix}} = 1$ and $\xi = 1$, which recovers the memoryless sphere multiplier. $\xi$ is a structure-only descriptor. In a device, the late survivors are shaped by $\mathrm{diag}\,\mathbf{R}\,B_R$; replacing $\rho(B_R)$ by $\rho(\mathrm{diag}\,\mathbf{R}\,B_R)$ for the Al device at 550 nm leaves the ranking of the structures nearly unchanged (rank correlation 0.96).

### 2.5 Band-resolved radiance-ordered bound

The power that a passive structure can send into a set of output directions in one encounter is bounded by filling the corresponding output étendue with the brightest internal radiance first [44–46]. For a uniform internal radiance the bound reduces to the étendue ratio, $1/n^2$ for the full hemisphere [21] and $(\sin^2\theta_{\mathrm{hi}}-\sin^2\theta_{\mathrm{lo}})/n^2$ for a polar band in air. An OLED microcavity, however, launches an anisotropic radiance into the substrate. We resolve the bound into the four polar bands used below, with the actual first-pass radiance of the device, and use it to rate how much of the available phase space each structure uses. In discrete form, the substrate bins are sorted by their radiance $v^{(0)}_k/w_k$; the étendue of an air band, mapped into the substrate by Snell's law, is $(\sin^2\theta_{\mathrm{hi}}-\sin^2\theta_{\mathrm{lo}})/n^2$ in units of $w$; and this étendue is filled with the brightest bins first, the last bin fractionally. The bound assumes lossless, polarization-independent coupling (no Fresnel loss) and is integrated over the full azimuth. Each band's bound is computed on its own, so different bands may claim the same bright substrate directions, and the per-band bounds do not add up to the bound for the whole hemisphere, which fills the étendue $1/n^2$ once.

---

## 3. Methods

### 3.1 Devices

Dipole emission and the device reflectance are computed with the Chance–Prock–Silbey model [53] and transfer matrices, using measured optical constants [15,18]. Birefringent layers are treated as uniaxial with the optic axis along the normal, using the anisotropic extension of the model [54,55] as implemented for OLEDs [56]. The substrate is treated incoherently, as is standard in OLED outcoupling models [57]. Three bottom-emitting stacks are used. The first is the platform of the freeform benchmark at its optimized cavity: Al (100 nm) / B3PyMPM (61.3 nm) / TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd (25 nm) / TCTA (10 nm) / TAPC (97.5 nm) / ITO (150 nm) / glass ($n = 1.51$). The second replaces Al by Ag. The third is a low-loss stack, Ag (100 nm) / B3PyMPM (200 nm) / EML (25 nm) / TAPC (180 nm) / ITO (50 nm) / glass [18]. In all three the emitters form a single dipole sheet at the centre of the birefringent EML, with a horizontal dipole ratio of 0.865. At 550 nm the uniform-radiance round-trip loss is $A'$ = 0.167 (Al), 0.044 (Ag), and 0.026 (low loss). The full layer list and $R(\theta)$ and $\mathbf{v}^{(0)}(\theta)$ at 550 nm are given in Table S14 and Fig. S11. The angular distribution of the light delivered to the substrate and the device reflectance are computed per wavelength.

### 3.2 Structures and their return matrices

All structures are index-matched to the substrate, so the structure/air boundary sets the escape condition. Their transmission and return matrices are computed on 1° polar bins (microlenses and wavy surfaces) or 0.5° bins (planar interface and scattering layers); the latter are merged into 1° bins, assuming uniform radiance within each 1° bin, so that all structures are evaluated on the same grid.

- **Planar interface.** Fresnel reflection, returned at the incident angle.
- **Scattering layers.** Index-matched, non-absorbing slabs with a Henyey–Greenstein phase function [58] (asymmetry $g$ = 0–0.9) and transport optical thickness $\tau^* = S(1-g)$ = 0.25–4, combined with the flat exit face by the adding method [59], the interface entering as a layer with its own Fresnel reflection and transmission operators [60]; Monte Carlo with 20,000 rays per incidence bin.
- **Hexagonal MLAs.** Close-packed axisymmetric lenslets of the spline family of the benchmark (Table 1), traced with a height-field Monte-Carlo tracer, as in OLED substrate-extraction [2] and silicon-texture [61] ray tracing, that reproduces the LightTools return matrix of the hemispherical MLA (escape probability 0.4275 against 0.4274; Note S11). The set includes the hemispherical and freeform optima of the benchmark.
- **Random MLAs.** A 6 × 6 supercell in which every lenslet carries its own profile (shape spread $s$ about the quarter circle), its own height (aspect spread), or a jittered position and radius; overlapping lenslets take the higher surface.
- **Wavy surfaces.** Periodic Gaussian-correlated random height fields with rms slope 0.2–1.5. Ray optics has no length scale, so the results are invariant under uniform scaling of the height field; for a fixed correlation shape they depend on the height and correlation length only through the rms slope. Such surfaces resemble the spontaneously formed buckles that have been used to extract light from OLEDs [62].

The ray treatment requires radii of curvature and correlation lengths well above the wavelength [63,64]; for microlenses, full-wave simulations converge to ray tracing for lens diameters above about 5 µm (≈10λ) [9]. The random structures (22 in total) are listed with their parameters and results in Table S12. All return matrices are traced with 20,000 rays per incidence bin (10,000 for wavy surfaces); the fraction of rays lost to the iteration limit is below $6\times10^{-4}$ for every structure. Monte Carlo uncertainty is small: re-tracing with other ray seeds changes $N_{\mathrm{mix}}$ by at most 0.01 and $\eta_{\mathrm{sta}}$ by at most 0.06 pp, and drawing other geometric realizations with the same statistics changes $N_{\mathrm{mix}}$ by at most 0.24, $\xi$ by at most 0.021, and $\eta_{\mathrm{sta}}$ by at most 0.23 pp (Table S13). The metrics are converged with respect to the discretization: halving or quadrupling the number of rays changes $N_{\mathrm{mix}}$ by at most 0.05 and $\xi$ by at most 0.005; coarsening the bins to 2° or 3° changes $\eta_{\mathrm{sta}}$ by at most 0.3 pp (planar interface) and by at most 0.06 pp for all other structures; and enlarging the periodic domain to a 10 × 10 supercell or a 24 × 24 wavy box changes $N_{\mathrm{mix}}$ by at most 0.06, $\xi$ by at most 0.025, and $\eta_{\mathrm{sta}}$ by at most 0.34 pp, comparable to the realization spread (Table S15). The traced return matrices satisfy detailed balance within Monte Carlo noise: the relative asymmetry of $w_i B_R(i\to j)$ is 5–13% on 1° bins and falls to 0.2–0.7% on 15° bins (Note S10).

### 3.3 Freeform benchmark

The case study of Section 4.4 compares an axisymmetric freeform MLA with an equally optimized hemispherical MLA on the Al platform. LightTools ray tracing is driven from MATLAB, with the CPS source on a 1 mm-radius disc and a 1.295 mm glass substrate. The pipeline has been validated against fabricated devices with nearly the same geometry, including attached and inverted MLAs [15,16]. Both lens classes share the material, pitch, height range, and cavity bounds (Table 1); the hemisphere is a point of the freeform design space. Because the ETL and HTL thicknesses are optimized together with the lens, Note S9 repeats the comparison in a simplified multi-pass model at a common cavity, with the stack fixed at the hemispherical total-EQE optimum and only the lens changed. Total EQE and the EQE in four polar bands (0–20°, 20–40°, 40–60°, 60–80°) are optimized as separate objectives by a radial-basis-function surrogate search (MATLAB surrogateopt [65,66]) followed by generalized pattern search (MATLAB patternsearch [67,68]). Here EQE is the optical outcoupling: the radiative efficiency (0.98) times the fraction of generated photons that reach air through the textured substrate, weighted by the emission spectrum over 151 wavelengths (453–753 nm) and integrated over the full air hemisphere (total EQE) or over a polar band and the full azimuth (band EQE). It thus includes the cavity's coupling into the substrate and is not directly comparable with $\eta_{\mathrm{sta}}$, which starts from the light already in the substrate. The searches use 10,000 rays and 31 wavelengths per evaluation, so their objective is noisy; every optimum is then re-evaluated three times, in independent Monte Carlo runs, with 50,000 rays and 151 wavelengths, which gives a spread of about $10^{-4}$ in absolute EQE. A restart control repeats every search from the hemispherical optimum. The gain over the re-measured hemisphere is tested with a pooled two-sample $t$ statistic, $t = \Delta/(s_p\sqrt{2/n})$ with $n = 3$ re-evaluations per design and four degrees of freedom (one-sided 95% threshold 2.13), with Holm correction [69] for the five objectives (Note S2).

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

**Fig. 1 | Multi-pass extraction as a Markov chain on angle.** (a) Light delivered to the substrate meets the passive structure, which transmits a fraction $p(\theta)$ and returns the rest with the redistribution $B_R(\theta_i\to\theta_r)$; the device reflects the returned light with its angle-resolved reflectance $R(\theta)$, which carries all losses of the full model (the scalar $A'$ is used only in the memoryless comparison). (b–d) Column-normalized return matrices $M$, shown as $\log_{10}M$ per 1° bin with a floor at $10^{-3}$ (shared colour bar), of a planar interface, the hemispherical MLA, and a scattering layer ($g$ = 0.8, $\tau^*$ = 1), with $\lambda_2$ and $N_{\mathrm{mix}}$. *(`fig1_mixing_concept.png`)*

---

## 4. Results

### 4.1 Two parameters predict multi-pass extraction

Figure 2 compares four models with the full angle-resolved series for 342 combinations: three devices, three wavelengths (500, 550, and 600 nm), and 38 structures (planar, six scattering layers, nine hemispherical and freeform MLAs, 17 random MLAs, and five wavy surfaces). The memoryless formula fails for structures that retain the angle (Fig. 2a). It overestimates the planar interface by up to 60 pp, because it lets returned light reach the escape cone that a flat face never lets it reach, and the hemispherical MLA by 2–3 pp; over all cases its mean error is 3.9 pp. Assuming that every ray returns at its own angle fails in the opposite direction for scattering layers (Fig. 2b).

With $m = 1-\lambda_2$ and Lambertian redistribution ($q \propto w$), the reduced model is accurate for lens arrays but overestimates the shallow wavy surfaces by up to 10.5 pp (Fig. 2c). Their $\lambda_2$ is low, yet light trapped beyond the critical angle is returned mostly to other trapped angles: on the slope-0.2 surface, only 3–8% of the light returned from 70–80° lands at angles with $p(\theta) > 0.3$, against 50% for Lambertian redistribution. The reciprocal redistribution $q \propto w(1-p)$ of Section 2.3 builds this in, since angles that return much light also receive much of the mixed light. It reduces the mean error to 0.30 pp and the maximum to 3.0 pp over all 342 cases, with no additional parameter (Fig. 2d, Table S11). The mean error is 0.29 pp for scattering layers, 0.41 pp for hemispherical and freeform MLAs, 0.16 pp for random MLAs, 0.19 pp for wavy surfaces of rms slope ≥ 0.7, and 1.4 pp for the two shallowest wavy surfaces, which remain the hardest case. Stratified by device, the mean (maximum) error is 0.41 (2.7) pp for Al, 0.33 (3.0) pp for Ag, and 0.16 (1.8) pp for the low-loss stack; by wavelength it is 0.38 (2.5), 0.29 (3.0), and 0.23 (2.6) pp at 500, 550, and 600 nm (Fig. 2e, Table S11). The reduced model tends to slightly underestimate $\eta_{\mathrm{sta}}$, and no device or wavelength concentrates the error.

A structure is therefore described, for the purpose of multi-pass extraction, by its angle-resolved escape probability and one number for how fast it forgets the angle, provided that the mixing respects reciprocity.

**Fig. 2 | Validation of the reduced model.** Extraction efficiency $\eta_{\mathrm{sta}}$ predicted by (a) the memoryless formula $\eta = p/[p+(1-p)A']$, (b) the angle-preserving limit ($m = 0$), and the reduced model with $m = 1-\lambda_2$ and (c) Lambertian redistribution $q \propto w$ or (d) reciprocal redistribution $q \propto w(1-p)$, against the full angle-resolved series. Symbols: planar (squares), scattering layers (circles), hemispherical (diamonds), freeform (triangles), and random (inverted triangles) MLAs, and wavy surfaces (crosses); colour and symbol encode the structure class; each class appears for the Al, Ag, and low-loss devices at 500, 550, and 600 nm (342 cases). The outliers in (c) are the wavy surfaces of rms slope 0.2 and 0.4. (e) Error of the reciprocal reduced model for each device and wavelength; bars mark the mean absolute error. *(`fig2_reduced_model.png`)*

### 4.2 The structure map

Figure 3a places all structures on the $(p, N_{\mathrm{mix}})$ plane, with contours of $\eta_{\mathrm{sta}}$ for the Al device. The classes separate:

- **Planar interface:** $p$ = 0.40 and $N_{\mathrm{mix}} = \infty$. Light outside the escape cone never reaches it, and $\eta_{\mathrm{sta}}$ = 0.38.
- **Scattering layers:** $N_{\mathrm{mix}}$ = 1.1–2.3, but $p$ = 0.19–0.39, falling with thickness because more light is back-scattered before it reaches the exit face. $\eta_{\mathrm{sta}}$ = 0.58–0.75.
- **Hemispherical MLAs:** the highest escape probability, $p$ = 0.42, but $N_{\mathrm{mix}}$ = 6–8. $\eta_{\mathrm{sta}}$ = 0.78–0.79.
- **Freeform MLAs:** $p$ = 0.42–0.43, like the hemisphere, with $N_{\mathrm{mix}}$ = 3–7.5.

Randomness lowers $N_{\mathrm{mix}}$ without lowering $p$ (Fig. 3b, Table S12). Drawing every lenslet profile independently reduces $N_{\mathrm{mix}}$ from 7.0 to 1.9 ± 0.1 (spread over geometric realizations) while $p$ stays at 0.42–0.43. Jittering the lenslet positions by up to 30% of the pitch reduces it to 3.5 ± 0.1, and adding a 20% radius spread to 3.2. Varying only the lenslet heights is the least effective (7.1 → 4.1 ± 0.2 at ±60%). A wavy surface of rms slope 0.4–1.5 has $\lambda_2$ = 0.12–0.18 ($N_{\mathrm{mix}}$ = 1.1–1.2), as low as a scattering layer, and its escape probability rises with slope from 0.41 to 0.43, reaching that of the hemispherical MLA at slope 1.0. At slope 1.5 the wavy surface gives the highest $\eta_{\mathrm{sta}}$ of all structures on the Al device, 0.801, against 0.792 for the hemispherical MLA.

The late-pass escape ratio adds what $N_{\mathrm{mix}}$ cannot see (Fig. 3c). For lens arrays the two are tied: $\xi$ falls from 0.80 to 0.39 as $N_{\mathrm{mix}}$ rises from 2 to 8, because light that keeps its angle for many passes ends up at the angles the lenses transmit worst. The hemispherical MLA has $\xi$ = 0.47: light that survives many passes escapes at less than half the rate of fresh Lambertian light. Randomizing the lenslets raises $\xi$ to 0.7–0.8, wavy surfaces of rms slope ≥ 0.7 reach 0.81–0.88, and thick scattering layers 0.93–0.99. Two structures fall well below the lens-array trend at low $N_{\mathrm{mix}}$: the wavy surface of slope 0.2 ($N_{\mathrm{mix}}$ = 2.2, $\xi$ = 0.35) and the thin isotropic scattering layer ($N_{\mathrm{mix}}$ = 2.2, $\xi$ = 0.54). Both mix quickly but locally, so trapped light is recycled among trapped angles. These are the structures with the largest residual errors in Fig. 2d. $\xi$ thus flags the cases in which a single mixing rate is least reliable. Its absolute value depends on how the near-grazing states are treated: for the hemispherical MLA it is 0.65, 0.47, and 0.31 for cutoffs of 80°, 85°, and 90°, whereas the ordering of the structures is preserved (Section 2.2).

The contours of Fig. 3a show what lower $N_{\mathrm{mix}}$ and higher $\xi$ are worth, and the answer is: little, once $p$ is high. The contours run almost vertically, so $\eta_{\mathrm{sta}}$ is set mainly by $p$. They are computed with the angular shape of $p(\theta)$ of the hemispherical MLA. Repeating the calculation with the $p(\theta)$ shapes of a scattering layer and of a wavy surface gives the same picture: at fixed $p$, varying $N_{\mathrm{mix}}$ from 1 to 16 changes $\eta_{\mathrm{sta}}$ by 2.1 pp (lens template), 1.2–1.4 pp (scattering layer), and 3.2 pp (wavy surface), whereas raising $p$ by 0.01 changes it by about 0.6 pp for all three templates (Note S12). Over the whole range of $N_{\mathrm{mix}}$, mixing is therefore worth a change of $p$ by 0.02–0.05. At fixed $p(\theta)$, the span available to mixing lies between full angle retention ($m = 0$) and complete reciprocal mixing ($m = 1$). For the hemispherical MLA on the Al device this span is 0.780–0.806, and the real structure sits at 0.792. Randomizing the lenslet profiles moves $N_{\mathrm{mix}}$ from 7.0 to 1.9 but changes $\eta_{\mathrm{sta}}$ by less than 1 pp (0.794 → 0.785), because the same randomness also reshapes $p(\theta)$. A counterfactual decomposition in the reduced model separates the two effects (Note S12). For $s$ = 1 on the Al device, changing only the mixing ($m$ from the randomized structure, $p(\theta)$ from the hemisphere) gains 0.91 pp, changing only $p(\theta)$ loses 2.20 pp, and their interaction adds 0.30 pp, for a net −0.98 pp against −0.91 pp in the full series. On the Ag device the two effects cancel (+0.91 and −0.94 pp), and on the low-loss device both are below 0.12 pp. The planar interface makes the reciprocity constraint explicit. Even complete mixing at its own $p(\theta)$ raises $\eta_{\mathrm{sta}}$ only from 0.38 to 0.49, because a face that transmits nearly all light inside the escape cone also returns almost none from there, and reciprocity sends the mixed light back in the same proportion. Mixing cannot open angles at which a structure does not transmit; only a change of $p(\theta)$ can (Section 5).

**Fig. 3 | The structure map.** (a) Escape probability for Lambertian incidence and angular mixing number of planar, scattering, hemispherical, freeform, random-lens, and wavy structures, with contours of $\eta_{\mathrm{sta}}$ from the reciprocal reduced model for the Al device at 550 nm, computed with the angular shape of $p(\theta)$ of the hemispherical MLA scaled to each $p$. The vertical axis is logarithmic. The hemispherical and freeform MLAs are the benchmark optima of Section 4.4 (Table 3). (b) $N_{\mathrm{mix}}$ (logarithmic axis) against each randomness parameter normalized by its largest value: lenslet shape spread $s$ (0–1), position jitter as a fraction of the pitch (0–0.3), lenslet height spread (0–±0.6), and rms slope of a wavy surface (0–1.5); the zero point of the jitter and height series is the hemispherical array; the blue dotted line marks the hemispherical MLA and the green dotted line a scattering layer ($g$ = 0.8, $\tau^*$ = 1). (c) Late-pass escape ratio $\xi = p_\infty/\tilde p$ against $N_{\mathrm{mix}}$ (logarithmic axis); the planar interface ($N_{\mathrm{mix}} = \infty$, $\xi$ = 0) is placed off-scale at the right edge. *(`fig3_structure_map.png`)*

### 4.3 How much of the radiance bound passive structures use

On the first pass, the microcavity delivers light to the substrate with a radiance that is far from uniform. For the Al platform at 550 nm, filling each air band with the brightest internal radiance allows 0.085, 0.205, 0.230, and 0.156 of the substrate power into the 0–20°, 20–40°, 40–60°, and 60–80° bands (Fig. 4a). These per-band bounds are separate maxima: each band may draw on the same bright substrate directions, so they cannot be reached simultaneously and their sum (0.676) exceeds the bound for the whole air hemisphere, 0.540. The uniform-radiance bound is additive: 0.051, 0.130, 0.148, and 0.096 in the four bands, plus 0.014 above 80°, give $1/n^2$ = 0.439. Every passive structure we examined stays well below the band-resolved bound. The MLAs reach 76–80% of it in total and 55–70% per band, a wavy surface of rms slope 1.0 77%, the planar interface 65%, and the scattering layers 35–68%. The lenses transmit almost exactly as much as a structure facing uniform radiance would (0.43 against 0.44), so the anisotropy that the cavity provides is not used.

The output composition, finally, is nearly independent of the structure type in the lossy devices (Fig. 4b). We define the selectivity $S_{40\text{–}60} = P_{\mathrm{air}}(40\text{–}60°)/P_{\mathrm{air}}(0\text{–}90°)$ of the multi-pass output, both powers in air and including all Fresnel losses. On the Al device at 550 nm, over planar, scattering, hemispherical, freeform, random, and wavy structures, it stays between 0.327 and 0.360, close to the Lambertian partition of 0.337; at 500 and 600 nm the window is 0.320–0.364 and 0.339–0.361, and on the Ag device 0.322–0.361. The low-loss device is different. Its window widens to 0.274–0.366, and the lower end is set by the hemispherical MLAs, whose long angle retention ($\xi \approx 0.45$) lets late passes, shaped by the structure, dominate the output when little is lost per round trip. Within the class of angle-insensitive passive structures, the angular composition of the output is therefore set by the platform in lossy devices, while in low-loss devices angle-retaining structures can shift it.

**Fig. 4 | Radiance bound and composition window.** (a) First-pass power into each polar band for the Al platform at 550 nm: uniform-radiance bound, radiance-ordered bound for the actual cavity emission, hemispherical MLA, wavy surface (rms slope 1.0), scattering layer ($g$ = 0.8, $\tau^*$ = 1), and planar interface. First-pass power is $\sum_k v^{(0)}_k T_k(\text{band})$, the power that reaches the band in air on the first encounter, as a fraction of the power launched into the substrate. (b) Selectivity $S_{40\text{–}60} = P_{\mathrm{air}}(40\text{–}60°)/P_{\mathrm{air}}(0\text{–}90°)$ of the multi-pass output against $\eta_{\mathrm{sta}}$ for all structures (Al, 550 nm); the band marks 0.33–0.36 and the dotted line the Lambertian partition. Windows for the other devices and wavelengths are given in Section 4.3 and Note S12. *(`fig4_bound_window.png`)*

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

The benchmark optima can be placed on the map directly, because the hemispherical and freeform MLAs of Fig. 3 are these designs, rebuilt from their optimized control points and traced with the height-field tracer (Table 3). Their escape probabilities lie within $p$ = 0.417–0.429, a range worth about 0.7 pp of $\eta_{\mathrm{sta}}$, while $N_{\mathrm{mix}}$ spans 3.0–8.0 and $\xi$ 0.39–0.79. Every optimized freeform lens thus has the same escape probability as the hemisphere and differs from it only in $N_{\mathrm{mix}}$ and in the shape of $p(\theta)$.

**Table 3 | Benchmark optima on the structure map (550 nm).** $\eta_{\mathrm{sta}}$ is for the Al device in the multi-pass model.

| Design | $p$ | $N_{\mathrm{mix}}$ | $\xi$ | $\eta_{\mathrm{sta}}$ |
|---|---|---|---|---|
| hemisphere, 0–20° optimum | 0.423 | 7.3 | 0.44 | 0.790 |
| hemisphere, 20–40° optimum | 0.422 | 7.0 | 0.47 | 0.793 |
| hemisphere, 40–60° optimum | 0.425 | 8.0 | 0.42 | 0.783 |
| hemisphere, 60–80° optimum | 0.420 | 6.2 | 0.50 | 0.792 |
| hemisphere, total-EQE optimum | 0.422 | 7.1 | 0.47 | 0.792 |
| freeform, 20–40° optimum | 0.429 | 3.8 | 0.58 | 0.788 |
| freeform, 40–60° optimum | 0.425 | 7.5 | 0.39 | 0.789 |
| freeform, 60–80° optimum | 0.417 | 3.0 | 0.79 | 0.785 |
| freeform, total-EQE optimum | 0.424 | 4.6 | 0.46 | 0.771 |

The 0–20° freeform optimum ends in a vertical sidewall that a height-field tracer cannot represent and is omitted. Shape freedom thus moves a structure vertically on the map, where the $\eta_{\mathrm{sta}}$ contours are flat, and cannot move it to the right, where they are steep. The band selectivities of all designs follow the same picture: across 606 sampled designs the 40–60° band EQE is proportional to the total EQE ($R^2$ = 0.968), and band-dedicated optimization raises the band power by at most a factor of 1.22 net of the loss in total EQE (Note S4). The same collapse holds for inverted and randomly assembled arrays (Note S7). Because the benchmark optimizes the cavity together with the lens, part of the band steering could come from the stack. With the stack fixed at the hemispherical total-EQE optimum, a simplified version of the multi-pass model (non-dispersive organic index 1.8; Note S9) gives a spread of only 2.5% (Al) and 1.7% (Ag) in $\eta_{\mathrm{sta}}$ among the nine optimized lenses, and almost all of the band steering of the dedicated designs comes from the cavity: at the common cavity the 60–80° share of the hemispherical 60–80° design falls from 0.303 to 0.249, close to the 0.247 of the total-EQE design (Note S9). A fixed-cavity LightTools control at full fidelity remains to be done.

**Fig. 5 | Freeform versus an equally optimized hemisphere.** (a) Band-optimized freeform profiles against the hemispherical reference at its optimized height. (b) Best EQE of each class at each objective, labelled with $G_j$. (c) Restart control: gain over the hemisphere when each search restarts from the hemispherical optimum, labelled with $t$; green marks gains above the one-sided 95% threshold. *(`fig5_hemisphere_benchmark.png`)*

---

## 5. Discussion

**Why lens shape saturates.** The escape probability of a lens array for uniform radiance is close to the étendue bound $1/n^2$ already for a hemisphere, so profile changes cannot raise $p$ appreciably. What they can change is how the returned light is redistributed, that is, $N_{\mathrm{mix}}$, and the shape of $p(\theta)$. Along that direction the efficiency contours are flat. The anisotropy of the cavity emission would in principle allow more. Being coextensive with the source, a refractive film cannot enlarge the étendue available to a band, but it can select, through its local slopes, which internal directions are refracted into it, and a surface with binary angular transmission can reach the radiance-ordered bound [46]; skewness conservation does not forbid this for coextensive, rotationally symmetric transfer [70]. The rotationally symmetric microlenses and textures studied here, however, fill each band from a broad range of internal angles, so in our simulations they extract close to what a uniform internal radiance would give. Approaching the band-resolved bound thus appears to call for angle-selective transmission rather than further freedom in the lens profile.

**Why mixing matters little.** Reciprocity ties mixing to escape. A structure that sends trapped light into the escape cone must, by detailed balance, return light from the cone to trapped angles in the same étendue-weighted proportion, which lowers its escape probability for fresh light. This is the trade-off seen in scattering layers: they reach $\xi \approx 1$ but lose $p$ with thickness, so that extraction peaks at an intermediate scatterance or haze [24,27]. For a structure with a given $p(\theta)$, the span between full retention and complete reciprocal mixing is 1.8–4.6 pp on the Al device, 1.7–7.0 pp on the Ag device, and 0.1–1.7 pp on the low-loss device for the lens arrays we examined, and real structures already capture part of it. The random lens arrays show that the remainder is hard to collect: the randomness that lowers $N_{\mathrm{mix}}$ and raises $\xi$ also flattens $p(\theta)$, and the effects nearly cancel. Lowering $N_{\mathrm{mix}}$ further therefore buys at most a few percentage points, and only in lossy devices. In low-loss devices the choice of outcoupling structure, and its shape in particular, matters even less, in line with ref. [18].

**Design implications.** Three numbers are worth reporting for any new outcoupling structure: its escape probability for uniform radiance $p$, its angular mixing number $N_{\mathrm{mix}}$, and its late-pass escape ratio $\xi$. The first sets the efficiency, the second the accuracy of any memoryless description, and the third whether light that survives many passes is still extractable. A structure that combines a high $p$ with a low $N_{\mathrm{mix}}$ and a high $\xi$, such as a random MLA or a wavy surface of rms slope one or more, avoids both the low escape probability of thick scattering layers and the angle retention of periodic lenses, without introducing absorption. Its advantage in efficiency over a periodic MLA is small (≤ 1 pp); randomly placed or irregular microlenses have been used on OLEDs mainly to avoid diffraction and angular colour shift [71]. Optimizers, including machine-learning ones [12,13], can then be benchmarked against the $(p, N_{\mathrm{mix}})$ map and the band-resolved bound rather than against each other.

**Scope.** The analysis covers passive, angle-insensitive external structures index-matched to the substrate, in the ray-optical regime. Elements that select the angle of transmission, such as angular filters, are also passive and also obey the radiance bound per encounter, but in combination with recycling they can change the output composition; they are outside the present scope. Internal structures that act on waveguided modes, absorbing scatterers, and wavelength-converting layers are not covered either. As in earlier substrate models [2,24], the device is treated as an angle-preserving reflector; diffuse reflection from planar metal cathodes is small, and including it would add a non-diagonal factor to the return operator, whose product with $B_R$ need not keep a real spectrum [41]. The reduced model further assumes that $p(\theta)$ is azimuthally symmetric. Its largest residual errors (up to 3 pp) occur for shallow surfaces with local mixing, which the late-pass escape ratio identifies. The three wavelengths (500, 550, and 600 nm) bracket the green emission band; the conclusions hold at each of them (Table S11).

---

## 6. Conclusion

We characterized passive light-extraction structures by two quantities from the spectrum of their return matrix, an angular mixing number, $N_{\mathrm{mix}} = 1/(1-|\lambda_2|)$, and a late-pass escape ratio, $\xi = p_\infty/\tilde p$ with $p_\infty = 1-\rho(B_R)$. Reciprocity, a known property of passive scattering, makes the spectrum real, bounds the late-pass escape probability by the escape probability for uniform radiance, and fixes how a reduced model must redistribute the returned light. With this reciprocal redistribution, the angle-resolved escape probability and $N_{\mathrm{mix}}$ predict multi-pass extraction to within 0.30 pp on average and 3.0 pp at most across 38 planar, scattering, microlens, random-lens, and wavy structures, where the memoryless formula deviates by up to 60 pp. On the resulting map, periodic microlenses combine high escape probability with long angle retention and a low late-pass escape ratio, scattering layers the opposite, and random lens arrays and wavy surfaces combine the escape probability of the former with the mixing of the latter. Because reciprocity prevents mixing from opening angles that a structure does not transmit, the efficiency is set mainly by the escape probability, and this combination is worth less than 1 pp. Every passive structure uses only part of the band-resolved radiance bound, and in lossy devices angle-insensitive structures share one output composition, which angle-retaining lenses can shift only when the round-trip loss is small. Freeform optimization of microlens shape mainly changes $N_{\mathrm{mix}}$, consistent with its gain of at most 7% in any band over an equally optimized hemisphere on an extended OLED emitter; going further likely calls for angle-selective rather than shape-optimized structures.

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

15. Kim, M.; Kim, J.; Yoo, S. **Near-Planar Light Outcoupling Structures with Finite Lateral Dimensions for Ultra-Efficient and Optical Crosstalk-Free OLED Displays.** *Nat. Commun.* **2025**, *16*, 11606. https://doi.org/10.1038/s41467-025-66538-6.

16. Kim, J.; Kim, E.; Park, J.; Song, J.; Kim, S.; Moon, H.; Yoo, S. **Toward Near-Foldable Surface Light Sources with Ultimate Efficiency: Ultrathin Substrates Embedded with Micron-Scale Inverted Lens Arrays.** *ACS Photonics* **2023**, *10*, 1775–1782. https://doi.org/10.1021/acsphotonics.3c00017.

17. Li, W.; Li, Y.; Shen, Y.; Zhang, Y.; Jin, T.; Chen, J.; Zhang, X.; Tang, J. **Releasing the Trapped Light for Efficient Silver Nanowires‐Based White Flexible Organic Light‐Emitting Diodes.** *Adv. Opt. Mater.* **2019**, *7*, 1900985. https://doi.org/10.1002/adom.201900985.

18. Kim, J.; Lee, S.; Yoo, S.; Gather, M. C. **Minimizing Round-Trip Loss for Near-Unity Quantum Efficiency in Organic Light-Emitting Diodes.** *arXiv* **2026**, 2610.12065.

19. Kang, K.; Byeon, I.; Kim, Y. G.; Choi, J.; Kim, D. **Nanostructures in Organic Light‐Emitting Diodes: Principles and Recent Advances in the Light Extraction Strategy.** *Laser Photonics Rev.* **2024**, *18*, 2400547. https://doi.org/10.1002/lpor.202400547.

20. Shinar, R.; Shinar, J. **Light extraction from organic light emitting diodes (OLEDs).** *J. Phys. Photonics* **2022**, *4*, 032002. https://doi.org/10.1088/2515-7647/ac6ea4.

21. Yablonovitch, E. **Statistical Ray Optics.** *J. Opt. Soc. Am.* **1982**, *72*, 899–907.

22. Schnitzer, I.; Yablonovitch, E.; Caneau, C.; Gmitter, T. J.; Scherer, A. **30% External Quantum Efficiency from Surface Textured, Thin-Film Light-Emitting Diodes.** *Appl. Phys. Lett.* **1993**, *63*, 2174–2176. https://doi.org/10.1063/1.110575.

23. Green, M. A. **Lambertian light trapping in textured solar cells and light‐emitting diodes: analytical solutions.** *Prog. Photovolt.* **2002**, *10*, 235–241. https://doi.org/10.1002/pip.404.

24. Shiang, J. J.; Duggal, A. R. **Application of radiative transport theory to light extraction from organic light emitting diodes.** *J. Appl. Phys.* **2004**, *95*, 2880–2888. https://doi.org/10.1063/1.1644037.

25. Eisenlohr, J.; Tucher, N.; Höhn, O.; Hauser, H.; Peters, M.; Kiefel, P.; Goldschmidt, J. C.; Bläsi, B. **Matrix Formalism for Light Propagation and Absorption in Thick Textured Optical Sheets.** *Opt. Express* **2015**, *23*, A502–A518. https://doi.org/10.1364/OE.23.00A502.

26. David, A. **Surface-Roughened Light-Emitting Diodes: An Accurate Model.** *J. Display Technol.* **2013**, *9*, 301–316. https://doi.org/10.1109/JDT.2013.2240373.

27. Altazin, S.; Reynaud, C.; Mayer, U. M.; Lanz, T.; Lapagna, K.; Knaack, R.; Penninck, L.; Kirsch, C.; et al. **Simulations, Measurements, and Optimization of OLEDs with Scattering Layer.** *SID Symp. Dig. Tech. Pap.* **2015**, *46*, 564–567. https://doi.org/10.1002/sdtp.10478.

28. Kovačič, M.; Samigullina, D.; Bouchard, F.; Krč, J.; Lipovšek, B.; Soldera, M.; Lasagni, A. F.; Reineke, S.; et al. **Analysis and optimization of light outcoupling in OLEDs with external hierarchical textures.** *Opt. Express* **2021**, *29*, 23701. https://doi.org/10.1364/oe.428021.

29. Deckman, H. W.; Roxlo, C. B.; Yablonovitch, E. **Maximum Statistical Increase of Optical Absorption in Textured Semiconductor Films.** *Opt. Lett.* **1983**, *8*, 491–493. https://doi.org/10.1364/OL.8.000491.

30. Manzoor, S.; Filipič, M.; Onno, A.; Topič, M.; Holman, Z. C. **Visualizing Light Trapping within Textured Silicon Solar Cells.** *J. Appl. Phys.* **2020**, *127*, 063104. https://doi.org/10.1063/1.5131173.

31. Battaglia, C.; Boccard, M.; Haug, F.-J.; Ballif, C. **Light Trapping in Solar Cells: When Does a Lambertian Scatterer Scatter Lambertianly?** *J. Appl. Phys.* **2012**, *112*, 094504. https://doi.org/10.1063/1.4761988.

32. van Rossum, M. C. W.; Nieuwenhuizen, T. M. **Multiple Scattering of Classical Waves: Microscopy, Mesoscopy, and Diffusion.** *Rev. Mod. Phys.* **1999**, *71*, 313–371. https://doi.org/10.1103/RevModPhys.71.313.

33. Levin, D. A.; Peres, Y. **Markov Chains and Mixing Times**, 2nd ed.; American Mathematical Society: Providence, RI, **2017**. https://doi.org/10.1090/mbk/107.

34. Darroch, J. N.; Seneta, E. **On Quasi-Stationary Distributions in Absorbing Discrete-Time Finite Markov Chains.** *J. Appl. Probab.* **1965**, *2*, 88–100. https://doi.org/10.2307/3211876.

35. Collet, P.; Martínez, S.; San Martín, J. **Quasi-Stationary Distributions: Markov Chains, Diffusions and Dynamical Systems**; Springer: Berlin, Heidelberg, **2013**. https://doi.org/10.1007/978-3-642-33131-2.

36. Altmann, E. G.; Portela, J. S. E.; Tél, T. **Leaking Chaotic Systems.** *Rev. Mod. Phys.* **2013**, *85*, 869–918. https://doi.org/10.1103/RevModPhys.85.869.

37. Ryu, J.-W.; Lee, S.-Y.; Kim, C.-M.; Park, Y.-J. **Survival Probability Time Distribution in Dielectric Cavities.** *Phys. Rev. E* **2006**, *73*, 036207. https://doi.org/10.1103/PhysRevE.73.036207.

38. Altmann, E. G.; Portela, J. S. E.; Tél, T. **Chaotic Systems with Absorption.** *Phys. Rev. Lett.* **2013**, *111*, 144101. https://doi.org/10.1103/PhysRevLett.111.144101.

39. Clauß, K.; Altmann, E. G.; Bäcker, A.; Ketzmerick, R. **Structure of Resonance Eigenfunctions for Chaotic Systems with Partial Escape.** *Phys. Rev. E* **2019**, *100*, 052205. https://doi.org/10.1103/PhysRevE.100.052205.

40. Veach, E. **Robust Monte Carlo Methods for Light Transport Simulation.** Ph.D. Thesis, Stanford University, Stanford, CA, **1997**.

41. Soler, C.; Subr, K. **Spectral Theory of Light Transport Operators.** *ACM Trans. Graph.* **2026**, *45*, 1–18. https://doi.org/10.1145/3774756.

42. Ashikhmin, M.; Shirley, P. **An Anisotropic Phong BRDF Model.** *J. Graph. Tools* **2000**, *5*, 25–32. https://doi.org/10.1080/10867651.2000.10487522.

43. Kelemen, C.; Szirmay-Kalos, L. **A Microfacet Based Coupled Specular-Matte BRDF Model with Importance Sampling.** In *Eurographics 2001 Short Presentations*; Eurographics Association, **2001**.

44. Bortz, J. C.; Shatz, N. E.; Ries, H.; Winston, R. **Consequences of Étendue and Skewness Conservation for Nonimaging Devices with Inhomogeneous Sources and Targets.** *Proc. SPIE* **1997**, *3139*, 59–75. https://doi.org/10.1117/12.290222.

45. Muschaweck, J.; Wiesmann, C. **Extracting Light out of LEDs.** *Adv. Opt. Technol.* **2013**, *2*, 291–300. https://doi.org/10.1515/aot-2013-0027.

46. Lenef, A.; Piquette, A.; Kelso, J. **Thermodynamics of Light Extraction from Luminescent Materials.** *ECS J. Solid State Sci. Technol.* **2018**, *7*, R3211–R3226. https://doi.org/10.1149/2.0271801jss.

47. Clarke, F. J. J.; Parry, D. J. **Helmholtz Reciprocity: Its Validity and Application to Reflectometry.** *Light. Res. Technol.* **1985**, *17*, 1–11. https://doi.org/10.1177/14771535850170010301.

48. Koenderink, J. J.; van Doorn, A. J. **Geometrical modes as a general method to treat diffuse interreflections in radiometry.** *J. Opt. Soc. Am.* **1983**, *73*, 843. https://doi.org/10.1364/josa.73.000843.

49. Neyts, K.; Ullan Nieto, A. **Importance of scattering and absorption for the outcoupling efficiency in organic light-emitting devices.** *J. Opt. Soc. Am. A* **2006**, *23*, 1201. https://doi.org/10.1364/josaa.23.001201.

50. Alves-Pereira, A. R.; Nunes-Pereira, E. J.; Martinho, J. M. G.; Berberan-Santos, M. N. **Radiation trapping in 1D using the Markov chain formalism: a computational physics project.** *Eur. J. Phys.* **2007**, *28*, 1105–1124. https://doi.org/10.1088/0143-0807/28/6/008.

51. Xu, H. **Photon escape from slab Thomson media: a scattering-order-resolved recursive formalism for comptonization applications.** *Mon. Not. R. Astron. Soc.* **2026**, *550*, stag1392. https://doi.org/10.1093/mnras/stag1392.

52. Langer, M. S. **When Shadows Become Interreflections.** *Int. J. Comput. Vis.* **1999**, *34*, 193–204. https://doi.org/10.1023/A:1008131719047.

53. Chance, R. R.; Prock, A.; Silbey, R. **Molecular Fluorescence and Energy Transfer Near Interfaces.** *Adv. Chem. Phys.* **1978**, *37*, 1–65. https://doi.org/10.1002/9780470142561.ch1.

54. Wasey, J. A. E.; Safonov, A.; Samuel, I. D. W.; Barnes, W. L. **Effects of dipole orientation and birefringence on the optical emission from thin films.** *Opt. Commun.* **2000**, *183*, 109–121. https://doi.org/10.1016/s0030-4018(00)00806-3.

55. Penninck, L.; De Visschere, P.; Beeckman, J.; Neyts, K. **Dipole radiation within one-dimensional anisotropic microcavities: a simulation method.** *Opt. Express* **2011**, *19*, 18558. https://doi.org/10.1364/oe.19.018558.

56. Moon, C.-K.; Kim, S.-Y.; Lee, J.-H.; Kim, J.-J. **Luminescence from oriented emitting dipoles in a birefringent medium.** *Opt. Express* **2015**, *23*, A279. https://doi.org/10.1364/oe.23.00a279.

57. Mladenovski, S.; Hofmann, S.; Reineke, S.; Penninck, L.; Verschueren, T.; Neyts, K. **Integrated optical model for organic light-emitting devices.** *J. Appl. Phys.* **2011**, *109*, 083114. https://doi.org/10.1063/1.3576114.

58. Henyey, L. G.; Greenstein, J. L. **Diffuse Radiation in the Galaxy.** *Astrophys. J.* **1941**, *93*, 70–83. https://doi.org/10.1086/144246.

59. van de Hulst, H. C. **Multiple Light Scattering: Tables, Formulas, and Applications**, Vols. 1–2; Academic Press: New York, **1980**.

60. Prahl, S. A. **The Adding-Doubling Method.** In *Optical-Thermal Response of Laser-Irradiated Tissue*; Plenum: New York, **1995**; pp 101–129. https://doi.org/10.1007/978-1-4757-6092-7_5.

61. McIntosh, K. R.; Baker-Finch, S. C. **OPAL 2: Rapid optical simulation of silicon solar cells.** In *Proc. 2012 38th IEEE Photovoltaic Specialists Conference*; **2012**; pp 265–271. https://doi.org/10.1109/pvsc.2012.6317616.

62. Koo, W. H.; Jeong, S. M.; Araoka, F.; Ishikawa, K.; Nishimura, S.; Toyooka, T.; Takezoe, H. **Light Extraction from Organic Light-Emitting Diodes Enhanced by Spontaneously Formed Buckles.** *Nat. Photonics* **2010**, *4*, 222–226. https://doi.org/10.1038/nphoton.2010.7.

63. Ogilvy, J. A. **Wave scattering from rough surfaces.** *Rep. Prog. Phys.* **1987**, *50*, 1553–1608. https://doi.org/10.1088/0034-4885/50/12/001.

64. Thorsos, E. I. **The validity of the Kirchhoff approximation for rough surface scattering using a Gaussian roughness spectrum.** *J. Acoust. Soc. Am.* **1988**, *83*, 78–92. https://doi.org/10.1121/1.396188.

65. Gutmann, H.-M. **A Radial Basis Function Method for Global Optimization.** *J. Glob. Optim.* **2001**, *19*, 201–227. https://doi.org/10.1023/A:1011255519438.

66. Regis, R. G.; Shoemaker, C. A. **A Stochastic Radial Basis Function Method for the Global Optimization of Expensive Functions.** *INFORMS J. Comput.* **2007**, *19*, 497–509. https://doi.org/10.1287/ijoc.1060.0182.

67. Torczon, V. **On the Convergence of Pattern Search Algorithms.** *SIAM J. Optim.* **1997**, *7*, 1–25. https://doi.org/10.1137/S1052623493250780.

68. Audet, C.; Dennis, J. E. **Analysis of Generalized Pattern Searches.** *SIAM J. Optim.* **2002**, *13*, 889–903. https://doi.org/10.1137/S1052623400378742.

69. Holm, S. **A Simple Sequentially Rejective Multiple Test Procedure.** *Scand. J. Stat.* **1979**, *6*, 65–70.

70. Ries, H.; Shatz, N.; Bortz, J.; Spirkl, W. **Performance Limitations of Rotationally Symmetric Nonimaging Devices.** *J. Opt. Soc. Am. A* **1997**, *14*, 2855–2862. https://doi.org/10.1364/JOSAA.14.002855.

71. Lim, J.; Oh, S. S.; Kim, D. Y.; Cho, S. H.; Kim, I. T.; Han, S. H.; Takezoe, H.; Choi, E. H.; et al. **Enhanced out-coupling factor of microcavity organic light-emitting devices with irregular microlens array.** *Opt. Express* **2006**, *14*, 6564. https://doi.org/10.1364/oe.14.006564.
