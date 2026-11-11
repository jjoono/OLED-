# Reciprocal Angular Transport Explains the Limited Benefit of Mixing in OLED Light Extraction

**Authors** — [TBD]

---

## Abstract

How much can angular mixing improve light extraction from organic light-emitting diodes? We describe repeated encounters between an OLED and its extraction structure as angle-resolved transport constrained by optical reciprocity. In a fully mixing, reciprocal return model, recirculating light escapes with the étendue-weighted mean escape probability reduced by a term proportional to its angular variance. For 37 microlens, scattering, randomized-lens, and textured structures, complete mixing would add no more than about two percentage points to the spectrally weighted extraction, more only where the cavity sends light into poorly transmitted angles. A reduced model with one spectral mixing parameter reproduces full transport calculations with a mean absolute error of 0.30 points over 342 cases; ray tracing and out-of-sample tests support it except for cone arrays.

---

## 1. Introduction

A planar OLED loses much of its generated light to waveguided modes, electrode absorption, and total internal reflection at the substrate/air interface [1,2]. External structures extract the substrate modes without touching the electrical stack: microlens arrays (MLAs) [3,4], scattering layers [5,6], and textured films [7]. Their extraction has been mapped against lens profile, aspect ratio, and fill factor [3,8], and their shapes are increasingly optimized by inverse design and machine learning [9,10]. With such structures the highest reported external quantum efficiencies of single-unit OLEDs lie near 55–62% [5,11–13], a ceiling recently traced to the loss that returned light suffers on each round trip through the device [14].

What is missing is a compact way to compare the structures themselves. Statistical ray optics assumes that light in the substrate is fully randomized on every pass [15–17], so that $\eta = p/[p+(1-p)A']$ depends on one escape probability $p$ and one round-trip loss $A'$ [14]. Angle-resolved models propagate the power distribution from pass to pass with redistribution matrices [2,14,18–22]. They are accurate but numerical and do not separate escape from redistribution; where randomization has been quantified, it has been done empirically [23–25]. The distinction matters: a hemispherical MLA returns unextracted light largely at its original angle [14], so later passes see mostly light at angles it transmits poorly, and the memoryless formula overestimates extraction by a structure-dependent amount. This leaves open a basic question: how much can better angular mixing gain at all?

Here we answer it with the spectrum of the return matrix. Reciprocity makes this matrix self-adjoint in the étendue inner product, and within a rank-one return model it uniquely determines the angular redistribution of returned light: into each angle in proportion to how much that angle itself returns. Completely mixed light then escapes with $p - \sigma^2/(1-p)$ rather than $p$, so that what mixing can gain is tied to the angular variance of a structure's escape probability. Reciprocity alone does not bound the efficiency of every return operator, but it imposes nontrivial constraints on angular redistribution that explain why faster mixing gives only limited gains in common passive structures. The second eigenvalue gives an angular mixing number, the per-pass analogue of the transport factor $1/(1-g)$ of radiative transfer [26] and the relaxation time of a reversible chain [27]; the leading eigenvalue gives the escape probability of the quasi-stationary distribution [28], the escape rate of open dynamical systems [29–31]. These objects, and the reciprocity that makes them real and bounded, are standard in Markov-chain, open-system, and light-transport theory [32–34], and reciprocal coupling of a diffuse lobe to angle-dependent reflection is standard in reflectance modelling [35,36]. What is new is their quantitative consequence for passive light extraction: a closed-form complete-mixing reference for what mixing can gain, a reciprocal reduced model tested against direct ray tracing and against structures not used in its development, including where it fails, and a map of OLED outcoupling structures, which we connect to the radiance-ordered extraction bound [37–39] and test on a controlled freeform-versus-hemisphere benchmark [11,12].

---

## 2. Theory

### 2.1 Multi-pass extraction as a Markov chain on angle

The substrate angle is divided into 1° polar bins $k$; $v_k$ is the power travelling towards the structure in bin $k$, integrated over azimuth, and $w_k = \sin^2\theta_{k+1}-\sin^2\theta_k$ is the bin's étendue, so uniform radiance is $v_k \propto w_k$. On each encounter a fraction $p_k$ escapes and $B_R(i\to r)$ returns from bin $i$ to bin $r$, so that

$$
B_R = M\,\mathrm{diag}(1-\mathbf{p}),
$$

where the column-stochastic $M$ is the distribution of the returned angle conditioned on return, the transition matrix of a Markov chain on angle. The escape probability for uniform radiance is $p = \sum_k w_k p_k/\sum_k w_k$. The device reflects the returned light with its angle-resolved reflectance $R(\theta)$, and the substrate-to-air extraction efficiency is

$$
\eta_{\mathrm{sta}} = \mathbf{p}^{\mathsf T}\left(I - \mathrm{diag}\,\mathbf{R}\; B_R\right)^{-1}\mathbf{v}^{(0)} ,
$$

with $\mathbf{v}^{(0)}$ the power launched into the substrate [14]. All losses enter through $R(\theta)$; the scalar $A' = 1-\sum_k w_k R_k/\sum_k w_k$ is used only in the memoryless comparison. The model is scalar and averaged over position and azimuth (Note S1): it is an azimuth- and position-marginalized chain, accurate when correlations of azimuth and impact position between successive encounters are short-lived, which direct $(\theta,\phi)$ tracing confirms to within 0.6 pp (Sec. 4.2). Note S10 gives the estimators of $p_k$, $B_R$, and $M$ as pseudocode, with the bookkeeping of escape, return, and lost rays. Resolving s and p polarization without conversion between them, the opposite limit, raises $\eta_{\mathrm{sta}}$ by 0.6–1.5 pp for a hemispherical MLA and a wavy surface but changes the complete-mixing gain by at most 0.4 pp and the reduced-model error by at most 0.1 pp (Table S18).

### 2.2 Reciprocity and the angular mixing number

For a passive, reciprocal structure, $w_i B_R(i\to j) = w_j B_R(j\to i)$ [40], as noted for the redistribution matrices of textured sheets [19]; it is the discrete form of the self-adjointness of reciprocal scattering operators in the étendue inner product [32]. $M$ is therefore similar to a symmetric matrix, its spectrum is real [33], and on the full angular domain its stationary distribution is $\propto w(1-p)$; for the truncated chain used below it is $\propto w s$, with $s_i$ the fraction returned into the retained bins (Note S1). The eigenvalue of second-largest modulus, $\lambda_2$, sets how fast the angle is forgotten, and we define

$$
N_{\mathrm{mix}} = \frac{1}{1-|\lambda_2|},
$$

a spectral measure of angular relaxation in units of encounters ($\infty$ for a planar interface, 1 for a Lambertian scatterer); it is a relaxation time, not an exact count of the passes a photon needs. For all structures here $\lambda_2$ is real and positive. Near-grazing states carry almost no étendue but can retain angle, so the spectral metrics are evaluated on $\theta < 85°$, which excludes 0.76% of the étendue, while $\eta_{\mathrm{sta}}$ is always computed on 0–90°. The cutoff controls the influence of these physically minor states rather than tuning accuracy: the reduced-model error is 0.32, 0.30, 0.40, and 0.55 pp on average for cutoffs of 80°, 85°, 88°, and 90°. It changes the absolute values of the metrics (for the 40–60° freeform MLA, $N_{\mathrm{mix}}$ is 7.5 at 85° and 17 at 90°) but largely preserves the relative ranking of the structures investigated (rank correlation 0.90–0.99; Note S1, Table S1).

### 2.3 Reciprocal reduced model

The reduced model keeps $p(\theta)$ and lets a fraction $1-m$ of the returned light keep its angle while a fraction $m$ is redistributed with a fixed distribution $q$:

$$
B_R(\theta_i\to\theta_r) = \left[1-p(\theta_i)\right]\left[(1-m)\,\delta_{ir} + m\,q(\theta_r)\right].
$$

Detailed balance holds only if $q \propto w(1-p)$: light is redistributed into each angle in proportion to how much that angle itself returns, as in reciprocal reflectance models that couple a diffuse lobe to angle-dependent specular reflection [35,36]. To our knowledge this has not been imposed in multi-pass extraction models, which redistribute isotropically [15–17,41]. Isotropic redistribution is reciprocal inside a separate scattering layer, but applied to the combined return of a structure whose escape depends on angle it sends returned light into an escape cone that returns almost nothing. The normalized kernel has eigenvalues 1 and $1-m$, so $m = 1-\lambda_2$.

### 2.4 Complete reciprocal mixing as a reference

For $m = 1$ the reciprocal kernel has rank one, $B_R(i\to r) = (1-p_i)\,q_r$ with $q \propto w(1-p)$. Mixed light then escapes with probability

$$
p_{\mathrm{mix}} = \sum_r q_r p_r = p - \frac{\sigma_p^2}{1-p}, \qquad \sigma_p^2 = \frac{\sum_k w_k (p_k - p)^2}{\sum_k w_k},
$$

which is also its late-pass escape probability, so that even a complete mixer has $\xi = 1 - \sigma_p^2/[p(1-p)] < 1$ unless $p(\theta)$ is constant. The series then sums in closed form,

$$
\eta_{\mathrm{mix}} = \mathbf{p}\cdot\mathbf{v}^{(0)} + \frac{(\mathbf{p}\cdot\mathbf{u})\,[(\mathbf{1}-\mathbf{p})\cdot\mathbf{v}^{(0)}]}{1-(\mathbf{1}-\mathbf{p})\cdot\mathbf{u}}, \qquad u_k = R_k\,q_k ,
$$

and for an angle-independent reflectance $1-A'$ it becomes $\eta_{\mathrm{mix}} = p_0 + (1-p_0)\,p_{\mathrm{mix}}(1-A')/[p_{\mathrm{mix}} + (1-p_{\mathrm{mix}})A']$, with $p_0$ the first-pass escape probability. This reduces to the memoryless formula only when $p(\theta)$ is constant. A refractive surface transmits inside the critical cone and returns everything beyond it, so $\sigma_p^2$ is large: a planar interface has $p$ = 0.40 but $p_{\mathrm{mix}}$ = 0.06. By letting mixed light escape with $p$, the memoryless formula overstates what mixing can achieve by an amount set by the angular variance of escape. $\eta_{\mathrm{mix}}$ is a physically constrained reference for angular redistribution, not a universal upper bound: other reciprocal return operators with the same $p(\theta)$ can extract more (Note S1 gives a two-state example in which a non-rank-one reciprocal kernel exceeds $\eta_{\mathrm{mix}}$ by 5 pp).

### 2.5 Late-pass escape ratio

$\lambda_2$ tells how fast the angle is forgotten, not where the light ends up. Light that circulates many times approaches the leading eigenvector of $B_R$, the quasi-stationary distribution [28], and escapes with probability $p_\infty = 1-\rho(B_R)$ per encounter, the escape rate of an open system [29–31]. Self-adjointness and the Rayleigh quotient with the uniform-radiance vector, together with the column sums [42], give

$$
\min_\theta p(\theta) \le p_\infty \le \tilde p ,
$$

where $\tilde p$ is the uniform-radiance probability of leaving the truncated chain (escape, or return beyond the cutoff; $\tilde p - p < 0.008$). We call $\xi = p_\infty/\tilde p$ the late-pass escape ratio: 1 when long-circulating light escapes as easily as Lambertian light, 0 when some angles never escape and never mix. All structures satisfy the bounds numerically (Note S1). $p_\infty$ is defined in the lossless device-return limit and describes the structure alone. In an OLED with angle-dependent reflectance the relevant survival operator is $\mathrm{diag}\,\mathbf{R}\,B_R$, whose leading eigenvalue accounts for both extraction and parasitic loss; it ranks the structures like $p_\infty$ (Spearman 0.96, Al, 550 nm).

### 2.6 Band-resolved radiance-ordered bound

The power that a passive structure can send into an output band in one encounter is bounded by filling the band's étendue with the brightest internal radiance first [37–39]; for uniform radiance this is the étendue ratio, $1/n^2$ for the hemisphere [15]. We evaluate it per polar band with the actual first-pass radiance of the microcavity: substrate bins are sorted by radiance per étendue element, $v^{(0)}_k/w_k$, and each air band's étendue (polar bins, full azimuth, scalar) is filled with the brightest bins first (Note S1). Bands are always defined on the air-side far field; for real structures the power reaching a band is accumulated through the traced transmission matrix from substrate bins into air bins.

---

## 3. Methods

**Devices.** Dipole emission and stack reflectance are computed with the Chance–Prock–Silbey model [43] extended to birefringent layers [44,45] and transfer matrices, with measured optical constants and an incoherent substrate [46]. The Al platform of the benchmark is Al (100 nm) / B3PyMPM (61.3 nm) / TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd (25 nm) / TCTA (10 nm) / TAPC (97.5 nm) / ITO (150 nm) / glass ($n = 1.51$), with a horizontal dipole ratio of 0.865; a second stack replaces Al by Ag, and a low-loss stack follows ref. [14]. At 550 nm, $A'$ = 0.167, 0.044, and 0.026 (Table S2, Fig. S1). Device optics are computed per wavelength with dispersive constants; the structure matrices are wavelength independent (ray optics, non-dispersive index 1.51). Spectrally weighted quantities use 400–700 nm in 10 nm steps, weighted with the emitter spectrum times $\eta_{\mathrm{sub}}(\lambda)$ (Note S10).

**Structures.** All structures are index-matched to the substrate. The library comprises a planar interface; six non-absorbing Henyey–Greenstein [47] scattering slabs ($g$ = 0–0.9, $\tau^*$ = 0.25–4) combined with the exit face by the adding method [48], which agrees with a direct Monte Carlo slab trace to 0.0002 in $p$ and 0.05 pp in $\eta_{\mathrm{sta}}$ (Note S10); hemispherical and freeform hexagonal MLAs, the optima of the benchmark; 17 random MLAs on a 6 × 6 supercell with independent lenslet profiles (spread $s$), heights, or positions and radii; and five Gaussian-correlated wavy surfaces of rms slope 0.2–1.5 (correlation length $\ell$, rms height 0.13–0.96 $\ell$), which resemble the buckles used on OLEDs [49]. Ray optics is scale invariant, so these results hold for any $\ell$ well above the wavelength, e.g. $\ell \gtrsim$ 10 µm; only the shape of the height spectrum matters, which the ring-spectrum surfaces of the blind test vary. MLAs and textures are traced with a height-field Monte-Carlo tracer that reproduces the LightTools return matrix of the hemispherical MLA (escape probability 0.4275 against 0.4274). Ray optics requires features well above the wavelength [50]; full-wave simulations of microlenses converge to ray tracing above about 5 µm [8]. Ray noise, realization, ray count, bin width, domain size, and the residual reciprocity asymmetry change $N_{\mathrm{mix}}$ by at most 0.24, $\xi$ by 0.025, and $\eta_{\mathrm{sta}}$ by 0.34 pp (Tables S5, S6). The traced matrices obey detailed balance within Monte Carlo noise (relative residual $\lVert F-F^{\mathsf T}\rVert/\lVert F\rVert$ = 0.013–0.15, median 0.047, with $F_{ji} = w_i B_R(i\to j)$), and projecting them onto exactly reciprocal matrices changes $N_{\mathrm{mix}}$ by at most 2.7% and $\xi$ by 0.005 (Note S10). Two independent checks are made. Direct Monte Carlo tracing of the whole multi-pass process in $(\theta, \phi)$, with specular device reflection $R(\theta)$ and each ray's azimuth kept between passes, tests the reduction to a chain on polar angle; and nine structures not used in developing the model (spherical caps, paraboloids, cones, Gaussian bumps, a random mix of cones and hemispheres, and wavy surfaces with a ring-shaped spectrum) form a blind test set.

**Freeform benchmark.** An axisymmetric freeform MLA (13 variables: spline profile, height, ETL and HTL thicknesses) is compared with an equally optimized hemispherical MLA (3 variables) on the Al platform with LightTools ray tracing validated against fabricated devices [11,12]; the hemisphere is a point of the freeform design space (Table S10). Total EQE and the EQE in four polar bands (0–20°, 20–40°, 40–60°, 60–80°) are optimized separately by surrogate and pattern search (MATLAB surrogateopt, patternsearch). Each optimum is re-evaluated three times at full fidelity, and a restart control repeats every search from the hemispherical optimum; gains are tested with a pooled $t$ statistic and Holm correction (Note S6).

**Fig. 1 | Multi-pass extraction as a Markov chain on angle.** (a) The passive structure transmits $p(\theta)$ and returns the rest with $B_R(\theta_i\to\theta_r)$; the device reflects it with $R(\theta)$. (b–d) Return matrices $M$ ($\log_{10}$, floor $10^{-3}$) of a planar interface, the hemispherical MLA, and a scattering layer ($g$ = 0.8, $\tau^*$ = 1). The off-diagonal bands of (c) are deterministic pathways of rays totally reflected at the lens flank (e.g. 70° → 34°); they are unchanged in a 6 × 6 supercell and with twice the rays (Note S10). *(`fig1_mixing_concept.png`)*

---

## 4. Results

### 4.1 What complete mixing can gain

Figure 2 asks the central question directly: how much would complete reciprocal mixing, at each structure's own $p(\theta)$, add to its efficiency? Completely mixed light escapes with $p_{\mathrm{mix}}$ well below $p$ whenever escape depends strongly on angle (Fig. 2a): $p_{\mathrm{mix}}/p$ is 0.87–0.91 for hemispherical MLAs, 0.81–0.93 for freeform and random MLAs, 0.37–0.90 for wavy surfaces, and 0.14 for the planar interface. At 550 nm the resulting gain for all 37 non-planar structures is −0.4 to +1.8 pp on Al, at most 2.0 pp on Ag, and at most 0.4 pp in the low-loss device (Fig. 2b); hemispherical MLAs, the slowest mixers with $N_{\mathrm{mix}}$ = 6–8, would gain only 1.3–1.7 pp on Al. Weighted with the emission spectrum (peak 572 nm), the gain is at most 1.2, 1.5, and 2.2 pp for the three devices. It is larger where the cavity launches light into angles the structure transmits poorly. At 500 nm the first-pass escape probability of the hemispherical MLA drops to 0.36 (Al) and 0.33 (Ag), against $p$ = 0.42, and complete mixing would add 5.7–6.4 and 8.4–9.4 pp to hemispherical MLAs; in the low-loss device at 600 nm it would add 3.2–4.0 pp. There, random lenses and wavy surfaces realize most of it: the wavy surface of slope 1.5 gains 5.8, 8.4, and 4.0 pp over the hemispherical MLA (Table S19). The reference thus also tells when mixing is worth pursuing: when the first pass is mismatched to $p(\theta)$. The memoryless formula predicts gains larger by 0.6–2.9 pp for lens arrays, by 17 pp for the shallowest wavy surface, and by 31 pp for the planar interface, whose complete mixing raises $\eta_{\mathrm{sta}}$ only from 0.38 to 0.49 ($p_{\mathrm{mix}}$ = 0.06). Mixing cannot open angles at which a structure does not transmit; only a change of $p(\theta)$ can. Two wavy surfaces, of rms slope 0.4 and 0.7, exceed the reference by 0.42 and 0.01 pp. They are the only structures whose escape probability is anticorrelated between successive encounters (one-step autocorrelation of $p$ under the stationary chain, $C$ = −0.09 and −0.08): light that has just failed to escape is sent preferentially to angles of high escape, which complete mixing ($C$ = 0) does not do. Across the library $C$ and $\eta_{\mathrm{sta}}-\eta_{\mathrm{mix}}$ are correlated with coefficient −0.76 (Note S10), so exceeding the reference is a predictable property of the operator, consistent with $\eta_{\mathrm{mix}}$ being a reference rather than a bound. Here $p$, $\sigma_p^2$, and $p_{\mathrm{mix}}$ are computed on 0–90°; on 0–80° or 0–85°, $p_{\mathrm{mix}}/p$ changes by at most 0.02.

**Fig. 2 | What complete reciprocal mixing can gain.** (a) Escape probability of completely mixed light, $p_{\mathrm{mix}} = p - \sigma_p^2/(1-p)$, against the uniform-radiance escape probability $p$ for all structures; dotted line, $p_{\mathrm{mix}} = p$ as assumed by the memoryless formula. (b) Gain over the actual $\eta_{\mathrm{sta}}$ predicted for complete reciprocal mixing at each structure's own $p(\theta)$ (filled) and by the memoryless formula (open), against $N_{\mathrm{mix}}$ (Al, 550 nm). *(`fig2_complete_mixing.png`)*

### 4.2 Two parameters predict multi-pass extraction

Figure 3 compares four models with the full series for 342 cases: three devices, three wavelengths (500, 550, 600 nm), and 38 structures. The memoryless formula overestimates the planar interface by up to 60 pp, because a flat face never lets returned light into the escape cone, and the hemispherical MLA by 2–3 pp (Fig. 3a). Assuming that every ray keeps its angle fails in the opposite direction for scattering layers (Fig. 3b). With $m = 1-\lambda_2$ and isotropic redistribution the reduced model fits the lens arrays but overestimates shallow wavy surfaces by up to 10.5 pp (Fig. 3c): their $\lambda_2$ is low, yet trapped light is returned mostly to other trapped angles (from 70–80°, only 3–8% reaches angles with $p > 0.3$, against 50% for isotropic redistribution). The reciprocal redistribution builds this in and reduces the mean absolute error to 0.30 pp and the maximum to 3.0 pp, without an additional parameter (Fig. 3d); the error is not concentrated in any device or wavelength (Fig. 3e, Table S4). The model slightly underestimates on average (mean signed error −0.24 pp); the largest error, −3.0 pp, is the shallowest wavy surface on Ag at 550 nm (Note S10). For the structures in the library, angle-dependent escape and a single spectral mixing parameter thus give an accurate reduced description of multi-pass extraction, provided the redistribution respects reciprocity. Reciprocity alone does not guarantee that a single rate suffices, as the cone arrays below show.

Two independent checks test this hierarchy (Note S9, Table S16). Direct ray tracing in $(\theta, \phi)$ agrees with the full series within 0.6 pp for eight structures; the largest differences (−0.4 pp for the hemispherical MLA, −0.6 pp for 60° cones) vanish when the azimuth is re-randomized at each pass, so they measure the azimuthal memory that a chain on polar angle neglects. On the blind set the reduced model errs by 0.48 pp on average and 3.3 pp at most for spherical caps, paraboloids, Gaussian bumps, a cone–hemisphere mix, and ring-spectrum wavy surfaces (63 cases), as on the library, but fails for the two cone arrays (mean 4–6 pp, up to 20 pp). The modes explain why. For 60° cones the slow modes ($\lambda_2$–$\lambda_4$ = 0.985, 0.91, 0.85) carry only 0.4% of the returned light's deviation from the stationary distribution, so $\lambda_2$ describes states the light hardly occupies; the relaxation rate weighted by the actual returned light, $\lambda_{\mathrm{eff}}$ = 0.43, reduces the error to 0.23 pp (maximum 0.8 pp). The 45° cones retroreflect near-normal light into itself, and this closed set holds 98% of that deviation, a genuinely multi-rate case. $\lambda_{\mathrm{eff}}$ is not better on the library (mean error 0.48 against 0.31 pp), so we keep $\lambda_2$ and use the slow-mode occupancy as a post hoc diagnostic of when one rate suffices (Note S9). The full matrix and the bounds are unaffected. In practice, one rate is accurate when the light returned by the structure has appreciable overlap with the slow eigenmodes, so that $\lambda_2$ describes states the light actually occupies; a large $\lambda_2$ with negligible slow-mode occupancy signals that a single rate will misstate the mixing.

**Fig. 3 | Validation of the reduced model.** $\eta_{\mathrm{sta}}$ from (a) the memoryless formula, (b) the angle-preserving limit ($m = 0$), and the reduced model with $m = 1-\lambda_2$ and (c) isotropic ($q \propto w$) or (d) reciprocal ($q \propto w(1-p)$) redistribution, against the full series; symbols encode the structure class (342 cases: three devices, three wavelengths). (e) Error of (d) by device and wavelength; bars, mean absolute error. *(`fig3_reduced_model.png`)*

### 4.3 The structure map

On the $(p, N_{\mathrm{mix}})$ plane the classes separate (Fig. 4a). The planar interface has $p$ = 0.40 and $N_{\mathrm{mix}} = \infty$ ($\eta_{\mathrm{sta}}$ = 0.38 on Al). Scattering layers mix within 1.1–2.3 passes but escape less, $p$ = 0.19–0.39 ($\eta_{\mathrm{sta}}$ = 0.58–0.75). Hemispherical MLAs have the highest escape probability, $p$ = 0.42, but $N_{\mathrm{mix}}$ = 6–8 ($\eta_{\mathrm{sta}}$ = 0.78–0.79); freeform MLAs share $p$ = 0.42–0.43 with $N_{\mathrm{mix}}$ = 3–7.5.

Randomness lowers $N_{\mathrm{mix}}$ without lowering $p$ (Fig. 4b, Table S3). Each randomness parameter is traced at 3–6 levels with one realization per level; ± is the standard deviation over three geometric realizations at the end point (Table S5). Independent lenslet profiles reduce it from 7.0 to 1.9 ± 0.1, position jitter of 30% of the pitch to 3.5 ± 0.1, and height spread of ±60% only to 4.1 ± 0.2. Wavy surfaces of rms slope 0.4–1.5 reach $N_{\mathrm{mix}}$ = 1.1–1.2, as low as a scattering layer, with $p$ rising from 0.41 to 0.43; at slope 1.5 they give the highest $\eta_{\mathrm{sta}}$ of all structures, 0.801 against 0.792 for the hemispherical MLA.

The late-pass escape ratio adds what $N_{\mathrm{mix}}$ cannot see (Fig. 4c). For lens arrays the two are tied: $\xi$ falls from 0.80 to 0.39 as $N_{\mathrm{mix}}$ rises from 2 to 8, because light that keeps its angle ends up where the lenses transmit worst; the hemispherical MLA has $\xi$ = 0.47. Random lenses reach 0.7–0.8, wavy surfaces of slope ≥ 0.7 0.81–0.88, and thick scattering layers 0.93–0.99. The shallowest wavy surface ($\xi$ = 0.35) and a thin isotropic scatterer (0.54) fall below the lens trend at $N_{\mathrm{mix}} \approx 2$: they mix quickly but locally, recycling trapped light among trapped angles, and they carry the largest residual errors in Fig. 3d. $\xi$ thus flags where a single mixing rate is least reliable. Its absolute value depends on the cutoff (hemispherical MLA 0.64 at 80°, 0.31 at 90°), but these two structures stay below random lenses of similar $N_{\mathrm{mix}}$ (0.31–0.56 against 0.72–0.83) for every cutoff from 80° to 90° (Note S10).

The contours show what this is worth: they run almost vertically, so $\eta_{\mathrm{sta}}$ is set mainly by $p$. For lens, scattering-layer, and wavy templates of $p(\theta)$ alike, the whole range of $N_{\mathrm{mix}}$ is worth a change of $p$ by only 0.02–0.05 (Table S7); the contour field is a device- and template-dependent illustration, whereas the points are structure descriptors. Randomizing the lenslet profiles moves $N_{\mathrm{mix}}$ from 7.0 to 1.9 but changes $\eta_{\mathrm{sta}}$ by less than 1 pp on Al at 550 nm, because the faster mixing (+0.9 pp on Al) is offset by the reshaped $p(\theta)$ (−2.2 pp; Table S8).

**Fig. 4 | The structure map.** (a) Escape probability and angular mixing number (log axis) of all structures, with contours of $\eta_{\mathrm{sta}}$ (reciprocal reduced model, Al, 550 nm, hemispherical $p(\theta)$ shape scaled to each $p$); the hemispherical and freeform MLAs are the benchmark optima (Table 3). (b) $N_{\mathrm{mix}}$ against each randomness parameter normalized by its maximum; dotted lines, hemispherical MLA (blue) and scattering layer (green). (c) $\xi$ against $N_{\mathrm{mix}}$; the planar interface ($N_{\mathrm{mix}} = \infty$, $\xi$ = 0) is placed off-scale. *(`fig4_structure_map.png`)*

### 4.4 Radiance bound and output composition

On the first pass the microcavity launches a strongly anisotropic radiance (Fig. S1). For the Al platform at 550 nm, the radiance-ordered bound allows 0.085, 0.205, 0.230, and 0.156 of the substrate power into the 0–20°, 20–40°, 40–60°, and 60–80° bands; these are separate maxima that draw on the same bright directions, and the bound for the whole hemisphere is 0.540, against $1/n^2$ = 0.439 for uniform radiance (Fig. 5a). Every passive structure stays well below it: MLAs reach 76–80% in total and 55–70% per band, the wavy surface of slope 1.0 77%, the planar interface 65%, and scattering layers 35–68%. The lenses transmit almost exactly what they would for uniform radiance (0.43 against 0.44), so the cavity's anisotropy is not used.

The output composition is nearly independent of the structure in lossy devices (Fig. 5b). The selectivity $S_{40\text{–}60} = P_{\mathrm{air}}(40\text{–}60°)/P_{\mathrm{air}}(0\text{–}90°)$ of the multi-pass output stays within 0.33–0.36 on Al at 550 nm, close to the Lambertian partition of 0.337, and within 0.32–0.36 at 500 and 600 nm and on Ag (Table S9, Fig. S2); its Monte Carlo uncertainty is at most 0.003 over ray seeds and realizations, an order of magnitude below the window width. In the low-loss device it widens to 0.27–0.37, because the long angle retention of hemispherical MLAs lets late passes, shaped by the structure, dominate when little is lost per round trip.

**Fig. 5 | Radiance bound and composition window.** (a) First-pass power into each polar band (Al, 550 nm): uniform-radiance and radiance-ordered bounds, and four structures. (b) $S_{40\text{–}60}$ of the multi-pass output against $\eta_{\mathrm{sta}}$ for all structures; band, 0.33–0.36; dotted line, Lambertian partition. *(`fig5_bound_window.png`)*

### 4.5 Case study: freeform versus an equally optimized hemisphere

The map explains a controlled shape-optimization benchmark. Band-optimized freeform profiles differ visibly from the hemisphere (Fig. 6a), but the gains are small: $G_j$ = 1.003, 1.002, 1.070, and 1.036 for the four bands and 1.013 for total EQE, relative ratios that correspond to +1.3 and +0.71 pp of absolute EQE for 40–60° and total (Fig. 6b, Table 2). Their reproducibility is judged from independent searches, not from repeated evaluation of one optimum. For each band the original search (random seeds plus the best designs of a separate weighted-sum campaign) and a restart from the hemispherical optimum give 1.046–1.070 (40–60°) and 1.020–1.036 (60–80°); three random starts and a restart give 1.005–1.013 for total EQE; and no search gains more than 0.3% at 0–40° (Fig. 6c). The original search even returned 0.946 and 0.966 at 0–40°, a search failure, since the hemisphere lies in the freeform design space. The hemisphere received more evaluations per variable (about 15 against 4–11), which favours it. Re-evaluating each optimum three times (standard deviation about $10^{-4}$) only shows that these differences exceed the Monte Carlo noise of one evaluation; best-so-far curves of all searches show the hemispherical searches converging within about 30 evaluations and the freeform searches still improving slowly at the end of their budgets (Fig. S5, Note S6). Evaluated through the freeform parameterization with quarter-circle control points, the hemisphere reproduces its own campaign within 0.09%.

**Table 2 | Freeform versus hemisphere (Al platform, LightTools).** Last column: freeform/hemisphere ratio over independent searches (number of searches in parentheses; data in Table S11).

| Objective | Hemisphere | Best freeform | $G_j$ | Independent searches |
|---|---|---|---|---|
| 0–20° | 0.06682 | 0.06699 | 1.003 | 0.946–1.003 (2) |
| 20–40° | 0.16591 | 0.16628 | 1.002 | 0.966–1.002 (2) |
| 40–60° | 0.18503 | 0.19802 | 1.070 | 1.046–1.070 (2) |
| 60–80° | 0.13589 | 0.14075 | 1.036 | 1.020–1.036 (2) |
| Total EQE | 0.54679 | 0.5539 | 1.013 | 1.005–1.013 (4) |

Placed on the map (Table 3), the nine optima span only $p$ = 0.417–0.429, worth about 0.7 pp of $\eta_{\mathrm{sta}}$, while $N_{\mathrm{mix}}$ spans 3.0–8.0 and $\xi$ 0.39–0.79. Shape freedom thus moves a design vertically on the map, where the contours are flat, and not to the right, where they are steep. Across 606 sampled designs the 40–60° EQE is proportional to the total EQE ($R^2$ = 0.968) (Note S7). Because the cavity is optimized together with the lens, we separated the two in the full device model, spectrally weighted, with each design's own cavity (Note S8, Table S13). Lens and cavity interact. Spectrally weighted, the freeform 40–60° lens beats the hemispherical one in its band by 0.4% at the hemisphere's band-optimized cavity, 1.8% at the total-EQE cavity, and 4.1% at its own; at 550 nm alone, by 5.6–9.2%. A lens-only optimization at a fixed stack in the same model (Al, 550 nm) gains 9.0% in the 40–60° band and 1.25% in $\eta_{\mathrm{sta}}$ over an equally optimized hemisphere. The second gain is the one the map predicts: at unchanged $p$ (0.421 against 0.422) the optimized lens mixes faster ($N_{\mathrm{mix}}$ 7.1 → 2.1) and has a flatter $p(\theta)$, recovering 1.0 pp of the 1.4 pp that complete mixing would add to the hemispherical MLA (Sec. 4.1). Lens shape thus steers light into a band by several percent at a fixed cavity but raises the total only through mixing. The model reproduces the LightTools EQE of the hemispheres within 4% but rates the freeform optima 1–5% lower relative to them, so the two tools agree that shape gains are a few percent at most but not on their sign below about 5%. The low $\eta_{\mathrm{sta}}$ of the total-EQE freeform optimum in Table 3 (0.771) reflects this: at the hemisphere's cavity its lens extracts less, and its LightTools advantage rests on its own cavity (ETL 67.5, HTL 120.6 nm) and on a profile with an annular trough that our spline reconstruction captures least well. A LightTools control at fixed cavity, which would settle the size of the lens-shape gains in the benchmark itself, is provided as a script (Note S8).

**Table 3 | Benchmark optima on the structure map** (550 nm; $\eta_{\mathrm{sta}}$ for Al in the multi-pass model; the 0–20° freeform optimum has a vertical sidewall and is omitted).

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

**Fig. 6 | Freeform versus an equally optimized hemisphere.** (a) Band-optimized freeform profiles and the hemispherical reference. (b) Best EQE of each class per objective, labelled with $G_j$. (c) Restart control: relative gain over the re-measured hemisphere, labelled with $t$; green, significant relative to Monte Carlo evaluation noise (one-sided 95%), which is not a measure of uncertainty in the global optimum. *(`fig6_hemisphere_benchmark.png`)*

---

## 5. Discussion

**Why lens shape saturates.** A hemispherical array already escapes close to the étendue bound $1/n^2$ for uniform radiance, so profile changes cannot raise $p$ appreciably; they change $N_{\mathrm{mix}}$ and the shape of $p(\theta)$, along which the efficiency contours are flat. The cavity's anisotropy would in principle allow more: a coextensive surface cannot enlarge a band's étendue, but one with binary angular transmission can reach the radiance-ordered bound [39], and skewness conservation does not forbid this [51]. The rotationally symmetric lenses and textures studied here, however, fill each band from a broad range of internal angles and extract close to what uniform radiance would give. Approaching the band-resolved bound thus appears to call for angle-selective transmission rather than further freedom in the lens profile.

**Why angular mixing gives limited gains.** Optical reciprocity couples angular redistribution to escape. Within a reciprocal rank-one mixing model, light is redistributed towards each angle in proportion to the fraction of incident light that the structure returns from that angle, so the escape probability after complete mixing is reduced by the étendue-weighted variance of the angle-dependent escape probability. This explains why faster angular mixing does not necessarily translate into more efficient extraction. Across the structures examined, faster mixing often comes with a reshaped or reduced escape profile, as in randomized lenses (Table S8) and scattering layers, whose extraction peaks at intermediate scatterance [18,21]. For the structures examined, improving the escape profile is generally more effective than increasing the rate of angular mixing alone.

**Design implications.** Three numbers are worth reporting for a new outcoupling structure: $p$, which sets the efficiency; $N_{\mathrm{mix}}$, which sets the accuracy of any memoryless description; and $\xi$, which tells whether long-circulating light is still extractable. Random MLAs and wavy surfaces of rms slope one or more combine a high $p$ with low $N_{\mathrm{mix}}$ and high $\xi$, but their spectrally weighted advantage over periodic MLAs is 0.2–0.5 pp in the Al and Ag devices and 2 pp in the low-loss device, larger only at wavelengths where the first pass is mismatched to $p(\theta)$; irregular microlenses have been used on OLEDs mainly to avoid diffraction and colour shift [52]. Optimizers, including machine-learning ones [9,10], can be benchmarked against the $(p, N_{\mathrm{mix}})$ map and the band-resolved bound.

**Scope.** The analysis covers passive external structures index-matched to the substrate, in the ray-optical regime; deliberately angle-selective filters, internal structures, absorbing scatterers, and wavelength converters are outside it. As in earlier substrate models [2,18], the device is an angle-preserving reflector; diffuse reflection would add a non-diagonal factor whose product with $B_R$ need not keep a real spectrum [33].

---

## 6. Conclusion

We asked how much better angular mixing can raise the outcoupling of a passive structure. Reciprocal transport imposes nontrivial constraints on angular redistribution: in a fully mixing, reciprocal return model, recirculating light escapes with $p - \sigma_p^2/(1-p)$, so for 37 scattering, microlens, random-lens, and wavy structures complete mixing would add no more than about 2 pp to the spectrally weighted extraction (up to 9 pp at single wavelengths where the first pass is mismatched to $p(\theta)$), whereas the memoryless formula predicts up to 19 pp (41 pp for a planar interface, where complete mixing adds 10 pp). This reference is not a universal bound, since other reciprocal operators can extract more, but it explains why faster mixing gives limited gains in common passive structures. The same reciprocity makes the spectrum of the return matrix real, bounds the late-pass escape probability, and fixes the redistribution of a rank-one reduced model, which predicts multi-pass extraction with a mean absolute error of 0.30 pp, agrees with direct ray tracing, and predicts seven of nine out-of-sample structure families; it fails for cone arrays, whose slowest modes are not occupied by the light that escapes. Lens-shape optimization changes mainly $N_{\mathrm{mix}}$ and the angular composition, consistent with its gains of a few percent at most. Going further likely calls for angle-selective rather than shape-optimized structures.

---

## Funding

[TBD]

## Acknowledgments

[TBD]

## Disclosures

The authors declare no conflicts of interest.

## Data availability

Return matrices, model inputs, analysis and plotting scripts, and the data underlying all figures will be deposited in a version-controlled public repository [DOI to be added]; Supplement 1 lists the scripts and their inputs.

## Supplemental document

See Supplement 1 for supporting content.

---

## References

1. Brütting, W.; Frischeisen, J.; Schmidt, T. D.; Scholz, B. J.; Mayr, C. **Device Efficiency of Organic Light-Emitting Diodes: Progress by Improved Light Outcoupling.** *Phys. Status Solidi A* **2013**, *210*, 44–65.

2. Greiner, H. **Light Extraction from Organic Light Emitting Diode Substrates: Simulation and Experiment.** *Jpn. J. Appl. Phys.* **2007**, *46*, 4125–4137. https://doi.org/10.1143/JJAP.46.4125.

3. Möller, S.; Forrest, S. R. **Improved Light Out-Coupling in Organic Light Emitting Diodes Employing Ordered Microlens Arrays.** *J. Appl. Phys.* **2002**, *91*, 3324–3327.

4. Qu, Y.; Kim, J.; Coburn, C.; Forrest, S. R. **Efficient, Nonintrusive Outcoupling in Organic Light Emitting Devices Using Embedded Microlens Arrays.** *ACS Photonics* **2018**, *5*, 2453–2458. https://doi.org/10.1021/acsphotonics.8b00255.

5. Song, J.; et al. **Lensfree OLEDs with over 50% External Quantum Efficiency via External Scattering and Horizontally Oriented Emitters.** *Nat. Commun.* **2018**, *9*, 3207. https://doi.org/10.1038/s41467-018-05671-x.

6. Levell, J. W.; Harkema, S.; Pendyala, R. K.; Rensing, P.; Senes, A.; Bollen, D.; MacKerron, D.; Wilson, J. S. **Device Reflectivity as a Simple Rule for Predicting the Suitability of Scattering Foils for Improved OLED Light Extraction.** *Proc. SPIE* **2013**, *8829*, 88291L. https://doi.org/10.1117/12.2023021.

7. Kovačič, M. **Design of External Microtextures for Efficient Light Outcoupling in OLEDs with Different Preferential Orientation of Emission Dipoles.** *Opto-Electron. Rev.* **2022**, *30*, e141542. https://doi.org/10.24425/opelre.2022.141542.

8. Sun, Y.; Forrest, S. R. **Organic Light Emitting Devices with Enhanced Outcoupling via Microlenses Fabricated by Imprint Lithography.** *J. Appl. Phys.* **2006**, *100*, 073106. https://doi.org/10.1063/1.2356904.

9. Kim, S.; Shin, J. M.; Lee, J.; Park, C.; Lee, S.; Park, J.; Seo, D.; Park, S.; Park, C. Y.; Jang, M. S. **Inverse Design of Organic Light-Emitting Diode Structure Based on Deep Neural Networks.** *Nanophotonics* **2021**, *10*, 4533–4541. https://doi.org/10.1515/nanoph-2021-0434.

10. Zhao, Z.-R.; et al. **Deep Q-Network-Guided TiO2/NOA61 Microstructures for Light Outcoupling and Angular Emission Redistribution in Blue Top-Emitting OLEDs.** *Adv. Funct. Mater.* **2026**, e78939. https://doi.org/10.1002/adfm.78939.

11. Kim, M.; Kim, J.; Yoo, S. **Near-Planar Light Outcoupling Structures with Finite Lateral Dimensions for Ultra-Efficient and Optical Crosstalk-Free OLED Displays.** *Nat. Commun.* **2025**, *16*, 11606. https://doi.org/10.1038/s41467-025-66538-6.

12. Kim, J.; Kim, E.; Park, J.; Song, J.; Kim, S.; Moon, H.; Yoo, S. **Toward Near-Foldable Surface Light Sources with Ultimate Efficiency: Ultrathin Substrates Embedded with Micron-Scale Inverted Lens Arrays.** *ACS Photonics* **2023**, *10*, 1775–1782. https://doi.org/10.1021/acsphotonics.3c00017.

13. Li, W.; Li, Y.; Shen, Y.; Zhang, Y.; Jin, T.; Chen, J.; Zhang, X.; Tang, J. **Releasing the Trapped Light for Efficient Silver Nanowires‐Based White Flexible Organic Light‐Emitting Diodes.** *Adv. Opt. Mater.* **2019**, *7*, 1900985. https://doi.org/10.1002/adom.201900985.

14. Kim, J.; Lee, S.; Yoo, S.; Gather, M. C. **Minimizing Round-Trip Loss for Near-Unity Quantum Efficiency in Organic Light-Emitting Diodes.** *arXiv* **2026**, 2610.12065.

15. Yablonovitch, E. **Statistical Ray Optics.** *J. Opt. Soc. Am.* **1982**, *72*, 899–907.

16. Schnitzer, I.; Yablonovitch, E.; Caneau, C.; Gmitter, T. J.; Scherer, A. **30% External Quantum Efficiency from Surface Textured, Thin-Film Light-Emitting Diodes.** *Appl. Phys. Lett.* **1993**, *63*, 2174–2176. https://doi.org/10.1063/1.110575.

17. Green, M. A. **Lambertian light trapping in textured solar cells and light‐emitting diodes: analytical solutions.** *Prog. Photovolt.* **2002**, *10*, 235–241. https://doi.org/10.1002/pip.404.

18. Shiang, J. J.; Duggal, A. R. **Application of radiative transport theory to light extraction from organic light emitting diodes.** *J. Appl. Phys.* **2004**, *95*, 2880–2888. https://doi.org/10.1063/1.1644037.

19. Eisenlohr, J.; Tucher, N.; Höhn, O.; Hauser, H.; Peters, M.; Kiefel, P.; Goldschmidt, J. C.; Bläsi, B. **Matrix Formalism for Light Propagation and Absorption in Thick Textured Optical Sheets.** *Opt. Express* **2015**, *23*, A502–A518. https://doi.org/10.1364/OE.23.00A502.

20. David, A. **Surface-Roughened Light-Emitting Diodes: An Accurate Model.** *J. Display Technol.* **2013**, *9*, 301–316. https://doi.org/10.1109/JDT.2013.2240373.

21. Altazin, S.; Reynaud, C.; Mayer, U. M.; Lanz, T.; Lapagna, K.; Knaack, R.; Penninck, L.; Kirsch, C.; et al. **Simulations, Measurements, and Optimization of OLEDs with Scattering Layer.** *SID Symp. Dig. Tech. Pap.* **2015**, *46*, 564–567. https://doi.org/10.1002/sdtp.10478.

22. Kovačič, M.; Samigullina, D.; Bouchard, F.; Krč, J.; Lipovšek, B.; Soldera, M.; Lasagni, A. F.; Reineke, S.; et al. **Analysis and optimization of light outcoupling in OLEDs with external hierarchical textures.** *Opt. Express* **2021**, *29*, 23701. https://doi.org/10.1364/oe.428021.

23. Deckman, H. W.; Roxlo, C. B.; Yablonovitch, E. **Maximum Statistical Increase of Optical Absorption in Textured Semiconductor Films.** *Opt. Lett.* **1983**, *8*, 491–493. https://doi.org/10.1364/OL.8.000491.

24. Manzoor, S.; Filipič, M.; Onno, A.; Topič, M.; Holman, Z. C. **Visualizing Light Trapping within Textured Silicon Solar Cells.** *J. Appl. Phys.* **2020**, *127*, 063104. https://doi.org/10.1063/1.5131173.

25. Battaglia, C.; Boccard, M.; Haug, F.-J.; Ballif, C. **Light Trapping in Solar Cells: When Does a Lambertian Scatterer Scatter Lambertianly?** *J. Appl. Phys.* **2012**, *112*, 094504. https://doi.org/10.1063/1.4761988.

26. van Rossum, M. C. W.; Nieuwenhuizen, T. M. **Multiple Scattering of Classical Waves: Microscopy, Mesoscopy, and Diffusion.** *Rev. Mod. Phys.* **1999**, *71*, 313–371. https://doi.org/10.1103/RevModPhys.71.313.

27. Levin, D. A.; Peres, Y. **Markov Chains and Mixing Times**, 2nd ed.; American Mathematical Society: Providence, RI, **2017**. https://doi.org/10.1090/mbk/107.

28. Darroch, J. N.; Seneta, E. **On Quasi-Stationary Distributions in Absorbing Discrete-Time Finite Markov Chains.** *J. Appl. Probab.* **1965**, *2*, 88–100. https://doi.org/10.2307/3211876.

29. Altmann, E. G.; Portela, J. S. E.; Tél, T. **Leaking Chaotic Systems.** *Rev. Mod. Phys.* **2013**, *85*, 869–918. https://doi.org/10.1103/RevModPhys.85.869.

30. Ryu, J.-W.; Lee, S.-Y.; Kim, C.-M.; Park, Y.-J. **Survival Probability Time Distribution in Dielectric Cavities.** *Phys. Rev. E* **2006**, *73*, 036207. https://doi.org/10.1103/PhysRevE.73.036207.

31. Clauß, K.; Altmann, E. G.; Bäcker, A.; Ketzmerick, R. **Structure of Resonance Eigenfunctions for Chaotic Systems with Partial Escape.** *Phys. Rev. E* **2019**, *100*, 052205. https://doi.org/10.1103/PhysRevE.100.052205.

32. Veach, E. **Robust Monte Carlo Methods for Light Transport Simulation.** Ph.D. Thesis, Stanford University, Stanford, CA, **1997**.

33. Soler, C.; Subr, K. **Spectral Theory of Light Transport Operators.** *ACM Trans. Graph.* **2026**, *45*, 1–18. https://doi.org/10.1145/3774756.

34. Koenderink, J. J.; van Doorn, A. J. **Geometrical modes as a general method to treat diffuse interreflections in radiometry.** *J. Opt. Soc. Am.* **1983**, *73*, 843. https://doi.org/10.1364/josa.73.000843.

35. Ashikhmin, M.; Shirley, P. **An Anisotropic Phong BRDF Model.** *J. Graph. Tools* **2000**, *5*, 25–32. https://doi.org/10.1080/10867651.2000.10487522.

36. Kelemen, C.; Szirmay-Kalos, L. **A Microfacet Based Coupled Specular-Matte BRDF Model with Importance Sampling.** In *Eurographics 2001 Short Presentations*; Eurographics Association, **2001**.

37. Bortz, J. C.; Shatz, N. E.; Ries, H.; Winston, R. **Consequences of Étendue and Skewness Conservation for Nonimaging Devices with Inhomogeneous Sources and Targets.** *Proc. SPIE* **1997**, *3139*, 59–75. https://doi.org/10.1117/12.290222.

38. Muschaweck, J.; Wiesmann, C. **Extracting Light out of LEDs.** *Adv. Opt. Technol.* **2013**, *2*, 291–300. https://doi.org/10.1515/aot-2013-0027.

39. Lenef, A.; Piquette, A.; Kelso, J. **Thermodynamics of Light Extraction from Luminescent Materials.** *ECS J. Solid State Sci. Technol.* **2018**, *7*, R3211–R3226. https://doi.org/10.1149/2.0271801jss.

40. Clarke, F. J. J.; Parry, D. J. **Helmholtz Reciprocity: Its Validity and Application to Reflectometry.** *Light. Res. Technol.* **1985**, *17*, 1–11. https://doi.org/10.1177/14771535850170010301.

41. Neyts, K.; Ullan Nieto, A. **Importance of scattering and absorption for the outcoupling efficiency in organic light-emitting devices.** *J. Opt. Soc. Am. A* **2006**, *23*, 1201. https://doi.org/10.1364/josaa.23.001201.

42. Langer, M. S. **When Shadows Become Interreflections.** *Int. J. Comput. Vis.* **1999**, *34*, 193–204. https://doi.org/10.1023/A:1008131719047.

43. Chance, R. R.; Prock, A.; Silbey, R. **Molecular Fluorescence and Energy Transfer Near Interfaces.** *Adv. Chem. Phys.* **1978**, *37*, 1–65. https://doi.org/10.1002/9780470142561.ch1.

44. Wasey, J. A. E.; Safonov, A.; Samuel, I. D. W.; Barnes, W. L. **Effects of dipole orientation and birefringence on the optical emission from thin films.** *Opt. Commun.* **2000**, *183*, 109–121. https://doi.org/10.1016/s0030-4018(00)00806-3.

45. Moon, C.-K.; Kim, S.-Y.; Lee, J.-H.; Kim, J.-J. **Luminescence from oriented emitting dipoles in a birefringent medium.** *Opt. Express* **2015**, *23*, A279. https://doi.org/10.1364/oe.23.00a279.

46. Mladenovski, S.; Hofmann, S.; Reineke, S.; Penninck, L.; Verschueren, T.; Neyts, K. **Integrated optical model for organic light-emitting devices.** *J. Appl. Phys.* **2011**, *109*, 083114. https://doi.org/10.1063/1.3576114.

47. Henyey, L. G.; Greenstein, J. L. **Diffuse Radiation in the Galaxy.** *Astrophys. J.* **1941**, *93*, 70–83. https://doi.org/10.1086/144246.

48. van de Hulst, H. C. **Multiple Light Scattering: Tables, Formulas, and Applications**, Vols. 1–2; Academic Press: New York, **1980**.

49. Koo, W. H.; Jeong, S. M.; Araoka, F.; Ishikawa, K.; Nishimura, S.; Toyooka, T.; Takezoe, H. **Light Extraction from Organic Light-Emitting Diodes Enhanced by Spontaneously Formed Buckles.** *Nat. Photonics* **2010**, *4*, 222–226. https://doi.org/10.1038/nphoton.2010.7.

50. Thorsos, E. I. **The validity of the Kirchhoff approximation for rough surface scattering using a Gaussian roughness spectrum.** *J. Acoust. Soc. Am.* **1988**, *83*, 78–92. https://doi.org/10.1121/1.396188.

51. Ries, H.; Shatz, N.; Bortz, J.; Spirkl, W. **Performance Limitations of Rotationally Symmetric Nonimaging Devices.** *J. Opt. Soc. Am. A* **1997**, *14*, 2855–2862. https://doi.org/10.1364/JOSAA.14.002855.

52. Lim, J.; Oh, S. S.; Kim, D. Y.; Cho, S. H.; Kim, I. T.; Han, S. H.; Takezoe, H.; Choi, E. H.; et al. **Enhanced out-coupling factor of microcavity organic light-emitting devices with irregular microlens array.** *Opt. Express* **2006**, *14*, 6564. https://doi.org/10.1364/oe.14.006564.
