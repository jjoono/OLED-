# Reciprocal Angular Transport Explains the Limited Benefit of Mixing in OLED Light Extraction

**Authors** — [TBD]

---

## Abstract

We quantify how much angular mixing can improve light extraction from organic light-emitting diodes, using angle-resolved transport constrained by optical reciprocity. With complete reciprocal mixing, recirculating light escapes with the étendue-weighted mean escape probability minus a term proportional to its angular variance. This reduction arises because the escape interface redistributes the light that it reflects. It disappears when the light is mixed in transmission through a lossless volume. For 37 microlens, scattering, randomized-lens, and textured structures, complete mixing by the structure increases the spectrally weighted extraction by at most about two percentage points. A reduced model with one mixing parameter reproduces full transport calculations over 342 cases with a mean absolute error of 0.30 percentage points and fails only for cone arrays in out-of-sample tests.

---

## 1. Introduction

A planar OLED loses much of its generated light to waveguided modes, electrode absorption, and total internal reflection at the substrate/air interface [1,2]. External structures extract the substrate modes without modifying the electrical stack. Examples are microlens arrays (MLAs) [3,4], scattering layers [5,6], and textured films [7]. Their extraction efficiency has been mapped against lens profile, aspect ratio, and fill factor [3,8], and their shapes are increasingly optimized by inverse design and machine learning [9,10]. With such structures, the highest reported external quantum efficiencies of single-unit OLEDs are 55–62% [5,11–13]. This limit has recently been attributed to the loss of returned light on each round trip through the device [14].

However, no compact method exists for comparing the structures themselves. Statistical ray optics assumes that light in the substrate is fully randomized on every pass [15–17]. The substrate-to-air efficiency $\eta = p/[p+(1-p)A']$ then depends only on the escape probability $p$ and the round-trip loss $A'$ [14]. Angle-resolved models propagate the power distribution from pass to pass with redistribution matrices [2,14,18–22]. They are accurate but numerical, and the degree of randomization has been quantified only empirically [23–25]. Because a hemispherical MLA returns unextracted light largely at its original angle [14], light on later passes is concentrated at angles that it transmits poorly. The fully randomized model therefore overestimates extraction by a structure-dependent amount, and it is not known how much improved angular mixing can increase extraction.

We address this question by analyzing the eigenvalues of the return matrix. Reciprocity makes this matrix self-adjoint in the étendue inner product. Within a rank-one return model, reciprocity also determines the angular redistribution of returned light uniquely. The light redistributed into each angle is proportional to the fraction that this angle returns. As a result, completely mixed light escapes with probability $p - \sigma^2/(1-p)$ instead of $p$, and the possible gain from mixing is limited by the angular variance $\sigma^2$ of the escape probability.

The second eigenvalue defines an angular mixing number, the per-pass analogue of the transport factor $1/(1-g)$ in radiative transfer [26] and the relaxation time of a reversible Markov chain [27]. The leading eigenvalue gives the escape probability of the quasi-stationary distribution [28], the escape rate of open dynamical systems [29–31]. These tools are standard in light-transport and reflectance theory [32–36]. Applied to passive light extraction, they give three results. First, a closed-form reference for the gain from complete reciprocal mixing, whose reduction occurs only when the light is mixed in reflection at the escape interface. Second, a reciprocal reduced model, tested against direct ray tracing and unseen structures. Third, a plane of escape probability and mixing number, related to the radiance-ordered bound [37–39] and tested on a freeform-versus-hemisphere benchmark [11,12].

---

## 2. Theory

### 2.1 Multi-pass extraction as a Markov chain on angle

The substrate angle is divided into 1° polar bins $k$. The quantity $v_k$ is the azimuth-integrated power travelling towards the structure in bin $k$, and $w_k = \sin^2\theta_{k+1}-\sin^2\theta_k$ is the étendue of the bin, so that uniform radiance corresponds to $v_k \propto w_k$. At each encounter with the structure, a fraction $p_k$ of the light in bin $k$ escapes, and a fraction $B_R(i\to r)$ returns from bin $i$ to bin $r$. The return matrix is therefore

$$
B_R = M\,\mathrm{diag}(1-\mathbf{p}),
$$

where the column-stochastic matrix $M$ gives the distribution of the returned angle for light that is returned. $M$ is thus the transition matrix of a Markov chain on angle. The escape probability for uniform radiance is $p = \sum_k w_k p_k/\sum_k w_k$. The device reflects the returned light with its angle-resolved reflectance $R(\theta)$, and the substrate-to-air extraction efficiency is

$$
\eta_{\mathrm{sta}} = \mathbf{p}^{\mathsf T}\left(I - \mathrm{diag}\,\mathbf{R}\; B_R\right)^{-1}\mathbf{v}^{(0)} ,
$$

where $\mathbf{v}^{(0)}$ is the power launched into the substrate [14]. All losses enter through $R(\theta)$. The scalar $A' = 1-\sum_k w_k R_k/\sum_k w_k$ is used only for comparison with the fully randomized model.

The model is scalar and averages over azimuth and impact position. Direct $(\theta,\phi)$ tracing agrees with it within 0.6 pp (Sec. 4.3), and separate s and p channels change the conclusions by at most 0.4 pp (Table S7). Note S10 gives the estimators as pseudocode.

### 2.2 Reciprocity and the angular mixing number

For a passive, reciprocal structure [40],

$$
w_i B_R(i\to j) = w_j B_R(j\to i).
$$

This relation has been noted for the redistribution matrices of textured sheets [19] and is the discrete form of the self-adjointness of reciprocal scattering operators in the étendue inner product [32]. $M$ is therefore similar to a symmetric matrix and has a real spectrum [33]. On the full angular domain, its stationary distribution is $\propto w(1-p)$. For the truncated chain used below, the stationary distribution is $\propto w s$, where $s_i$ is the fraction returned into the retained bins (Note S1). The eigenvalue of second-largest modulus, $\lambda_2$, determines how quickly the angular distribution relaxes. We define

$$
N_{\mathrm{mix}} = \frac{1}{1-|\lambda_2|}
$$

as the relaxation time of the angular distribution, expressed in encounters. $N_{\mathrm{mix}}$ is infinite for a planar interface and equals 1 for a Lambertian scatterer. For all structures studied here, $\lambda_2$ is real and positive.

Near-grazing states carry almost no étendue but can retain their angle for many encounters. The spectral metrics are therefore evaluated for $\theta < 85°$ (0.76% of the étendue is excluded), whereas $\eta_{\mathrm{sta}}$ is always computed over 0–90°. Cutoffs between 80° and 90° change the mean error by less than 0.3 pp (Note S1).

### 2.3 Reciprocal reduced model

The reduced model keeps the $p(\theta)$ of the structure. Of the returned light, a fraction $1-m$ keeps its angle and a fraction $m$ is redistributed with a fixed distribution $q$:

$$
B_R(\theta_i\to\theta_r) = \left[1-p(\theta_i)\right]\left[(1-m)\,\delta_{ir} + m\,q(\theta_r)\right].
$$

Detailed balance holds only if $q \propto w(1-p)$. In this case, the light redistributed into each angle is proportional to the fraction that this angle returns. Reciprocal reflectance models that couple a diffuse lobe to angle-dependent specular reflection have the same form [35,36]. To our knowledge, this condition has not been imposed in multi-pass extraction models, which redistribute returned light isotropically [15–17,41]. For a structure whose escape probability depends on angle, isotropic redistribution places part of the returned light in the escape cone, although the escape cone itself returns almost no light. The normalized kernel has eigenvalues 1 and $1-m$, so that $m = 1-\lambda_2$.

### 2.4 Complete reciprocal mixing as a reference

For $m = 1$, the reciprocal kernel has rank one, $B_R(i\to r) = (1-p_i)\,q_r$ with $q \propto w(1-p)$. The escape probability of the mixed light is then

$$
p_{\mathrm{mix}} = \sum_r q_r p_r = p - \frac{\sigma_p^2}{1-p}, \qquad \sigma_p^2 = \frac{\sum_k w_k (p_k - p)^2}{\sum_k w_k}.
$$

This quantity is also the late-pass escape probability (Sec. 2.6). Consequently, $\xi = 1 - \sigma_p^2/[p(1-p)] < 1$ even for complete mixing unless $p(\theta)$ is constant. The multi-pass series can then be summed in closed form,

$$
\eta_{\mathrm{mix}} = \mathbf{p}\cdot\mathbf{v}^{(0)} + \frac{(\mathbf{p}\cdot\mathbf{u})\,[(\mathbf{1}-\mathbf{p})\cdot\mathbf{v}^{(0)}]}{1-(\mathbf{1}-\mathbf{p})\cdot\mathbf{u}}, \qquad u_k = R_k\,q_k .
$$

For an angle-independent reflectance $1-A'$, this expression becomes $\eta_{\mathrm{mix}} = p_0 + (1-p_0)\,p_{\mathrm{mix}}(1-A')/[p_{\mathrm{mix}} + (1-p_{\mathrm{mix}})A']$, where $p_0$ is the first-pass escape probability. It reduces to the fully randomized model only when $p(\theta)$ is constant.

A refractive surface transmits light inside the critical cone and returns all light outside it, so $\sigma_p^2$ is large. For a planar interface, $p$ = 0.40 but $p_{\mathrm{mix}}$ = 0.06. The fully randomized model assumes that mixed light escapes with probability $p$. It therefore overestimates the gain from mixing by an amount determined by the angular variance of escape.

$\eta_{\mathrm{mix}}$ is a reference for reciprocal angular redistribution and is not a universal upper bound. Other reciprocal return operators with the same $p(\theta)$ can extract more light. Note S1 gives a two-state example in which a reciprocal kernel of higher rank exceeds $\eta_{\mathrm{mix}}$ by 5 pp.

### 2.5 Reflective and transmissive mixing

The reduction derived in Sec. 2.4 arises because the escape interface itself redistributes the returned light. Reciprocity weights redistribution in reflection by the return probability of the reflecting element, which is $1-p$ for the escape interface. Two other ways of mixing the light behave differently.

First, the light can be mixed in transmission by a lossless layer with negligible reflection, such as an ideal volume randomizer between the device and the escape interface. Each pass through such a layer distributes the light in proportion to $w$, independently of $p$ and $R$. The mixed light therefore escapes with probability $p$, and complete mixing of this kind reproduces the fully randomized model exactly. In a lossless, index-matched substrate, the position of the layer between the device and the escape interface does not affect the result (Note S11). Real scattering layers also reflect part of the light back towards the device, which lowers their own escape probability (Sec. 4.4).

Second, the device can mix the light in reflection, for example at a diffusely reflecting electrode. The structure then keeps its own return matrix $B_R$, and the device reflection becomes a redistribution kernel $D$. Reciprocity of the device requires $w_i D(i\to j) = w_j D(j\to i)$, and the column sums of $D$ equal the reflectance $R_i$. The rank-one kernel that satisfies both conditions is

$$
D(i\to j) = \frac{R_i\,R_j\,w_j}{\sum_k w_k R_k} .
$$

Light mixed by the device is therefore distributed in proportion to $wR$, and its escape probability at the next encounter is

$$
p_D = \frac{\sum_k w_k R_k\,p_k}{\sum_k w_k R_k} = p + \frac{\mathrm{Cov}_w(R,p)}{\langle R\rangle_w},
$$

where the covariance and the mean are étendue weighted. The variance term of Sec. 2.4 does not appear. If $R$ does not depend on angle, $p_D = p$. When the device absorbs most strongly inside the escape cone, as in strong microcavities, the covariance is negative and mixing by the device is also reduced (Note S11).

### 2.6 Late-pass escape ratio

$\lambda_2$ describes how quickly the angular distribution relaxes, but not the distribution that the light finally reaches. Light that circulates many times approaches the leading eigenvector of $B_R$, the quasi-stationary distribution [28]. This light escapes with probability $p_\infty = 1-\rho(B_R)$ per encounter, which corresponds to the escape rate of an open system [29–31]. The self-adjointness of $B_R$, the Rayleigh quotient with the uniform-radiance vector, and the column sums [42] give

$$
\min_\theta p(\theta) \le p_\infty \le \tilde p ,
$$

where $\tilde p$ is the probability that uniform-radiance light leaves the truncated chain by escape or by return beyond the cutoff ($\tilde p - p < 0.008$). We define the late-pass escape ratio $\xi = p_\infty/\tilde p$. $\xi$ equals 1 when long-circulating light escapes as easily as Lambertian light and approaches 0 when some angles neither escape nor mix. All structures satisfy these bounds numerically (Note S1). $p_\infty$ characterizes the structure alone, in the limit of lossless return by the device. In an OLED, the survival operator $\mathrm{diag}\,\mathbf{R}\,B_R$ also contains parasitic loss, but it ranks the structures in the same order (Spearman correlation 0.96).

---

## 3. Methods

**Devices.** Dipole emission and stack reflectance are computed with the Chance–Prock–Silbey model [43] extended to birefringent layers [44,45] and with transfer matrices. The calculations use measured optical constants and an incoherent substrate [46]. The Al platform of the benchmark is Al (100 nm) / B3PyMPM (61.3 nm) / TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd (25 nm) / TCTA (10 nm) / TAPC (97.5 nm) / ITO (150 nm) / glass ($n = 1.51$), with a horizontal dipole ratio of 0.865. A second stack replaces Al by Ag, and a low-loss stack follows ref. [14]. At 550 nm, $A'$ is 0.167, 0.044, and 0.026 for the three stacks (Table S2, Fig. S1). Device optics are computed for each wavelength with dispersive optical constants. The structure matrices are independent of wavelength (ray optics, non-dispersive index 1.51). Spectrally weighted quantities cover 400–700 nm in 10 nm steps and are weighted with the emitter spectrum multiplied by $\eta_{\mathrm{sub}}(\lambda)$ (Note S10).

**Structures.** All structures are index-matched to the substrate. The 38 structures are: a planar interface; six non-absorbing Henyey–Greenstein [47] scattering slabs ($g$ = 0–0.9, $\tau^*$ = 0.25–4), combined with the exit face by the adding method [48]; the hemispherical and freeform hexagonal MLAs of the benchmark; 17 random MLAs on a 6 × 6 supercell with independent lenslet profiles (spread $s$), heights, or positions and radii; and five Gaussian-correlated wavy surfaces with rms slope 0.2–1.5, which resemble the buckles used on OLEDs [49]. MLAs and textures are traced with a height-field Monte Carlo ray tracer that reproduces the LightTools return matrix of the hemispherical MLA (escape probability 0.4275 compared with 0.4274). Ray optics applies to features much larger than the wavelength [50], which for microlenses means above about 5 µm [8]. Numerical parameters and the residual reciprocity asymmetry change $N_{\mathrm{mix}}$ by at most 0.24, $\xi$ by at most 0.025, and $\eta_{\mathrm{sta}}$ by at most 0.34 pp (Tables S5, S6, Note S10).

**Freeform benchmark.** An axisymmetric freeform MLA with 13 variables (spline profile, height, and ETL and HTL thicknesses) is compared with an equally optimized hemispherical MLA with 3 variables on the Al platform, using LightTools ray tracing validated against fabricated devices [11,12]. The hemisphere is contained in the freeform design space (Table S12). Total EQE and the EQE in four polar bands (0–20°, 20–40°, 40–60°, and 60–80°) are optimized separately by surrogate optimization and pattern search in MATLAB. Each optimum is re-evaluated three times at full fidelity, and a restart control repeats every search from the hemispherical optimum (Note S6).

**Fig. 1 | Multi-pass extraction as a Markov chain on angle.** (a) The passive structure transmits a fraction $p(\theta)$ of the incident light and returns the remainder with $B_R(\theta_i\to\theta_r)$. The device reflects the returned light with $R(\theta)$. (b–d) Return matrices $M$ ($\log_{10}$, floor $10^{-3}$) of a planar interface, the hemispherical MLA, and a scattering layer ($g$ = 0.8, $\tau^*$ = 1). The off-diagonal bands in (c) correspond to deterministic paths of rays that are totally reflected at the lens flank (e.g., 70° → 34°). These bands are unchanged in a 6 × 6 supercell and with twice the number of rays (Note S10). *(`fig1_mixing_concept.png`)*

---

## 4. Results

### 4.1 Efficiency gain from complete reciprocal mixing

Figure 2 shows how much complete reciprocal mixing by the structure itself would increase the efficiency of each structure at its own $p(\theta)$. When the escape probability depends strongly on angle, $p_{\mathrm{mix}}$ is well below $p$ (Fig. 2a): $p_{\mathrm{mix}}/p$ is 0.81–0.93 for lens arrays, 0.37–0.90 for wavy surfaces, and 0.14 for the planar interface.

At 550 nm, complete mixing changes $\eta_{\mathrm{sta}}$ of the 37 non-planar structures by −0.4 to +1.8 pp in the Al device. The maximum increase is 2.0 pp in the Ag device and 0.4 pp in the low-loss device (Fig. 2b). With spectral weighting by the emission spectrum (peak 572 nm), the maximum increases are 1.1, 1.4, and 2.2 pp for the Al, Ag, and low-loss devices. The fully randomized model predicts larger gains, by 0.6–2.9 pp for lens arrays, by 17 pp for the shallowest wavy surface, and by 31 pp for the planar interface. Mixing cannot raise transmission at angles where the structure does not transmit.

The gain from mixing is larger when the cavity emits into angles that the structure transmits poorly. At 500 nm, the first-pass escape probability of the hemispherical MLA decreases to 0.36 in the Al device and 0.33 in the Ag device, compared with $p$ = 0.42. Under these conditions, complete mixing would increase $\eta_{\mathrm{sta}}$ by 5.7–6.4 pp (Al) and 8.4–9.4 pp (Ag), and random lenses and wavy surfaces achieve most of this increase (Table S10). Mixing is therefore most useful when the first-pass angular distribution does not match $p(\theta)$.

Two wavy surfaces (rms slope 0.4 and 0.7) exceed $\eta_{\mathrm{mix}}$ by 0.42 and 0.01 pp. Only in these structures are the escape probabilities at successive encounters anticorrelated (Note S10).

**Fig. 2 | Efficiency gain from complete reciprocal mixing and the role of how the light is mixed.** (a) Escape probability of completely mixed light, $p_{\mathrm{mix}} = p - \sigma_p^2/(1-p)$, versus the uniform-radiance escape probability $p$ for all structures. The dotted line, $p_{\mathrm{mix}} = p$, corresponds to the fully randomized model. (b) Increase of $\eta_{\mathrm{sta}}$ predicted by complete reciprocal mixing at the $p(\theta)$ of each structure (filled symbols) and by the fully randomized model (open symbols), versus $N_{\mathrm{mix}}$ (Al, 550 nm). (c) Difference between the $\eta_{\mathrm{sta}}$ predicted for complete reciprocal mixing in reflection at the device and at the escape interface, versus $p - p_{\mathrm{mix}} = \sigma_p^2/(1-p)$, for all structures (Al, 550 nm). (d) Extraction efficiency of a slab with an internal isotropic source and an angle-independent back reflectance $R$ (solid, 0.90; dashed, 0.98) versus its refractive index. Lines: planar exit face without mixing, with complete reciprocal mixing at the exit face, and with complete mixing at the back reflector, which for constant $R$ equals complete transmissive mixing and the fully randomized model. Symbols: textured exit faces without mixing ($R$ = 0.90, Table S28). *(`fig2_complete_mixing.png`)*

### 4.2 Reflective mixing at the escape interface versus other forms of mixing

Figure 2c compares complete mixing in reflection at the escape interface with complete mixing in reflection at the device (Al, 550 nm). The difference between the two efficiencies increases with $p - p_{\mathrm{mix}} = \sigma_p^2/(1-p)$ (correlation coefficient 0.97). For the planar interface, $\eta_{\mathrm{sta}}$ is 0.384 without mixing, 0.488 with mixing at the escape interface, 0.770 with mixing at the device, and 0.798 with a lossless transmissive randomizer, which equals the fully randomized model. For lens arrays, where $\sigma_p^2$ is small, mixing by the device increases $\eta_{\mathrm{sta}}$ by 1.5–4.1 pp, compared with at most 1.8 pp at the escape interface. Mixing by the device comes within 2.2 pp of the fully randomized model for all non-planar structures in the Al device and within 1.2 pp in the Ag and low-loss devices, because their reflectance is nearly uncorrelated with $p(\theta)$ ($|p_D - p| < 0.008$, Note S11). With spectral weighting, mixing by the device adds 1.1–3.3 pp for lens arrays in the Al device, compared with at most 0.6 pp at the escape interface (Table S29).

The same distinction applies to any slab with an exit face and a back reflector (Fig. 2d). We consider a slab with an internal isotropic source and an angle-independent back reflectance $R$ = 0.90, for which complete mixing at the back reflector equals complete transmissive mixing. For a planar exit face, the extraction efficiency is 0.44 without mixing, 0.60 with complete mixing at the exit face, and 0.87 with complete mixing at the back reflector at $n$ = 1.5, and 0.08, 0.19, and 0.37 at $n$ = 3.5. A textured exit face transmits a similar fraction from every internal angle, so its $\sigma_p^2$ is small (Table S28). Hemispherical-MLA and steep wavy faces therefore come within 1.1 pp of the fully randomized result at $n \ge 2.4$, and complete mixing at the exit face changes the efficiency by at most 2.0 pp. A shallow wavy face (rms slope 0.4) stays up to 5.8 pp below this result. Mixing at its exit face recovers none of the difference, whereas mixing at the back reflector recovers all of it.

### 4.3 Validation of the reduced model

We compare four models with the full multi-pass series for 342 cases (three devices, three wavelengths, and 38 structures). The fully randomized model overestimates $\eta_{\mathrm{sta}}$ of the planar interface by up to 60 pp and that of the hemispherical MLA by 2–3 pp (Fig. S2a). The opposite assumption, that every ray keeps its angle, underestimates $\eta_{\mathrm{sta}}$ of scattering layers (Fig. S2b). The reduced model with $m = 1-\lambda_2$ and isotropic redistribution fits the lens arrays but overestimates shallow wavy surfaces by up to 10.5 pp (Fig. S2c). These surfaces have a small $\lambda_2$, but they return trapped light mostly to other trapped angles. Reciprocal redistribution accounts for this behaviour without an additional parameter. It reduces the mean absolute error to 0.30 pp (mean signed error −0.24 pp) and the maximum error to 3.0 pp (Fig. 3a). The error is not concentrated in any device or wavelength (Fig. 3b, Table S4).

Direct ray tracing in $(\theta, \phi)$ agrees with the full series within 0.6 pp for eight structures (Note S9, Table S19). The remaining differences disappear when the azimuth is randomized at each pass, so they result from the azimuthal correlation that the chain neglects. For test structures not used in developing the model (spherical caps, paraboloids, Gaussian bumps, a cone–hemisphere mixture, and wavy surfaces with a ring spectrum), the mean error is 0.48 pp and the maximum error is 3.3 pp. For two cone arrays, the error reaches 20 pp, because the slow eigenmodes carry little of the returned light (60° cones) or the light is retroreflected within a closed group of angles (45° cones), which a single rate cannot describe (Note S9).

**Fig. 3 | Validation of the reciprocal reduced model.** (a) $\eta_{\mathrm{sta}}$ from the reduced model with $m = 1-\lambda_2$ and reciprocal redistribution ($q \propto w(1-p)$) versus the full series for 342 cases (three devices and three wavelengths). Symbols indicate the structure class. (b) Error of (a) for each device and wavelength. Points are offset horizontally for visibility, and bars show the mean absolute error. The three simpler models are compared in Fig. S2. *(`fig3_reduced_model.png`)*

### 4.4 Escape probability and mixing number of the structures

The structure classes occupy separate regions of the $(p, N_{\mathrm{mix}})$ plane (Fig. 4a). The planar interface has $p$ = 0.40 and $N_{\mathrm{mix}} = \infty$. Scattering layers relax within 1.1–2.3 encounters but have lower escape probabilities, $p$ = 0.19–0.39. Hemispherical MLAs have the highest escape probability, $p$ = 0.42, and $N_{\mathrm{mix}}$ = 6–8. Freeform MLAs have similar escape probabilities, $p$ = 0.42–0.43, and $N_{\mathrm{mix}}$ = 3–7.5.

Randomization lowers $N_{\mathrm{mix}}$ without lowering $p$ (Fig. 4b, Table S3). Independent lenslet profiles reduce $N_{\mathrm{mix}}$ from 7.0 to 1.9, position jitter of 30% of the pitch to 3.5, and a height spread of ±60% only to 4.1 (Table S5). Wavy surfaces with rms slope 0.4–1.5 reach $N_{\mathrm{mix}}$ = 1.1–1.2, as low as a scattering layer.

The late-pass escape ratio $\xi$ provides information that $N_{\mathrm{mix}}$ does not contain (Fig. 4c). For lens arrays, $\xi$ decreases from 0.80 to 0.39 as $N_{\mathrm{mix}}$ increases from 2 to 8, because light that keeps its angle accumulates where the lenses transmit least. The shallowest wavy surface ($\xi$ = 0.35) and a thin isotropic scattering layer ($\xi$ = 0.54) lie below this trend at $N_{\mathrm{mix}} \approx 2$. They relax quickly but return trapped light mainly to other trapped angles, and they have the largest residual errors in Fig. 3a.

The $\eta_{\mathrm{sta}}$ contours in Fig. 4a are nearly vertical. The efficiency is therefore determined mainly by $p$. For the $p(\theta)$ shapes of lenses, scattering layers, and wavy surfaces, the full range of $N_{\mathrm{mix}}$ changes $\eta_{\mathrm{sta}}$ only as much as a change of $p$ by 0.02–0.05 (Table S8). Randomizing the lenslet profiles reduces $N_{\mathrm{mix}}$ from 7.0 to 1.9 but changes $\eta_{\mathrm{sta}}$ by less than 1 pp for the Al device at 550 nm. The gain from faster mixing (+0.9 pp) is offset by the change of $p(\theta)$ (−2.2 pp; Table S9).

**Fig. 4 | Escape probability and mixing number of all structures.** (a) Escape probability and angular mixing number (logarithmic axis) of all structures, with contours of $\eta_{\mathrm{sta}}$ (reciprocal reduced model, Al, 550 nm, hemispherical $p(\theta)$ shape scaled to each $p$). The hemispherical and freeform MLAs are the optima of the benchmark (Table S15). (b) $N_{\mathrm{mix}}$ versus each randomization parameter normalized by its maximum. Dotted lines show the hemispherical MLA (blue) and the scattering layer (green). (c) $\xi$ versus $N_{\mathrm{mix}}$. The planar interface ($N_{\mathrm{mix}} = \infty$, $\xi$ = 0) is shown outside the axis range. *(`fig4_structure_map.png`)*

### 4.5 Case study: freeform versus an equally optimized hemisphere

The results above also apply to a controlled shape-optimization benchmark. The band-optimized freeform profiles differ visibly from the hemisphere (Fig. 5a), but the efficiency gains are small. The ratio $G_j$ of the best freeform EQE to the best hemispherical EQE is 1.003, 1.002, 1.070, and 1.036 for the four bands and 1.013 for total EQE (Fig. 5b, Table 1). For the 40–60° band and for total EQE, these ratios correspond to absolute EQE gains of 1.3 and 0.71 pp.

Independent searches give similar ratios: 1.046–1.070 for 40–60°, 1.020–1.036 for 60–80°, and 1.005–1.013 for total EQE. At 0–40°, no search improves the EQE by more than 0.3%; the lower original values at 0–40° (0.946, 0.966) are search failures, since the hemisphere lies in the freeform design space. We use the spread over independent searches as the uncertainty (Fig. 5c, Fig. S6).

**Table 1 | Freeform versus hemisphere (Al platform, LightTools).** The last column gives the freeform/hemisphere ratio over independent searches, with the number of searches in parentheses (data in Table S13).

| Objective | Hemisphere | Best freeform | $G_j$ | Independent searches |
|---|---|---|---|---|
| 0–20° | 0.06682 | 0.06699 | 1.003 | 0.946–1.003 (2) |
| 20–40° | 0.16591 | 0.16628 | 1.002 | 0.966–1.002 (2) |
| 40–60° | 0.18503 | 0.19802 | 1.070 | 1.046–1.070 (2) |
| 60–80° | 0.13589 | 0.14075 | 1.036 | 1.020–1.036 (2) |
| Total EQE | 0.54679 | 0.5539 | 1.013 | 1.005–1.013 (4) |

On the $(p, N_{\mathrm{mix}})$ plane, the nine optima span only $p$ = 0.417–0.429, which corresponds to about 0.7 pp of $\eta_{\mathrm{sta}}$, whereas $N_{\mathrm{mix}}$ spans 3.0–8.0 and $\xi$ spans 0.39–0.79 (Table S15). Shape freedom therefore changes a design mainly along the $N_{\mathrm{mix}}$ axis, along which the efficiency contours are flat.

Because the benchmark optimizes the cavity together with the lens, we separated the two effects (Note S8). With a fixed stack, a lens-only optimization increases the 40–60° band EQE by 9.0% but $\eta_{\mathrm{sta}}$ by only 1.25%. At unchanged $p$, the optimized lens relaxes faster ($N_{\mathrm{mix}}$ 7.1 → 2.1) and recovers 1.0 of the 1.4 pp that complete mixing would add. At a fixed cavity, the lens shape thus redistributes power among the bands but raises the total only through faster mixing.

**Fig. 5 | Freeform versus an equally optimized hemisphere.** (a) Band-optimized freeform profiles and the hemispherical reference. Profiles are drawn from the optimized control points with a monotone spline. (b) Best EQE of each lens class for each objective, labelled with $G_j$. (c) Freeform/hemisphere ratio for each independent search: the original per-band search (circles), restarts from the hemispherical optimum (squares), and three random starts for total EQE (triangles). Ratios below one indicate search failures, because the hemisphere is contained in the freeform design space. *(`fig6_hemisphere_benchmark.png`)*

---

## 5. Discussion

**Saturation of lens-shape optimization.** A hemispherical array already reaches an escape probability close to the uniform-radiance étendue limit $1/n^2$, so changes of the lens profile mainly change $N_{\mathrm{mix}}$ and the shape of $p(\theta)$, along which the efficiency contours are flat. The anisotropic first-pass emission of the cavity would in principle allow more: for the Al platform at 550 nm, the radiance-ordered bound is 0.540, compared with $1/n^2$ = 0.439 (Fig. S4, Note S4). The lenses, however, transmit almost the same power as for uniform radiance (0.43), because they fill each band from a broad range of internal angles. A surface with binary angular transmission can reach the bound [39], and skewness conservation does not exclude this [51], so reaching it appears to require angle-selective transmission. The angular composition of the output is also nearly independent of the structure: the fraction in 40–60° is 0.32–0.36, close to the Lambertian 0.337 (Fig. S3, Table S11).

**Limited gain from angular mixing.** Optical reciprocity couples angular redistribution to escape: light mixed in reflection by the escape interface is redistributed in proportion to the fraction that each angle returns, so its escape probability is reduced by the variance of $p(\theta)$. Faster mixing therefore does not necessarily increase extraction, and in the structures examined it is often accompanied by a reshaped or lower escape probability. Examples are randomized lenses (Table S9) and scattering layers, whose extraction is highest at intermediate scattering strength [18,21]. For these structures, improving the escape profile is generally more effective than increasing the mixing rate. Transmissive mixing in a lossless volume, or mixing by a device whose reflectance does not depend on angle, avoids the reduction (Sec. 2.5, Sec. 4.2). A scattering volume with little backscattering can therefore approach the fully randomized model, whereas a refractive structure that mixes the light through its own escape interface cannot.

**Design implications.** We suggest reporting three quantities for a new outcoupling structure: $p$, which determines the efficiency; $N_{\mathrm{mix}}$, which indicates how accurately the fully randomized model describes the structure; and $\xi$, which indicates whether long-circulating light can still be extracted. Random MLAs and wavy surfaces with rms slope of one or more combine high $p$, low $N_{\mathrm{mix}}$, and high $\xi$. Their spectrally weighted advantage over periodic MLAs is nevertheless only 0.2–0.5 pp in the Al and Ag devices and 2 pp in the low-loss device. Irregular microlenses have been used on OLEDs mainly to avoid diffraction and colour shift [52]. Where mixing is desired, it is more effective in transmission or at the device than in reflection at the escape interface. Optimizers, including machine-learning methods [9,10], can be evaluated against the $(p, N_{\mathrm{mix}})$ plane.

**Scope.** The analysis applies to passive external structures index-matched to the substrate in the ray-optical regime. The extraction of waveguided modes by internal structures, and the change of the cavity that such structures cause, are not included. Angle-selective filters, absorbing scatterers, and wavelength converters are also not included. Except in Sec. 2.5 and Sec. 4.2, the device is treated as an angle-preserving reflector, as in earlier substrate models [2,18]. For a partly diffuse electrode, the gain increases nearly linearly with the diffuse fraction (Table S30), but the product of the device kernel with $B_R$ need not have a real spectrum [33]. Light trapping in solar cells is the reciprocal problem, in which a low escape probability is desirable, and is not analysed here.

---

## 6. Conclusion

We have quantified how much improved angular mixing can increase the outcoupling efficiency of passive structures. In a reciprocal return model with complete mixing at the escape interface, recirculating light escapes with probability $p - \sigma_p^2/(1-p)$. For 37 scattering, microlens, random-lens, and wavy structures, complete mixing by the structure therefore increases the spectrally weighted extraction by at most about 2 pp, whereas the fully randomized model predicts up to 19 pp. This reference is not a universal upper bound, but it explains the limited gain from faster mixing in common passive structures. Transmissive mixing in a lossless volume and mixing by a reflector with angle-independent reflectance avoid the reduction, also in slabs of higher refractive index. A reduced model based on the same reciprocity condition predicts multi-pass extraction with a mean absolute error of 0.30 pp, agrees with direct ray tracing, and fails only for cone arrays, which a single relaxation rate does not describe. Lens-shape optimization mainly changes $N_{\mathrm{mix}}$ and the angular composition of the output, consistent with gains of a few percent at most. Further improvement will likely require angle-selective structures.

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
