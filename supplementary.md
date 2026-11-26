# Supplement 1

Reflective and transmissive angular mixing affect substrate light extraction differently in organic light-emitting diodes

Notes S1–S4 give the model details, the devices and the 38 structures, the accuracy and convergence checks, and the separation of mixing from escape. Notes S5–S8 document the freeform-versus-hemisphere benchmark. Note S9 reports direct ray tracing, the out-of-sample test, and the prospective test of the reduced model. Note S10 gives the estimators and numerical diagnostics, Note S11 compares reflective mixing at the escape interface with transmissive mixing and with mixing at the device, Note S12 summarizes the freeform benchmark, and Note S13 describes the fabricated devices and the finite-geometry model. The last section lists the scripts that reproduce every figure and table. References of the main text are cited as "main ref. n".

---

## Note S1. Model details and numerical checks

**Polarization and azimuth.** The transport model is scalar. $R(\theta)$ is the polarization average $(|r_s|^2+|r_p|^2)/2$ of the stack reflectance, $\mathbf{v}^{(0)}$ is the sum of the s- and p-polarized power launched into the substrate, and the tracers apply the unpolarized Fresnel coefficients at every encounter. Textures that are not axisymmetric (hexagonal arrays, supercells, wavy fields) are reduced to a polar-angle kernel by launching rays at the center angle of each 1° bin, uniformly over the cell and in azimuth, and binning the exit directions by polar angle only; $B_R$ is the position- and azimuth-averaged return operator.

**Full series.** $\eta_{\mathrm{sta}}$ is obtained by a direct linear solve on the full 0–90° grid. Over all cases the spectral radius of $\mathrm{diag}\,\mathbf{R}\,B_R$ stays below 0.9991, the residual below $10^{-16}$, and a pass-by-pass sum agrees with the solve to $10^{-14}$.

**Angular cutoff.** Bins near 90° carry almost no étendue but can return light at their own angle, and a few such states can dominate $\lambda_2$ without affecting extraction. The spectral metrics are therefore evaluated on $\theta < 85°$, which excludes 0.76% of the étendue and less than 0.1% of $\mathbf{v}^{(0)}$: $B_R$ is restricted to these bins, $M$ is formed by renormalizing the columns of the sub-block, and $\rho$ is taken from the sub-block. $\eta_{\mathrm{sta}}$, in the full series and in the reduced model, always uses 0–90°. On the full domain detailed balance gives the stationary distribution of $M$ as $\pi \propto w(1-p)$. The renormalized sub-block has stationary distribution $\pi^{(S)}_i \propto w_i s_i$ instead, with $s_i = \sum_{j<85°} B_R(i\to j)$ the fraction returned into the retained bins, because $\sum_i [B_R(i\to j)/s_i]\,w_i s_i = \sum_i w_j B_R(j\to i) = w_j s_j$; $s_i$ differs from $1-p_i$ by the light returned beyond the cutoff (up to 0.25 in the bins next to 85° for the hemispherical MLA). The reduced model and the closed-form reference use $q \propto w(1-p)$ on the full domain; only the spectral metrics use the truncated chain. For every structure $\lambda_2$ is real and positive (0.09–0.88), so $|\lambda_2| = \lambda_2$. $N_{\mathrm{mix}}$ is the relaxation time, not the e-folding pass number $-1/\ln\lambda_2$ (2 against 1.44 at $\lambda_2$ = 0.5).

**Table S1.** Sensitivity of $|\lambda_2|$ to the angular range.

| Structure | 0–90° | 0–85° | 0–80° | 0–70° | weight of second eigenvector above 80° |
|---|---|---|---|---|---|
| hemispherical MLA (total-EQE height) | 0.876 | 0.859 | 0.846 | 0.840 | 0.23 |
| freeform MLA (40–60° optimum) | 0.941 | 0.867 | 0.689 | 0.532 | 0.49 |
| random lenslet shapes, s = 0.2 | 0.679 | 0.628 | 0.514 | 0.419 | 0.30 |
| random lenslet heights, spread 0.2 | 0.935 | 0.859 | 0.826 | 0.787 | 0.46 |

For the freeform 40–60° lens and height-randomized hemispheres the full-range value is set by grazing states ($N_{\mathrm{mix}}$ of the 40–60° freeform MLA: 7.5 at 85°, 17 at 90°). The reduced-model error depends little on the cutoff between 80° and 85° and grows above it, when grazing states start to dominate $\lambda_2$:

| Cutoff | mean error (pp) | max error (pp) | hemispherical MLA | freeform MLA | random MLA | wavy | scattering |
|---|---|---|---|---|---|---|---|
| 80° | 0.32 | 3.02 | 0.21 | 0.64 | 0.20 | 0.66 | 0.29 |
| 85° | 0.30 | 2.99 | 0.26 | 0.60 | 0.16 | 0.67 | 0.29 |
| 88° | 0.40 | 3.55 | 0.38 | 0.95 | 0.25 | 0.67 | 0.30 |
| 90° | 0.55 | 4.72 | 0.45 | 1.16 | 0.52 | 0.68 | 0.30 |

(class columns: mean error in pp over three devices and three wavelengths). The cutoff therefore controls the influence of physically minor states (0.76% of the étendue, less than 0.1% of $\mathbf{v}^{(0)}$) rather than tuning the accuracy, which is within 0.02 pp of its best value anywhere between 80° and 85°. The absolute values of $N_{\mathrm{mix}}$ and $\xi$ depend on the cutoff ($\xi$ of the hemispherical MLA is 0.65, 0.47, and 0.31 at 80°, 85°, and 90°), whereas the relative ranking is largely preserved: the rank correlation with the 85° values is 0.96–0.99 for $N_{\mathrm{mix}}$ and 0.90–0.98 for $\xi$ at 80°, 88°, and 90° (Table S6).

**Detailed balance of the traced matrices.** The relative asymmetry $\sum|F-F^{\mathsf T}|/\sum|F|$ of $F_{ji} = w_i B_R(i\to j)$ is 0.05–0.13 on 1° bins, 0.009–0.024 on 5° bins, and 0.002–0.007 on 15° bins, falling as expected for Monte Carlo noise; the planar matrix is exactly symmetric. Symmetrizing $F$ while keeping each column sum changes $N_{\mathrm{mix}}$ by at most 0.20 (2.7%), $\xi$ by 0.005, $\eta_{\mathrm{sta}}$ by 0.04 pp, and the reduced-model error by 0.04 pp over all 37 non-planar structures and three devices.

**Bin width of $p$ and $p_{\mathrm{mix}}$.** Merging the 1° bins (étendue-weighted) into 2°, 3°, and 5° bins leaves $p$ unchanged to $10^{-4}$. $p_{\mathrm{mix}}$ changes by at most 0.0008 for the hemispherical MLA (0.3782 → 0.3790 at 5°), 0.0005 for a random MLA, and 0.0017 for the wavy surface of rms slope 0.2; $\eta_{\mathrm{sta}}$ (Al, 550 nm) changes by at most 0.0006 for these structures up to 3° bins. The planar interface, whose $p(\theta)$ steps at the critical angle, is the most sensitive case: $p_{\mathrm{mix}}$ = 0.054 (native 0.5° Fresnel), 0.056 (1°), 0.058 (2°), and 0.061 (3°), and $\eta_{\mathrm{sta}}$ = 0.3841 (1°) and 0.3864 (3°). Values are in `selectivity/data/rev3_checks.json` (`selectivity/rev3_checks.py`).

**Correlation in main-text Fig. 2(b).** Over the 38 structures (Al, 550 nm), the Pearson correlation between $p - p_{\mathrm{mix}}$ and the device-minus-interface difference of $\eta_{\mathrm{sta}}$ is 0.970 (0.976 without the planar interface, N = 37); the Spearman correlation is 0.949 (0.945 without the planar interface).

**Reciprocal redistribution.** For $B_R(i\to r) = (1-p_i)[(1-m)\delta_{ir} + m q_r]$, detailed balance of the mixing term requires $w_i(1-p_i)q_r = w_r(1-p_r)q_i$, i.e. $q \propto w(1-p)$. The normalized kernel $(1-m)I + m\,q\mathbf{1}^{\mathsf T}$ has eigenvalues 1 and $1-m$, so $m = 1-\lambda_2$ for any $q$.

**Complete reciprocal mixing.** With $q_r = w_r(1-p_r)/\sum_k w_k(1-p_k)$, $p_{\mathrm{mix}} = \sum_r q_r p_r = (p - \langle p^2\rangle)/(1-p) = p - \sigma_p^2/(1-p)$, where $\langle\cdot\rangle$ is the étendue-weighted mean; $p_{\mathrm{mix}} \le p$ is the Cauchy–Schwarz inequality $p^2 \le \langle p^2\rangle$. The rank-one kernel $B_R = q(\mathbf{1}-\mathbf{p})^{\mathsf T}$ with $\mathrm{diag}\,\mathbf{R}$ gives $(I-\mathbf{u}\mathbf{v}^{\mathsf T})^{-1} = I + \mathbf{u}\mathbf{v}^{\mathsf T}/(1-\mathbf{v}^{\mathsf T}\mathbf{u})$ with $u_k = R_k q_k$, $\mathbf{v} = \mathbf{1}-\mathbf{p}$, hence the closed form of main-text Sec. 2.2; it agrees with the numerical series to $10^{-9}$ for all 342 cases. Values for every structure are in `selectivity/data/reciprocal_limit.json`; $p_{\mathrm{mix}}/p$ is 0.87–0.91 for hemispherical MLAs, 0.81–0.87 for freeform MLAs, 0.84–0.93 for random MLAs, 0.37–0.90 for wavy surfaces, 0.59–0.99 for scattering layers, and 0.14 for the planar interface.

**Reciprocity alone does not bound the efficiency.** The complete-mixing kernel is the only reciprocal rank-one kernel with given $p(\theta)$, but other reciprocal kernels exist. Two states with $w$ = (1, 1), $p$ = (0.4, 0), $R$ = 0.8, and all light launched into the second state give, for the rank-one kernel $B$ = [[0.225, 0.375], [0.375, 0.625]], $\eta$ = 0.375; the kernel [[0, 0.6], [0.6, 0.4]] has the same column sums and is also reciprocal, yet sends all light returned from the trapped state into the escaping one and gives $\eta$ = 0.427. $\eta_{\mathrm{mix}}$ is therefore a reference, not a bound; among the structures studied, only two wavy surfaces exceed it, by at most 0.4 pp.

**Bounds on $p_\infty$.** With $D = \mathrm{diag}(w)$, $D^{-1/2}B_R D^{1/2}$ is symmetric, so $\rho(B_R)$ is the maximum of its Rayleigh quotient. With $x = D^{1/2}\mathbf{1}$ the quotient is $\sum_{i,r} w_i B_R(i\to r)/\sum_i w_i = 1-\tilde p$, both sums running over the bins below the cutoff; hence $p_\infty \le \tilde p$. $\tilde p$ is the probability that uniform-radiance light leaves the truncated chain in one encounter, by escape or by return beyond the cutoff, and exceeds $p$ by less than 0.008. The column sums are $1-p(\theta_i)$, so $p_\infty \ge \min_\theta p(\theta)$. All 37 non-planar structures satisfy both bounds, with smallest margins of 0.0014 (upper) and 0.106 (lower); the planar interface meets the lower bound with equality. With $p$ in place of $\tilde p$ the upper bound would fail for the thick scattering layer ($\tau^*$ = 4) by 0.005. $\xi$ is a structure-only descriptor; the device-aware variant $1-\rho(\mathrm{diag}\,\mathbf{R}\,B_R)$ ranks the structures like $p_\infty$ (Spearman 0.956, Al, 550 nm). A two-parameter reciprocal kernel fitted to both $\lambda_2$ and $\rho$, and a two-rate scalar formula, were not more accurate overall than the one-parameter model (0.32 pp mean) and are not used.

**Band-resolved bound.** The substrate bins are sorted by radiance $v^{(0)}_k/w_k$; the étendue of an air band, mapped into the substrate by Snell's law, is $(\sin^2\theta_{\mathrm{hi}}-\sin^2\theta_{\mathrm{lo}})/n^2$ in units of $w$, and it is filled with the brightest bins first, the last one fractionally. The bound is lossless, polarization-independent, and integrated over azimuth. Each band's bound is a separate maximum, so the per-band bounds (sum 0.676) do not add up to the bound for the whole hemisphere (0.540), which fills the étendue $1/n^2$ once; the uniform-radiance bound is additive (0.439, of which 0.014 lies above 80°).

---

## Note S2. Devices and the 38 structures

**Devices.** The emitters form a single dipole sheet at the center of the 25 nm birefringent EML (measured ordinary and extraordinary indices), with a horizontal dipole ratio of 0.865 (Table S2, Fig. S1).

**Table S2.** Device stacks (from the cathode to the substrate).

| Device | Cathode | ETL | EML | HTL-side layers | Anode | Substrate | $A'$ at 550 nm |
|---|---|---|---|---|---|---|---|
| Al | Al 100 nm | B3PyMPM 61.3 nm | TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd 25 nm | TCTA 10 nm / TAPC 97.5 nm | ITO 150 nm | glass, $n$ = 1.51 | 0.167 |
| Ag | Ag 100 nm | B3PyMPM 61.3 nm | as Al | TCTA 10 nm / TAPC 97.5 nm | ITO 150 nm | glass, $n$ = 1.51 | 0.044 |
| Low loss | Ag 100 nm | B3PyMPM 200 nm | as Al | TAPC 180 nm | ITO 50 nm | glass, $n$ = 1.51 | 0.026 |

**Fig. S1.** Device inputs at 550 nm. (a) Radiance launched into the substrate on the first pass, $v^{(0)}_k/w_k$, relative to its mean over 0–20°. (b) Angle-resolved reflectance $R(\theta)$ seen from the substrate, with the uniform-radiance round-trip loss $A'$. Dotted line, critical angle of the glass/air interface. *(`figS11_device_inputs.png`)*

**Height-field tracer.** The surface $z = h(x, y)$ separates the index-matched medium ($n = 1.51$) from air; rays are marched with a Lipschitz-safe step, refined by bisection, and reflected or refracted with the unpolarized Fresnel coefficients. 20,000 rays per 1° incidence bin are used (10,000 for wavy surfaces). For a hemispherical lenslet the tracer reproduces the LightTools return matrix of main ref. 14: escape probability 0.4275 against 0.4274, transmittance within 0.008 at every angle, and fewer than 0.1% of rays undetermined. The freeform profiles are rebuilt from the optimized control points with a PCHIP spline.

**Random lens arrays** tile a 6 × 6 supercell periodically. Shape spread $s$ perturbs the five control points of the quarter circle by Gaussian noise of standard deviation $0.25s$ (radial) and $0.4s$ (height), so that $s$ = 1 approaches independent profiles; height spread scales each lenslet by $1\pm s$; position jitter displaces each center by up to a fraction of the pitch, optionally with ±20% radius spread. Where lenslets overlap the higher surface is taken; profiles with near-vertical walls are rejected. **Wavy surfaces** are periodic Gaussian-correlated height fields on a 12 × 12 box (correlation length one), scaled to a prescribed rms slope; only the slope statistics matter in ray optics. **Scattering layers** are non-absorbing, index-matched Henyey–Greenstein slabs (20,000 rays per 0.5° bin) combined with the flat exit face by the adding method.

**Table S3.** The 38 structures at 550 nm. $p$: escape probability for uniform radiance; $N_{\mathrm{mix}}$ and $p_\infty$ on $\theta < 85°$; $\xi = p_\infty/\tilde p$; $\eta$: full series for the Al, Ag, and low-loss (LL) devices; $S_{40\text{–}60}$: 40–60° share of the multi-pass output (Al); last two columns: reciprocal and isotropic ($q \propto w$) reduced model minus full series (Al). A, lenslet shape spread $s$; B, height spread ±$s$; C, position jitter (fraction of the pitch); C2, as C with ±20% radius spread; W, wavy surface of the given rms slope. Scattering layers are labeled by $g$ and $\tau^*$, MLAs by the benchmark objective of the optimum they reproduce.

| Structure | p | N_mix | p∞ | ξ | η Al | η Ag | η LL | S40–60 (Al) | err. recip. Al (pp) | err. q∝w Al (pp) |
|---|---|---|---|---|---|---|---|---|---|---|
| planar | 0.398 | ∞ | 0.000 | 0.00 | 0.384 | 0.356 | 0.664 | 0.334 | -0.00 | -0.00 |
| scatter g0.0_t0.25 | 0.381 | 2.2 | 0.212 | 0.54 | 0.675 | 0.874 | 0.945 | 0.349 | -2.07 | +4.21 |
| scatter g0.0_t1.0 | 0.324 | 1.1 | 0.309 | 0.93 | 0.725 | 0.907 | 0.949 | 0.353 | -0.10 | +1.24 |
| scatter g0.8_t0.25 | 0.386 | 1.2 | 0.320 | 0.81 | 0.749 | 0.914 | 0.959 | 0.359 | -0.19 | +3.04 |
| scatter g0.8_t1.0 | 0.326 | 1.2 | 0.326 | 0.98 | 0.737 | 0.913 | 0.951 | 0.354 | -0.08 | +0.34 |
| scatter g0.9_t1.0 | 0.326 | 1.2 | 0.326 | 0.98 | 0.738 | 0.913 | 0.951 | 0.354 | -0.10 | +0.29 |
| scatter g0.8_t4.0 | 0.188 | 1.2 | 0.194 | 0.99 | 0.579 | 0.838 | 0.901 | 0.353 | -0.06 | +0.14 |
| MLA free_20-40 | 0.429 | 3.8 | 0.251 | 0.58 | 0.788 | 0.925 | 0.967 | 0.352 | -0.69 | +0.32 |
| MLA free_40-60 | 0.425 | 7.5 | 0.167 | 0.39 | 0.789 | 0.929 | 0.961 | 0.360 | -1.11 | -0.44 |
| MLA free_60-80 | 0.417 | 3.0 | 0.337 | 0.79 | 0.785 | 0.929 | 0.966 | 0.349 | -1.46 | -0.21 |
| MLA free_total | 0.424 | 4.6 | 0.199 | 0.46 | 0.771 | 0.914 | 0.961 | 0.341 | -0.89 | +0.44 |
| MLA hemi_0-20 | 0.423 | 7.3 | 0.190 | 0.44 | 0.790 | 0.924 | 0.964 | 0.328 | -0.38 | +0.13 |
| MLA hemi_20-40 | 0.422 | 7.0 | 0.203 | 0.47 | 0.793 | 0.926 | 0.964 | 0.327 | -0.21 | +0.23 |
| MLA hemi_40-60 | 0.425 | 8.0 | 0.181 | 0.42 | 0.783 | 0.920 | 0.965 | 0.339 | -0.64 | -0.03 |
| MLA hemi_60-80 | 0.420 | 6.2 | 0.214 | 0.50 | 0.792 | 0.926 | 0.965 | 0.331 | -0.24 | +0.20 |
| MLA hemi_total | 0.422 | 7.1 | 0.199 | 0.47 | 0.792 | 0.925 | 0.964 | 0.327 | -0.27 | +0.19 |
| A_shape0.0 | 0.422 | 7.0 | 0.214 | 0.50 | 0.794 | 0.927 | 0.965 | 0.328 | -0.15 | +0.24 |
| A_shape0.1 | 0.421 | 3.8 | 0.279 | 0.65 | 0.797 | 0.931 | 0.966 | 0.336 | -0.30 | +0.27 |
| A_shape0.2 | 0.421 | 2.7 | 0.334 | 0.78 | 0.799 | 0.934 | 0.966 | 0.349 | -0.28 | +0.39 |
| A_shape0.4 | 0.422 | 2.5 | 0.343 | 0.80 | 0.794 | 0.932 | 0.966 | 0.358 | -0.40 | +0.57 |
| A_shape0.7 | 0.426 | 1.9 | 0.344 | 0.79 | 0.788 | 0.929 | 0.966 | 0.358 | -0.27 | +1.31 |
| A_shape1.0 | 0.428 | 1.9 | 0.335 | 0.77 | 0.785 | 0.927 | 0.966 | 0.354 | -0.23 | +1.47 |
| B_aspect0.2 | 0.422 | 7.1 | 0.184 | 0.43 | 0.795 | 0.928 | 0.965 | 0.328 | -0.20 | +0.18 |
| B_aspect0.4 | 0.422 | 6.0 | 0.215 | 0.50 | 0.797 | 0.930 | 0.966 | 0.332 | -0.18 | +0.20 |
| B_aspect0.6 | 0.423 | 4.1 | 0.294 | 0.69 | 0.799 | 0.933 | 0.968 | 0.341 | -0.19 | +0.24 |
| C2_jitter0.05_r20 | 0.421 | 4.2 | 0.281 | 0.66 | 0.798 | 0.931 | 0.966 | 0.331 | -0.12 | +0.34 |
| C2_jitter0.1_r20 | 0.421 | 3.5 | 0.300 | 0.70 | 0.796 | 0.930 | 0.966 | 0.329 | -0.15 | +0.45 |
| C2_jitter0.2_r20 | 0.420 | 3.7 | 0.294 | 0.69 | 0.799 | 0.932 | 0.966 | 0.335 | -0.13 | +0.35 |
| C2_jitter0.3_r20 | 0.420 | 3.2 | 0.314 | 0.74 | 0.796 | 0.932 | 0.966 | 0.341 | -0.21 | +0.40 |
| C_jitter0.05 | 0.421 | 6.2 | 0.215 | 0.50 | 0.796 | 0.928 | 0.965 | 0.328 | -0.15 | +0.25 |
| C_jitter0.1 | 0.421 | 5.0 | 0.245 | 0.57 | 0.797 | 0.930 | 0.966 | 0.331 | -0.15 | +0.29 |
| C_jitter0.2 | 0.421 | 3.9 | 0.281 | 0.66 | 0.798 | 0.932 | 0.966 | 0.336 | -0.18 | +0.32 |
| C_jitter0.3 | 0.420 | 3.5 | 0.298 | 0.70 | 0.798 | 0.932 | 0.966 | 0.340 | -0.20 | +0.34 |
| W_slope0.2 | 0.402 | 2.2 | 0.142 | 0.35 | 0.608 | 0.815 | 0.933 | 0.334 | -2.71 | +9.57 |
| W_slope0.4 | 0.410 | 1.2 | 0.264 | 0.63 | 0.719 | 0.896 | 0.957 | 0.348 | -1.27 | +6.47 |
| W_slope0.7 | 0.417 | 1.2 | 0.348 | 0.81 | 0.773 | 0.923 | 0.965 | 0.335 | -0.37 | +2.78 |
| W_slope1.0 | 0.422 | 1.1 | 0.377 | 0.87 | 0.792 | 0.932 | 0.968 | 0.329 | -0.07 | +1.70 |
| W_slope1.5 | 0.427 | 1.2 | 0.387 | 0.88 | 0.801 | 0.935 | 0.968 | 0.337 | -0.09 | +1.21 |

---

## Note S3. Accuracy and convergence

**Table S4.** Error of the reduced and fully randomized models by structure class (pp; 342 cases, three devices, 500/550/600 nm).

| Class | cases | reciprocal mean | reciprocal max | q ∝ w mean | q ∝ w max | fully randomized mean | fully randomized max |
|---|---|---|---|---|---|---|---|
| Planar | 9 | 0.00 | 0.00 | 0.00 | 0.00 | 44.05 | 59.89 |
| Scattering layers | 54 | 0.29 | 2.49 | 0.96 | 4.67 | 1.83 | 11.14 |
| Hemispherical and freeform MLAs | 81 | 0.41 | 2.09 | 0.38 | 2.22 | 3.45 | 12.33 |
| Random lens arrays (A, B, C) | 153 | 0.16 | 0.51 | 0.44 | 1.93 | 2.41 | 10.89 |
| Wavy, rms slope ≥ 0.7 | 27 | 0.19 | 1.06 | 1.11 | 3.07 | 1.61 | 5.57 |
| Wavy, rms slope ≤ 0.4 | 18 | 1.40 | 2.99 | 5.40 | 10.51 | 8.66 | 21.73 |
| All | 342 | 0.30 | 2.99 | 0.81 | 10.51 | 3.93 | 59.89 |

Fig. S2 compares the three simpler models of Table S4 with the full series.

**Fig. S2.** Simpler models against the full series. $\eta_{\mathrm{sta}}$ from (a) the fully randomized model, (b) the angle-preserving limit ($m = 0$), and (c) the reduced model with $m = 1-\lambda_2$ and isotropic redistribution ($q \propto w$), versus the full series (342 cases: three devices and three wavelengths). Symbols indicate the structure class. Main-text Fig. 3(a) shows the reciprocal reduced model. *(`figS_model_comparison.png`)*

By device, the mean (maximum) reciprocal-model error is 0.41 (2.7) pp for Al, 0.33 (3.0) pp for Ag, and 0.16 (1.8) pp for the low-loss stack; by wavelength 0.38 (2.5), 0.29 (3.0), and 0.23 (2.6) pp at 500, 550, and 600 nm. The model tends to slightly underestimate $\eta_{\mathrm{sta}}$.

**Monte Carlo uncertainty.** Ray noise was estimated by re-tracing with two further ray seeds, and realization spread by drawing two further geometries with the same statistics (Table S5). Ray noise is negligible; realization spread is at most 0.24 in $N_{\mathrm{mix}}$, 0.021 in $\xi$, and 0.23 pp in $\eta_{\mathrm{sta}}$, and changes the reduced-model error by at most 0.03 pp.

**Table S5.** Monte Carlo uncertainty (550 nm). ±: standard deviation over three ray seeds; in parentheses: over three geometric realizations. Last column: reciprocal reduced model minus full series (Al).

| Structure | N_mix | ξ | η Al | err. Al (pp) |
|---|---|---|---|---|
| A_shape0.0 | 7.05 ± 0.01 (0.00) | 0.499 ± 0.003 (0.000) | 0.7944 ± 0.0001 (0.0000) | -0.15 ± 0.01 (0.00) |
| A_shape1.0 | 1.95 ± 0.00 (0.12) | 0.769 ± 0.000 (0.009) | 0.7853 ± 0.0003 (0.0023) | -0.23 ± 0.00 (0.03) |
| B_aspect0.6 | 4.08 ± 0.00 (0.24) | 0.686 ± 0.002 (0.021) | 0.7994 ± 0.0002 (0.0001) | -0.19 ± 0.01 (0.01) |
| C_jitter0.3 | 3.50 ± 0.01 (0.10) | 0.699 ± 0.001 (0.006) | 0.7979 ± 0.0001 (0.0005) | -0.20 ± 0.00 (0.01) |
| W_slope0.2 | 2.17 ± 0.01 (0.01) | 0.345 ± 0.001 (0.001) | 0.6084 ± 0.0006 (0.0008) | -2.71 ± 0.01 (0.02) |
| W_slope1.0 | 1.13 ± 0.00 (0.01) | 0.871 ± 0.001 (0.007) | 0.7924 ± 0.0002 (0.0023) | -0.07 ± 0.02 (0.01) |
| scatter g0.8_t1.0 | 1.17 ± 0.00 | 0.977 ± 0.000 | 0.7374 ± 0.0000 | -0.08 ± 0.00 |
| scatter g0.0_t0.25 | 2.25 ± 0.01 | 0.545 ± 0.000 | 0.6754 ± 0.0001 | -2.07 ± 0.01 |

**Table S6.** Convergence with bin width, cutoff, ray count, and domain size (Al, 550 nm). Wider bins are obtained by étendue-weighted merging.

| Structure | Variant | $N_{\mathrm{mix}}$ | $\xi$ | $\eta_{\mathrm{sta}}$ | reduced − full (pp) |
|---|---|---|---|---|---|
| planar | 1° bins, cutoff 85° (reference) | ∞ | 0.000 | 0.3841 | +0.00 |
| planar | 3° bins, cutoff 85° | ∞ | 0.000 | 0.3864 | +0.00 |
| hemispherical MLA | 1° bins, cutoff 85° (reference) | 7.10 | 0.468 | 0.7918 | -0.27 |
| hemispherical MLA | 3° bins, cutoff 85° | 6.81 | 0.519 | 0.7920 | -0.25 |
| hemispherical MLA | 1° bins, cutoff 80° | 6.49 | 0.651 | 0.7918 | -0.22 |
| hemispherical MLA | 1° bins, no cutoff | 8.08 | 0.306 | 0.7918 | -0.35 |
| scattering layer ($g$ 0.8, $\tau^*$ 1) | 0.5° bins (native), cutoff 85° | 1.17 | 0.995 | 0.7374 | -0.08 |
| scattering layer ($g$ 0.8, $\tau^*$ 1) | 1° bins, cutoff 85° (reference) | 1.17 | 0.995 | 0.7374 | -0.08 |
| scattering layer ($g$ 0.8, $\tau^*$ 1) | 1° bins, no cutoff | 1.23 | 0.971 | 0.7374 | -0.11 |
| random lenslet shapes, $s$ = 1 | 1° bins, cutoff 85° (reference) | 1.95 | 0.777 | 0.7853 | -0.23 |
| random lenslet shapes, $s$ = 1 | 3° bins, cutoff 85° | 1.89 | 0.795 | 0.7854 | -0.19 |
| random lenslet shapes, $s$ = 1 | 1° bins, no cutoff | 2.13 | 0.720 | 0.7853 | -0.34 |
| wavy, rms slope 0.2 | 1° bins, cutoff 85° (reference) | 2.17 | 0.351 | 0.6084 | -2.71 |
| wavy, rms slope 0.2 | 3° bins, cutoff 85° | 2.16 | 0.369 | 0.6090 | -2.67 |
| wavy, rms slope 0.2 | 1° bins, no cutoff | 2.18 | 0.307 | 0.6084 | -2.75 |
| wavy, rms slope 1.0 | 1° bins, cutoff 85° (reference) | 1.13 | 0.886 | 0.7924 | -0.07 |
| wavy, rms slope 1.0 | 1° bins, no cutoff | 1.14 | 0.868 | 0.7924 | -0.08 |
| hemispherical array ($s$ = 0) | 5,000 rays per bin | 7.08 | 0.502 | 0.7944 | -0.17 |
| hemispherical array ($s$ = 0) | 20,000 rays per bin (main set) | 7.05 | 0.499 | 0.7944 | -0.15 |
| hemispherical array ($s$ = 0) | 40,000 rays per bin | 7.02 | 0.497 | 0.7943 | -0.14 |
| wavy, rms slope 0.2 | 2,500 rays per bin | 2.17 | 0.345 | 0.6094 | -2.76 |
| wavy, rms slope 0.2 | 10,000 rays per bin (main set) | 2.17 | 0.345 | 0.6084 | -2.71 |
| random lenslet shapes, $s$ = 1 | 6 × 6 supercell (main set) | 1.95 | 0.769 | 0.7853 | -0.23 |
| random lenslet shapes, $s$ = 1 | 10 × 10 supercell, new realization | 1.89 | 0.794 | 0.7887 | -0.17 |
| position jitter 0.3 | 6 × 6 supercell (main set) | 3.50 | 0.699 | 0.7979 | -0.20 |
| position jitter 0.3 | 10 × 10 supercell, new realization | 3.46 | 0.703 | 0.7962 | -0.21 |
| wavy, rms slope 0.2 | 12 × 12 box (main set) | 2.17 | 0.345 | 0.6084 | -2.71 |
| wavy, rms slope 0.2 | 24 × 24 box, new realization | 2.17 | 0.347 | 0.6098 | -2.76 |
| wavy, rms slope 1.0 | 12 × 12 box (main set) | 1.13 | 0.871 | 0.7924 | -0.07 |
| wavy, rms slope 1.0 | 24 × 24 box, new realization | 1.13 | 0.867 | 0.7912 | -0.04 |

$\eta_{\mathrm{sta}}$ changes by at most 0.3 pp (planar, 3° bins) and 0.06 pp for all other structures. The planar interface is the one structure whose returned light never changes angle, so its $\eta_{\mathrm{sta}}$ depends on how the Fresnel edge at the critical angle (41.47°) is resolved within a bin: with 1° bins and uniform radiance within each bin it is 0.384 (Al, 550 nm), whereas an exact per-angle sum at 0.01° gives 0.380 (0.380 with s and p separate). The 0.4 pp difference comes from the bin containing the critical angle (41.47°) and does not affect any statement about the planar interface; for all other structures, which redistribute angle, bin effects are at most 0.06 pp; a 10 × 10 supercell or 24 × 24 wavy box changes $N_{\mathrm{mix}}$ by at most 0.06, $\xi$ by 0.025, and $\eta_{\mathrm{sta}}$ by 0.34 pp, comparable to the realization spread.

**Polarization.** The scalar model averages $R$ over s and p and sums $\mathbf{v}^{(0)}$, and the tracers use unpolarized Fresnel coefficients. The opposite limit treats s and p as separate channels that keep their label at every encounter: $\eta = \sum_{s,p} \mathbf{p}_{\sigma}^{\mathsf T}(I-\mathrm{diag}\,\mathbf{R}_{\sigma}B_{R,\sigma})^{-1}\mathbf{v}^{(0)}_{\sigma}$, with return matrices traced with $|r_s|^2$ or $|r_p|^2$ only and polarization-resolved $\mathbf{v}^{(0)}_\sigma$ and $R_\sigma$ (s carries 65–70% of $\mathbf{v}^{(0)}$). The limit is exact for a planar interface; for curved surfaces, which rotate the local plane of incidence, the true value lies between the two (Table S7). The polarization treatment shifts the absolute $\eta_{\mathrm{sta}}$ of the hemispherical MLA and the wavy surface by 0.6–1.5 pp, but the complete-mixing gain by at most 0.4 pp and the reduced-model error by at most 0.1 pp, so the comparisons in the main text are unaffected. For the planar interface the separate-channel result is exact and agrees with the scalar model within 0.11 pp; only its hypothetical complete-mixing reference is 2.6–3.1 pp higher when each polarization is mixed separately, which does not change any conclusion.

**Table S7.** Scalar model versus polarization-resolved limit without s/p conversion (Al). Last two columns: scalar / polarization-resolved, pp. Planar matrices are integrated within each 1° bin (0.384 at 550 nm, as in Table S3; 0.380 for an exact per-angle sum, Note S3).

| Structure | λ (nm) | $\eta_{\mathrm{sta}}$, scalar | $\eta_{\mathrm{sta}}$, s and p separate | difference (pp) | complete-mixing gain | reduced − full |
|---|---|---|---|---|---|---|
| planar interface | 500 | 0.3658 | 0.3649 | -0.09 | 10.09 / 12.72 | +0.00 / +0.00 |
| planar interface | 550 | 0.3837 | 0.3834 | -0.03 | 10.37 / 13.16 | +0.00 / +0.00 |
| planar interface | 600 | 0.4165 | 0.4154 | -0.11 | 8.27 / 11.42 | +0.00 / +0.00 |
| hemispherical MLA | 500 | 0.7019 | 0.7081 | +0.62 | 5.86 / 5.99 | -0.61 / -0.60 |
| hemispherical MLA | 550 | 0.7913 | 0.8038 | +1.25 | 1.41 / 1.05 | -0.28 / -0.30 |
| hemispherical MLA | 600 | 0.8032 | 0.8156 | +1.25 | 0.13 / -0.02 | -0.04 / -0.04 |
| wavy, rms slope 1.0 | 500 | 0.7546 | 0.7637 | +0.91 | 0.12 / 0.04 | -0.31 / -0.42 |
| wavy, rms slope 1.0 | 550 | 0.7916 | 0.8015 | +0.99 | 0.12 / 0.06 | -0.07 / -0.10 |
| wavy, rms slope 1.0 | 600 | 0.7906 | 0.8052 | +1.46 | 0.08 / 0.03 | +0.01 / -0.03 |

---

## Note S4. Mixing versus escape and the angular composition of the output

**Contour templates.** Fig. 4(a) of the main text scales the $p(\theta)$ shape of the hemispherical MLA. Table S8 repeats the calculation with the shapes of a scattering layer and a wavy surface: at fixed $p$, the whole range of $N_{\mathrm{mix}}$ is worth 1.2–3.2 pp, whereas 0.01 in $p$ is worth about 0.6 pp for every template.

**Table S8.** Sensitivity of $\eta_{\mathrm{sta}}$ to mixing and escape (Al, 550 nm). Span: change of $\eta_{\mathrm{sta}}$ for $N_{\mathrm{mix}}$ from 1 to 16 at fixed $p$. Slope: change per 0.01 in $p$ at $N_{\mathrm{mix}}$ = 4.

| Template | span at $p$ = 0.30 (pp) | at 0.38 | at 0.42 | slope per 0.01 in $p$ (pp) |
|---|---|---|---|---|
| hemispherical MLA | 2.05 | 2.10 | 2.09 | 0.61 |
| scattering layer ($g$ 0.8, $\tau^*$ 1) | 1.40 | 1.30 | 1.23 | 0.62 |
| wavy surface, rms slope 1.0 | 3.22 | 3.23 | 3.17 | 0.59 |

**Decomposition.** Along the shape sweep, the change of $\eta_{\mathrm{sta}}$ relative to the hemispherical array ($s$ = 0) is split in the reciprocal reduced model into a mixing part ($m$ of the randomized array, $p(\theta)$ of the hemisphere), an escape part (the reverse), and their interaction (Table S9). On Al the escape part dominates at large spread ($-2.20$ against $+0.91$ pp at $s$ = 1), on Ag the two cancel, and on the low-loss stack both stay below 0.15 pp.

**Table S9.** Decomposition of the change in $\eta_{\mathrm{sta}}$ along the shape sweep (550 nm, pp).

| $s$ | Device | mixing only | $p(\theta)$ only | interaction | total (reduced) | total (full series) |
|---|---|---|---|---|---|---|
| 0.2 | Al | +0.65 | -0.22 | -0.08 | +0.35 | +0.48 |
| 0.2 | Ag | +0.68 | +0.15 | -0.15 | +0.68 | +0.71 |
| 0.2 | low loss | +0.09 | +0.13 | -0.02 | +0.20 | +0.15 |
| 0.4 | Al | +0.72 | -0.99 | -0.07 | -0.34 | -0.09 |
| 0.4 | Ag | +0.74 | -0.11 | -0.19 | +0.44 | +0.52 |
| 0.4 | low loss | +0.09 | +0.05 | -0.02 | +0.13 | +0.10 |
| 1.0 | Al | +0.91 | -2.20 | +0.30 | -0.98 | -0.91 |
| 1.0 | Ag | +0.91 | -0.94 | +0.13 | +0.10 | +0.07 |
| 1.0 | low loss | +0.11 | +0.07 | -0.01 | +0.18 | +0.13 |

**Mixing gain by device and wavelength.** The complete-mixing reference depends on how the first-pass distribution $\mathbf{v}^{(0)}$ overlaps $p(\theta)$ (Table S10). When the cavity launches light into angles that the hemispherical MLA transmits poorly (first-pass escape probability 0.33–0.38 against $p$ = 0.42: Al and Ag at 500 nm, low loss at 600 nm), complete mixing would add 3–9 pp, and random lenses and wavy surfaces realize most of it. Elsewhere the reference is below about 2 pp. Weighted with the emission spectrum of the emitter (Ir(dmppy-ph)$_2$tmd, peak 572 nm; relative weight 0.01 at 500 nm, 0.52 at 550 nm, 0.71 at 600 nm) and with $\eta_{\mathrm{sub}}(\lambda)$, the reference is at most 1.1, 1.4, and 2.2 pp for the Al, Ag, and low-loss devices.

**Table S10.** Complete-mixing gain and realized randomization gain (pp). First pass: first-pass escape probability $\mathbf{p}^{\mathsf T}\mathbf{v}^{(0)}$ of the hemispherical MLA ($p$ = 0.422). Random − hemisphere: best random MLA minus the hemispherical array; wavy − hemisphere: wavy surface of slope 1.5 minus the hemispherical MLA.

| Device, wavelength | first pass | gain, hemispherical MLAs | gain, all non-planar | best random − hemisphere | wavy (1.5) − hemisphere |
|---|---|---|---|---|---|
| Al, 500 nm | 0.363 | 5.7 to 6.4 | -0.6 to 6.4 | +3.6 | +5.8 |
| Al, 550 nm | 0.433 | 1.3 to 1.7 | -0.4 to 1.8 | +0.5 | +0.9 |
| Al, 600 nm | 0.462 | 0.1 to 0.3 | -0.2 to 1.2 | +0.0 | -0.3 |
| Al, spectrally weighted | — | 0.4 to 0.6 | -0.2 to 1.1 | +0.2 (best random or wavy) | — |
| Ag, 500 nm | 0.326 | 8.4 to 9.4 | -0.4 to 9.4 | +5.7 | +8.4 |
| Ag, 550 nm | 0.413 | 1.1 to 1.5 | -0.3 to 2.0 | +0.7 | +1.0 |
| Ag, 600 nm | 0.462 | 0.3 to 0.4 | -0.1 to 1.2 | +0.0 | +0.1 |
| Ag, spectrally weighted | — | 0.6 to 0.9 | -0.2 to 1.4 | +0.5 (best random or wavy) | — |
| low loss, 500 nm | 0.467 | 0.2 to 0.3 | -0.1 to 0.8 | +0.1 | +0.1 |
| low loss, 550 nm | 0.458 | 0.1 to 0.1 | -0.0 to 0.4 | +0.3 | +0.4 |
| low loss, 600 nm | 0.382 | 3.2 to 4.0 | -0.1 to 4.3 | +3.0 | +4.0 |
| low loss, spectrally weighted | — | 1.8 to 2.2 | -0.1 to 2.2 | +2.0 (best random or wavy) | — |

**Band-resolved radiance-ordered bound.** In one encounter, the power that a passive structure can transmit into an output band is bounded by the power obtained when the étendue of the band is filled with the brightest internal radiance first (main refs. 48 and 49). For uniform radiance, this bound equals the étendue ratio, which is $1/n^2$ for the hemisphere. We evaluate it for each air-side polar band with the actual first-pass radiance of the microcavity, filling the étendue of each band with the substrate bins of highest $v^{(0)}_k/w_k$ first. For the structures, the band power is accumulated through the traced transmission matrix into air bins. For the Al platform at 550 nm, the bound allows 0.085, 0.205, 0.230, and 0.156 of the substrate power into the 0–20°, 20–40°, 40–60°, and 60–80° bands, and 0.540 for the whole hemisphere. MLAs reach 76–80% of the hemispherical bound, wavy surfaces 77%, the planar interface 65%, and scattering layers 35–68% (Fig. S4(a)).

**Angular composition of the output.** Fig. S4 shows the band-resolved radiance bound and the composition of the output for the Al device at 550 nm. $S_{40\text{–}60} = P_{\mathrm{air}}(40\text{–}60°)/P_{\mathrm{air}}(0\text{–}90°)$, both powers in air including Fresnel losses, over all 38 structures (Table S11). In the low-loss device the lower end (0.274–0.293) belongs to the hemispherical MLAs; all other structures stay within 0.31–0.37.

**Table S11.** Range of $S_{40\text{–}60}$ over all structures.

| Device, wavelength | min | max |
|---|---|---|
| Al, 500 nm | 0.320 | 0.364 |
| Al, 550 nm | 0.327 | 0.360 |
| Al, 600 nm | 0.339 | 0.361 |
| Ag, 550 nm | 0.322 | 0.361 |
| low loss, 550 nm | 0.274 | 0.366 |

**Fig. S3.** Angular composition of the output for all devices and wavelengths. $S_{40\text{–}60}$ of the multi-pass output against $\eta_{\mathrm{sta}}$ for all 38 structures; rows, Al, Ag, and low-loss devices; columns, 500, 550, and 600 nm; dotted line, Lambertian partition (0.337). Monte Carlo uncertainty of $S_{40\text{–}60}$ is at most 0.003 (ray seeds and realizations, six random and wavy structures). *(`figS2_composition_window.png`)*

**Fig. S4.** Radiance bound and angular composition of the output (main text, Sec. 3.2). (a) First-pass power into each polar band (Al, 550 nm) for the uniform-radiance bound, the radiance-ordered bound, and four structures. (b) $S_{40\text{–}60}$ of the multi-pass output versus $\eta_{\mathrm{sta}}$ for all structures. The shaded band marks 0.33–0.36, and the dotted line marks the Lambertian value. *(`fig5_bound_window.png`)*

---

## Note S5. Benchmark platform and design space

All campaigns share the platform of Table S12; the hemisphere is a point of the freeform design space. Bands are polar bins of the far field over the full azimuth, with the Lambertian partition 0.117 / 0.296 / 0.337 / 0.220.

**Table S12.** Common platform and design variables.

| Quantity | Freeform | Hemispherical reference |
|---|---|---|
| Source | disk, radius 1 mm, CPS dipole distribution, $\eta_{\mathrm{rad}}$ = 0.98, emission 453–753 nm | same |
| Substrate, textured patch | glass, 1.295 mm, $n$ = 1.51; 25 × 25 mm | same |
| Lenslets | hexagonal, ~10 µm radius, N-BK7 index-matched | same |
| Profile | 7-point spline, end points (0, 1) and (1, 0); $x_2$–$x_6$ in 0–1, ascending; $y_2$–$y_6$ in 0–1.5 | quarter circle |
| Height scale (stretch$_Z$) | 0.1–3 | 0.1–3 |
| ETL and HTL thickness | 10–150 nm each | 10–150 nm each |
| Feasibility | no self-intersection or collinear triples; point spacing 0.05–1.0; turning angle ≤ 120° | always satisfied |
| Free variables | 13 | 3 |

**Searches and fidelity.** Each arm is optimized by MATLAB `surrogateopt` (radial-basis surrogate) followed by `patternsearch` at 10,000 rays and 31 wavelengths; every optimum is re-evaluated three times at 50,000 rays and 151 wavelengths. EQE is the radiative efficiency times the spectrally weighted fraction of generated photons reaching air, over the hemisphere or a band; it includes the cavity's coupling into the substrate and is not directly comparable with $\eta_{\mathrm{sta}}$. Infeasible candidates are returned to the optimizer as constraint violations and not traced. The 25 mm patch underestimates the total of an infinite film by about 4%, but no band share changes by more than 0.24 pp between 15 and 100 mm, and all comparisons are made at equal patch.

**Fig. S5.** Platform and benchmark design. (a) OLED–substrate–MLA architecture. (b) The hemispherical reference and an axisymmetric freeform profile (40–60° optimum; open circles, free control points) in units of the lens radius. (c) The four polar bands. (d) Optimization workflow. *(`fig1_platform.png`)*

---

## Note S6. Restart control and statistics

In the original per-band campaign (random seeds plus the best designs of a separate weighted-sum campaign) the 0–20° and 20–40° gains were below unity (0.946, 0.966). Because the hemisphere lies in the freeform design space, this was a search failure. The restart control repeats every search from that arm's hemispherical optimum: the hemisphere and eight perturbations within 8% of each range seed `surrogateopt`, and `patternsearch` runs both from the surrogate winner and directly from the hemisphere. The hemisphere is re-measured in the same session.

**Table S13.** Restart control. Mean ± standard deviation over three re-evaluations.

| Arm | hemisphere (archived) | hemisphere (re-measured) | dev. | restarted | gain | $t$ | winning branch |
|---|---|---|---|---|---|---|---|
| 0–20° | 0.06682 | 0.06677 ± 0.00017 | 0.07% | 0.06699 ± 0.00013 | +0.32% | 1.7 | polish from hemisphere |
| 20–40° | 0.16591 | 0.16595 ± 0.00017 | 0.02% | 0.16628 ± 0.00020 | +0.20% ‡ | 2.1 ‡ | surrogateopt |
| 40–60° | 0.18503 | 0.18486 ± 0.00008 | 0.09% | 0.19355 ± 0.00010 | **+4.70%** | **114.6** | polish from surrogate |
| 60–80° | 0.13589 | 0.13581 ± 0.00014 | 0.06% | 0.14075 ± 0.00025 | **+3.64%** | **29.6** | surrogateopt |
| total EQE | 0.54679 | 0.54679 ± 0.00013 | 0.00% | 0.55231 ± 0.00004 | **+1.01%** | **68.7** | polish from surrogate |

‡ Settled by re-measuring both stored designs five times each: 0.16583 ± 0.00015 against 0.16631 ± 0.00021, +0.29% at $t$ = 4.1 (threshold 1.86). In Table S31, the 40–60° entry is the per-band optimum (0.19802), which the restart did not surpass, and the total is the best of three dedicated starts (0.5539).

**Independent searches.** Freeform/hemisphere ratios: 0–20° 0.946 (original campaign, 60 + 15 evaluations) and 1.003 (restart, about 47); 20–40° 0.966 and 1.002; 40–60° 1.070 and 1.046; 60–80° 1.020 and 1.036; total EQE 1.005, 1.011, and 1.013 (three random starts, 140 evaluations each) and 1.010 (restart). The hemisphere used about 45 evaluations for 3 variables (15 per variable), the freeform searches 47–140 for 13 (4–11 per variable), so the comparison favors the hemisphere. Inverted (concave) and randomly assembled MLA families optimized on the same platform reached total EQEs of 0.517 and 0.522, below the hemisphere.

**What the $t$ values test.** $t = \Delta/(s_p\sqrt{2/n})$ with $n$ = 3 and pooled standard deviation $s_p$ compares two fixed designs against the Monte Carlo noise of one evaluation (about $10^{-4}$ in EQE); with a Holm correction over five objectives [S. Holm, Scand. J. Stat. 6, 65 (1979)] the thresholds are 3.75, 3.50, 3.19, 2.78, and 2.13. The large values ($t$ = 30–115) reflect this small noise. They show that a reported difference is not an evaluation artifact but say nothing about the uncertainty of the optimizer, which only the spread over independent searches above measures. The restart winner is also the best of about seventy noisy evaluations, so its +0.3% residues at 0–40° are upper limits.

**Fig. S6.** Best-so-far search curves. Best objective found so far against the evaluation count, at search fidelity (10,000 rays, 31 wavelengths), for the hemispherical campaign (3 variables), the original freeform per-band campaign, and the restart control, for (a)–(d) the four angular bands and (e) the total EQE. Dotted lines mark the final re-evaluated values (50,000 rays, 151 wavelengths), and dots in (e) mark the final values of three independent random starts (140 evaluations each). Search-fidelity maxima are biased upward by noise and can exceed the re-evaluated values. *(`figS5_optimizer_convergence.png`)*

**Fig. S7.** Restart control. (a) Re-measured hemisphere and restarted optimum for every arm, with standard deviations and $t$. (b) Five-repeat confirmation of the 20–40° arm. *(`figS1_restart_control.png`)*

---

## Note S7. Total efficiency fixes the angular composition

Over the 606 usable evaluations of a weighted-sum campaign (114 random feasible designs, 492 optimizer-visited), each band EQE is nearly proportional to the total (Table S14, Fig. S8); random designs alone reproduce the collapse. $S_{40\text{–}60}$ spans only 0.336–0.375 (5–95%), and the 40–60° band carries the largest share in every design with total EQE ≥ 0.40. Band-dedicated optima raise their band's share by 1–27% but lose 1–11% of total EQE, for net band-power gains of 1.00–1.22 against the natural composition of the twenty most efficient designs.

**Table S14.** Band EQE against total EQE (least squares, 95% intervals).

| Band | slope | intercept | $R^2$ (all 606) | $R^2$ (114 random only) |
|---|---|---|---|---|
| 0–20° | 0.125 ± 0.004 | −0.014 ± 0.002 | 0.848 | 0.828 |
| 20–40° | 0.318 ± 0.007 | −0.020 ± 0.004 | 0.931 | 0.925 |
| 40–60° | 0.348 ± 0.005 | 0.003 ± 0.003 | 0.968 | 0.967 |
| 60–80° | 0.189 ± 0.008 | 0.026 ± 0.004 | 0.769 | 0.724 |

**Fig. S8.** Total efficiency fixes the angular composition. (a) Total EQE against 40–60° EQE for random (gray) and optimizer-visited (green) designs; line, $R^2$ = 0.968. (b) Selectivity of the band-dedicated optima (red) against the natural composition (gray, 10–90% range) and the Lambertian partition (dashed), labeled with the net band-power gain. (c) Selectivity against total EQE for all bands. *(`fig3_efficiency_composition.png`)*

---

## Note S8. Lens versus cavity

**Table S15.** Benchmark optima on the $(p, N_{\mathrm{mix}})$ plane (Note S12). Values are for 550 nm, and $\eta_{\mathrm{sta}}$ is for the Al device in the multi-pass model. The 0–20° freeform optimum has a vertical sidewall and is omitted.

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

**Full device model.** The benchmark optimizes the cavity together with the lens. To separate the two, the eight optimized lenses of Table S15 were evaluated in the full device model of the main text (birefringent CPS + TMM, dispersive constants), weighted with the emission spectrum (400–700 nm, 10 nm steps), with $\eta_{\mathrm{rad}}$ = 0.98 and the polar bands of the benchmark (Table S16). "Lens only" places both lenses at the hemisphere's cavity; "cavity only" places the hemisphere at the freeform's cavity.

**Table S16.** Freeform/hemisphere EQE ratio: LightTools and the multi-pass model.

| Objective | cavity ETL/HTL (nm), hemisphere / freeform | LightTools | model, own cavities | model, lens only | model, cavity only |
|---|---|---|---|---|---|
| 20–40° | 54.9/52.1 / 61.2/75.8 | 0.966 | 0.959 | 0.946 | 0.974 |
| 40–60° | 74.6/81.7 / 68.5/114.4 | 1.070 | 1.026 | 1.004 | 0.986 |
| 60–80° | 99.8/71.1 / 90.5/76.3 | 1.020 | 0.972 | 0.964 | 0.996 |
| total EQE | 61.3/97.5 / 67.5/120.6 | 1.013 | 0.974 | 0.971 | 1.003 |

The freeform entries are the designs of the original per-band campaign and of the total-EQE campaign (for 60–80° the restart optimum, ratio 1.036, has no traced return matrix). For the hemispheres the model reproduces the LightTools EQE within +0.7% (total) to +3.8% (60–80°); for the freeform optima it lies between −3.2% (total) and +0.5% (20–40°). The two tools therefore agree on the size of the differences (a few percent) but not on their sign below about 5%; the main reason is the profile: the freeform curves are rebuilt with a PCHIP spline, whereas LightTools uses its own spline, and the total-EQE optimum has a non-monotonic profile with an annular trough (control-point heights 1.05, 0.15, 0.19, 0.53, 0.50), for which the reconstruction is least certain. Lens and cavity interact (Table S17): the freeform 40–60° lens gains 0.4%, 1.8%, and 4.1% over the hemispherical one at the hemisphere's band cavity, the total-EQE cavity, and its own cavity when spectrally weighted, and 5.6–9.2% at 550 nm alone, because the cavity resonance moves with wavelength and the spectral average smears the angular distribution that the lens acts on.

**Table S17.** 40–60° band EQE, freeform/hemisphere (40–60° designs of the benchmark), at three cavities.

| Cavity ETL/HTL (nm) | 550 nm | spectrally weighted |
|---|---|---|
| 74.6/81.7 (hemisphere, 40–60° optimum) | 1.056 | 1.004 |
| 61.3/97.5 (hemisphere, total-EQE optimum) | 1.070 | 1.018 |
| 68.5/114.4 (freeform, 40–60° optimum) | 1.092 | 1.041 | This also explains Table S15, where the total-EQE freeform lens, at the hemisphere's cavity and 550 nm, has $\eta_{\mathrm{sta}}$ = 0.771 against 0.792.

**Lens-only optimization at a fixed stack.** With the stack fixed (Al, ETL 61.3 / HTL 97.5 nm, 550 nm), the hemisphere's height and an 11-variable freeform profile (5 radial and 5 height control points, height scale; same bounds and feasibility rules as the benchmark) were optimized for the same objective in the multi-pass model with the same tracer: a grid search for the hemisphere and a (1+8) evolution strategy with step-size adaptation, started at the hemispherical optimum, for the freeform (80 evaluations at 4,000 rays per bin; winners re-evaluated at 20,000). Results (Table S18): the freeform gains 1.25% in $\eta_{\mathrm{sta}}$ and 9.0% in the 40–60° band. The total-objective winner keeps $p$ (0.422 against 0.421) but lowers $N_{\mathrm{mix}}$ from 7.1 to 2.1 and raises $p_{\mathrm{mix}}/p$ from 0.91 to 0.93, i.e. it recovers 1.0 of the 1.4 pp that complete reciprocal mixing would add to the hemisphere. The LightTools freeform optima, evaluated at the same stack, reach 0.771 ($\eta_{\mathrm{sta}}$) and 0.284 (40–60°).

**Table S18.** Lens-only optimization at a fixed stack (multi-pass model, Al, 550 nm).

| Objective | lens | $\eta_{\mathrm{sta}}$ | 40–60° power | $p$ | $N_{\mathrm{mix}}$ | $\xi$ |
|---|---|---|---|---|---|---|
| $\eta_{\mathrm{sta}}$ | hemisphere (height 0.93) | 0.7936 | 0.2609 | 0.421 | 7.1 | 0.47 |
| $\eta_{\mathrm{sta}}$ | freeform, optimized at fixed stack | 0.8035 | 0.2728 | 0.422 | 2.1 | 0.85 |
| 40–60° power | hemisphere (height 1.73) | 0.7773 | 0.2692 | 0.427 | 8.8 | 0.38 |
| 40–60° power | freeform, optimized at fixed stack | 0.7946 | 0.2935 | 0.428 | 2.7 | 0.79 |

**LightTools control.** `opt_fixedcavity_freeform.m` repeats the per-band freeform optimization of `opt_4band_freeform.m` with ETL and HTL fixed at each arm's hemispherical optimum and the hemisphere as the first seed, so that $G_j \ge 1$ holds by construction and any gain is due to the lens shape; it requires the LightTools pipeline and has not been run here.

## Note S9. Direct ray tracing and out-of-sample test

**Direct ray tracing.** Rays are launched with the first-pass distribution $\mathbf{v}^{(0)}$ (Al, 550 nm; polar angle drawn uniformly in $\sin^2\theta$ within each 1° bin, azimuth and position uniform), traced through the height field, and, when returned, absorbed with probability $1-R(\theta)$ or reflected specularly, keeping their azimuth, and sent to a new uniformly drawn position (the substrate is thick compared with the lattice). $10^6$ rays per case; the standard error is 0.04–0.05 pp. A second run re-randomizes the azimuth at every pass, which is what the chain on polar angle assumes.

**Out-of-sample structures.** Nine structures were chosen after the reduced model and its reciprocal redistribution were fixed and were not used to adjust anything: hexagonal arrays of spherical caps (contact angle 60° and 35°), paraboloids $h = 1-\rho^2$, cones of 45° and 60° side angle, truncated Gaussian bumps (peak height 1.5), a 6 × 6 supercell randomly filled with hemispheres and 60° cones (radius spread ±20%, position jitter 15% of the pitch), and quasi-periodic wavy surfaces whose spectrum is concentrated on a ring $|\mathbf{k}| = 2$ (rms slope 0.6 and 1.2). Return matrices: 20,000 rays per bin (10,000 for wavy surfaces); errors over three devices and three wavelengths.

**Table S19.** Direct ray tracing and out-of-sample test (Al, 550 nm; reduced and fully randomized models minus full series).

| Structure | $p$ | $N_{\mathrm{mix}}$ | $\xi$ | $\eta_{\mathrm{sta}}$ full series | direct, $(\theta,\phi)$ | direct, $\phi$ re-randomized | reduced − full (pp) | fully randomized − full (pp) |
|---|---|---|---|---|---|---|---|---|
| hemispherical MLA (main set) | 0.422 | 7.1 | 0.47 | 0.7918 | 0.7878 ± 0.0004 | 0.7913 | -0.27 | +2.24 |
| freeform MLA, 40–60° (main set) | 0.425 | 7.5 | 0.39 | 0.7887 | 0.7882 ± 0.0004 | 0.7886 | -1.11 | +2.69 |
| wavy, rms slope 0.2 (main set) | 0.402 | 2.2 | 0.35 | 0.6084 | 0.6087 ± 0.0005 | 0.6093 | -2.71 | +19.29 |
| wavy, rms slope 1.0 (main set) | 0.422 | 1.1 | 0.87 | 0.7924 | 0.7916 ± 0.0004 | 0.7917 | -0.07 | +2.15 |
| spherical caps, 60° contact | 0.418 | 2.0 | 0.74 | 0.7761 | — | — | -0.24 | +3.51 |
| spherical caps, 35° contact | 0.410 | 1.4 | 0.66 | 0.7273 | 0.7258 ± 0.0004 | 0.7260 | -1.82 | +7.90 |
| paraboloids | 0.425 | 2.7 | 0.71 | 0.8010 | — | — | -0.27 | +1.50 |
| cones, 45° side | 0.423 | ∞ | 0.23 | 0.7838 | — | — | -3.53 | +3.06 |
| cones, 60° side | 0.431 | 68.2 | 0.68 | 0.7866 | 0.7808 ± 0.0004 | 0.7865 | -2.55 | +3.26 |
| Gaussian bumps | 0.424 | 1.6 | 0.70 | 0.7573 | — | — | -0.34 | +5.76 |
| random mix of cones and hemispheres | 0.427 | 2.7 | 0.76 | 0.7951 | 0.7937 ± 0.0004 | 0.7947 | -0.11 | +2.17 |
| ring-spectrum wavy, slope 0.6 | 0.413 | 1.3 | 0.76 | 0.7565 | 0.7551 ± 0.0004 | 0.7567 | -0.72 | +5.18 |
| ring-spectrum wavy, slope 1.2 | 0.424 | 1.2 | 0.86 | 0.7934 | — | — | -0.11 | +2.18 |

Over all 81 out-of-sample cases the reciprocal reduced model has a mean error of 1.50 pp; without the two cone arrays (63 cases) 0.48 pp and at most 3.3 pp (35° caps), against 2.62 and 8.8 pp for the fully randomized model and 1.30 and 5.2 pp for isotropic redistribution. The cone arrays fail (mean 6.0 and 4.2 pp, maximum 20 and 14 pp, low-loss device). The 45° cones retroreflect light within 0–4° (0.5% of the étendue) into itself, giving $\lambda_2 = 1$ like the planar interface; the 60° cones have a ladder of slow modes ($\lambda_2$–$\lambda_4$ = 0.985, 0.91, 0.85) spread over most angles. In both cases $1-\lambda_2$ underestimates how fast the light that escapes is redistributed. The direct traces confirm the full series for the 60° cones within 0.6 pp, so the failure is that of the one-rate reduction, not of the angle-resolved model.

**Mode occupancy (post hoc).** Because the truncated $M$ is self-adjoint in the inner product weighted by its stationary distribution $\pi^{(S)} \propto w s$ (Note S1), the light returned after the first pass, $u_0 \propto B_R\mathbf{v}^{(0)}$, can be expanded in its orthogonal eigenmodes. Table S20 lists the share of the deviation $u_0/\pi^{(S)} - 1$ carried by slow modes ($\lambda > 0.8$) and the relaxation rate weighted by the actual returned light, $\lambda_{\mathrm{eff}} = \sum_{\ell \ge 2} a_\ell^2 \lambda_\ell/\sum_{\ell \ge 2} a_\ell^2$ (defined in `selectivity/mixing.py` before the out-of-sample test, but examined only after it). For the 60° cones the slow modes carry 0.4% of the deviation, $\lambda_{\mathrm{eff}}$ = 0.43 against $\lambda_2$ = 0.985, and the reduced model with $m = 1-\lambda_{\mathrm{eff}}$ is accurate; for the 45° cones the closed retroreflecting set carries 98% of it, and no single rate works. Over the 38 structures $\lambda_{\mathrm{eff}}$ is less accurate than $\lambda_2$ (mean 0.48 against 0.31 pp, maximum 3.7 against 3.0 pp), mainly for hemispherical MLAs, so $\lambda_2$ is kept and a large $\lambda_2$ with a negligible slow-mode share serves as a warning that one rate misstates the mixing. Main-text Fig. 3(c) uses the gap $\lambda_2 - \lambda_{\mathrm{eff}}$ as the measure: all irreducible structures with a gap of at most 0.35 (45 of 47) have a maximum error of 3.3 pp; the 60° cones (gap 0.555) and the reducible 45° cones are the only failures.

**Prospective test.** Because the two conditions were found after the out-of-sample test, twelve further structures were traced afterwards with the same tracer and settings ($2\times10^4$ rays per bin; `multipass_reset/prospective.py`): hexagonal cone arrays with side angles of 30°, 40°, 50°, 55° and 70°, square pyramids (30°, 45°, 54.7°), inverted square pyramids (54.7°), and one-dimensional prism gratings (30°, 45°, 60°). Each structure was classified from $M$ and the Al first-pass distribution at 550 nm as valid ($\lambda_2 \le 0.999$ and gap $\le 0.35$) or outside, and the classification was written to file before the multi-pass series was evaluated (`selectivity/prospective_test.py`). Results (Table S20): the 45° pyramids and 45° prisms are reducible ($\lambda_2$ = 1, $\xi$ = 0.01–0.02), were predicted to fail, and fail (maximum 21 and 23 pp). The gap condition misclassified in both directions: the 55° cones and 54.7° pyramids (gaps 0.37 and 0.52) were flagged but have maximum errors of 2.5 and 2.0 pp, and the eight structures classified as valid have a mean error of 1.4 pp and a maximum of 7.5 pp (40° cones; 5.7 pp for 30° pyramids, 4.9 pp for 30° cones, 4.2 pp for 60° prisms). The fully randomized model has mean errors of 1.8–4.5 pp for these structures. Irreducibility is therefore a necessary condition that predicts the large failures, whereas the gap is not a reliable predictor; faceted structures carry deterministic paths between flat facets that one rate describes less well than the curved structures of the library.

**Table S20.** Slow-mode occupancy of the returned light and reduced-model error (mean / max over three devices and three wavelengths, pp). Prospective rows: predictions recorded before the multi-pass series was evaluated (`selectivity/data/prospective_predictions.json`).

| Structure | $\lambda_2$ | $\lambda_{\mathrm{eff}}$ | slow-mode share | error with $\lambda_2$ | error with $\lambda_{\mathrm{eff}}$ |
|---|---|---|---|---|---|
| hemispherical MLA | 0.859 | 0.734 | 0.243 | 0.24 / 0.62 | 0.67 / 2.93 |
| freeform MLA, 40–60° | 0.867 | 0.723 | 0.572 | 0.87 / 2.09 | 0.83 / 2.47 |
| wavy, rms slope 0.2 | 0.540 | 0.195 | 0.000 | 1.98 / 2.99 | 0.72 / 2.33 |
| random lenslet shapes, $s$ = 1 | 0.487 | 0.376 | 0.000 | 0.12 / 0.23 | 0.19 / 0.82 |
| spherical caps, 35° (out-of-sample) | 0.296 | 0.237 | 0.000 | 1.42 / 3.30 | 1.18 / 2.83 |
| paraboloids (out-of-sample) | 0.636 | 0.379 | 0.000 | 0.33 / 0.83 | 0.37 / 1.51 |
| cone–hemisphere mix (out-of-sample) | 0.634 | 0.397 | 0.000 | 0.27 / 0.82 | 0.38 / 1.27 |
| cones, 60° (out-of-sample) | 0.985 | 0.430 | 0.004 | 4.15 / 14.36 | 0.23 / 0.77 |
| cones, 45° (out-of-sample) | 1.000 | 0.988 | 0.978 | 5.99 / 20.15 | 5.50 / 18.61 |
| cones, 30° side (prospective, predicted valid) | 0.511 | 0.438 | 0.000 | 1.94 / 4.89 | — |
| cones, 40° side (prospective, predicted valid) | 0.761 | 0.539 | 0.000 | 2.34 / 7.49 | — |
| cones, 50° side (prospective, predicted valid) | 0.608 | 0.421 | 0.000 | 0.72 / 2.66 | — |
| cones, 55° side (prospective, predicted outside) | 0.641 | 0.269 | 0.000 | 0.77 / 2.46 | — |
| cones, 70° side (prospective, predicted valid) | 0.689 | 0.476 | 0.000 | 0.30 / 0.56 | — |
| inverted pyramids, 54.7° (prospective, predicted valid) | 0.578 | 0.243 | 0.000 | 0.89 / 2.30 | — |
| prism grating, 30° (prospective, predicted valid) | 0.341 | 0.268 | 0.000 | 1.11 / 2.90 | — |
| prism grating, 45° (prospective, predicted outside) | 1.000 | 0.992 | 0.981 | 8.19 / 23.40 | — |
| prism grating, 60° (prospective, predicted valid) | 0.792 | 0.461 | 0.000 | 1.73 / 4.19 | — |
| square pyramids, 30° (prospective, predicted valid) | 0.519 | 0.319 | 0.000 | 2.30 / 5.68 | — |
| square pyramids, 45° (prospective, predicted outside) | 1.000 | 0.986 | 0.988 | 6.61 / 21.09 | — |
| square pyramids, 54.7° (prospective, predicted outside) | 0.546 | 0.024 | 0.000 | 0.70 / 1.96 | — |

## Note S10. Estimators and numerical diagnostics

**Estimators.** The traced matrices are built as follows (height-field tracer; the scattering layers use the slab Monte Carlo of `selectivity/slab.py` plus the adding method with the flat Fresnel face).

```
for each incidence bin i = 0..89 (1 deg):
    for each of N rays (N = 20,000; 10,000 for wavy surfaces):
        theta = bin center i + 0.5 deg; phi ~ U(0, 2 pi); (x, y) ~ U(unit cell or supercell)
        start just below the base plane, moving up into the structure
        repeat up to 60-80 interface events:
            march to the next crossing of z = h(x, y) (Lipschitz-safe step, bisection)
            reflect with probability (R_s + R_p)/2 at the local normal, else refract
            if the ray is in air moving up above max h:   escape; record air polar angle
            if the ray is in glass moving down below z=0: return; record substrate polar angle
        rays that reach the event limit are counted as lost (not renormalized)
    BT[i, j_air] = escaped count / N;  BR[i, j_sub] = returned count / N
p_i = sum_j BT[i, j];   M[:, i] = BR[i, :]^T / sum_j BR[i, j]   (full domain)
spectral metrics: principal sub-block of BR^T below 85 deg; M_S[:, i] = B_S[:, i] / s_i, s_i = sum over retained bins
```

An encounter starts when a ray leaves the substrate upwards and ends when it escapes into air or re-enters the substrate below the base plane; Fresnel reflection inside the escape cone is part of the return. The lost fraction is at most 0.3% in any bin and is not renormalized, so $p_i + \sum_j B_R(i\to j) = 1 -$ lost. $\lambda_2$ is the eigenvalue of second-largest modulus of $M_S$ from `numpy.linalg.eigvals` (LAPACK geev, double precision) on the non-symmetrized matrix; $p_\infty = 1-\rho(B_S)$ uses the un-renormalized sub-block, so return beyond the cutoff counts as leaving the chain, consistent with $\tilde p$. $\eta_{\mathrm{sta}}$ uses `numpy.linalg.solve` (LAPACK gesv) on 0–90°, with $\rho(\mathrm{diag}\,\mathbf{R}\,B_R) \le 0.9991$ verified for every case; $\eta_{\mathrm{mix}}$ uses the closed form on the same grid, and the first pass is $\mathbf{p}^{\mathsf T}\mathbf{v}^{(0)}$ for every structure, since the slab matrices already contain the multiple scattering inside the layer and at its face. Spectral weighting: $\bar\eta = \sum_\lambda s(\lambda)\,\eta_{\mathrm{sub}}(\lambda)\,\eta(\lambda)/\sum_\lambda s(\lambda)\,\eta_{\mathrm{sub}}(\lambda)$, 400–700 nm in 10 nm steps, with the device recomputed at each wavelength and the structure matrices fixed.

**Detailed balance.** With $F_{ji} = w_i B_R(i\to j)$ below 85°, the relative residual $\lVert F-F^{\mathsf T}\rVert/\lVert F\rVert$ is 0.013–0.154 (median 0.047) over the 37 non-planar structures; entry-wise, $|F_{ji}-F_{ij}|/(F_{ji}+F_{ij})$ has a median of 0.06 and a 90th percentile of 0.25, dominated by weak entries. Projecting onto the nearest reciprocal matrix with the same column sums changes $N_{\mathrm{mix}}$ by at most 2.7%, $\xi$ by 0.005, and $\eta_{\mathrm{sta}}$ by 0.04 pp (Note S1), so the physical conclusions do not depend on Monte Carlo asymmetry.

**Adding method versus direct slab Monte Carlo.** A single trace through the HG slab with the Fresnel face (Lambertian within each 1° input bin) reproduces the adding-method matrices (Table S21; adding / direct). $S_{40\text{–}60}$ differs by up to 0.006 because the main calculation maps substrate bins to air bins after the slab, whereas the direct trace bins the air angle exactly.

**Table S21.** Scattering layers: adding method / direct Monte Carlo (Al, 550 nm). Last column: largest difference in $p(\theta)$ below 85°.

| Slab | $p$ | $N_{\mathrm{mix}}$ | $\xi$ | $\eta_{\mathrm{sta}}$ | $S_{40\text{–}60}$ | max $|\Delta p(\theta)|$ |
|---|---|---|---|---|---|---|
| $g$ 0.8, $\tau^*$ 1 | 0.3258 / 0.3257 | 1.172 / 1.170 | 0.977 / 0.977 | 0.7374 / 0.7373 | 0.354 / 0.348 | 0.011 |
| $g$ 0, $\tau^*$ 0.25 | 0.3812 / 0.3814 | 2.246 / 2.255 | 0.545 / 0.546 | 0.6754 / 0.6759 | 0.349 / 0.343 | 0.009 |

**Bands in the hemispherical return matrix.** The share of returned light that lands more than 10° from its incident angle (incidence 45–80°) is 0.130 for the single-lenslet unit cell, 0.127 for a 6 × 6 supercell of identical hemispheres traced with the other tracer, and 0.127 with twice the rays; the supercell matrices at 20,000 and 40,000 rays differ by 2.4% (relative Frobenius norm). The bands are deterministic pathways of rays totally reflected at the lens flank, e.g. from 70° to about 34° and from 50° to about 39°.

**Signed errors of the reduced model.** The reciprocal reduced model underestimates on average (Table S22); the largest error is −2.99 pp for the wavy surface of rms slope 0.2 on Ag at 550 nm, the structure with the lowest $\xi$ among fast mixers.

**Table S22.** Reciprocal reduced model minus full series (pp).

| Subset | mean signed | mean absolute | max absolute |
|---|---|---|---|
| all 342 cases | -0.24 | 0.30 | 2.99 |
| scattering layers | -0.28 | 0.29 | 2.49 |
| hemispherical MLAs | -0.24 | 0.26 | 1.12 |
| freeform MLAs | -0.40 | 0.60 | 2.09 |
| random MLAs | -0.08 | 0.16 | 0.51 |
| wavy surfaces | -0.67 | 0.67 | 2.99 |
| Al | -0.41 | 0.41 | 2.71 |
| Ag | -0.28 | 0.33 | 2.99 |
| low loss | -0.04 | 0.16 | 1.78 |
| 500 nm | -0.32 | 0.38 | 2.49 |
| 550 nm | -0.29 | 0.29 | 2.99 |
| 600 nm | -0.11 | 0.23 | 2.58 |

**Cutoff dependence of $\xi$.** Table S23 gives $\xi$ for representative structures. Hemispherical and random lenses, whose second eigenvector reaches into near-grazing states, change most; thick scattering layers barely change. The shallow wavy surface and the thin isotropic scatterer stay below random lenses of similar $N_{\mathrm{mix}}$ at every cutoff.

**Table S23.** $\xi$ versus cutoff.

| Structure | 80° | 85° | 88° | 90° |
|---|---|---|---|---|
| planar | 0.000 | 0.000 | 0.000 | 0.000 |
| hemispherical MLA | 0.639 | 0.465 | 0.381 | 0.306 |
| random lenslet shapes, $s$ = 1 | 0.828 | 0.769 | 0.730 | 0.720 |
| wavy, rms slope 0.2 | 0.426 | 0.345 | 0.313 | 0.307 |
| wavy, rms slope 1.0 | 0.890 | 0.871 | 0.868 | 0.868 |
| scattering layer, $g$ 0, $\tau^*$ 0.25 | 0.563 | 0.545 | 0.538 | 0.537 |
| scattering layer, $g$ 0.8, $\tau^*$ 4 | 0.996 | 0.993 | 0.990 | 0.989 |

**Exceeding the complete-mixing reference.** For the chain started in its stationary distribution $q$, the one-step autocorrelation of the escape probability, $C = \sum_i q_i f_i \sum_j M_{ji} f_j / \sum_i q_i f_i^2$ with $f = p - p_{\mathrm{mix}}$, is zero for complete mixing and positive for angle retention. Table S24 lists the three structures with the lowest $C$ and three typical ones; $C$ and $\eta_{\mathrm{sta}}-\eta_{\mathrm{mix}}$ (Al, 550 nm) are correlated with coefficient -0.76 over the 37 non-planar structures.

**Table S24.** Escape autocorrelation and the gap to the reference (Al, 550 nm).

| Structure | $C$ | $\eta_{\mathrm{sta}}-\eta_{\mathrm{mix}}$ (pp) |
|---|---|---|
| W_slope0.4 | -0.094 | +0.42 |
| W_slope0.7 | -0.084 | +0.01 |
| scatter g0.0_t1.0 | +0.026 | -0.02 |
| scatter g0.8_t1.0 | +0.145 | -0.03 |
| A_shape1.0 | +0.430 | -0.67 |
| MLA hemi_total | +0.814 | -1.40 |

**LightTools settings of the benchmark.** Disc source (radius 1 mm) with the angular distribution of the CPS substrate radiance (1° polar mesh, uniform in azimuth), monochromatic runs weighted by the CPS power per wavelength; search at 10 nm steps (31 wavelengths, 453–753 nm) and 10,000 rays, final evaluation at 2 nm steps (151 wavelengths) and 50,000 rays, three repeats; device reflectance as a polarization-averaged angle- and wavelength-resolved user coating on the bottom face; textured patch 25 × 25 mm on a 1.295 mm substrate; far-field intensity mesh in 1° polar bins integrated over azimuth, bands summed with $\sin\theta$ weighting. Freeform profiles are passed to LightTools as spline control points, and each constructed geometry is checked by a control-point round trip with tolerance $10^{-4}$ (lens radius units); geometries failing it are rejected. Evaluated through the freeform pipeline with quarter-circle control points, the hemisphere reproduces its own campaign within 0.00–0.09% (Table S13).

---

## Note S11. Reflective and transmissive mixing

**Device-side kernel.** Let the device reflect returned light with a kernel $D(i\to j)$ instead of the specular reflectance $\mathrm{diag}\,\mathbf{R}$. Energy conservation at each incidence angle fixes the column sums, $\sum_j D(i\to j) = R_i$. Reciprocity of the device requires $w_i D(i\to j) = w_j D(j\to i)$. A rank-one kernel $D(i\to j) = a_i b_j$ satisfies the second condition only if $b_j \propto w_j a_j$, and the first condition then gives $a_i \propto R_i$. The result is $D(i\to j) = R_i R_j w_j/\sum_k w_k R_k$ (main text, Sec. 2.3). The multi-pass efficiency is $\eta = \mathbf{p}^{\mathsf T}(I - D\,B_R)^{-1}\mathbf{v}^{(0)}$, where $B_R$ is the return matrix of the structure itself. Light that has been reflected by $D$ is distributed in proportion to $w_j R_j$, independently of its incidence angle. Its escape probability at the next encounter is therefore $p_D = \sum_k w_k R_k p_k/\sum_k w_k R_k$.

**Comparison with reflective mixing at the escape interface.** At the escape interface, the same construction with column sums $1-p_i$ gives the kernel $(1-p_i)\,w_r(1-p_r)/\sum_k w_k(1-p_k)$ of Sec. 2.2. The redistribution weight $1-p$ is then perfectly anticorrelated with the escape probability, which produces the term $-\sigma_p^2/(1-p)$. At the device, the weight is $R$, and the correction is $\mathrm{Cov}_w(R,p)/\langle R\rangle_w$. For the three bottom-emitting stacks of this work, $|p_D - p| < 0.008$ for all structures at 500, 550, and 600 nm.

**Transmissive mixing.** Main-text Fig. 2(d) compares the ideal case below with non-absorbing Henyey–Greenstein volumes under the planar face (`selectivity/scatter_volume.py`, data in `selectivity/data/scatter_volume.json`): at $\tau^*$ = 0.5 and $g$ = 0.9 they reach 0.764, 0.923 and 0.959 in the Al, Ag and low-loss devices, against 0.798, 0.937 and 0.962 for the ideal randomizer; isotropic volumes ($g$ = 0) reach 0.725, 0.907 and 0.952 at $\tau^*$ = 0.5–1. An ideal lossless randomizer with Lambertian transmission and no reflection, placed between the OLED and a planar glass/air face, gives $\eta_{\mathrm{sta}}$ = 0.798, 0.937, and 0.962 for the Al, Ag, and low-loss stacks at 550 nm. These values equal the fully randomized model $p/[p+(1-p)A']$ to three digits. Every pass through the randomizer distributes the light in proportion to $w$, so the escape probability of the mixed light is $p$, and neither $\sigma_p^2$ nor $\mathrm{Cov}_w(R,p)$ enters. For the six Henyey–Greenstein layers of the main text, placing the layer directly on the OLED or directly under the planar face gives identical $\eta_{\mathrm{sta}}$ (Table S25), because the index-matched glass between them neither absorbs nor changes the angle of the light. Under the hemispherical MLA, the same layers lower $\eta_{\mathrm{sta}}$ from 0.792 to 0.734 for $\tau^*$ = 1. They reflect part of the light back to the device and replace the high escape probability of the lenses by that of the layer. Values in Table S25 use 0.5° bins; the planar face without a layer gives 0.381, 0.351, and 0.654, compared with 0.384, 0.356, and 0.664 at the 1° resolution of the main text.

**Table S38.** Henyey–Greenstein volumes under a planar exit face (Al device, 550 nm): diffuse reflectance of the bare volume for Lambertian incidence, escape probability $p$ of the volume together with its exit face, and $\eta_{\mathrm{sta}}$. The roll-off of main-text Fig. 2(d) at large $\tau^*$ follows the rise of the backscatter: more light returns to the device on each pass.

| $g$ | $\tau^*$ | Volume reflectance | $p$ (volume + face) | $\eta_{\mathrm{sta}}$ |
|---|---|---|---|---|
| 0 | 0.1 | 0.084 | 0.392 | 0.578 |
| 0 | 0.25 | 0.180 | 0.381 | 0.675 |
| 0 | 0.5 | 0.296 | 0.362 | 0.722 |
| 0 | 1 | 0.446 | 0.324 | 0.725 |
| 0 | 2 | 0.610 | 0.263 | 0.676 |
| 0 | 4 | 0.754 | 0.189 | 0.580 |
| 0.9 | 0.1 | 0.100 | 0.395 | 0.694 |
| 0.9 | 0.25 | 0.194 | 0.387 | 0.755 |
| 0.9 | 0.5 | 0.306 | 0.368 | 0.764 |
| 0.9 | 1 | 0.451 | 0.326 | 0.738 |
| 0.9 | 2 | 0.612 | 0.262 | 0.676 |
| 0.9 | 4 | 0.754 | 0.188 | 0.579 |

**Table S25.** Scattering layers on the OLED or under the exit face (550 nm, $\eta_{\mathrm{sta}}$). "Exit side": layer under the planar glass/air face. "Device side": layer on the OLED, planar face. The hemispherical-MLA column places the layer on the OLED under the MLA.

| Layer ($g$, $\tau^*$) | Al: exit side / device side | Al: hemispherical MLA alone / with layer | Ag: exit side / device side | Low loss: exit side / device side |
|---|---|---|---|---|
| 0.0, 0.25 | 0.675 / 0.675 | 0.792 / 0.791 | 0.874 / 0.874 | 0.944 / 0.944 |
| 0.0, 1.0 | 0.725 / 0.725 | 0.792 / 0.737 | 0.907 / 0.907 | 0.949 / 0.949 |
| 0.8, 0.25 | 0.749 / 0.749 | 0.792 / 0.792 | 0.914 / 0.914 | 0.959 / 0.959 |
| 0.8, 1.0 | 0.737 / 0.737 | 0.792 / 0.734 | 0.913 / 0.913 | 0.951 / 0.951 |
| 0.9, 1.0 | 0.738 / 0.738 | 0.792 / 0.734 | 0.913 / 0.913 | 0.951 / 0.951 |
| 0.8, 4.0 | 0.579 / 0.579 | 0.792 / 0.575 | 0.838 / 0.838 | 0.901 / 0.901 |

**Strong microcavities.** For a top-emitting microcavity with a 20 nm Ag cathode, an optical clear resin (n = 1.59) as the exit medium, and a planar exit face, the modeled device reflectance seen from the resin is 0.63 inside the escape cone and 0.92 outside it. The device absorbs most strongly at the angles at which it emits. The covariance is then negative, and $p_D$ = 0.316 compared with $p$ = 0.355. Mixing by the device is therefore also reduced in strong microcavities, although less than mixing at the escape interface. This example uses the stack model of `selectivity/advsci_model.py` and is illustrative only.

**Table S26.** Gain over the actual $\eta_{\mathrm{sta}}$ for complete reciprocal mixing in reflection at the escape interface and at the device (37 non-planar structures). The last column gives the difference between mixing by the device and the fully randomized model.

| Device | λ (nm) | Mixing at the escape interface (pp) | Mixing by the device (pp) | Mixing by the device − fully randomized (pp) |
|---|---|---|---|---|
| Al | 500 | −0.6 to +6.4 (median 3.5) | +0.3 to +18.8 (median 4.5) | −3.3 to −0.7 |
| Al | 550 | −0.4 to +1.8 (median 0.8) | +0.2 to +17.1 (median 2.0) | −2.1 to +0.3 |
| Al | 600 | −0.2 to +1.2 (median 0.1) | +0.2 to +15.7 (median 1.3) | −1.4 to +1.5 |
| Ag | 500 | −0.4 to +9.4 (median 4.7) | +0.3 to +15.8 (median 5.0) | −2.7 to −0.7 |
| Ag | 550 | −0.3 to +2.0 (median 0.6) | +0.1 to +11.2 (median 1.2) | −1.1 to −0.1 |
| Ag | 600 | −0.1 to +1.2 (median 0.2) | +0.1 to +9.2 (median 0.6) | −0.7 to +0.5 |
| low loss | 500 | −0.1 to +0.8 (median 0.1) | +0.1 to +6.7 (median 0.4) | −0.1 to +0.4 |
| low loss | 550 | 0.0 to +0.4 (median 0.0) | +0.1 to +4.1 (median 0.3) | +0.2 to +1.1 |
| low loss | 600 | −0.1 to +4.3 (median 2.1) | +0.1 to +4.8 (median 2.1) | −0.9 to +1.2 |

**Slab model.** The slab has refractive index $n$, an internal isotropic source ($\mathbf{v}^{(0)} \propto w$), a planar exit face with unpolarized Fresnel transmission integrated within each 0.25° bin, and an angle-independent back reflectance $R$. No mixing means specular reflection at both faces. Mixing at the exit face uses the kernel of Sec. 2.2 with the Fresnel $p(\theta)$. Mixing at the back reflector uses $D$ with constant $R$, which gives the fully randomized result $p/[p+(1-p)(1-R)]$. A lossless transmissive randomizer in the slab gives the same result.

**Table S27.** Slab with a planar exit face and a back reflector.

| $n$ | $R$ | $p$ | $p_{\mathrm{mix}}$ | No mixing | Mixing at exit face | Mixing at back reflector |
|---|---|---|---|---|---|---|
| 1.5 | 0.90 | 0.404 | 0.054 | 0.440 | 0.599 | 0.871 |
| 1.5 | 0.98 | 0.404 | 0.054 | 0.445 | 0.837 | 0.971 |
| 2.0 | 0.90 | 0.210 | 0.040 | 0.245 | 0.419 | 0.726 |
| 2.0 | 0.98 | 0.210 | 0.040 | 0.249 | 0.733 | 0.930 |
| 2.4 | 0.90 | 0.137 | 0.032 | 0.170 | 0.331 | 0.614 |
| 2.4 | 0.98 | 0.137 | 0.032 | 0.174 | 0.665 | 0.888 |
| 3.0 | 0.90 | 0.080 | 0.024 | 0.107 | 0.243 | 0.467 |
| 3.0 | 0.98 | 0.080 | 0.024 | 0.111 | 0.576 | 0.814 |
| 3.5 | 0.90 | 0.055 | 0.019 | 0.079 | 0.192 | 0.369 |
| 3.5 | 0.98 | 0.055 | 0.019 | 0.082 | 0.509 | 0.745 |

**Textured exit faces at higher index.** The same slab was traced with an index-matched hemispherical MLA and with two Gaussian wavy surfaces (rms slope 0.4 and 1.0) as the exit face, at $n$ = 1.51, 2.4, and 3.5 (height-field tracer of Note S2; $10^4$–$2\times10^4$ rays per 1° bin, up to 400 surface events, unaccounted power below $6\times10^{-4}$). A textured face transmits a comparable fraction from every internal angle, so $\sigma_p^2$ is small and $p_{\mathrm{mix}}$ is close to $p$ (Table S28). For the hemispherical MLA and the steeper wavy surface, the slab without mixing comes within 2.9 pp of the fully randomized result at $n$ = 1.51 and within 1.1 pp at $n \ge 2.4$, and complete mixing at the exit face changes it by at most 2.0 pp. The shallow wavy surface ($p_{\mathrm{mix}}/p$ = 0.60–0.88) stays 2.8–5.8 pp below the fully randomized result at $R$ = 0.90; complete mixing at its exit face recovers none of this difference, whereas mixing at the back reflector or in transmission recovers all of it. Mixing at the exit face is therefore ineffective in both cases: when $\sigma_p^2$ is small, there is little to gain, and when $\sigma_p^2$ is large, the reduction cancels the gain.

**Table S28.** Slab with a textured exit face and a back reflector. Isotropic internal source, angle-independent specular back reflectance $R$.

| Exit face | $n$ | $R$ | $p$ | $p_{\mathrm{mix}}$ | No mixing | Mixing at exit face | Mixing at back reflector |
|---|---|---|---|---|---|---|---|
| Hemispherical MLA | 1.51 | 0.90 | 0.421 | 0.384 | 0.850 | 0.870 | 0.879 |
| Hemispherical MLA | 1.51 | 0.98 | 0.421 | 0.384 | 0.962 | 0.971 | 0.973 |
| Hemispherical MLA | 2.4 | 0.90 | 0.152 | 0.150 | 0.632 | 0.640 | 0.643 |
| Hemispherical MLA | 2.4 | 0.98 | 0.152 | 0.150 | 0.894 | 0.899 | 0.900 |
| Hemispherical MLA | 3.5 | 0.90 | 0.064 | 0.064 | 0.405 | 0.406 | 0.407 |
| Hemispherical MLA | 3.5 | 0.98 | 0.064 | 0.064 | 0.772 | 0.774 | 0.774 |
| Wavy, rms slope 0.4 | 1.51 | 0.90 | 0.410 | 0.244 | 0.818 | 0.815 | 0.874 |
| Wavy, rms slope 0.4 | 1.51 | 0.98 | 0.410 | 0.244 | 0.955 | 0.954 | 0.972 |
| Wavy, rms slope 0.4 | 2.4 | 0.90 | 0.141 | 0.108 | 0.564 | 0.563 | 0.622 |
| Wavy, rms slope 0.4 | 2.4 | 0.98 | 0.141 | 0.108 | 0.863 | 0.863 | 0.891 |
| Wavy, rms slope 0.4 | 3.5 | 0.90 | 0.057 | 0.050 | 0.349 | 0.350 | 0.377 |
| Wavy, rms slope 0.4 | 3.5 | 0.98 | 0.057 | 0.050 | 0.726 | 0.727 | 0.751 |
| Wavy, rms slope 1.0 | 1.51 | 0.90 | 0.422 | 0.369 | 0.866 | 0.866 | 0.879 |
| Wavy, rms slope 1.0 | 1.51 | 0.98 | 0.422 | 0.369 | 0.970 | 0.970 | 0.973 |
| Wavy, rms slope 1.0 | 2.4 | 0.90 | 0.150 | 0.145 | 0.631 | 0.631 | 0.639 |
| Wavy, rms slope 1.0 | 2.4 | 0.98 | 0.150 | 0.145 | 0.895 | 0.895 | 0.898 |
| Wavy, rms slope 1.0 | 3.5 | 0.90 | 0.062 | 0.061 | 0.396 | 0.395 | 0.399 |
| Wavy, rms slope 1.0 | 3.5 | 0.98 | 0.062 | 0.061 | 0.766 | 0.765 | 0.768 |

**Spectrally weighted mixing modes.** Table S29 repeats the comparison of Table S26 with the spectral weighting of the main text (400–700 nm in 10 nm steps, emitter spectrum multiplied by $\eta_{\mathrm{sub}}(\lambda)$; device recomputed at each wavelength). For lens arrays in the Al device, mixing by the device adds 1.1–3.3 pp, compared with at most 0.6 pp at the escape interface. In the low-loss device, the three modes give similar gains for lens arrays (0.4–2.9 pp), because there the gain comes mainly from the mismatch between the first-pass distribution and $p(\theta)$ at long wavelengths (Table S10).

**Table S29.** Spectrally weighted gain over the actual $\eta_{\mathrm{sta}}$ for the three mixing modes. Ranges over each structure class; the number of structures is given in parentheses.

| Device | Structures | $\eta_{\mathrm{sta}}$ | Mixing at the escape interface (pp) | Mixing by the device (pp) | Transmissive mixing (pp) |
|---|---|---|---|---|---|
| Al | Hemispherical MLAs (5) | 0.790–0.794 | +0.4 to +0.6 | +1.5 to +2.3 | +0.7 to +1.3 |
| Al | Freeform MLAs (4) | 0.771–0.790 | 0.0 to +0.6 | +2.0 to +3.3 | +1.5 to +3.1 |
| Al | Random MLAs (17) | 0.779–0.796 | +0.2 to +0.4 | +1.1 to +2.6 | +0.5 to +2.6 |
| Al | Wavy surfaces (5) | 0.619–0.796 | −0.2 to +1.1 | +1.4 to +15.8 | +0.8 to +16.8 |
| Al | Scattering layers (6) | 0.562–0.740 | 0.0 to +1.1 | +0.2 to +9.2 | 0.0 to +10.3 |
| Al | Planar interface | 0.419 | +8.4 | +35.0 | +36.5 |
| Ag | Hemispherical MLAs (5) | 0.935–0.938 | +0.6 to +0.9 | +1.0 to +1.4 | +0.8 to +1.2 |
| Ag | Freeform MLAs (4) | 0.926–0.936 | +0.1 to +0.9 | +0.9 to +1.9 | +1.0 to +2.0 |
| Ag | Random MLAs (17) | 0.935–0.941 | +0.2 to +0.6 | +0.6 to +1.1 | +0.6 to +1.3 |
| Ag | Wavy surfaces (5) | 0.838–0.943 | −0.2 to +1.4 | +0.6 to +9.8 | +0.5 to +10.4 |
| Ag | Scattering layers (6) | 0.847–0.921 | 0.0 to +0.8 | +0.1 to +4.7 | +0.1 to +5.4 |
| Ag | Planar interface | 0.396 | +31.3 | +53.6 | +54.5 |
| low loss | Hemispherical MLAs (5) | 0.936–0.943 | +1.8 to +2.2 | +1.9 to +2.5 | +2.2 to +2.9 |
| low loss | Freeform MLAs (4) | 0.940–0.955 | +0.8 to +2.1 | +1.2 to +2.5 | +1.0 to +2.5 |
| low loss | Random MLAs (17) | 0.941–0.960 | +0.4 to +1.9 | +0.8 to +2.1 | +0.6 to +2.4 |
| low loss | Wavy surfaces (5) | 0.923–0.962 | −0.1 to +1.0 | +0.4 to +4.9 | +0.3 to +4.0 |
| low loss | Scattering layers (6) | 0.898–0.958 | 0.0 to +0.2 | +0.1 to +2.2 | −0.1 to +1.3 |
| low loss | Planar interface | 0.634 | +23.0 | +33.7 | +32.8 |

**Partially diffuse electrode.** A device that reflects a fraction $m$ of the returned light with the reciprocal kernel $D$ and the rest specularly, $D_m = (1-m)\,\mathrm{diag}\,\mathbf{R} + m D$, is reciprocal for every $m$. Table S30 gives the spectrally weighted $\eta_{\mathrm{sta}}$ for six structures. The gain increases nearly linearly with $m$ for lens arrays and steep wavy surfaces and is at most 2.3 pp at $m$ = 1 in the Al device. For the planar interface and the shallowest wavy surface, a weakly diffuse electrode ($m$ = 0.1) already adds 1.2–29 pp, because these structures retain trapped light at its angle. The product $D_m B_R$ need not have a real spectrum, so the mixing number of Sec. 2.1 is not defined for $m > 0$.

**Table S30.** Spectrally weighted $\eta_{\mathrm{sta}}$ with a partially diffuse electrode (diffuse fraction $m$).

| Device | Structure | $m$ = 0 | 0.1 | 0.25 | 0.5 | 1 |
|---|---|---|---|---|---|---|
| Al | Planar interface | 0.419 | 0.504 | 0.589 | 0.678 | 0.769 |
| Al | Scattering layer ($g$ = 0.8, $\tau^*$ = 1) | 0.724 | 0.725 | 0.725 | 0.726 | 0.728 |
| Al | Hemispherical MLA | 0.794 | 0.797 | 0.801 | 0.806 | 0.811 |
| Al | Random MLA (profile spread 0.7) | 0.781 | 0.784 | 0.788 | 0.794 | 0.804 |
| Al | Wavy, rms slope 0.2 | 0.619 | 0.647 | 0.680 | 0.722 | 0.778 |
| Al | Wavy, rms slope 1.0 | 0.787 | 0.789 | 0.792 | 0.796 | 0.805 |
| Ag | Planar interface | 0.396 | 0.683 | 0.810 | 0.884 | 0.932 |
| Ag | Scattering layer ($g$ = 0.8, $\tau^*$ = 1) | 0.919 | 0.919 | 0.919 | 0.920 | 0.921 |
| Ag | Hemispherical MLA | 0.937 | 0.941 | 0.943 | 0.946 | 0.949 |
| Ag | Random MLA (profile spread 0.7) | 0.936 | 0.937 | 0.939 | 0.942 | 0.946 |
| Ag | Wavy, rms slope 0.2 | 0.838 | 0.861 | 0.885 | 0.910 | 0.936 |
| Ag | Wavy, rms slope 1.0 | 0.939 | 0.940 | 0.941 | 0.943 | 0.947 |
| low loss | Planar interface | 0.634 | 0.850 | 0.917 | 0.951 | 0.971 |
| low loss | Scattering layer ($g$ = 0.8, $\tau^*$ = 1) | 0.948 | 0.948 | 0.949 | 0.949 | 0.949 |
| low loss | Hemispherical MLA | 0.940 | 0.947 | 0.952 | 0.957 | 0.961 |
| low loss | Random MLA (profile spread 0.7) | 0.960 | 0.961 | 0.963 | 0.965 | 0.968 |
| low loss | Wavy, rms slope 0.2 | 0.923 | 0.935 | 0.947 | 0.960 | 0.972 |
| low loss | Wavy, rms slope 1.0 | 0.962 | 0.963 | 0.964 | 0.965 | 0.967 |

**Light trapping.** Light trapping in a solar cell is the reciprocal problem. To test whether the placement of mixing carries over, we used the face matrices of Table S28 as the front face of a weakly absorbing slab with a lossless back reflector (`multipass_reset/solar_trapping.py`). Light arrives from air within 10° of the normal. The power coupled into each internal angle follows from the transmission matrix by reciprocity. A one-way crossing at internal angle $\theta$ absorbs $1-\exp(-\alpha d/\cos\theta)$. The four cases share the same front-face transmission $p(\theta)$. A: the traced face and a specular back reflector. B: the face mixes the light it reflects completely, with the kernel of Sec. 2.2, and the back reflector is specular. C: the traced face and a Lambertian back reflector, which mixes the light away from the escape face. D: both B and C. Table S39 gives the path enhancement $F$ = absorptance/($\alpha d$) at $\alpha d = 10^{-3}$. $F$ is given per unit incident power, so it includes the front-surface reflection. In all 12 cases, B gives a lower or equal $F$ than C, and D is within 0.1 of C. Near-normal light enters at angles with a high escape probability. On its first return to the front face, much of it leaves before mixing in reflection can act on it. A Lambertian back reflector moves this light into trapped angles before it reaches the front face. This is the same mechanism that lowers $p_{\mathrm{mix}}$ below $p$ in Sec. 2.2. The difference is large for faces that couple light in near the normal (planar, shallow wavy surface: 3.4 vs 9.0 at $n$ = 1.51). It is small for lens arrays and the steep wavy surface, because these faces already send the incoming light into oblique angles. Values of $F$ slightly above $4n^2$ occur because the front faces reflect part of the light even inside the escape cone. The model uses ray optics and a single wavelength. It is a check of the mechanism, not a solar-cell design.

**Table S39.** Weakly absorbing slab with a textured front face and a lossless back reflector, illuminated within 10° of the normal. Coupled: fraction of the incident power that enters the slab. A–D: path enhancement $F$ at $\alpha d = 10^{-3}$ (see text). Last column: absorptance at $\alpha d$ = 0.05 for B / C.

| $n$ | Front face | $p$ | $p_{\mathrm{mix}}$ | Coupled | A | B: mixing in reflection at the front | C: Lambertian back | D: B and C | $4n^2$ | Absorptance, $\alpha d$ = 0.05 (B / C) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.51 | Planar | 0.398 | 0.055 | 0.96 | 2.0 | 5.3 | 9.4 | 9.4 | 9.1 | 0.120 / 0.303 |
| 1.51 | Hemispherical MLA | 0.421 | 0.384 | 0.99 | 5.9 | 8.8 | 8.9 | 9.0 | 9.1 | 0.281 / 0.301 |
| 1.51 | Wavy, rms slope 0.4 | 0.410 | 0.244 | 0.95 | 3.1 | 3.4 | 9.0 | 9.0 | 9.1 | 0.124 / 0.299 |
| 1.51 | Wavy, rms slope 1.0 | 0.422 | 0.369 | 0.97 | 7.3 | 7.4 | 8.8 | 8.8 | 9.1 | 0.239 / 0.300 |
| 2.4 | Planar | 0.137 | 0.033 | 0.83 | 2.0 | 17.8 | 23.4 | 23.5 | 23.0 | 0.188 / 0.477 |
| 2.4 | Hemispherical MLA | 0.152 | 0.150 | 0.91 | 20.5 | 22.7 | 23.0 | 23.0 | 23.0 | 0.494 / 0.508 |
| 2.4 | Wavy, rms slope 0.4 | 0.141 | 0.108 | 0.83 | 17.1 | 17.1 | 22.6 | 22.6 | 23.0 | 0.317 / 0.474 |
| 2.4 | Wavy, rms slope 1.0 | 0.150 | 0.145 | 0.91 | 23.1 | 23.0 | 23.5 | 23.5 | 23.0 | 0.483 / 0.515 |
| 3.5 | Planar | 0.055 | 0.020 | 0.69 | 2.0 | 38.4 | 46.6 | 46.7 | 49.0 | 0.241 / 0.533 |
| 3.5 | Hemispherical MLA | 0.064 | 0.064 | 0.82 | 45.0 | 47.6 | 48.1 | 48.1 | 49.0 | 0.613 / 0.621 |
| 3.5 | Wavy, rms slope 0.4 | 0.057 | 0.050 | 0.72 | 43.4 | 43.5 | 47.3 | 47.3 | 49.0 | 0.469 / 0.556 |
| 3.5 | Wavy, rms slope 1.0 | 0.062 | 0.061 | 0.88 | 53.1 | 52.9 | 53.4 | 53.3 | 49.0 | 0.652 / 0.672 |

**Summary of robustness checks.** Table S40 lists the conditions under which the central result of the main text was tested. In every case, reflective mixing by the escape interface adds little, while mixing in transmission or at an angle-independent reflector recovers the fully randomized model. The exception is a first-pass distribution that does not match $p(\theta)$, as at 500 nm.

**Table S40.** Robustness of the central result. Gains are changes of $\eta_{\mathrm{sta}}$ in pp for complete mixing, for the 37 non-planar structures unless stated otherwise.

| Condition | Result | Source |
|---|---|---|
| Al, Ag, and low-loss devices, 550 nm | Reflective mixing by the structure: −0.4 to +1.8 (Al), at most 2.0 (Ag) and 0.4 (low loss). Mixing at the device: +1.5 to +4.1 for lens arrays (Al) | Tables S26 and S27 |
| Spectral weighting, 400–700 nm | Reflective mixing: at most 1.1, 1.4, and 2.2 (Al, Ag, low loss). Mixing at the device: 1.1–3.3 for lens arrays (Al) | Table S29 |
| First pass not matched to $p(\theta)$, 500 nm | Reflective mixing: 5.7–6.4 (Al) and 8.4–9.4 (Ag) for the hemispherical MLA | Table S10 |
| Exit medium with $n$ = 1.5–3.5 (slab, $R$ = 0.90) | Textured faces within 1.1 of the fully randomized model at $n \ge 2.4$; mixing at the exit face at most 2.0 | Table S28 |
| Finite substrate (25 mm, absorbing edges) | Reflective mixing: −5.2 to +0.8. Transmissive mixing: −1.0 to +22.7 (six structures, three devices) | Fig. S10 |
| Partly diffuse electrode (diffuse fraction 0–1) | Gain nearly linear in the diffuse fraction, at most 2.3 | Table S30 |
| Numerical parameters, bin width, and reciprocity residual | $\eta_{\mathrm{sta}}$ changes by at most 0.34, $N_{\mathrm{mix}}$ by at most 0.24 | Note S1, Tables S5 and S6 |
| Separate s and p channels | At most 0.4 | Table S7 |
| Direct $(\theta,\phi)$ ray tracing | Within 0.6 of the full series | Note S9 |
| Angular cutoff of the spectral quantities, 80–90° | Mean error of the reduced model changes by less than 0.3 | Note S1 |


## Note S12. Freeform benchmark

This note collects the LightTools benchmark of the Al platform (Notes S5, S6, S8). An axisymmetric freeform MLA with 13 variables (spline profile, height, and ETL and HTL thicknesses) is compared with an equally optimized hemispherical MLA with 3 variables on the Al platform, using LightTools ray tracing validated against fabricated devices (main refs. 11 and 12). The hemisphere is contained in the freeform design space (Table S12). Total EQE and the EQE in four polar bands (0–20°, 20–40°, 40–60°, and 60–80°) are optimized separately by surrogate optimization and pattern search in MATLAB. Each optimum is re-evaluated three times at full fidelity, and a restart control repeats every search from the hemispherical optimum (Note S6).

The band-optimized freeform profiles differ visibly from the hemisphere (Fig. S9(a)), but the efficiency gains are small. The ratio $G_j$ of the best freeform EQE to the best hemispherical EQE is 1.003, 1.002, 1.070, and 1.036 for the four bands and 1.013 for total EQE (Fig. S9(b), Table S31). For the 40–60° band and for total EQE, these ratios correspond to absolute EQE gains of 1.3 and 0.71 pp.

Independent searches give similar ratios: 1.046–1.070 for 40–60°, 1.020–1.036 for 60–80°, and 1.005–1.013 for total EQE. At 0–40°, no search improves the EQE by more than 0.3%; the lower original values at 0–40° (0.946, 0.966) are search failures, since the hemisphere lies in the freeform design space. We use the spread over independent searches as the uncertainty (Fig. S9(c), Fig. S6).

**Table S31.** Freeform versus hemisphere (Al platform, LightTools). The last column gives the freeform/hemisphere ratio over independent searches, with the number of searches in parentheses (data in Table S13).

| Objective | Hemisphere | Best freeform | $G_j$ | Independent searches |
|---|---|---|---|---|
| 0–20° | 0.06682 | 0.06699 | 1.003 | 0.946–1.003 (2) |
| 20–40° | 0.16591 | 0.16628 | 1.002 | 0.966–1.002 (2) |
| 40–60° | 0.18503 | 0.19802 | 1.070 | 1.046–1.070 (2) |
| 60–80° | 0.13589 | 0.14075 | 1.036 | 1.020–1.036 (2) |
| Total EQE | 0.54679 | 0.5539 | 1.013 | 1.005–1.013 (4) |

On the $(p, N_{\mathrm{mix}})$ plane, the nine optima span only $p$ = 0.417–0.429, which corresponds to about 0.7 pp of $\eta_{\mathrm{sta}}$, whereas $N_{\mathrm{mix}}$ spans 3.0–8.0 and $\xi$ spans 0.39–0.79 (Table S15). Shape freedom therefore changes a design mainly along the $N_{\mathrm{mix}}$ axis, along which the efficiency contours are flat.

Because the benchmark optimizes the cavity together with the lens, we separated the two effects (Note S8). With a fixed stack, a lens-only optimization increases the 40–60° band EQE by 9.0% but $\eta_{\mathrm{sta}}$ by only 1.25%. At unchanged $p$, the optimized lens relaxes faster ($N_{\mathrm{mix}}$ 7.1 → 2.1) and recovers 1.0 of the 1.4 pp that complete mixing would add. At a fixed cavity, the lens shape thus redistributes power among the bands but raises the total only through faster mixing.

**Fig. S9.** Freeform versus an equally optimized hemisphere. (a) Band-optimized freeform profiles and the hemispherical reference. Profiles are drawn from the optimized control points with a monotone spline. (b) Best EQE of each lens class for each objective, labeled with $G_j$. (c) Freeform/hemisphere ratio for each independent search: the original per-band search (circles), restarts from the hemispherical optimum (squares), and three random starts for total EQE (triangles). Ratios below one indicate search failures, because the hemisphere is contained in the freeform design space. *(`fig6_hemisphere_benchmark.png`)*

## Note S13. Fabricated devices and the finite-geometry model

**Devices and layout.** The devices of main-text Sec. 3.3 share one substrate design (main-text Fig. 5(a)): 25.4 × 25.4 mm glass of 0.7 mm thickness, ITO on the upper half, five Al fingers (width 2.78 mm, pitch 3.99 mm, straight length 9.93 mm, semicircular ends) and an Al rectangle near the upper edge. A pixel-defining layer opens a circle of 1 mm radius centered on the upper end-cap of the middle finger. The stack is ITO (150 nm) / TAPC (50 nm) / TCTA (10 nm) / TCTA:B3PyMPM (1:1):Ir(dmppy-ph)$_2$tmd (4 wt%, 25 nm) / B3PyMPM (45 nm) / LiF (1 nm) / Al (100 nm), with the emission zone at the center of the EML and a horizontal dipole ratio of 0.865. Lens devices carry a second 0.7 mm glass, coupled with index-matching oil, with a hexagonal close-packed array of hemispherical or freeform lenslets (radius 10 µm, index-matched) on a 3 × 3 mm square centered on the pixel. The freeform profile was digitized from its seven-point spline (height 2.28 lens radii at the axis, maximum 2.50 at 0.46 radii, vertical sidewall at the rim). The substrate sides were blackened. The green devices use the same layout and lens arrays with ITO (150 nm) / TAPC (105 nm) / TCTA (10 nm) / TCTA:B3PyMPM:Ir(ppy)$_2$acac (4 wt%, 25 nm) / B3PyMPM (42 nm) / LiF (1 nm) / Al (100 nm), a horizontal dipole ratio of 0.77, and a freeform lenslet of a different profile (cup with a central spike and an undercut, radius 10.2 µm, height 14 µm, taken from the STL model).

**Measurement.** Current–voltage–luminance curves and angle-resolved spectra were recorded with a photodiode goniometer (photodiode at 59 mm from the rotation axis), from +40° to −80° in 10° steps with four readings per angle. The EQE is the on-axis EQE, inferred for a Lambertian emitter, multiplied by the angular factor $2\int_0^{\pi/2}[I(\theta)/I(0)]\sin\theta\,d\theta$ (natural cubic spline of the negative-angle side, $I(90°) = 0$). Table S32 lists the maximum EQE above turn-on and the maximum EQE at 3.4 V and above, which avoids the low-luminance range near turn-on, where the photodiode dark current contributes 3–5% of the signal. The main text uses the latter. Device S4D3 (freeform) is excluded: its intensity maximum lies at −40°, only one angular file exists, and its EQE exceeds the optical ceiling of the simulation.

The green devices (Table S33) were measured in the same way, with two to seven angular scans per device. Several scans were taken with the sample tilted. For each scan, the tilt $\delta$ is the shift that minimizes the difference between the $+\alpha$ and $-\alpha$ intensities over 5–30°, and the corrected profile is $I(\delta-\alpha)/I(\delta)$; for the bare device G-S1D4, the corrected scans with $\delta = -13°$ agree with the scans at $\delta \approx 0°$ within 0.02 up to 50° and within 0.06 at 60°. Angles beyond the measured range after the shift (70–80° for $|\delta| > 10°$) are taken from the model when the angular factor is computed. Scan 7 of the hemispherical device is excluded (outliers at 40° and 60°). All freeform devices show a second intensity maximum at −70° that the model does not contain; we attribute it to stray light from the edges of the lens glass and exclude 70–80° for these devices (with these angles, the angular factors are 0.91–1.00). Device G-S5D3 (freeform) is excluded: its intensity at 40–60° exceeds that of the other freeform devices by up to 0.15. The EQE is the Lambertian-equivalent on-axis EQE of each device times its mean angular factor.

**Model.** Device optics follow main-text Sec. 2.4 (CPS, PLQY 0.98 with the Purcell factor of the stack, 400–700 nm in 10 nm steps). The back-face reflectance seen from the glass is computed for four regions: pixel (full stack), Al-covered with or without ITO, and uncovered (organics/air) with or without ITO; the pixel-defining layer is neglected. Each encounter with the front face inside the lens patch uses the traced 1° return and transmission kernels of the lens array, including the change of azimuth (height-field tracer of Note S2, $2\times10^4$ rays per bin); outside the patch, the planar glass/air face acts specularly with Fresnel transmission. Between encounters, light moves laterally by $2t\tan\theta$ with $t$ = 1.4 mm (0.7 mm for bare devices), and light that reaches the blackened sides is lost. Light leaving the front face is binned in air angle. With mirror sides, the estimator reproduces $\eta_{\mathrm{sta}}$ of the infinite substrate to $2\times10^{-4}$.

**Lattice orientation and goniometer.** Because the hexagonal lattice is not axially symmetric, the intensity in a goniometer scan plane depends on the azimuth of the plane relative to the lattice. The lens kernels are therefore also traced with the incident azimuth resolved relative to the lattice (revolved-profile tracer; $6\times10^4$ rays per incident angle; incident azimuth reduced to 0–30° by the mirror symmetry of the lattice in six bins; outgoing azimuth in 5° bins), and the lateral Monte Carlo tracks the absolute azimuth of each ray ($3\times10^6$ rays per wavelength). The goniometer is modeled as a square photodiode of side $s$ at $D$ = 59 mm from the rotation axis: at goniometer angle $\theta_g$, light leaving the front face at lateral offset $x_s$ (in the scan plane) and $y_\perp$ is collected if its in-plane angle lies within $\theta_g - b/D \pm s/2D$, with $b = x_s\cos\theta_g - (t - z_c)\sin\theta_g$, and its out-of-plane angle within $-y_\perp/D \pm s/2D$; the angular factor is evaluated from the 10° points as in the measurement. Because the kernels are traced at the bin-center incidence angles, light refracted by the flat glass between the lenslets forms a comb with teeth every 1.5° in air, which is removed by a power-weighted three-bin average. Table S34 lists the sensitivity to the lattice orientation $\beta$, scan plane $\phi_s$, rotation axis $z_c$ (device plane, mid-substrate, lens surface) and photodiode size (5 and 10 mm). The lattice orientation dominates: for the hemispherical arrays it changes the angular factor by up to 0.05 and $I(10°)/I(0)$ from 0.955 to 0.983. The rotation axis changes the angular factor by at most 0.03, and the photodiode size by at most 0.02. The ranges in main-text Fig. 5 cover all these geometries.

**Cavity of the green devices.** With the nominal thicknesses, the model gives a bare angular factor of 0.88–0.90 and an rms deviation of 0.08 from the measured bare profile, compared with a measured angular factor of 0.99–1.03. A grid over the B3PyMPM (35–80 nm) and TAPC (80–130 nm) thicknesses fitted to the bare profile alone gives 65–70 nm for B3PyMPM, with the TAPC thickness weakly constrained (80–120 nm; rms deviation ≤ 0.03). Over this set, the predicted angular factors are 0.68–0.72 (hemispherical) and 0.86–0.90 (freeform) with the azimuthally averaged definition. The main text uses 65 and 80 nm. No thickness pair reproduces the lens-device EQEs together with the bare devices: the predicted hemispherical/bare on-axis ratio is at least 1.31 over the whole grid, compared with 1.21 measured. Fig. S10 shows the same model for a generic substrate.

**Table S32.** Measured devices. EQE$_{\max}$: maximum above turn-on (voltage, luminance). EQE$_{\geq 3.4\,\mathrm{V}}$: maximum at 3.4 V and above (main text). Angular factor at EQE$_{\max}$.

| Device | Structure | EQE$_{\max}$ (%) | at V (cd/m²) | EQE$_{\geq 3.4\,\mathrm{V}}$ (%) | Angular factor |
|---|---|---|---|---|---|
| S1D1 | bare | 31.3 | 3.0 (24) | 29.7 | 1.01 |
| S1D2 | bare | 30.6 | 2.2 (1.6) | 27.7 | 1.08 |
| S2D2 | bare | 33.3 | 2.2 (2.2) | 30.7 | 0.95 |
| S4D1 | bare | 34.1 | 2.4 (5.6) | 31.4 | 1.02 |
| S5D2 | bare | 31.8 | 2.6 (70) | 30.8 | 0.97 |
| S1D3 | freeform | 25.9 | 2.3 (1.1) | 22.9 | 1.13 |
| S2D3 | freeform | 31.4 | 2.2 (2.3) | 26.7 | 0.99 |
| S3D3 | freeform | 30.4 | 2.2 (2.6) | 23.7 | 1.11 |
| S5D3 | freeform | 33.7 | 2.2 (3.3) | 28.7 | 1.03 |
| S1D4 | hemispherical | 23.8 | 2.6 (38) | 22.9 | 0.63 |
| S2D4 | hemispherical | 27.7 | 2.6 (84) | 27.1 | 0.69 |
| S3D2 | hemispherical | 24.4 | 2.2 (2.0) | 22.5 | 0.63 |
| S4D2 | hemispherical | 25.9 | 2.4 (9.0) | 25.3 | 0.71 |
| S5D4 | hemispherical | 24.7 | 2.6 (92) | 24.1 | 0.66 |

On the EQE$_{\max}$ basis, the measured means are 32.2% (bare), 30.3% (freeform) and 25.3% (hemispherical); the order of the structures is the same on both bases. Over eleven angular measurements of the five hemispherical devices, the angular factor is 0.62–0.73. The rms deviation of $I/I(0)$ from the model for the 3 × 3 mm patch is 0.026–0.056 per device for the best-matching geometry and 0.035–0.070 for the worst (10–80°); above 60°, the measured intensity of four of the five devices exceeds the upper end of the predicted range by up to 0.09. With the azimuthally averaged far-field model, the deviation is 0.055; a 0.5 mm offset of the patch raises it to 0.067, and patches of 5 and 10 mm give 0.094 and 0.118.

**Table S33.** Green devices (prefix G- distinguishes them from the orange devices of Table S32). Tilt: shift fitted to the ± asymmetry. Angular factor per scan (lab definition; freeform without 70–80°). EQE$_{0°}$: Lambertian-equivalent on-axis EQE. EQE: EQE$_{0°}$ × mean angular factor. Model (bare-calibrated cavity, internal efficiency 0.95): bare 22.5%, hemispherical 22.6%, freeform 20.9%.

| Device | Structure | Scans used | Tilt (°) | Angular factor | EQE$_{0°}$ (%) | EQE (%) |
|---|---|---|---|---|---|---|
| G-S1D4 | bare | 4 | −13.3 to 1.3 | 0.99–1.03 | 22.4 | 22.7 |
| G-S2D2 | bare | 3 | −1.6 to 0.6 | 1.01–1.03 | 21.9 | 22.3 |
| G-S3D2 | hemispherical | 6 (of 7) | −13.4 to −8.2 | 0.69–0.73 | 26.9 | 19.1 |
| G-S1D3 | freeform | 4 | −0.9 to −0.1 | 0.82–0.83 | 21.1 | 17.4 |
| G-S2D3 | freeform | 5 | −1.0 to 1.0 | 0.85–0.90 | 22.4 | 19.7 |
| G-S3D3 | freeform | 3 | −2.0 to −1.2 | 0.89–0.90 | 22.7 | 20.3 |

**Table S34.** Sensitivity of the predicted angular factor to the goniometer geometry (3 × 3 mm patch). Lattice/scan: range over lattice orientation $\beta$ = 0°, 30° and scan plane $\phi_s$ = 0°, 45°, 90° (mid-substrate axis, 5 mm photodiode). Axis: largest change between rotation axes at the device plane, mid-substrate and lens surface. Photodiode: change from 5 to 10 mm. ±40°: ratio of the intensities at +40° and −40°, which the layout makes asymmetric in the scan plane along the fingers.

| Emitter | Structure | All geometries | Lattice/scan | Axis | Photodiode | ±40° |
|---|---|---|---|---|---|---|
| orange | bare | 0.99–1.02 | 1.00–1.01 | 0.01 | −0.00 to +0.01 | 0.99–1.00 |
| orange | hemispherical | 0.62–0.68 | 0.63–0.67 | 0.02 | −0.01 to +0.01 | 0.95–1.00 |
| orange | freeform | 0.92–1.04 | 0.93–1.01 | 0.03 | +0.01 to +0.02 | 0.89–1.00 |
| green | bare | 1.00–1.02 | 1.01 | 0.01 | −0.00 to +0.00 | 0.98–1.02 |
| green | hemispherical | 0.65–0.71 | 0.66–0.71 | 0.02 | −0.01 to +0.01 | 0.94–1.00 |
| green | freeform | 0.83–0.91 | 0.84–0.89 | 0.02 | +0.00 to +0.01 | 0.93–1.00 |

**Loss in the attached lens glass.** The measured EQE of the lens devices is 84–91% of the prediction, and the on-axis ratio is lower by a similar factor, whereas the angular factors agree. Three angle-independent or weakly angle-dependent losses were added to the model (azimuthally averaged, orange devices) and scaled to lower the hemispherical EQE by 14% (`labdev/lens_loss_test.py`): a loss of the light transmitted by the lens surface; a transmission $\tau$ = 0.925 per crossing of the oil-coupled interface between the two glasses; and bulk attenuation $\exp(-\alpha\,0.7\,\mathrm{mm}/\cos\theta)$ per crossing of the lens glass with $\alpha$ = 0.094 mm$^{-1}$. All three change the angular factor by at most 0.02 (hemispherical 0.662 → 0.659–0.669, freeform 0.994 → 0.971–0.979), because the lens devices act mainly on first-pass light; the measurement therefore cannot distinguish them. The bulk coefficient is about two orders of magnitude above that of display glass, so the loss is more likely at the coupling interface or in the lens film. With $\tau$ = 0.925 fixed on the orange hemispherical devices, the model gives 0.89 (orange freeform), 0.87 (green hemispherical) and 0.90 (green freeform) of the lossless EQE, compared with 0.90, 0.84 and 0.91 measured (Table S37). Main-text Fig. 5(e) shows the lens-to-bare ratios for the lossless model and for the model with this interface transmission.

**Table S37.** Candidate causes of the lens-device EQE deficit. Changes relative to the lossless model of the 3 × 3 mm patch (orange devices unless stated).

| Candidate | Lens-device EQE | Angular factor | Consistent with the measurement |
|---|---|---|---|
| cavity of the bare devices | fixed by the bare devices; common to bare and lens devices | — | not a cause of a lens/bare difference |
| lens shape: caps, $a$ = 0.7 and 0.5 | −3% and −7% | +0.04 and +0.13 | no: angular factor above the measured 0.63–0.71 |
| lens shape: gaps, $r_b$ = 0.95–0.7 | 0 to −2% | +0.03 to +0.20 | no: angular factor and a maximum near 10° |
| goniometer geometry (lattice, scan plane, axis, photodiode) | no change of the EQE | ±0.03 | included in the predicted ranges |
| patch offset of 0.5 mm | +1% (on-axis EQE −5%) | +0.05; profile rms 0.067 instead of 0.055 | no: raises the EQE |
| loss at the lens surface, interface or in the lens glass | −14% (fitted), −11% to −13% for the other groups | ≤ 0.02 | yes: one parameter reproduces all four lens-device groups within 0.03 |

**Table S35.** Lens-array variants in the device layout (3 × 3 mm patch, unscaled simulation). Caps: spherical caps with height/base-radius ratio $a$. Gaps: hemispheres of base radius $r_b$ (in units of the cell radius) with flat glass between them. Measured (EQE$_{\geq 3.4\,\mathrm{V}}$ / 0.89): 27.3%, angular factor 0.66 ± 0.04 over five devices (0.63–0.71; 0.62–0.73 over eleven measurements).

| Lenslet | EQE (%) | EQE$_{0°}$ (%) | Angular factor | $I/I(0)$ at 20°, 40°, 60° |
|---|---|---|---|---|
| hemisphere | 31.7 | 47.8 | 0.66 | 0.93, 0.53, 0.23 |
| cap, $a$ = 0.7 | 30.6 | 44.0 | 0.70 | 0.94, 0.57, 0.24 |
| cap, $a$ = 0.5 | 29.4 | 37.4 | 0.79 | 0.97, 0.66, 0.30 |
| gaps, $r_b$ = 0.95 | 31.4 | 45.7 | 0.69 | 1.01, 0.55, 0.22 |
| gaps, $r_b$ = 0.9 | 31.2 | 44.0 | 0.71 | 1.06, 0.57, 0.22 |
| gaps, $r_b$ = 0.85 | 31.3 | 42.1 | 0.74 | 1.13, 0.59, 0.22 |
| gaps, $r_b$ = 0.8 | 31.4 | 40.4 | 0.78 | 1.20, 0.62, 0.21 |
| gaps, $r_b$ = 0.7 | 31.7 | 37.0 | 0.86 | 1.31, 0.67, 0.20 |

For the green devices (B3PyMPM 65 nm, TAPC 80 nm), the hemispherical/bare on-axis ratio is 1.43 for close-packed hemispheres, 1.36, 1.30, 1.19 and 1.09 for base radii of 0.95, 0.9, 0.8 and 0.7, and 1.31 and 1.11 for caps with $a$ = 0.7 and 0.5. The variants that approach the measured 1.21 have angular factors of 0.82–0.83 (azimuthally averaged), above the measured 0.69–0.73, and the gaps also produce a maximum near 10° that the measurement does not show.

**Table S36.** Simulated EQE and angular factor against the width of the lens patch (device layout, unscaled). Bare device: 33.7%, angular factor 1.00.

| Patch width (mm) | Hemispherical EQE (%) | Angular factor | Freeform EQE (%) | Angular factor |
|---|---|---|---|---|
| 2 | 30.9 | 0.87 | 31.0 | 1.40 |
| 3 | 31.7 | 0.66 | 31.8 | 0.99 |
| 4 | 35.6 | 0.70 | 35.0 | 0.91 |
| 5 | 39.3 | 0.76 | 37.7 | 0.93 |
| 7 | 42.8 | 0.81 | 40.4 | 0.96 |
| 10 | 45.4 | 0.85 | 43.1 | 0.98 |
| 14 | 47.1 | 0.88 | 45.0 | 1.01 |
| 20 | 48.5 | 0.90 | 46.6 | 1.03 |

**Fig. S10.** Finite lateral size for a generic substrate (Al device, 550 nm; disk source of 1 mm radius, 1.295 mm glass). (a) $\eta_{\mathrm{sta}}$ against the half-width of the substrate (glass/air sides) without additional mixing (solid) and with complete reciprocal mixing at the escape interface (dashed); dotted lines mark the infinite substrate. (b) $\eta_{\mathrm{sta}}$ against the pixel radius when the stack reflectance is present only under the pixel and the back face outside it is glass/air (solid) or absorbing (dash-dotted). (c) Gain from complete mixing at the escape interface for the infinite substrate, for the benchmark size, and for the benchmark size with an absorbing back face outside the pixel. *(`figS_finite_size.png`)*

---

## Reproducing the figures and data

All figures share the style defined in `figstyle.py`. The numbered reference list of the main text is generated from `tools/cite.py`. Main-text Figs. 1–4 and Fig. S4 are produced by `selectivity/oe_figs.py`; panels (b) and (c) of Fig. 2 and Tables S26–S27 use `selectivity/placement.py`, Fig. 3(c) uses `selectivity/data/cutoff_modes.json` and the prospective test (`multipass_reset/prospective.py`, `selectivity/prospective_test.py`), Fig. 2(d) uses `selectivity/scatter_volume.py`, Table S25 uses `selectivity/device_side_slabs.py`, Table S28 uses `multipass_reset/run_highindex.py` and `selectivity/highindex_slabs.py`, Table S39 uses `multipass_reset/solar_trapping.py`, and Tables S29–S30 use `selectivity/spectral_placement.py`. `oe_figs.py` evaluates the 38 structures with the full and reduced models, takes the return matrices from `multipass_reset/run_random_sweep.py`, and writes PNG, PDF, and an Excel workbook of the plotted arrays to `figures_oe/`. Fig. S9 (former benchmark figure) and Figs. S5, S7, and S8 come from the campaign archives through `make_figures.py` and `make_supp_figures_and_data.py`; main-text Fig. 5 and Tables S32–S36 come from `labdev/` (`device.py`, `layout_stacks.py`, `trace_lenses.py`, `layout.py`, `layout_sweep.py`, `trace_variants.py`, `layout_variants.py`, `device_green.py`, `trace_green_freeform.py`, `green_joint_fit.py`, `green_patch.py`, `trace_abs.py`, `gonio_az.py`, `fig_experiment.py`); Fig. S10 from `selectivity/finite_size.py` and `selectivity/fig_finite_size.py`; Fig. S1 and Table S2 from `selectivity/fig_devices.py`; Tables S5–S11 from `selectivity/unc_analysis.py`, `selectivity/review_checks.py`, `selectivity/review_traces_analysis.py`, and `selectivity/make_note_s12.py`; the complete-mixing reference (Fig. 2, Note S1) from `selectivity/reciprocal_limit.py`; Tables S1 and S20 from `selectivity/cutoff_modes.py`; Tables S16–S17 from `selectivity/cavity_decomp.py` and `selectivity/cav_lens_cross.py`; Table S18 from `multipass_reset/fixedcav_opt.py`; Table S19 from `multipass_reset/direct_blind.py` and `blind_table.py`; Fig. S3 and Tables S22–S24 from `selectivity/review2_diag.py`; Table S21 from `selectivity/slab_direct_check.py`; Fig. S6 from `make_optimizer_convergence.py`; Table S7 from `multipass_reset/polcheck.py`; Table S10 from `selectivity/data/mixing_gain_by_condition.json` and `mixing_gain_spectral.json`.

```
cd multipass_reset && python3 run_random_sweep.py && cd ..
cd selectivity && python3 extract_inputs.py && python3 run_devices.py && python3 oe_figs.py && cd ..
python3 make_figures.py && python3 make_supp_figures_and_data.py
```
