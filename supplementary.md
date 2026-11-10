# Supplement 1

Reciprocity Constrains Angular Mixing in Passive Light-Extraction Structures:
Two Spectral Metrics for OLED Outcoupling

Notes S1–S4 support the spectral metrics: model details, the devices and the structure library, accuracy and convergence, and the separation of mixing from escape. Notes S5–S8 support the freeform-versus-hemisphere case study. Main-text references are cited as "main ref. n".

---

## Note S1 | Model details and numerical checks

**Polarization and azimuth.** The transport model is scalar. $R(\theta)$ is the polarization average $(|r_s|^2+|r_p|^2)/2$ of the stack reflectance, $\mathbf{v}^{(0)}$ is the sum of the s- and p-polarized power launched into the substrate, and the tracers apply the unpolarized Fresnel coefficients at every encounter. Textures that are not axisymmetric (hexagonal arrays, supercells, wavy fields) are reduced to a polar-angle kernel by launching rays at the centre angle of each 1° bin, uniformly over the cell and in azimuth, and binning the exit directions by polar angle only; $B_R$ is the position- and azimuth-averaged return operator.

**Full series.** $\eta_{\mathrm{sta}}$ is obtained by a direct linear solve on the full 0–90° grid. Over all cases the spectral radius of $\mathrm{diag}\,\mathbf{R}\,B_R$ stays below 0.9991, the residual below $10^{-16}$, and a pass-by-pass sum agrees with the solve to $10^{-14}$.

**Angular cutoff.** Bins near 90° carry almost no étendue but can return light at their own angle, and a few such states can dominate $\lambda_2$ without affecting extraction. The spectral metrics are therefore evaluated on $\theta < 85°$, which excludes 0.76% of the étendue and less than 0.1% of $\mathbf{v}^{(0)}$: $B_R$ is restricted to these bins, $M$ is formed by renormalizing the columns of the sub-block, and $\rho$ is taken from the sub-block. $\eta_{\mathrm{sta}}$, in the full series and in the reduced model, always uses 0–90°. For every structure $\lambda_2$ is real and positive (0.09–0.88), so $|\lambda_2| = \lambda_2$. $N_{\mathrm{mix}}$ is the relaxation time, not the e-folding pass number $-1/\ln\lambda_2$ (2 against 1.44 at $\lambda_2$ = 0.5).

**Table S1 | Sensitivity of $|\lambda_2|$ to the angular range.**

| Structure | 0–90° | 0–85° | 0–80° | 0–70° | weight of second eigenvector above 80° |
|---|---|---|---|---|---|
| hemispherical MLA (total-EQE height) | 0.876 | 0.859 | 0.846 | 0.840 | 0.23 |
| freeform MLA (40–60° optimum) | 0.941 | 0.867 | 0.689 | 0.532 | 0.49 |
| random lenslet shapes, s = 0.2 | 0.679 | 0.628 | 0.514 | 0.419 | 0.30 |
| random lenslet heights, spread 0.2 | 0.935 | 0.859 | 0.826 | 0.787 | 0.46 |

For the freeform 40–60° lens and height-randomized hemispheres the full-range value is set by grazing states. With the 85° cutoff the reciprocal reduced model has a mean error of 0.30 pp and a largest error of 3.0 pp over the 342 cases, against 0.55 and 4.7 pp without it. The absolute values of $N_{\mathrm{mix}}$ and $\xi$ depend on the cutoff ($\xi$ of the hemispherical MLA is 0.65, 0.47, and 0.31 at 80°, 85°, and 90°), but the ranking does not: the rank correlation with the 85° values is 0.96–0.99 for $N_{\mathrm{mix}}$ and 0.90–0.98 for $\xi$ at 80°, 88°, and 90° (Table S6).

**Detailed balance of the traced matrices.** The relative asymmetry $\sum|F-F^{\mathsf T}|/\sum|F|$ of $F_{ji} = w_i B_R(i\to j)$ is 0.05–0.13 on 1° bins, 0.009–0.024 on 5° bins, and 0.002–0.007 on 15° bins, falling as expected for Monte Carlo noise; the planar matrix is exactly symmetric. Symmetrizing $F$ while keeping each column sum changes $N_{\mathrm{mix}}$ by at most 0.20 (2.7%), $\xi$ by 0.005, $\eta_{\mathrm{sta}}$ by 0.04 pp, and the reduced-model error by 0.04 pp over all 37 non-planar structures and three devices.

**Reciprocal redistribution.** For $B_R(i\to r) = (1-p_i)[(1-m)\delta_{ir} + m q_r]$, detailed balance of the mixing term requires $w_i(1-p_i)q_r = w_r(1-p_r)q_i$, i.e. $q \propto w(1-p)$. The normalized kernel $(1-m)I + m\,q\mathbf{1}^{\mathsf T}$ has eigenvalues 1 and $1-m$, so $m = 1-\lambda_2$ for any $q$.

**Complete reciprocal mixing.** With $q_r = w_r(1-p_r)/\sum_k w_k(1-p_k)$, $p_{\mathrm{mix}} = \sum_r q_r p_r = (p - \langle p^2\rangle)/(1-p) = p - \sigma_p^2/(1-p)$, where $\langle\cdot\rangle$ is the étendue-weighted mean; $p_{\mathrm{mix}} \le p$ is the Cauchy–Schwarz inequality $p^2 \le \langle p^2\rangle$. The rank-one kernel $B_R = q(\mathbf{1}-\mathbf{p})^{\mathsf T}$ with $\mathrm{diag}\,\mathbf{R}$ gives $(I-\mathbf{u}\mathbf{v}^{\mathsf T})^{-1} = I + \mathbf{u}\mathbf{v}^{\mathsf T}/(1-\mathbf{v}^{\mathsf T}\mathbf{u})$ with $\mathbf{u} = \mathbf{R}\circ q$, $\mathbf{v} = \mathbf{1}-\mathbf{p}$, hence the closed form of main-text Sec. 2.4; it agrees with the numerical series to $10^{-9}$ for all 342 cases. Values for every structure are in `selectivity/data/reciprocal_limit.json`; $p_{\mathrm{mix}}/p$ is 0.87–0.91 for hemispherical MLAs, 0.81–0.87 for freeform MLAs, 0.84–0.93 for random MLAs, 0.37–0.90 for wavy surfaces, 0.59–0.99 for scattering layers, and 0.14 for the planar interface.

**Bounds on $p_\infty$.** With $D = \mathrm{diag}(w)$, $D^{-1/2}B_R D^{1/2}$ is symmetric, so $\rho(B_R)$ is the maximum of its Rayleigh quotient. With $x = D^{1/2}\mathbf{1}$ the quotient is $\sum_{i,r} w_i B_R(i\to r)/\sum_i w_i = 1-\tilde p$, both sums running over the bins below the cutoff; hence $p_\infty \le \tilde p$. $\tilde p$ is the probability that uniform-radiance light leaves the truncated chain in one encounter, by escape or by return beyond the cutoff, and exceeds $p$ by less than 0.008. The column sums are $1-p(\theta_i)$, so $p_\infty \ge \min_\theta p(\theta)$. All 37 non-planar structures satisfy both bounds, with smallest margins of 0.0014 (upper) and 0.106 (lower); the planar interface meets the lower bound with equality. With $p$ in place of $\tilde p$ the upper bound would fail for the thick scattering layer ($\tau^*$ = 4) by 0.005. $\xi$ is a structure-only descriptor; the device-aware variant $1-\rho(\mathrm{diag}\,\mathbf{R}\,B_R)$ ranks the structures like $p_\infty$ (Spearman 0.956, Al, 550 nm). A two-parameter reciprocal kernel fitted to both $\lambda_2$ and $\rho$, and a two-rate scalar formula, were not more accurate overall than the one-parameter model (0.32 pp mean) and are not used.

**Band-resolved bound.** The substrate bins are sorted by radiance $v^{(0)}_k/w_k$; the étendue of an air band, mapped into the substrate by Snell's law, is $(\sin^2\theta_{\mathrm{hi}}-\sin^2\theta_{\mathrm{lo}})/n^2$ in units of $w$, and it is filled with the brightest bins first, the last one fractionally. The bound is lossless, polarization-independent, and integrated over azimuth. Each band's bound is a separate maximum, so the per-band bounds (sum 0.676) do not add up to the bound for the whole hemisphere (0.540), which fills the étendue $1/n^2$ once; the uniform-radiance bound is additive (0.439, of which 0.014 lies above 80°).

---

## Note S2 | Devices and structure library

**Devices.** The emitters form a single dipole sheet at the centre of the 25 nm birefringent EML (measured ordinary and extraordinary indices), with a horizontal dipole ratio of 0.865 (Table S2, Fig. S1).

**Table S2 | Device stacks (from the cathode to the substrate).**

| Device | Cathode | ETL | EML | HTL-side layers | Anode | Substrate | $A'$ at 550 nm |
|---|---|---|---|---|---|---|---|
| Al | Al 100 nm | B3PyMPM 61.3 nm | TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd 25 nm | TCTA 10 nm / TAPC 97.5 nm | ITO 150 nm | glass, $n$ = 1.51 | 0.167 |
| Ag | Ag 100 nm | B3PyMPM 61.3 nm | as Al | TCTA 10 nm / TAPC 97.5 nm | ITO 150 nm | glass, $n$ = 1.51 | 0.044 |
| Low loss | Ag 100 nm | B3PyMPM 200 nm | as Al | TAPC 180 nm | ITO 50 nm | glass, $n$ = 1.51 | 0.026 |

**Fig. S1 | Device inputs at 550 nm.** (a) Radiance launched into the substrate on the first pass, $v^{(0)}_k/w_k$, relative to its mean over 0–20°. (b) Angle-resolved reflectance $R(\theta)$ seen from the substrate, with the uniform-radiance round-trip loss $A'$. Dotted line, critical angle of the glass/air interface. *(`figS11_device_inputs.png`)*

**Height-field tracer.** The surface $z = h(x, y)$ separates the index-matched medium ($n = 1.51$) from air; rays are marched with a Lipschitz-safe step, refined by bisection, and reflected or refracted with the unpolarized Fresnel coefficients. 20,000 rays per 1° incidence bin are used (10,000 for wavy surfaces). For a hemispherical lenslet the tracer reproduces the LightTools return matrix of main ref. 14: escape probability 0.4275 against 0.4274, transmittance within 0.008 at every angle, and fewer than 0.1% of rays undetermined. The freeform profiles are rebuilt from the optimized control points with a PCHIP spline.

**Random lens arrays** tile a 6 × 6 supercell periodically. Shape spread $s$ perturbs the five control points of the quarter circle by Gaussian noise of standard deviation $0.25s$ (radial) and $0.4s$ (height), so that $s$ = 1 approaches independent profiles; height spread scales each lenslet by $1\pm s$; position jitter displaces each centre by up to a fraction of the pitch, optionally with ±20% radius spread. Where lenslets overlap the higher surface is taken; profiles with near-vertical walls are rejected. **Wavy surfaces** are periodic Gaussian-correlated height fields on a 12 × 12 box (correlation length one), scaled to a prescribed rms slope; only the slope statistics matter in ray optics. **Scattering layers** are non-absorbing, index-matched Henyey–Greenstein slabs (20,000 rays per 0.5° bin) combined with the flat exit face by the adding method.

**Table S3 | Structure library at 550 nm.** $p$: escape probability for uniform radiance; $N_{\mathrm{mix}}$ and $p_\infty$ on $\theta < 85°$; $\xi = p_\infty/\tilde p$; $\eta$: full series for the Al, Ag, and low-loss (LL) devices; $S_{40\text{–}60}$: 40–60° share of the multi-pass output (Al); last two columns: reciprocal and isotropic ($q \propto w$) reduced model minus full series (Al). A, lenslet shape spread $s$; B, height spread ±$s$; C, position jitter (fraction of the pitch); C2, as C with ±20% radius spread; W, wavy surface of the given rms slope. Scattering layers are labelled by $g$ and $\tau^*$, MLAs by the benchmark objective of the optimum they reproduce.

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

## Note S3 | Accuracy and convergence

**Table S4 | Error of the reduced and memoryless models by structure class (pp; 342 cases, three devices, 500/550/600 nm).**

| Class | cases | reciprocal mean | reciprocal max | q ∝ w mean | q ∝ w max | memoryless mean | memoryless max |
|---|---|---|---|---|---|---|---|
| Planar | 9 | 0.00 | 0.00 | 0.00 | 0.00 | 44.05 | 59.89 |
| Scattering layers | 54 | 0.29 | 2.49 | 0.96 | 4.67 | 1.83 | 11.14 |
| Hemispherical and freeform MLAs | 81 | 0.41 | 2.09 | 0.38 | 2.22 | 3.45 | 12.33 |
| Random lens arrays (A, B, C) | 153 | 0.16 | 0.51 | 0.44 | 1.93 | 2.41 | 10.89 |
| Wavy, rms slope ≥ 0.7 | 27 | 0.19 | 1.06 | 1.11 | 3.07 | 1.61 | 5.57 |
| Wavy, rms slope ≤ 0.4 | 18 | 1.40 | 2.99 | 5.40 | 10.51 | 8.66 | 21.73 |
| All | 342 | 0.30 | 2.99 | 0.81 | 10.51 | 3.93 | 59.89 |

By device, the mean (maximum) reciprocal-model error is 0.41 (2.7) pp for Al, 0.33 (3.0) pp for Ag, and 0.16 (1.8) pp for the low-loss stack; by wavelength 0.38 (2.5), 0.29 (3.0), and 0.23 (2.6) pp at 500, 550, and 600 nm. The model tends to slightly underestimate $\eta_{\mathrm{sta}}$.

**Monte Carlo uncertainty.** Ray noise was estimated by re-tracing with two further ray seeds, and realization spread by drawing two further geometries with the same statistics (Table S5). Ray noise is negligible; realization spread is at most 0.24 in $N_{\mathrm{mix}}$, 0.021 in $\xi$, and 0.23 pp in $\eta_{\mathrm{sta}}$, and changes the reduced-model error by at most 0.03 pp.

**Table S5 | Monte Carlo uncertainty (550 nm).** ±: standard deviation over three ray seeds; in parentheses: over three geometric realizations. Last column: reciprocal reduced model minus full series (Al).

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

**Table S6 | Convergence with bin width, cutoff, ray count, and domain size (Al, 550 nm).** Wider bins are obtained by étendue-weighted merging.

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
| hemispherical array ($s$ = 0) | 20,000 rays per bin (library) | 7.05 | 0.499 | 0.7944 | -0.15 |
| hemispherical array ($s$ = 0) | 40,000 rays per bin | 7.02 | 0.497 | 0.7943 | -0.14 |
| wavy, rms slope 0.2 | 2,500 rays per bin | 2.17 | 0.345 | 0.6094 | -2.76 |
| wavy, rms slope 0.2 | 10,000 rays per bin (library) | 2.17 | 0.345 | 0.6084 | -2.71 |
| random lenslet shapes, $s$ = 1 | 6 × 6 supercell (library) | 1.95 | 0.769 | 0.7853 | -0.23 |
| random lenslet shapes, $s$ = 1 | 10 × 10 supercell, new realization | 1.89 | 0.794 | 0.7887 | -0.17 |
| position jitter 0.3 | 6 × 6 supercell (library) | 3.50 | 0.699 | 0.7979 | -0.20 |
| position jitter 0.3 | 10 × 10 supercell, new realization | 3.46 | 0.703 | 0.7962 | -0.21 |
| wavy, rms slope 0.2 | 12 × 12 box (library) | 2.17 | 0.345 | 0.6084 | -2.71 |
| wavy, rms slope 0.2 | 24 × 24 box, new realization | 2.17 | 0.347 | 0.6098 | -2.76 |
| wavy, rms slope 1.0 | 12 × 12 box (library) | 1.13 | 0.871 | 0.7924 | -0.07 |
| wavy, rms slope 1.0 | 24 × 24 box, new realization | 1.13 | 0.867 | 0.7912 | -0.04 |

$\eta_{\mathrm{sta}}$ changes by at most 0.3 pp (planar, 3° bins) and 0.06 pp for all other structures; a 10 × 10 supercell or 24 × 24 wavy box changes $N_{\mathrm{mix}}$ by at most 0.06, $\xi$ by 0.025, and $\eta_{\mathrm{sta}}$ by 0.34 pp, comparable to the realization spread.

---

## Note S4 | Mixing versus escape and the composition window

**Contour templates.** Fig. 3a of the main text scales the $p(\theta)$ shape of the hemispherical MLA. Table S7 repeats the calculation with the shapes of a scattering layer and a wavy surface: at fixed $p$, the whole range of $N_{\mathrm{mix}}$ is worth 1.2–3.2 pp, whereas 0.01 in $p$ is worth about 0.6 pp for every template.

**Table S7 | Sensitivity of $\eta_{\mathrm{sta}}$ to mixing and escape (Al, 550 nm).** Span: change of $\eta_{\mathrm{sta}}$ for $N_{\mathrm{mix}}$ from 1 to 16 at fixed $p$. Slope: change per 0.01 in $p$ at $N_{\mathrm{mix}}$ = 4.

| Template | span at $p$ = 0.30 (pp) | at 0.38 | at 0.42 | slope per 0.01 in $p$ (pp) |
|---|---|---|---|---|
| hemispherical MLA | 2.05 | 2.10 | 2.09 | 0.61 |
| scattering layer ($g$ 0.8, $\tau^*$ 1) | 1.40 | 1.30 | 1.23 | 0.62 |
| wavy surface, rms slope 1.0 | 3.22 | 3.23 | 3.17 | 0.59 |

**Decomposition.** Along the shape sweep, the change of $\eta_{\mathrm{sta}}$ relative to the hemispherical array ($s$ = 0) is split in the reciprocal reduced model into a mixing part ($m$ of the randomized array, $p(\theta)$ of the hemisphere), an escape part (the reverse), and their interaction (Table S8). On Al the escape part dominates at large spread ($-2.20$ against $+0.91$ pp at $s$ = 1), on Ag the two cancel, and on the low-loss stack both stay below 0.15 pp.

**Table S8 | Decomposition of the change in $\eta_{\mathrm{sta}}$ along the shape sweep (550 nm, pp).**

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

**Composition window.** $S_{40\text{–}60} = P_{\mathrm{air}}(40\text{–}60°)/P_{\mathrm{air}}(0\text{–}90°)$, both powers in air including Fresnel losses, over all 38 structures (Table S9). In the low-loss device the lower end (0.274–0.293) belongs to the hemispherical MLAs; all other structures stay within 0.31–0.37.

**Table S9 | Range of $S_{40\text{–}60}$ over all structures.**

| Device, wavelength | min | max |
|---|---|---|
| Al, 500 nm | 0.320 | 0.364 |
| Al, 550 nm | 0.327 | 0.360 |
| Al, 600 nm | 0.339 | 0.361 |
| Ag, 550 nm | 0.322 | 0.361 |
| low loss, 550 nm | 0.274 | 0.366 |

---

## Note S5 | Benchmark platform and design space

All campaigns share the platform of Table S10; the hemisphere is a point of the freeform design space. Bands are polar bins of the far field over the full azimuth, with the Lambertian partition 0.117 / 0.296 / 0.337 / 0.220.

**Table S10 | Common platform and design variables.**

| Quantity | Freeform | Hemispherical reference |
|---|---|---|
| Source | disc, radius 1 mm, CPS dipole distribution, $\eta_{\mathrm{rad}}$ = 0.98, emission 453–753 nm | same |
| Substrate, textured patch | glass, 1.295 mm, $n$ = 1.51; 25 × 25 mm | same |
| Lenslets | hexagonal, ~10 µm radius, N-BK7 index-matched | same |
| Profile | 7-point spline, end points (0, 1) and (1, 0); $x_2$–$x_6$ in 0–1, ascending; $y_2$–$y_6$ in 0–1.5 | quarter circle |
| Height scale (stretch$_Z$) | 0.1–3 | 0.1–3 |
| ETL and HTL thickness | 10–150 nm each | 10–150 nm each |
| Feasibility | no self-intersection or collinear triples; point spacing 0.05–1.0; turning angle ≤ 120° | always satisfied |
| Free variables | 13 | 3 |

**Searches and fidelity.** Each arm is optimized by MATLAB `surrogateopt` (radial-basis surrogate) followed by `patternsearch` at 10,000 rays and 31 wavelengths; every optimum is re-evaluated three times at 50,000 rays and 151 wavelengths. EQE is the radiative efficiency times the spectrally weighted fraction of generated photons reaching air, over the hemisphere or a band; it includes the cavity's coupling into the substrate and is not directly comparable with $\eta_{\mathrm{sta}}$. Infeasible candidates are returned to the optimizer as constraint violations and not traced. The 25 mm patch underestimates the total of an infinite film by about 4%, but no band share changes by more than 0.24 pp between 15 and 100 mm, and all comparisons are made at equal patch.

**Fig. S2 | Platform and benchmark design.** (a) OLED–substrate–MLA architecture. (b) The hemispherical reference and an axisymmetric freeform profile (40–60° optimum; open circles, free control points) in units of the lens radius. (c) The four polar bands. (d) Optimization workflow. *(`fig1_platform.png`)*

---

## Note S6 | Restart control and statistics

In the original per-band campaign (random seeds plus the best designs of a separate weighted-sum campaign) the 0–20° and 20–40° gains were below unity (0.946, 0.966). Because the hemisphere lies in the freeform design space, this was a search failure. The restart control repeats every search from that arm's hemispherical optimum: the hemisphere and eight perturbations within 8% of each range seed `surrogateopt`, and `patternsearch` runs both from the surrogate winner and directly from the hemisphere. The hemisphere is re-measured in the same session.

**Table S11 | Restart control.** Mean ± standard deviation over three re-evaluations.

| Arm | hemisphere (archived) | hemisphere (re-measured) | dev. | restarted | gain | $t$ | winning branch |
|---|---|---|---|---|---|---|---|
| 0–20° | 0.06682 | 0.06677 ± 0.00017 | 0.07% | 0.06699 ± 0.00013 | +0.32% | 1.7 | polish from hemisphere |
| 20–40° | 0.16591 | 0.16595 ± 0.00017 | 0.02% | 0.16628 ± 0.00020 | +0.20% ‡ | 2.1 ‡ | surrogateopt |
| 40–60° | 0.18503 | 0.18486 ± 0.00008 | 0.09% | 0.19355 ± 0.00010 | **+4.70%** | **114.6** | polish from surrogate |
| 60–80° | 0.13589 | 0.13581 ± 0.00014 | 0.06% | 0.14075 ± 0.00025 | **+3.64%** | **29.6** | surrogateopt |
| total EQE | 0.54679 | 0.54679 ± 0.00013 | 0.00% | 0.55231 ± 0.00004 | **+1.01%** | **68.7** | polish from surrogate |

‡ Settled by re-measuring both stored designs five times each: 0.16583 ± 0.00015 against 0.16631 ± 0.00021, +0.29% at $t$ = 4.1 (threshold 1.86). In Table 2 of the main text, the 40–60° entry is the per-band optimum (0.19802), which the restart did not surpass, and the total is the best of three dedicated starts (0.5539).

**Independent searches.** Freeform/hemisphere ratios: 0–20° 0.946 (original campaign, 60 + 15 evaluations) and 1.003 (restart, about 47); 20–40° 0.966 and 1.002; 40–60° 1.070 and 1.046; 60–80° 1.020 and 1.036; total EQE 1.005, 1.011, and 1.013 (three random starts, 140 evaluations each) and 1.010 (restart). The hemisphere used about 45 evaluations for 3 variables (15 per variable), the freeform searches 47–140 for 13 (4–11 per variable), so the comparison favours the hemisphere. Inverted (concave) and randomly assembled MLA families optimized on the same platform reached total EQEs of 0.517 and 0.522, below the hemisphere.

**What the $t$ values test.** $t = \Delta/(s_p\sqrt{2/n})$ with $n$ = 3 and pooled standard deviation $s_p$ compares two fixed designs against the Monte Carlo noise of one evaluation (about $10^{-4}$ in EQE); with a Holm correction over five objectives [S. Holm, Scand. J. Stat. 6, 65 (1979)] the thresholds are 3.75, 3.50, 3.19, 2.78, and 2.13. The large values ($t$ = 30–115) reflect this small noise. They show that a reported difference is not an evaluation artefact but say nothing about the uncertainty of the optimizer, which only the spread over independent searches above measures. The restart winner is also the best of about seventy noisy evaluations, so its +0.3% residues at 0–40° are upper limits.

**Fig. S3 | Restart control.** (a) Re-measured hemisphere and restarted optimum for every arm, with standard deviations and $t$. (b) Five-repeat confirmation of the 20–40° arm. *(`figS1_restart_control.png`)*

---

## Note S7 | Total efficiency fixes the angular composition

Over the 606 usable evaluations of a weighted-sum campaign (114 random feasible designs, 492 optimizer-visited), each band EQE is nearly proportional to the total (Table S12, Fig. S4); random designs alone reproduce the collapse. $S_{40\text{–}60}$ spans only 0.336–0.375 (5–95%), and the 40–60° band carries the largest share in every design with total EQE ≥ 0.40. Band-dedicated optima raise their band's share by 1–27% but lose 1–11% of total EQE, for net band-power gains of 1.00–1.22 against the natural composition of the twenty most efficient designs.

**Table S12 | Band EQE against total EQE (least squares, 95% intervals).**

| Band | slope | intercept | $R^2$ (all 606) | $R^2$ (114 random only) |
|---|---|---|---|---|
| 0–20° | 0.125 ± 0.004 | −0.014 ± 0.002 | 0.848 | 0.828 |
| 20–40° | 0.318 ± 0.007 | −0.020 ± 0.004 | 0.931 | 0.925 |
| 40–60° | 0.348 ± 0.005 | 0.003 ± 0.003 | 0.968 | 0.967 |
| 60–80° | 0.189 ± 0.008 | 0.026 ± 0.004 | 0.769 | 0.724 |

**Fig. S4 | Total efficiency fixes the angular composition.** (a) Total EQE against 40–60° EQE for random (grey) and optimizer-visited (green) designs; line, $R^2$ = 0.968. (b) Selectivity of the band-dedicated optima (red) against the natural composition (grey, 10–90% range) and the Lambertian partition (dashed), labelled with the net band-power gain. (c) Selectivity against total EQE for all bands. *(`fig3_efficiency_composition.png`)*

---

## Note S8 | Lens versus cavity

**Full device model.** The benchmark optimizes the cavity together with the lens. To separate the two, the eight optimized lenses of Table 3 were evaluated in the full device model of the main text (birefringent CPS + TMM, dispersive constants), weighted with the emission spectrum (400–700 nm, 10 nm steps), with $\eta_{\mathrm{rad}}$ = 0.98 and the polar bands of the benchmark (Table S13). "Lens only" places both lenses at the hemisphere's cavity; "cavity only" places the hemisphere at the freeform's cavity.

**Table S13 | Freeform/hemisphere EQE ratio: LightTools and the multi-pass model.**

| Objective | cavity ETL/HTL (nm), hemisphere / freeform | LightTools | model, own cavities | model, lens only | model, cavity only |
|---|---|---|---|---|---|
| 20–40° | 54.9/52.1 / 61.2/75.8 | 0.966 | 0.959 | 0.946 | 0.974 |
| 40–60° | 74.6/81.7 / 68.5/114.4 | 1.070 | 1.026 | 1.004 | 0.986 |
| 60–80° | 99.8/71.1 / 90.5/76.3 | 1.020 | 0.972 | 0.964 | 0.996 |
| total EQE | 61.3/97.5 / 67.5/120.6 | 1.013 | 0.974 | 0.971 | 1.003 |

The freeform entries are the designs of the original per-band campaign and of the total-EQE campaign (for 60–80° the restart optimum, ratio 1.036, has no traced return matrix). For the hemispheres the model reproduces the LightTools EQE within +0.7% (total) to +3.8% (60–80°); for the freeform optima it lies between −3.2% (total) and +0.5% (20–40°). The two tools therefore agree on the size of the differences (a few percent) but not on their sign below about 5%; the main reason is the profile: the freeform curves are rebuilt with a PCHIP spline, whereas LightTools uses its own spline, and the total-EQE optimum has a non-monotonic profile with an annular trough (control-point heights 1.05, 0.15, 0.19, 0.53, 0.50), for which the reconstruction is least certain. Within the model, the lens shape alone changes the band EQE by −5.4% to +0.4% and the cavity by −2.6% to +0.3%. This also explains Table 3, where the total-EQE freeform lens, at the hemisphere's cavity and 550 nm, has $\eta_{\mathrm{sta}}$ = 0.771 against 0.792.

**Lens-only optimization at a fixed stack.** {{ES}}

**LightTools control.** `opt_fixedcavity_freeform.m` repeats the per-band freeform optimization of `opt_4band_freeform.m` with ETL and HTL fixed at each arm's hemispherical optimum and the hemisphere as the first seed, so that $G_j \ge 1$ holds by construction and any gain is due to the lens shape; it requires the LightTools pipeline and has not been run here.

---

## Reproducing the figures and data

Main-text Figs. 1–4 are produced by `selectivity/oe_figs.py` (structure library, full and reduced models; PNG, PDF, and an Excel workbook of the plotted arrays in `figures_oe/`), with return matrices from `multipass_reset/run_random_sweep.py`. Fig. 5 and Figs. S2–S5 come from the campaign archives through `make_figures.py` and `make_supp_figures_and_data.py`; Fig. S1 and Table S2 from `selectivity/fig_devices.py`; Tables S5–S9 from `selectivity/unc_analysis.py`, `selectivity/review_checks.py`, `selectivity/review_traces_analysis.py`, and `selectivity/make_note_s12.py`.

```
cd multipass_reset && python3 run_random_sweep.py && cd ..
cd selectivity && python3 extract_inputs.py && python3 run_devices.py && python3 oe_figs.py && cd ..
python3 make_figures.py && python3 make_supp_figures_and_data.py
```
