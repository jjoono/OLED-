# Supplement 1

Reflective and transmissive angular mixing affect substrate light extraction differently in organic light-emitting diodes

Note S1 gives the model details and numerical checks, Note S2 the devices and structures, Note S3 the accuracy and limits of the reduced model, Note S4 the comparison of the mixing placements, Note S5 the étendue bound and the freeform benchmark, and Note S6 the fabricated devices and the finite-geometry model. References of the main text are cited as "main ref. n".

---

## Note S1. Model details and numerical checks

**Scalar model.** $R(\theta)$ is the polarization average of the stack reflectance, $\mathbf{v}^{(0)}$ is the sum of the s- and p-polarized power launched into the substrate, and the tracers apply the unpolarized Fresnel coefficients. Textures without axial symmetry are reduced to a polar-angle kernel by launching rays uniformly over the cell and in azimuth and binning the exit directions by polar angle. $\eta_{\mathrm{sta}}$ is obtained by a direct linear solve on the full 0–90° grid, and the spectral radius of $\mathrm{diag}\,\mathbf{R}\,B_R$ stays below 0.9991 in all cases.

**Angular cutoff.** Bins near 90° carry almost no étendue but can return light at their own angle, and a few such states can dominate $\lambda_2$ without affecting the extraction. The spectral quantities ($\lambda_2$, $N_{\mathrm{mix}}$, $p_\infty$) are therefore evaluated on $\theta < 85°$, which excludes 0.76% of the étendue and less than 0.1% of $\mathbf{v}^{(0)}$. $B_R$ is restricted to these bins, $M$ is formed by dividing each column of the sub-block by its sum $s_i$, and $p_\infty$ is taken from the sub-block without this renormalization. Everything else, including $\eta_{\mathrm{sta}}$, $p$, $\sigma_p^2$, and $p_{\mathrm{mix}}$, uses the full range. For every structure $\lambda_2$ is real and positive (0.09–0.88). The error of the reduced model depends little on the cutoff between 80° and 85° and grows above it (Table S1). The absolute values of $N_{\mathrm{mix}}$ and $\xi$ depend on the cutoff, but their ranking is largely preserved (rank correlation with the 85° values of 0.96–0.99 for $N_{\mathrm{mix}}$ and 0.90–0.98 for $\xi$).

**Table S1.** Error of the reciprocal reduced model against the angular cutoff of the spectral quantities (pp; class columns give the mean error over three devices and three wavelengths).

| Cutoff | mean error (pp) | max error (pp) | hemispherical MLA | freeform MLA | random MLA | wavy | scattering |
|---|---|---|---|---|---|---|---|
| 80° | 0.32 | 3.02 | 0.21 | 0.64 | 0.20 | 0.66 | 0.29 |
| 85° | 0.30 | 2.99 | 0.26 | 0.60 | 0.16 | 0.67 | 0.29 |
| 88° | 0.40 | 3.55 | 0.38 | 0.95 | 0.25 | 0.67 | 0.30 |
| 90° | 0.55 | 4.72 | 0.45 | 1.16 | 0.52 | 0.68 | 0.30 |

**Complete reciprocal mixing.** For $B_R(i\to r) = (1-p_i)[(1-m)\delta_{ir} + m q_r]$, detailed balance requires $w_i(1-p_i)q_r = w_r(1-p_r)q_i$, i.e., $q \propto w(1-p)$. With this $q$, $p_{\mathrm{mix}} = \sum_r q_r p_r = (p - \langle p^2\rangle)/(1-p) = p - \sigma_p^2/(1-p)$, where $\langle\cdot\rangle$ is the étendue-weighted mean, and $p_{\mathrm{mix}} \le p$ follows from the Cauchy–Schwarz inequality. The rank-one kernel gives $(I-\mathbf{u}\mathbf{v}^{\mathsf T})^{-1} = I + \mathbf{u}\mathbf{v}^{\mathsf T}/(1-\mathbf{v}^{\mathsf T}\mathbf{u})$ with $u_k = R_k q_k$ and $\mathbf{v} = \mathbf{1}-\mathbf{p}$, which is the closed form of main-text Eq. (7); it agrees with the numerical series to $10^{-9}$ for all 342 cases. The ratio $p_{\mathrm{mix}}/p$ is 0.87–0.91 for hemispherical MLAs, 0.81–0.87 for freeform MLAs, 0.84–0.93 for random MLAs, 0.37–0.90 for wavy surfaces, 0.59–0.99 for scattering layers, and 0.14 for the planar interface.

**Reference, not a bound.** The complete-mixing kernel is the only reciprocal rank-one kernel with a given $p(\theta)$, but other reciprocal kernels exist. Two states with $w$ = (1, 1), $p$ = (0.4, 0), $R$ = 0.8, and all light launched into the second state give $\eta$ = 0.375 for the rank-one kernel [[0.225, 0.375], [0.375, 0.625]]. The kernel [[0, 0.6], [0.6, 0.4]] has the same column sums and is also reciprocal, but it sends all light returned from the trapped state into the escaping one and gives $\eta$ = 0.427. Among the structures studied, only two wavy surfaces exceed $\eta_{\mathrm{mix}}$, by at most 0.4 pp.

**Bounds on $p_\infty$.** With $D = \mathrm{diag}(w)$, $D^{-1/2}B_R D^{1/2}$ is symmetric, so $\rho(B_R)$ is the maximum of its Rayleigh quotient. The uniform-radiance vector gives the quotient $1-\tilde p$, hence $p_\infty \le \tilde p$, where $\tilde p$ is the probability that uniform-radiance light leaves the truncated chain in one encounter, by escape or by return beyond the cutoff ($\tilde p - p < 0.008$). The column sums give $p_\infty \ge \min_\theta p(\theta)$. All structures satisfy both bounds. The device-aware variant $1-\rho(\mathrm{diag}\,\mathbf{R}\,B_R)$ ranks the structures like $p_\infty$ (Spearman coefficient 0.96).

**Numerical checks.** The relative asymmetry $\sum|F-F^{\mathsf T}|/\sum|F|$ of $F_{ji} = w_i B_R(i\to j)$ is 0.05–0.13 on 1° bins and 0.002–0.007 on 15° bins, as expected for Monte Carlo noise. Symmetrizing $F$ changes $N_{\mathrm{mix}}$ by at most 0.20 and $\eta_{\mathrm{sta}}$ by at most 0.04 pp. Merging the 1° bins into 2°, 3°, or 5° bins leaves $p$ unchanged to $10^{-4}$ and changes $\eta_{\mathrm{sta}}$ by at most 0.06 pp, except for the planar interface, whose $p(\theta)$ steps at the critical angle ($p_{\mathrm{mix}}$ = 0.054, 0.056, and 0.061 for 0.5°, 1°, and 3° bins; $\eta_{\mathrm{sta}}$ = 0.3841 and 0.3864 for 1° and 3° bins). Ray noise is negligible, and the spread over geometric realizations is at most 0.24 in $N_{\mathrm{mix}}$ and 0.23 pp in $\eta_{\mathrm{sta}}$ (Table S2). A larger supercell or wavy box changes $\eta_{\mathrm{sta}}$ by at most 0.34 pp. For the scattering layers, a direct Monte Carlo trace of the slab together with its exit face agrees with the adding method to within 0.0005 in $\eta_{\mathrm{sta}}$. The off-diagonal bands in the return matrix of the hemispherical MLA [main-text Fig. 1(c)] are fixed paths of rays that are totally reflected at the lens flank, for example from 70° to about 34°; they persist in a 6 × 6 supercell and with twice the number of rays.

**Table S2.** Monte Carlo uncertainty (550 nm). ±: standard deviation over three ray seeds; in parentheses: over three geometric realizations. Last column: reciprocal reduced model minus full series (Al device).

| Structure | N_mix | ξ | η Al | err. Al (pp) |
| --- | --- | --- | --- | --- |
| A_shape0.0 | 7.05 ± 0.01 (0.00) | 0.499 ± 0.003 (0.000) | 0.7944 ± 0.0001 (0.0000) | -0.15 ± 0.01 (0.00) |
| A_shape1.0 | 1.95 ± 0.00 (0.12) | 0.769 ± 0.000 (0.009) | 0.7853 ± 0.0003 (0.0023) | -0.23 ± 0.00 (0.03) |
| B_aspect0.6 | 4.08 ± 0.00 (0.24) | 0.686 ± 0.002 (0.021) | 0.7994 ± 0.0002 (0.0001) | -0.19 ± 0.01 (0.01) |
| C_jitter0.3 | 3.50 ± 0.01 (0.10) | 0.699 ± 0.001 (0.006) | 0.7979 ± 0.0001 (0.0005) | -0.20 ± 0.00 (0.01) |
| W_slope0.2 | 2.17 ± 0.01 (0.01) | 0.345 ± 0.001 (0.001) | 0.6084 ± 0.0006 (0.0008) | -2.71 ± 0.01 (0.02) |
| W_slope1.0 | 1.13 ± 0.00 (0.01) | 0.871 ± 0.001 (0.007) | 0.7924 ± 0.0002 (0.0023) | -0.07 ± 0.02 (0.01) |
| scatter g0.8_t1.0 | 1.17 ± 0.00 | 0.977 ± 0.000 | 0.7374 ± 0.0000 | -0.08 ± 0.00 |
| scatter g0.0_t0.25 | 2.25 ± 0.01 | 0.545 ± 0.000 | 0.6754 ± 0.0001 | -2.07 ± 0.01 |

**Polarization.** As the opposite limit to the scalar model, s and p were treated as separate channels that keep their label at every encounter. This limit is exact for a planar interface; for curved surfaces the true value lies between the two. It shifts $\eta_{\mathrm{sta}}$ of the hemispherical MLA and of a steep wavy surface by 0.6–1.5 pp, but the complete-mixing gain by at most 0.4 pp and the error of the reduced model by at most 0.1 pp.

---

## Note S2. Devices and structures

**Devices.** The emitters form a dipole sheet at the center of the 25 nm birefringent emission layer, with a horizontal dipole ratio of 0.865 (Table S3, Fig. S1).

**Table S3.** Device stacks (from the cathode to the substrate).

| Device | Cathode | ETL | EML | HTL-side layers | Anode | Substrate | $A'$ at 550 nm |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Al | Al 100 nm | B3PyMPM 61.3 nm | TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd 25 nm | TCTA 10 nm / TAPC 97.5 nm | ITO 150 nm | glass, $n$ = 1.51 | 0.167 |
| Ag | Ag 100 nm | B3PyMPM 61.3 nm | as Al | TCTA 10 nm / TAPC 97.5 nm | ITO 150 nm | glass, $n$ = 1.51 | 0.044 |
| Low loss | Ag 100 nm | B3PyMPM 200 nm | as Al | TAPC 180 nm | ITO 50 nm | glass, $n$ = 1.51 | 0.026 |

**Fig. S1.** Device inputs at 550 nm. (a) Radiance launched into the substrate on the first pass, $v^{(0)}_k/w_k$, relative to its mean over 0–20°. (b) Angle-resolved reflectance $R(\theta)$ seen from the substrate. The dotted line marks the critical angle of the glass/air interface. *(`figS11_device_inputs.png`)*

**Tracer and structures.** In the height-field tracer, the surface $z = h(x, y)$ separates the index-matched medium ($n$ = 1.51) from air, and rays are reflected or refracted with the unpolarized Fresnel coefficients (20,000 rays per 1° bin, 10,000 for wavy surfaces). For a hemispherical lenslet it reproduces the LightTools return matrix of main ref. 13 (escape probability 0.4275 against 0.4274). Random MLAs tile a 6 × 6 supercell: the shape spread $s$ perturbs the control points of the quarter-circle profile, the height spread scales each lenslet by $1\pm s$, and the position jitter displaces each center by a fraction of the pitch, optionally with a radius spread of ±20%. Wavy surfaces are periodic Gaussian-correlated height fields scaled to a prescribed rms slope. Scattering layers are non-absorbing, index-matched Henyey–Greenstein slabs combined with the flat exit face by the adding method. Table S4 lists all structures.

**Table S4.** The 38 structures at 550 nm. $N_{\mathrm{mix}}$ and $\xi$ are evaluated on $\theta < 85°$; $\eta$ is the full series for the Al, Ag, and low-loss (LL) devices; the last column is the reciprocal reduced model minus the full series (Al device). A, lenslet shape spread $s$; B, height spread; C, position jitter; C2, as C with ±20% radius spread; W, wavy surface of the given rms slope. Scattering layers are labeled by $g$ and $\tau^*$, and MLAs by the benchmark objective of the optimum they reproduce (Note S5).

| Structure | p | N_mix | ξ | η Al | η Ag | η LL | err. recip. Al (pp) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| planar | 0.398 | ∞ | 0.00 | 0.384 | 0.356 | 0.664 | -0.00 |
| scatter g0.0_t0.25 | 0.381 | 2.2 | 0.54 | 0.675 | 0.874 | 0.945 | -2.07 |
| scatter g0.0_t1.0 | 0.324 | 1.1 | 0.93 | 0.725 | 0.907 | 0.949 | -0.10 |
| scatter g0.8_t0.25 | 0.386 | 1.2 | 0.81 | 0.749 | 0.914 | 0.959 | -0.19 |
| scatter g0.8_t1.0 | 0.326 | 1.2 | 0.98 | 0.737 | 0.913 | 0.951 | -0.08 |
| scatter g0.9_t1.0 | 0.326 | 1.2 | 0.98 | 0.738 | 0.913 | 0.951 | -0.10 |
| scatter g0.8_t4.0 | 0.188 | 1.2 | 0.99 | 0.579 | 0.838 | 0.901 | -0.06 |
| MLA free_20-40 | 0.429 | 3.8 | 0.58 | 0.788 | 0.925 | 0.967 | -0.69 |
| MLA free_40-60 | 0.425 | 7.5 | 0.39 | 0.789 | 0.929 | 0.961 | -1.11 |
| MLA free_60-80 | 0.417 | 3.0 | 0.79 | 0.785 | 0.929 | 0.966 | -1.46 |
| MLA free_total | 0.424 | 4.6 | 0.46 | 0.771 | 0.914 | 0.961 | -0.89 |
| MLA hemi_0-20 | 0.423 | 7.3 | 0.44 | 0.790 | 0.924 | 0.964 | -0.38 |
| MLA hemi_20-40 | 0.422 | 7.0 | 0.47 | 0.793 | 0.926 | 0.964 | -0.21 |
| MLA hemi_40-60 | 0.425 | 8.0 | 0.42 | 0.783 | 0.920 | 0.965 | -0.64 |
| MLA hemi_60-80 | 0.420 | 6.2 | 0.50 | 0.792 | 0.926 | 0.965 | -0.24 |
| MLA hemi_total | 0.422 | 7.1 | 0.47 | 0.792 | 0.925 | 0.964 | -0.27 |
| A_shape0.0 | 0.422 | 7.0 | 0.50 | 0.794 | 0.927 | 0.965 | -0.15 |
| A_shape0.1 | 0.421 | 3.8 | 0.65 | 0.797 | 0.931 | 0.966 | -0.30 |
| A_shape0.2 | 0.421 | 2.7 | 0.78 | 0.799 | 0.934 | 0.966 | -0.28 |
| A_shape0.4 | 0.422 | 2.5 | 0.80 | 0.794 | 0.932 | 0.966 | -0.40 |
| A_shape0.7 | 0.426 | 1.9 | 0.79 | 0.788 | 0.929 | 0.966 | -0.27 |
| A_shape1.0 | 0.428 | 1.9 | 0.77 | 0.785 | 0.927 | 0.966 | -0.23 |
| B_aspect0.2 | 0.422 | 7.1 | 0.43 | 0.795 | 0.928 | 0.965 | -0.20 |
| B_aspect0.4 | 0.422 | 6.0 | 0.50 | 0.797 | 0.930 | 0.966 | -0.18 |
| B_aspect0.6 | 0.423 | 4.1 | 0.69 | 0.799 | 0.933 | 0.968 | -0.19 |
| C2_jitter0.05_r20 | 0.421 | 4.2 | 0.66 | 0.798 | 0.931 | 0.966 | -0.12 |
| C2_jitter0.1_r20 | 0.421 | 3.5 | 0.70 | 0.796 | 0.930 | 0.966 | -0.15 |
| C2_jitter0.2_r20 | 0.420 | 3.7 | 0.69 | 0.799 | 0.932 | 0.966 | -0.13 |
| C2_jitter0.3_r20 | 0.420 | 3.2 | 0.74 | 0.796 | 0.932 | 0.966 | -0.21 |
| C_jitter0.05 | 0.421 | 6.2 | 0.50 | 0.796 | 0.928 | 0.965 | -0.15 |
| C_jitter0.1 | 0.421 | 5.0 | 0.57 | 0.797 | 0.930 | 0.966 | -0.15 |
| C_jitter0.2 | 0.421 | 3.9 | 0.66 | 0.798 | 0.932 | 0.966 | -0.18 |
| C_jitter0.3 | 0.420 | 3.5 | 0.70 | 0.798 | 0.932 | 0.966 | -0.20 |
| W_slope0.2 | 0.402 | 2.2 | 0.35 | 0.608 | 0.815 | 0.933 | -2.71 |
| W_slope0.4 | 0.410 | 1.2 | 0.63 | 0.719 | 0.896 | 0.957 | -1.27 |
| W_slope0.7 | 0.417 | 1.2 | 0.81 | 0.773 | 0.923 | 0.965 | -0.37 |
| W_slope1.0 | 0.422 | 1.1 | 0.87 | 0.792 | 0.932 | 0.968 | -0.07 |
| W_slope1.5 | 0.427 | 1.2 | 0.88 | 0.801 | 0.935 | 0.968 | -0.09 |

---

## Note S3. Accuracy and limits of the reduced model

**Library structures.** Table S5 gives the errors by structure class, and Fig. S2 compares three simpler models with the full series. By device, the mean (maximum) error of the reciprocal reduced model is 0.41 (2.7) pp for Al, 0.33 (3.0) pp for Ag, and 0.16 (1.8) pp for the low-loss stack. The model slightly underestimates $\eta_{\mathrm{sta}}$ on average (mean signed error −0.24 pp).

**Table S5.** Error of the reduced and fully randomized models by structure class (pp; 342 cases, three devices, 500, 550, and 600 nm).

| Class | cases | reciprocal mean | reciprocal max | q ∝ w mean | q ∝ w max | fully randomized mean | fully randomized max |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Planar | 9 | 0.00 | 0.00 | 0.00 | 0.00 | 44.05 | 59.89 |
| Scattering layers | 54 | 0.29 | 2.49 | 0.96 | 4.67 | 1.83 | 11.14 |
| Hemispherical and freeform MLAs | 81 | 0.41 | 2.09 | 0.38 | 2.22 | 3.45 | 12.33 |
| Random lens arrays (A, B, C) | 153 | 0.16 | 0.51 | 0.44 | 1.93 | 2.41 | 10.89 |
| Wavy, rms slope ≥ 0.7 | 27 | 0.19 | 1.06 | 1.11 | 3.07 | 1.61 | 5.57 |
| Wavy, rms slope ≤ 0.4 | 18 | 1.40 | 2.99 | 5.40 | 10.51 | 8.66 | 21.73 |
| All | 342 | 0.30 | 2.99 | 0.81 | 10.51 | 3.93 | 59.89 |

**Fig. S2.** Simpler models against the full series for the 342 cases. $\eta_{\mathrm{sta}}$ from (a) the fully randomized model, (b) the angle-preserving limit ($m = 0$), and (c) the reduced model with isotropic redistribution ($q \propto w$). *(`figS_model_comparison.png`)*

**Direct ray tracing and out-of-sample structures.** Rays were launched with the first-pass distribution $\mathbf{v}^{(0)}$ (Al device, 550 nm), traced through the height field and, when returned, absorbed with probability $1-R(\theta)$ or reflected specularly at a new lateral position ($10^6$ rays per case, standard error 0.04–0.05 pp). A second run re-randomized the azimuth at every pass, as the chain on polar angle assumes. Nine structures were chosen after the reduced model had been fixed: spherical caps, paraboloids, cones with side angles of 45° and 60°, Gaussian bumps, a random mixture of hemispheres and cones, and wavy surfaces with a ring spectrum (Table S6). Without the two cone arrays, the reduced model has a mean error of 0.48 pp and a maximum of 3.3 pp over three devices and three wavelengths, against 2.62 and 8.8 pp for the fully randomized model. The cone arrays fail, with maximum errors of 20 and 14 pp. The direct traces confirm the full series for the 60° cones to within 0.6 pp, so the failure belongs to the one-rate reduction and not to the angle-resolved model.

**Table S6.** Direct ray tracing and out-of-sample test (Al device, 550 nm; reduced and fully randomized models minus full series).

| Structure | $p$ | $N_{\mathrm{mix}}$ | $\xi$ | $\eta_{\mathrm{sta}}$ full series | direct, $(\theta,\phi)$ | direct, $\phi$ re-randomized | reduced − full (pp) | fully randomized − full (pp) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
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

**Conditions of validity and prospective test.** The 45° cones retroreflect light within 0–4° into the same range, which gives a closed set and $\lambda_2 = 1$. The 60° cones have $\lambda_2$ = 0.985, but this slow mode carries little of the returned light. To quantify this, the deviation of the first returned light, $B_R\mathbf{v}^{(0)}$, from the stationary distribution is expanded in the eigenmodes of $M$ with amplitudes $a_\ell$, and the weighted relaxation rate is $\lambda_{\mathrm{eff}} = \sum_{\ell \ge 2} a_\ell^2 \lambda_\ell/\sum_{\ell \ge 2} a_\ell^2$. Over the library, $\lambda_{\mathrm{eff}}$ is less accurate than $\lambda_2$ as the model parameter, so $\lambda_2$ is kept and the gap $\lambda_2 - \lambda_{\mathrm{eff}}$ serves as a warning. Because these two conditions ($\lambda_2 \le 0.999$ and a gap of at most 0.35) were found after the out-of-sample test, twelve further structures were traced: cones with side angles of 30–70°, square and inverted pyramids, and prism gratings. Each structure was classified from $M$ and the first-pass distribution before the multi-pass series was evaluated (Table S7). The 45° pyramids and prisms are reducible, were predicted to fail, and fail with maximum errors of 21 and 23 pp. The gap condition did not discriminate: it flagged two structures that the model describes to within 2.5 pp, and the eight structures that passed both conditions have a mean error of 1.4 pp but a maximum of 7.5 pp (40° cones). Irreducibility is therefore a necessary condition that predicts the large failures, whereas the gap is not a reliable predictor.

**Table S7.** Relaxation rates and error of the reduced model (mean / maximum over three devices and three wavelengths, pp). For the prospective rows, the prediction was recorded before the multi-pass series was evaluated.

| Structure | $\lambda_2$ | $\lambda_{\mathrm{eff}}$ | error with $\lambda_2$ |
| --- | --- | --- | --- |
| hemispherical MLA | 0.859 | 0.734 | 0.24 / 0.62 |
| freeform MLA, 40–60° | 0.867 | 0.723 | 0.87 / 2.09 |
| wavy, rms slope 0.2 | 0.540 | 0.195 | 1.98 / 2.99 |
| random lenslet shapes, $s$ = 1 | 0.487 | 0.376 | 0.12 / 0.23 |
| spherical caps, 35° (out-of-sample) | 0.296 | 0.237 | 1.42 / 3.30 |
| paraboloids (out-of-sample) | 0.636 | 0.379 | 0.33 / 0.83 |
| cone–hemisphere mix (out-of-sample) | 0.634 | 0.397 | 0.27 / 0.82 |
| cones, 60° (out-of-sample) | 0.985 | 0.430 | 4.15 / 14.36 |
| cones, 45° (out-of-sample) | 1.000 | 0.988 | 5.99 / 20.15 |
| cones, 30° side (prospective, predicted valid) | 0.511 | 0.438 | 1.94 / 4.89 |
| cones, 40° side (prospective, predicted valid) | 0.761 | 0.539 | 2.34 / 7.49 |
| cones, 50° side (prospective, predicted valid) | 0.608 | 0.421 | 0.72 / 2.66 |
| cones, 55° side (prospective, predicted outside) | 0.641 | 0.269 | 0.77 / 2.46 |
| cones, 70° side (prospective, predicted valid) | 0.689 | 0.476 | 0.30 / 0.56 |
| inverted pyramids, 54.7° (prospective, predicted valid) | 0.578 | 0.243 | 0.89 / 2.30 |
| prism grating, 30° (prospective, predicted valid) | 0.341 | 0.268 | 1.11 / 2.90 |
| prism grating, 45° (prospective, predicted outside) | 1.000 | 0.992 | 8.19 / 23.40 |
| prism grating, 60° (prospective, predicted valid) | 0.792 | 0.461 | 1.73 / 4.19 |
| square pyramids, 30° (prospective, predicted valid) | 0.519 | 0.319 | 2.30 / 5.68 |
| square pyramids, 45° (prospective, predicted outside) | 1.000 | 0.986 | 6.61 / 21.09 |
| square pyramids, 54.7° (prospective, predicted outside) | 0.546 | 0.024 | 0.70 / 1.96 |

---

## Note S4. Reflective and transmissive mixing

**Device kernel.** Let the device reflect returned light with a kernel $D(i\to j)$. Energy conservation fixes the column sums, $\sum_j D(i\to j) = R_i$, and reciprocity requires $w_i D(i\to j) = w_j D(j\to i)$. A rank-one kernel $a_i b_j$ satisfies the second condition only if $b_j \propto w_j a_j$, and the first then gives $a_i \propto R_i$, which is main-text Eq. (8). The multi-pass efficiency is $\eta = \mathbf{p}^{\mathsf T}(I - D\,B_R)^{-1}\mathbf{v}^{(0)}$. For the three stacks of this work, $|p_D - p| < 0.008$ for all structures at 500, 550, and 600 nm, and complete mixing at the device adds 1.5–4.1 pp for lens arrays in the Al device at 550 nm, against at most 1.8 pp at the escape interface. As a counterexample, a top-emitting microcavity with a 20 nm Ag cathode and a planar exit face absorbs most strongly inside the escape cone (reflectance 0.63 inside and 0.92 outside), so the covariance is negative and $p_D$ = 0.316 against $p$ = 0.355. A device that reflects a fraction $m$ of the light with $D$ and the rest specularly is also reciprocal. Its gain grows nearly linearly with $m$ and is at most 2.3 pp at $m$ = 1 for lens arrays in the Al device (spectrally weighted).

**Transmissive mixing.** An ideal lossless randomizer with Lambertian transmission and no reflection, placed between the OLED and a planar exit face, gives $\eta_{\mathrm{sta}}$ = 0.798, 0.937, and 0.962 for the Al, Ag, and low-loss stacks at 550 nm, equal to the fully randomized model to three digits. Non-absorbing Henyey–Greenstein volumes with $g$ = 0.9 and $\tau^*$ = 0.5 reach 0.764, 0.923, and 0.959. At larger $\tau^*$ the efficiency falls because the backscatter of the volume rises (Table S8). Placing a layer directly on the OLED or directly under the planar face gives the same $\eta_{\mathrm{sta}}$, because the index-matched glass between them neither absorbs nor changes the angle of the light. Under the hemispherical MLA, a layer with $\tau^*$ = 1 lowers $\eta_{\mathrm{sta}}$ from 0.792 to 0.734, since it replaces the high escape probability of the lenses by its own.

**Table S8.** Henyey–Greenstein volumes under a planar exit face (Al device, 550 nm): diffuse reflectance of the bare volume for Lambertian incidence, escape probability $p$ of the volume with its exit face, and $\eta_{\mathrm{sta}}$.

| $g$ | $\tau^*$ | Volume reflectance | $p$ (volume + face) | $\eta_{\mathrm{sta}}$ |
| --- | --- | --- | --- | --- |
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

**Gain by device and wavelength.** The gain from complete mixing by the structure depends on how the first-pass distribution overlaps $p(\theta)$ (Table S9). When the cavity launches light into angles that the hemispherical MLA transmits poorly (Al and Ag at 500 nm, low loss at 600 nm), complete mixing would add 3–9 pp, and random MLAs and wavy surfaces realize most of it. Elsewhere the gain is below about 2 pp. Table S10 compares the three placements with spectral weighting (400–700 nm in 10 nm steps, emitter spectrum times the substrate-coupled fraction). Along the sweep of the lenslet shape spread, the reduced model separates the effect of faster mixing from that of the changed $p(\theta)$: at $s$ = 1 in the Al device they are +0.91 and −2.20 pp.

**Table S9.** Gain from complete mixing by the structure and realized gain from randomization (pp). First pass: first-pass escape probability of the hemispherical MLA ($p$ = 0.422). Last two columns: best random MLA and wavy surface of slope 1.5 minus the hemispherical MLA.

| Device, wavelength | first pass | gain, hemispherical MLAs | gain, all non-planar | best random − hemisphere | wavy (1.5) − hemisphere |
| --- | --- | --- | --- | --- | --- |
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

**Table S10.** Spectrally weighted gain over the actual $\eta_{\mathrm{sta}}$ for the three placements of complete mixing (ranges over each structure class; number of structures in parentheses).

| Device | Structures | $\eta_{\mathrm{sta}}$ | Mixing at the escape interface (pp) | Mixing by the device (pp) | Transmissive mixing (pp) |
| --- | --- | --- | --- | --- | --- |
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

**Slab with a back reflector.** The slab of main-text Fig. 2(c) has refractive index $n$, an internal isotropic source ($\mathbf{v}^{(0)} \propto w$), and an angle-independent back reflectance $R$. Mixing at the exit face uses the kernel of main-text Eq. (5) with $m$ = 1; mixing at the back reflector uses $D$ with constant $R$, which gives the fully randomized result $p/[p+(1-p)(1-R)]$. The slab was also traced with a hemispherical MLA and two wavy surfaces as the exit face (Table S11). For the MLA and the steeper wavy surface, the slab without mixing comes within 2.9 pp of the fully randomized result at $n$ = 1.51 and within 1.1 pp at $n \ge 2.4$. The shallow wavy surface stays 2.8–5.8 pp below it; complete mixing at its exit face recovers none of this difference, whereas mixing at the back reflector recovers all of it.

**Table S11.** Slab with a textured exit face and a back reflector (isotropic internal source, specular back reflectance $R$).

| Exit face | $n$ | $R$ | $p$ | $p_{\mathrm{mix}}$ | No mixing | Mixing at exit face | Mixing at back reflector |
| --- | --- | --- | --- | --- | --- | --- | --- |
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

**Light trapping.** The face matrices of Table S11 were used as the front face of a weakly absorbing slab with a lossless back reflector. Light arrives from air within 10° of the normal, the power coupled into each internal angle follows from the transmission matrix by reciprocity, and a one-way crossing at internal angle $\theta$ absorbs $1-\exp(-\alpha d/\cos\theta)$. Four cases share the same front-face transmission: (A) the traced face and a specular back reflector, (B) complete mixing of the light reflected by the front face, (C) a Lambertian back reflector, and (D) both B and C. In all 12 cases the path enhancement of B is lower than or equal to that of C, and D is within 0.1 of C (Table S12). Near-normal light enters at angles with a high escape probability, and much of it leaves on its first return before mixing in reflection can act on it. Values slightly above $4n^2$ occur because the faces reflect part of the light inside the escape cone.

**Table S12.** Weakly absorbing slab with a textured front face and a lossless back reflector. Coupled: fraction of the incident power that enters the slab. A–D: path enhancement, absorptance/($\alpha d$), at $\alpha d = 10^{-3}$. Last column: absorptance at $\alpha d$ = 0.05 for B / C.

| $n$ | Front face | $p$ | $p_{\mathrm{mix}}$ | Coupled | A | B: mixing in reflection at the front | C: Lambertian back | D: B and C | $4n^2$ | Absorptance, $\alpha d$ = 0.05 (B / C) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
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

**Finite substrate and summary.** Fig. S3 repeats the comparison for a substrate of finite lateral size, and Table S13 summarizes the conditions under which the central result was tested.

**Fig. S3.** Finite lateral size (Al device, 550 nm; disk source of 1 mm radius, 1.295 mm glass). (a) $\eta_{\mathrm{sta}}$ against the half-width of the substrate without mixing (solid) and with complete mixing at the escape interface (dashed). (b) $\eta_{\mathrm{sta}}$ against the pixel radius when the back face outside the pixel is glass/air (solid) or absorbing (dash-dotted). (c) Gain from complete mixing at the escape interface. *(`figS_finite_size.png`)*

**Table S13.** Robustness of the central result. Gains are changes of $\eta_{\mathrm{sta}}$ in pp for complete mixing, for the 37 non-planar structures unless stated otherwise.

| Condition | Result | Source |
| --- | --- | --- |
| Al, Ag, and low-loss devices, 550 nm | Reflective mixing by the structure: −0.4 to +1.8 (Al), at most 2.0 (Ag) and 0.4 (low loss). Mixing at the device: +1.5 to +4.1 for lens arrays (Al) | Table S9, Note S4 |
| Spectral weighting, 400–700 nm | Reflective mixing: at most 1.1, 1.4, and 2.2 (Al, Ag, low loss). Mixing at the device: 1.1–3.3 for lens arrays (Al) | Table S10 |
| First pass not matched to $p(\theta)$, 500 nm | Reflective mixing: 5.7–6.4 (Al) and 8.4–9.4 (Ag) for the hemispherical MLA | Table S9 |
| Exit medium with $n$ = 1.5–3.5 (slab, $R$ = 0.90) | Textured faces within 1.1 of the fully randomized model at $n \ge 2.4$; mixing at the exit face at most 2.0 | Table S11 |
| Finite substrate (25 mm, absorbing edges) | Reflective mixing: −5.2 to +0.8. Transmissive mixing: −1.0 to +22.7 (six structures, three devices) | Fig. S3 |
| Partly diffuse electrode (diffuse fraction 0–1) | Gain nearly linear in the diffuse fraction, at most 2.3 | Note S4 |
| Numerical parameters, bin width, and reciprocity residual | $\eta_{\mathrm{sta}}$ changes by at most 0.34, $N_{\mathrm{mix}}$ by at most 0.24 | Note S1, Table S2 |
| Separate s and p channels | Complete-mixing gain changes by at most 0.4 | Note S1 |
| Direct $(\theta,\phi)$ ray tracing | Within 0.6 of the full series | Table S6 |
| Angular cutoff of the spectral quantities, 80–90° | Mean error of the reduced model changes by less than 0.3 | Table S1 |

---

## Note S5. Étendue bound and freeform benchmark

**Radiance-ordered bound.** In one encounter, the power that a passive structure can transmit into air is bounded by the power obtained when the étendue of the escape cone, $1/n^2$ in units of $w$, is filled with the brightest internal radiance first (main ref. 32). For uniform radiance the bound is $1/n^2$ = 0.439. With the first-pass radiance of the Al device at 550 nm it is 0.540. The MLAs reach 76–80% of this bound, wavy surfaces 77%, and the planar interface 65% [Fig. S4(a)]. The angular composition of the output is also nearly fixed: in the Al device, the share of the 40–60° band is 0.32–0.36 for all structures, close to the Lambertian value of 0.337 [Fig. S4(b)].

**Fig. S4.** Radiance bound and angular composition (Al device, 550 nm). (a) First-pass power into each polar band for the uniform-radiance bound, the radiance-ordered bound, and four structures. (b) Share of the 40–60° band in the multi-pass output against $\eta_{\mathrm{sta}}$ for all structures; the dotted line marks the Lambertian value. *(`fig5_bound_window.png`)*

**Freeform benchmark.** An axisymmetric freeform MLA with 13 variables (a seven-point spline profile, the height scale, and the thicknesses of the two transport layers) was compared with a hemispherical MLA with 3 variables on the Al platform, using LightTools ray tracing validated against fabricated devices (main refs. 8 and 9). The hemisphere is contained in the freeform design space. The total EQE and the EQE in four polar bands were optimized separately with a surrogate optimizer followed by a pattern search (10,000 rays, 31 wavelengths), and each optimum was re-evaluated three times at 50,000 rays and 151 wavelengths. The search was repeated from the hemispherical optimum and, for the total EQE, from three random starts. The band-optimized profiles differ visibly from the hemisphere, but the gains are small (Fig. S5, Table S14): the freeform exceeds the hemisphere by at most 7.0% in one band and by 1.3% in total EQE. The spread over independent searches is used as the uncertainty. Ratios below unity at 0–40° are search failures, because the hemisphere lies in the freeform design space. In a lens-only optimization at a fixed stack with the multi-pass model, the freeform gained 9.0% in the 40–60° band but only 1.25% in $\eta_{\mathrm{sta}}$; at unchanged $p$, it lowered $N_{\mathrm{mix}}$ from 7.1 to 2.1. The lens shape thus redistributes power among the bands but raises the total only through faster mixing.

**Table S14.** Freeform versus hemisphere (Al platform, LightTools). $G_j$: ratio of the best freeform EQE to the best hemispherical EQE. Last column: range of this ratio over independent searches, with the number of searches in parentheses.

| Objective | Hemisphere | Best freeform | $G_j$ | Independent searches |
| --- | --- | --- | --- | --- |
| 0–20° | 0.06682 | 0.06699 | 1.003 | 0.946–1.003 (2) |
| 20–40° | 0.16591 | 0.16628 | 1.002 | 0.966–1.002 (2) |
| 40–60° | 0.18503 | 0.19802 | 1.070 | 1.046–1.070 (2) |
| 60–80° | 0.13589 | 0.14075 | 1.036 | 1.020–1.036 (2) |
| Total EQE | 0.54679 | 0.5539 | 1.013 | 1.005–1.013 (4) |

**Fig. S5.** Freeform MLA against an equally optimized hemispherical MLA. (a) Band-optimized freeform profiles and the hemispherical reference. (b) Best EQE of each lens class for each objective, labeled with $G_j$. (c) Freeform/hemisphere ratio for each independent search. *(`fig6_hemisphere_benchmark.png`)*

---

## Note S6. Fabricated devices and the finite-geometry model

**Devices.** All devices share one substrate layout [main-text Fig. 5(a)]: 25.4 × 25.4 mm glass of 0.7 mm thickness, ITO on the upper half, five Al fingers, and a pixel-defining layer that opens a circle of 1 mm radius. The orange stack is given in the main text (emitter concentration 4 wt%, horizontal dipole ratio 0.865). The green devices use ITO (150 nm) / TAPC (105 nm) / TCTA (10 nm) / TCTA:B3PyMPM:Ir(ppy)$_2$acac (4 wt%, 25 nm) / B3PyMPM (42 nm) / LiF (1 nm) / Al (100 nm), with a horizontal dipole ratio of 0.77. The lens glass carries hexagonal close-packed lenslets of 10 µm radius on a 3 × 3 mm square. The freeform profile differs between the two emitters.

**Measurement.** The angle-resolved intensity was recorded with a photodiode goniometer (photodiode 59 mm from the rotation axis) in 10° steps with four readings per angle. The angular factor uses a natural cubic spline through the measured points and $I(90°) = 0$. For the orange devices (Table S15), the main text uses the maximum EQE at 3.4 V and above, which avoids the range near turn-on where the dark current of the photodiode contributes 3–5% of the signal. For the green devices (Table S16), several scans were taken with a tilted sample; the tilt was obtained from the asymmetry between positive and negative angles and corrected. Three data sets were excluded: device S4D3 (intensity maximum at −40°, EQE above the optical ceiling of the model), device G-S5D3 (intensity at 40–60° above that of the other freeform devices by up to 0.15), and one scan of the green hemispherical device with outliers. The green freeform devices show a second maximum at −70° that we attribute to stray light from the edges of the lens glass; 70–80° is excluded for these devices.

**Table S15.** Orange devices. EQE$_{\max}$: maximum above turn-on (voltage, luminance). EQE$_{\geq 3.4\,\mathrm{V}}$: maximum at 3.4 V and above (main text). Angular factor at EQE$_{\max}$.

| Device | Structure | EQE$_{\max}$ (%) | at V (cd/m²) | EQE$_{\geq 3.4\,\mathrm{V}}$ (%) | Angular factor |
| --- | --- | --- | --- | --- | --- |
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

**Table S16.** Green devices. Tilt: shift fitted to the asymmetry between positive and negative angles. Angular factor per scan (freeform without 70–80°). EQE$_{0°}$: Lambertian-equivalent on-axis EQE. EQE: EQE$_{0°}$ times the mean angular factor.

| Device | Structure | Scans used | Tilt (°) | Angular factor | EQE$_{0°}$ (%) | EQE (%) |
| --- | --- | --- | --- | --- | --- | --- |
| G-S1D4 | bare | 4 | −13.3 to 1.3 | 0.99–1.03 | 22.4 | 22.7 |
| G-S2D2 | bare | 3 | −1.6 to 0.6 | 1.01–1.03 | 21.9 | 22.3 |
| G-S3D2 | hemispherical | 6 (of 7) | −13.4 to −8.2 | 0.69–0.73 | 26.9 | 19.1 |
| G-S1D3 | freeform | 4 | −0.9 to −0.1 | 0.82–0.83 | 21.1 | 17.4 |
| G-S2D3 | freeform | 5 | −1.0 to 1.0 | 0.85–0.90 | 22.4 | 19.7 |
| G-S3D3 | freeform | 3 | −2.0 to −1.2 | 0.89–0.90 | 22.7 | 20.3 |

**Model.** The device optics follow main-text Sec. 2.4. The back reflectance seen from the glass is computed for the pixel, for the Al-covered regions, and for the uncovered regions, each with or without ITO. Inside the lens patch, each encounter with the front face uses the traced return and transmission kernels of the lens array; outside it, the planar face acts specularly. Between encounters, the light moves laterally by $2t\tan\theta$ with $t$ = 1.4 mm (0.7 mm for bare devices), and light that reaches the blackened sides is lost. With mirror sides, the model reproduces $\eta_{\mathrm{sta}}$ of the infinite substrate to within $2\times10^{-4}$. Because the hexagonal lattice is not axially symmetric, the lens kernels are traced with the azimuth resolved relative to the lattice, and the goniometer is modeled with its scan plane, rotation axis, and photodiode size (5 or 10 mm). The lattice orientation changes the angular factor of the hemispherical devices by up to 0.05, the rotation axis by at most 0.03, and the photodiode size by at most 0.02; the ranges in main-text Fig. 5 cover all these geometries.

**Calibration on the bare devices.** One internal efficiency per emitter (0.89 for orange and 0.95 for green) is fixed on the bare devices. For the green stack, the nominal thicknesses give a bare angular factor of 0.88–0.90, against 0.99–1.03 measured. A fit of the B3PyMPM and TAPC thicknesses to the bare profile alone gives 65–70 nm for B3PyMPM with a weakly constrained TAPC thickness, and the model uses 65 and 80 nm.

**Loss in the attached lens glass.** The measured EQE of the lens devices is 84–91% of the lossless prediction, whereas the angular factors agree. Three candidate losses were scaled to lower the EQE of the orange hemispherical devices by 14%: a loss at the lens surface, a transmission of 0.925 per crossing of the oil-coupled interface, and bulk attenuation in the lens glass. All three change the angular factor by at most 0.02, so the measurement cannot distinguish them; the bulk coefficient would be about two orders of magnitude above that of display glass, which makes the interface or the lens film the more likely origin. With the interface transmission fixed on the orange hemispherical devices, the model gives 0.89, 0.87, and 0.90 of the lossless EQE for the orange freeform, green hemispherical, and green freeform devices, against 0.90, 0.84, and 0.91 measured. Variations of the lens shape do not explain the deficit, because those that lower the on-axis gain to the measured value raise the angular factor above its measured range (Table S17).

**Table S17.** Candidate causes of the EQE deficit of the lens devices. Changes are relative to the lossless model of the 3 × 3 mm patch (orange devices unless stated). Caps: spherical caps with height/base-radius ratio $a$. Gaps: hemispheres of base radius $r_b$, in units of the cell radius, with flat glass between them.

| Candidate | Lens-device EQE | Angular factor | Consistent with the measurement |
| --- | --- | --- | --- |
| cavity of the bare devices | fixed by the bare devices; common to bare and lens devices | — | not a cause of a lens/bare difference |
| lens shape: caps, $a$ = 0.7 and 0.5 | −3% and −7% | +0.04 and +0.13 | no: angular factor above the measured 0.63–0.71 |
| lens shape: gaps, $r_b$ = 0.95–0.7 | 0 to −2% | +0.03 to +0.20 | no: angular factor and a maximum near 10° |
| goniometer geometry (lattice, scan plane, axis, photodiode) | no change of the EQE | ±0.03 | included in the predicted ranges |
| patch offset of 0.5 mm | +1% (on-axis EQE −5%) | +0.05; profile rms 0.067 instead of 0.055 | no: raises the EQE |
| loss at the lens surface, interface or in the lens glass | −14% (fitted), −11% to −13% for the other groups | ≤ 0.02 | yes: one parameter reproduces all four lens-device groups within 0.03 |

**Patch width.** Table S18 gives the simulated EQE against the width of the lens patch. A 3 mm hemispherical patch lowers the EQE below that of the bare device, a patch of about 6 mm gives half the gain of a 20 mm patch, and 10 mm gives 80% of it.

**Table S18.** Simulated EQE and angular factor against the width of the lens patch (orange device layout, without the internal-efficiency scaling and without the lens-glass loss). Bare device: 33.7%, angular factor 1.00.

| Patch width (mm) | Hemispherical EQE (%) | Angular factor | Freeform EQE (%) | Angular factor |
| --- | --- | --- | --- | --- |
| 2 | 30.9 | 0.87 | 31.0 | 1.40 |
| 3 | 31.7 | 0.66 | 31.8 | 0.99 |
| 4 | 35.6 | 0.70 | 35.0 | 0.91 |
| 5 | 39.3 | 0.76 | 37.7 | 0.93 |
| 7 | 42.8 | 0.81 | 40.4 | 0.96 |
| 10 | 45.4 | 0.85 | 43.1 | 0.98 |
| 14 | 47.1 | 0.88 | 45.0 | 1.01 |
| 20 | 48.5 | 0.90 | 46.6 | 1.03 |
