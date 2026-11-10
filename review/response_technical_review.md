# Response to the technical review (Nature Research Assistant, 10 Nov 2026, 24 suggestions)

Each item lists what was changed and where. "Main" refers to `manuscript_draft.md`, "Supp." to `supplementary.md`. New computations are in `selectivity/review_checks.py`, `multipass_reset/run_review_traces.py`, `selectivity/review_traces_analysis.py` (data in `selectivity/data/review_checks.json`, `review_traces.json`).

| # | Suggestion | Action | Where |
|---|---|---|---|
| 1 | State-vector and binning definitions | Bin edges θ_k = k°, v_k = power per bin over the full azimuth, uniform radiance ↔ v ∝ w, p = Σw_k p_k/Σw_k; 0.5° slab matrices coarsened to 1° | Main 2.1 |
| 2 | Escape vs conditional return; meaning of "Markov chain" | B_R = M diag(1−p); M = column-stochastic chain conditioned on return; escape carried by the survival factor, not an absorbing state; stationary distribution ∝ w(1−p) | Main 2.1, 2.2 |
| 3 | A′ versus R(θ) | A′ defined as uniform-radiance average of 1−R, used only in the memoryless formula; all losses of the full model enter via R(θ); Fig. 1a label and caption changed | Main 2.1, Fig. 1 |
| 4 | Polarization reduction | R = (|r_s|²+|r_p|²)/2; v⁽⁰⁾ = s + p; tracers use unpolarized Fresnel coefficients at every encounter | Main 2.1, Supp. Note S12 |
| 5 | Azimuth handling | Rays launched at bin-centre angle, uniformly over unit cell/supercell and azimuth; exit binned by polar angle only; B_R is the position- and azimuth-averaged kernel | Main 2.1 |
| 6 | θ < 85° procedure; consistency of state spaces | Sub-block of B_R below 85°; M by column renormalization of the sub-block; ρ from the sub-block; η always on the full 0–90° grid; v⁽⁰⁾ above 85° < 0.1% | Main 2.2 |
| 7 | N_mix definition; negative λ₂ | N_mix = 1/(1−\|λ₂\|) everywhere (SLEM); all λ₂ real and positive (0.09–0.88) | Abstract, Main 1, 2.2, 6 |
| 8 | Robustness to MC asymmetry | Symmetrized w_iB_R(i→j) with column sums kept: ΔN_mix ≤ 2.7%, Δξ ≤ 0.005, Δη ≤ 0.04 pp, Δ(reduced error) ≤ 0.04 pp | Main 2.2, Supp. Note S12 |
| 9 | Convergence: bin size, cutoff, ray count | 0.5°/1°/2°/3° bins, cutoff 80/85/88/90°, 2.5k–40k rays: η unchanged; N_mix, ξ absolute values cutoff-dependent, rankings preserved (rank corr. 0.96–0.99 N_mix, 0.90–0.98 ξ) | Main 2.2, 3.2, 4.2; Supp. Table S15 |
| 10 | Discrete p and numerical check of the bounds | p̃ defined explicitly (escape plus leakage beyond cutoff, ≤ p + 0.008); bounds hold for all 37 non-planar structures (margins 0.0014 upper, 0.106 lower); with p instead of p̃ one thick scatterer would fail by 0.005 | Main 2.4, Supp. Note S12 |
| 11 | Device-aware late-pass metric | ξ framed as structure-only; 1−ρ(diag R B_R) ranks structures like p∞ (Spearman 0.956, Al 550 nm) | Main 2.4, Supp. Note S12 |
| 12 | Scaling conventions in Figs. 1 and 3 | Fig. 1: shared colour bar, log10 M with floor 10⁻³; Fig. 3a,c: log axes labelled; planar point marked off-scale | Figs. 1, 3, captions |
| 13 | Fig. 2 encoding; stratified errors | Caption states colour = class; new panel (e) by device × wavelength; stratified table (Al 0.41/2.7, Ag 0.33/3.0, low loss 0.16/1.8 pp; 500/550/600 nm 0.38/0.29/0.23 pp mean) | Fig. 2e, Main 4.1, Supp. Table S11 |
| 14 | Universality of Fig. 3a contours | Contours recomputed with scattering-layer and wavy p(θ) templates: N_mix 1→16 changes η by 1.2–3.2 pp vs 0.6 pp per 0.01 in p for all templates | Main 4.2, Supp. Table S16 |
| 15 | Disentangle mixing from p(θ) reshaping | Counterfactual decomposition along the shape sweep: s = 1, Al: mixing +0.91, p(θ) −2.20, interaction +0.30 pp (full series −0.91) | Main 4.2, Supp. Table S17 |
| 16 | Finite supercell / periodicity | 10 × 10 supercells and 24 × 24 wavy boxes (new realizations) compared with the library | Main 3.2, Supp. Table S15 |
| 17 | Band-bound arithmetic; algorithm | Per-band bounds are separate maxima (sum 0.676 > joint 0.540); uniform bound additive incl. 0.014 above 80°; discrete algorithm (radiance sort, Snell-mapped étendue, fractional last bin, lossless, azimuth-integrated) | Main 2.5, 4.3 |
| 18 | Definitions of first-pass power and selectivity; robustness | Definitions in caption and text; window per device/wavelength: Al 0.320–0.364, Ag 0.322–0.361, low loss 0.274–0.366 (hemispherical MLAs); claim restricted to lossy devices | Main 4.3, 6, Fig. 4 caption, Supp. Table S18 |
| 19 | Device inputs | Full layer list incl. ETL 61.3 nm/HTL 97.5 nm, dipole sheet at EML centre, A′ at 550 nm; R(θ) and v⁽⁰⁾(θ) plotted | Main 3.1, Supp. Table S14, Fig. S11 |
| 20 | Full-series numerics | Direct linear solve on the full grid; max ρ(diag R B_R) = 0.9991; residual < 1e-16; pass sum agrees to 1e-14 | Main 2.1, Supp. Note S12 |
| 21 | Fig. 3b axis | Each randomness parameter normalized by its maximum, ranges in the legend; hemispherical zero point added to jitter and height series | Fig. 3b, caption |
| 22 | EQE definition and statistics in the benchmark | EQE = η_rad × spectrally weighted outcoupling (151 λ, 453–753 nm) over hemisphere/band; search vs re-evaluation fidelity; n = 3, pooled t, df = 4, one-sided threshold 2.13; Holm | Main 3.3 |
| 23 | Control with the stack fixed | Simplified multi-pass model at a common cavity (Note S9): spread among nine lenses 2.5% (Al), 1.7% (Ag); band steering mostly from the cavity. A full-fidelity LightTools fixed-stack control is stated as not yet done | Main 3.3, 4.4 |
| 24 | Benchmark designs on the map | New Table 3: p, N_mix, ξ, η of the nine optima (p 0.417–0.429, N_mix 3.0–8.0, ξ 0.39–0.79); they are the MLA markers of Fig. 3 | Main 4.4, Table 3, Fig. 3 caption |

## Open items
- Item 23: the LightTools fixed-stack control at full fidelity requires the LightTools/MATLAB pipeline and is not run here.
