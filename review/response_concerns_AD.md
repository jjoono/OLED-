# Response to the four concerns (A–D) before submission

New computations: `selectivity/reciprocal_limit.py` (A), `multipass_reset/direct_blind.py` + `blind_table.py` (B),
`selectivity/cavity_decomp.py`, `multipass_reset/fixedcav_opt.py`, `opt_fixedcavity_freeform.m` (C).
Data: `selectivity/data/reciprocal_limit.json`, `multipass_reset/blind/results.json`, `selectivity/data/cavity_decomp.json`,
`multipass_reset/fixedcav/results.json`.

## A. New physical result or an application of known mathematics?

**Change.** The paper is now built around the question "how much can better angular mixing raise outcoupling?"
(new title, abstract, introduction, Sec. 2.4, Sec. 4.3, Fig. 3d, Discussion, Conclusion).

**New quantitative result.** For the only reciprocal rank-one return kernel, completely mixed light escapes with
p_mix = p − σ_p²/(1 − p) (σ_p²: étendue-weighted variance of p(θ)), not with p as the memoryless formula assumes;
the multi-pass series then has a closed form (Sec. 2.4, Note S1; verified against the numerical series to 1e-9).
Consequences:
- For the 37 non-planar structures, complete reciprocal mixing at a structure's own p(θ) would add −0.4 to +1.8 pp
  (Al), ≤ 2.0 pp (Ag), ≤ 0.4 pp (low loss). Hemispherical MLAs, the slowest mixers, would gain only 1.3–1.7 pp because
  p_mix = (0.87–0.91) p.
- The memoryless formula predicts gains larger by 0.6–2.9 pp (lens arrays), 17 pp (shallow wavy), 31 pp (planar,
  p_mix = 0.06 against p = 0.40).
- Even a complete mixer has ξ = 1 − σ_p²/[p(1 − p)] < 1 unless p(θ) is flat: the design lever is a higher p or a flatter
  p(θ), not faster mixing.

The Markov/reciprocity machinery is still stated as known; the novelty claim is restricted to these consequences.

## B. No independent validation

- **Direct ray tracing** of the whole multi-pass process in (θ, φ), azimuth kept between passes, specular R(θ)
  (1e6 rays per case): {{DIRECT}}
- **Blind test** on nine families not used in model development: {{BLIND}}
- **Not possible here:** experimental BSDF measurement, and an independent commercial ray tracer for structures other
  than the hemisphere. Recommended if the target is Nature Photonics / LSA: measure the return BSDF of two or three
  fabricated films (hemispherical MLA, random MLA, scattering film) and compute N_mix, ξ, p_mix from it.

## C. Cavity co-optimization in the freeform benchmark

- Full device model, spectrally weighted, own cavities (Note S8, Table S13): lens shape alone changes the band EQE by
  −5.4% to +0.4%, the cavity by −2.6% to +0.3%.
- The model reproduces the hemispheres' LightTools EQE within 4% but rates the freeform optima 1–5% lower relative to
  them; tools agree that shape gains are a few percent at most, not on their sign below ~5%.
- Table 3's 0.771 vs 0.792: the total-EQE freeform lens at the hemisphere's cavity extracts less; its LightTools
  advantage rests on its own cavity (ETL 67.5 / HTL 120.6 nm vs 61.3 / 97.5) and on an annular-trough profile.
- Lens-only optimization at a fixed stack in the model: {{ES}}
- `opt_fixedcavity_freeform.m`: LightTools control with ETL/HTL fixed at each arm's hemispherical optimum and the
  hemisphere as first seed. **To be run by the authors** (needs LightTools/MATLAB).

## D. Statistics

- Main text and Table 2 now judge reproducibility by independent searches: 40–60° 1.046–1.070, 60–80° 1.020–1.036,
  total 1.005–1.013 (4 searches), 0–40° never above 1.003; the hemisphere had more evaluations per variable
  (15 vs 4–11).
- The t values are described as tests of evaluation noise only (Note S6), not of optimizer uncertainty.
- Remaining recommendation: 2–3 further independent restarts per arm (e.g. `opt_fixedcavity_freeform.m` with different
  `RNG_SEED`) would give a proper spread.
