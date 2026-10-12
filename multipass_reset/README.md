# Angle-preserving vs Lambertian-reset rear (draft analysis)

Question: if the light returned to the device were re-randomized to a Lambertian
distribution on every pass (e.g. Ag / thin scatterer / ITO rear), how much do the
hemispherical and freeform MLAs of this paper gain, and does lens shape regain leverage?

- `hftrace.py` — numba Monte-Carlo BSDF of a hex close-packed array of axisymmetric
  lenslets with arbitrary profile h(rho), index-matched (n = 1.51). Validated against the
  LightTools hemisphere BSDF of the preprint branch: p = 0.4275 vs 0.4274, B_T(theta)
  within 0.01.
- `run_bsdfs.py`, `run_free.py` — BSDFs of the optimized hemisphere (5) and freeform (5)
  profiles from `opt_hemisphere_result.mat`, `opt_4band_result_25by25.mat`,
  `freeform_EQEtotal_result.mat` (`bsdfs.npz`). Profiles are rebuilt with PCHIP; the
  LightTools curve uses its own spline, so freeform BSDFs approximate the LightTools geometry.
- `multipass.py` — matrix multi-pass series (preprint Supplementary Note 1) with either an
  angle-preserving rear (diag R) or a Lambertian-reset rear; device from the preprint CPS
  code (`cps2`, branch `claude/oled-efficiency-paper-ag-955zsk`, `sim/fig3c`), 550 nm,
  Al or Ag cathode, ITO 150 nm, Theta = 0.865. Outputs `multipass_results.*` (own cavity),
  `mixing_sweep.csv` (partial randomization m), `common_cavity.csv` (all lenses at the
  hemisphere total-EQE cavity, ETL 61.3 / HTL 97.5 nm).
- `fig_reset.py` — draft Fig. S8 + Excel data (`figS8_angle_reset_draft.*`).

Excluded: the 0-20 deg freeform optimum ends in a vertical rim wall (two control points at
r = 1), which a height-field tracer cannot represent.
