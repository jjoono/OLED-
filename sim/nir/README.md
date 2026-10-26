# NIR (760 nm) bottom-emitting OLED: transparent electrode x DBR x substrate x MLA global optimisation

Stack: Ag 100 | ETL | EML 20 nm (n 1.75, dipole at centre, isotropic Θ = 2/3) | HTL | TE | (seed 5 nm for Ag/Au) | DBR (SiO2/TiO2 x 4, chirped, free thicknesses) | substrate n 1.5 / 1.8 | hemispherical MLA (BSDF `mla/lt_hemisphere_bsdf.mat`).
Emitter: Gaussian 760 nm, FWHM 80 nm. TE: ITO (measured, d ≥ 44 nm), IZO (measured, d ≥ 80 nm), IOH = H:In2O3 Drude model
(N 1.5e20 cm^-3, μ 100 cm²/Vs, d ≥ 84 nm), i.e. Rs ≤ 50 Ω/sq; Ag (HATCN-seeded, 8–12 nm), Au (8–15 nm).
Objective: EQE at PLQY q = 0.3 (Purcell-enhanced: EQE = q F η_sub η_ext / (1 − q + q F)); DE + Nelder–Mead, then cross-seeded refinement (`refine.py`).
Other q are recomputed from F and EQE(q=1) with the q = 0.3 structure (`summarize.py` -> `summary.csv/md`).

Two optimum families appear:
* 1st-order cavity (ETL ≈ 87–100 nm): F 1.4–1.9, SPP loss 25–34 %. Wins on n = 1.5 and for metal TEs.
* 2nd-order cavity (ETL ≈ 330 nm): SPP 4–6 %, F ≈ 1, η_sub ≈ 0.9–0.95. Only on n = 1.8 (≥ n_org, no waveguide), wins for TCOs.

Best at q = 0.3: IOH / n 1.8 / DBR4 27.0 % (q = 1: 87 %); Ag 12 nm / n 1.5 / DBR4 25.4 %; flat-glass ITO reference 8.4 %.
Caveats: isotropic dipoles, unity charge balance/spin factor, laterally infinite substrate, IOH optical model hypothetical,
structures optimised for q = 0.3 only, Ag/Au thicknesses at bound in several cases.
