# Re-running Extended Data Fig. 7d and Fig. 8 with MATLAB + LightTools

Files
- `planar_inputs.m`: dipole/TMM part (same equations as `Planar_sweep22_preprint_MLA.m`), returns w_sub, Psub_norm, R_OLED.
- `multipass_series.m`: round-trip series, summed to convergence (last term < 1e-6 of the sum).
- `check_against_table1.m`: reproduces `../table1_author/table1_author_model.csv` (20 terms) as a check.
- `ED7d_etl_absorption.m`: ED Fig. 7d (no LightTools needed; uses BSDF slice 11, n_MLA = 1.80). Writes `ED7d_etl_absorption.csv`.
- `ED8_eqe_from_bsdf.m` + `ED8_lighttools_geometries.csv`: ED Fig. 8. Needs one LightTools BSDF per row of the csv in `./bsdf/<id>.mat` (variable `BSDF`).
- `flat_bsdf.m`: Fresnel BSDF of a flat interface for the planar reference.
- `release_minimal/`: the only code proposed for public release (series + demo input + README + LICENSE).

LightTools set-up for ED8: `LT_ED8_bsdf_sweep.m` (user's COM script with the ED8 parameter lists)
- Lens/substrate index 1.80; hexagonal texture zone of hemispheres in the existing .lts.
- Aspect ratio = StretchZ (height / base radius), so every lens is a half-spheroid (a hemisphere at AR = 1).
- Fill factor: Scale = sqrt(FF/0.9069), assuming Scale shrinks the lens at fixed pitch.
- Writes `bsdf/<id>.mat` (variable `BSDF`, 180 x 90, columns normalized), then run `ED8_eqe_from_bsdf.m`.
- Check: AR1p00 must match slice 11 of the original BSDF; `ED8_eqe_from_bsdf.m` then gives EQE_Ag 0.838, EQE_Al 0.553.
