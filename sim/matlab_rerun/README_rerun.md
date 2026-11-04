# Re-running Extended Data Fig. 7d and Fig. 8 with MATLAB + LightTools

Files
- `planar_inputs.m`: dipole/TMM part (same equations as `Planar_sweep22_preprint_MLA.m`), returns w_sub, Psub_norm, R_OLED.
- `multipass_series.m`: round-trip series, summed to convergence (last term < 1e-6 of the sum).
- `check_against_table1.m`: reproduces `../table1_author/table1_author_model.csv` (20 terms) as a check.
- `ED7d_etl_absorption.m`: ED Fig. 7d (no LightTools needed; uses BSDF slice 11, n_MLA = 1.80). Writes `ED7d_etl_absorption.csv`.
- `ED8_eqe_from_bsdf.m` + `ED8_lighttools_geometries.csv`: ED Fig. 8. Needs one LightTools BSDF per row of the csv in `./bsdf/<id>.mat` (variable `BSDF`).
- `flat_bsdf.m`: Fresnel BSDF of a flat interface for the planar reference.
- `release_minimal/`: the only code proposed for public release (series + demo input + README + LICENSE).

LightTools set-up for ED8 (one run per row of `ED8_lighttools_geometries.csv`)
- Same procedure and output format as `hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat` (slice 11).
- Lens and substrate index 1.80, no absorption; hexagonal lattice with the listed pitch (centre-to-centre).
- `spherical_cap`: base radius r, height h, sphere radius (r^2+h^2)/2h, centre below the base plane.
- `half_ellipsoid` (AR > 1): (x^2+y^2)/r^2 + z^2/h^2 = 1.
- Only the ratios matter in ray optics; r = 10 um is used in the table.
- Incident angles 0-89 deg (1-deg steps) from the substrate side, same ray count per angle as the original BSDF.
- Store as 180 x 90: rows 1-90 transmitted, rows 91-180 reflected (row 180 = 0 deg).
- Check: the row `AR1p00` must reproduce slice 11 of the original file (EQE_Ag 0.838, EQE_Al 0.553).
