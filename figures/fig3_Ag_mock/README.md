# Fig. 3 mock-ups with an Ag reflector (generic model)

All panels: 550 nm, isotropic dipole, PLQY = 1, `Planar_sweep22_preprint_v2.m` bookkeeping
(five channels air + sub-confined + WG + SPP + abs = 1) with the branch-fixed `TMF_birefringence_whole*.m`.
Reflector: Ag, McPeak n,k (`material.l_Ag_McPeak(151)` = 0.044 + 3.819i), 100 nm.
Octave 8.4 was used; the scripts also run in MATLAB (they need `nk_JH_total.mat` and the three TMF files in the path).

| panel | sweep | Octave script | data | plot |
|---|---|---|---|---|
| (a) η_sub vs organic thickness for n_sub = 1.5 / 1.65 / 1.8 | d_org = 10…490 nm, n_sub | `fig3a_Ag.m` | `fig3a_Ag.csv` (n_sub, d_org, air, sub, wg, spp, abs) | `plot_fig3a_Ag.py` → `fig3a_Ag_mock.png` |
| (b) SPP vs ETL thickness for n_e,ETL = 1.8 / 1.7 / 1.6 / 1.5 | d_ETL = 5…400 nm, n_e,ETL | `fig3b_Ag.m` | `fig3b_Ag.csv` (n_e, d_ETL, …) | `plot_fig3b_Ag.py` → `fig3b_Ag_mock.png` |
| (c) 2-D map η_sub(d_ETL, n_sub) for n_e,ETL = 1.8 and 1.5 | d_ETL = 5…400 nm × n_sub = 1.40…2.00 | `fig3c_Ag.m` (`octave --eval "NE_ETL=1.8; OUTCSV='fig3c_Ag_ne18.csv'; fig3c_Ag"`) | `fig3c_Ag_ne18.csv`, `fig3c_Ag_ne15.csv` (n_sub, d_ETL, …) | `plot_fig3c_Ag.py` → `fig3c_Ag_mock.png` |

Stack for (b) and (c): Ag / ETL (n_o = 1.8, n_e varied, d_ETL) / EML 1.8 (20 nm, dipole at the centre) /
HTL 1.8 (50 nm) / TCO 1.8 + 0.02i (50 nm, index-matched so that no TCO-guided light lands in the u > 1 bin) / substrate.
Stack for (a): Ag / organic 1.8 (d_org, dipole 10 nm from the Ag side of the EML) / TCO 2 + 0.02i / substrate n_sub.

Key numbers (PLQY = 1):

* (a) η_sub at d_org = 490 nm: 0.449 (n_sub 1.5), 0.617 (1.65), 0.895 (1.8); absorption 2–5 %.
* (b) the u > 1 (SPP) power drops below 5 % at d_ETL ≈ 200 nm for n_e,ETL = 1.8/1.7/1.6 but at ≈ 80 nm for 1.5:
  the anisotropic SPP index n_SPP(n_e) = 2.04/1.92/1.80/1.69 (Ag, n_o = 1.8) falls below the organic index only for n_e = 1.5,
  so the SPP becomes leaky and its power re-appears as substrate light (a threshold, not a gradual barrier).
  For Al the same crossing happens already at n_e ≈ 1.7 (n_SPP = 1.86/1.76/1.66/1.55).
* (c) the map factorises into a horizontal transition at n_sub = n_EML (waveguide cut-off) and a vertical transition at the
  SPP threshold thickness; the top-right corner gives η_sub ≈ 0.91–0.95, limited only by absorption (TCO k = 0.02, Ag).
