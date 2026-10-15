# TEOLED optimisation: Ir(mphmq)2tmd, thin Ag on HATCN, Parylene-C TFE, n 1.65 MLA

Ag 100 | Liq 1 | Bebq2 d_E (l_bebq2_JH) | Bebq2:Ir(mphmq)2tmd 25 | TAPC d_H | HATCN | Ag d_Ag |
Al2O3 30 / ParyC 500 / Al2O3 30 / ParyC 500 / Al2O3 30 | lens 1.65 | hemisphere MLA (ray-traced BSDF).
CPS emission + full recycling series, 6 wavelengths (580-680 nm) weighted by l_I_Irmphmq2tmd_measure_JH,
Theta 0.82, eta_rad 0.96 (Purcell change neglected). `grid.py` coarse grid -> `grid.csv`; `refine.py` Nelder-Mead -> `refine.csv`.
Top Ag nk: 'meas' = l_Ag_15_onHATCN applied at every thickness, 'bulk' = McPeak (ideal film).

Best (EQE): Ag 6 / 8 / 10 nm = 0.631 / 0.605 / 0.581 (meas nk), 0.689 / 0.670 / 0.650 (bulk nk); d_E ~100-115, d_H ~50-85.
HATCN 5-20 nm changes EQE < 1 %p. Bottom: Ag/Liq 0.605, Ag/Al0.2/Liq 0.591, Al/Liq 0.485. Without TFE 0.626.
Index study (Ag 8, meas): low-index ETL (1.65) +5 %p; low-index EML alone -1 %p; with n_ETL <= n_lens the
2nd-order thick ETL (~290-310 nm) becomes as good or better (n 1.55: 0.687).
