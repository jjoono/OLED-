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

ZnS capping (l_ZnS, n 2.36 @550) between top Ag and TFE, d_E/d_H re-optimised (`zns.py` -> `zns.csv`):
meas Ag 8: ZnS 0/20/40/60/80 = 0.605/0.606/0.599/0.603/0.603; meas Ag 6: 0.631 -> 0.603-0.617; bulk Ag 8: 0.670 -> 0.683 (20 nm).
At 615 nm (Ag 8 meas): ZnS 20 raises emission into the lens 0.746 -> 0.771 but also the round-trip loss A' 9.2 -> 11.2 %;
ZnS 60 adds a guided mode (wg 0.04). Net effect within +-1.5 %p.
