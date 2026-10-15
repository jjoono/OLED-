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

Electrically plausible candidates (`cands.py` -> `cands.csv`, TFE + n1.65 MLA, 1st cavity order, d in nm):
A Al/Liq1/Bebq2 84/EML25/TAPC 55/HATCN50/Ag15 (Batch10 replica) 0.407
B Ag/Al0.2/Liq1/Bebq2 88/EML25/TAPC 38/HATCN50/Ag15 0.504
C Ag/Al0.2/Liq1/Bebq2 97/EML25/TAPC 76/HATCN10/Ag12 0.544
D Ag/Liq1/Bebq2 97/EML25/TAPC 76/HATCN10/Ag12 0.558
E Ag/Al0.2/Liq1/Bebq2:Ir 107 (RZ 5 nm from NPB)/NPB 81/HATCN10/Ag12 0.536
F Ag/Liq1/Bebq2 101/EML25/TAPC 71/HATCN5/Ag8 0.605
+-10 nm in d_ETL or d_HTL costs <= 0.6 %p; 2nd order (d_ETL ~240-270) is 4-8 %p lower in every case.

Push to >= 60 % (`push60.py` -> `push60.csv`; 'Al2O3->lens' = ALD Al2O3 30 nm, lens printed on it, ParyC dyads over the lens;
'film 50%' = top-Ag nk halfway between the 260703 #2 film and McPeak; 'test5' = library l_Ag_12nm_test5):
C Ag12 meas full TFE 0.544 | C1 Al2O3->lens 0.559 | C2 +film 50% 0.591 | Ag12 test5 full TFE 0.596, Al2O3->lens 0.616
G Ag8 meas full TFE 0.591 | G1 Al2O3->lens 0.610 | G2 film 50% 0.618 | G3 both 0.639 | Ag8 test5 0.635 / 0.659
H Ag/Liq, Ag8 meas, Al2O3->lens 0.625 | I Ag10 meas Al2O3->lens 0.584 | J Ag10 film 50% Al2O3->lens 0.615
