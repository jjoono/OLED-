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

2026-08 thin-Ag-on-HATCN nk (`joint_run.py` -> `joint_run.csv`; HATCN 5, top Ag 7/8 nm):
TRinv = T/R point inversion (nk_joint/*_TRinv.csv); jointA = TRinv rescaled to the SE+T/R joint n,k at 550 nm
(README_260910_joint.md: 7 nm 0.184+3.360i, 8 nm 0.153+3.450i) until the joint CSVs are available.
Ag 8 jointA: Ag/Al0.2 full TFE 0.584, Al2O3->lens 0.602; Ag/Liq 0.598 / 0.618.  Ag 7 jointA: 0.586 / 0.604; 0.600 / 0.619.
Ag 8 TRinv: 0.539 / 0.554; 0.552 / 0.567.

Joint SE+T/R nk (thin_Ag_nk_260910, `nk_joint/Ag_HATCN_*_nk.csv`, 'joint' columns; Ag at the SE-fitted thickness;
`joint_true.py` -> `joint_true.csv`). EQE, d_ETL ~100-105, d_HTL 57-82:
| Ag nom (fit) | Ag/Al0.2 full TFE | Ag/Al0.2 Al2O3->lens | Ag/Liq full TFE | Ag/Liq Al2O3->lens |
| 6 (6.18)  | 0.581 | 0.599 | 0.594 | 0.613 |
| 7 (6.84)  | 0.592 | 0.610 | 0.605 | 0.625 |
| 8 (7.87)  | 0.588 | 0.606 | 0.602 | 0.622 |
| 10 (8.97) | 0.581 | 0.599 | 0.596 | 0.615 |
| 12 (10.71)| 0.574 | 0.592 | 0.590 | 0.610 |

Bottom mirror IZO on the chirped ZnS/LiF DBR (`dbr_bottom.py` -> `dbr_bottom.csv`, top Ag 8 joint):
glass behind 0.26-0.45 (lens-side loss+leak 27.6 % vs 9.1 % for Ag at 615 nm); with IZO/Al1.5 0.34;
Ag behind the DBR 0.575 / 0.594 (DBR guided mode 0.11) vs Ag/Liq 0.602 / 0.622.

High-index 1.5-dyad TFE (`hi175.py` -> `hi175.csv`): Ag top | Al2O3 30 | n1.75 layer 3 um (incoherent, lossless) | Al2O3 30 |
n1.75 hemisphere MLA (BSDF `mla_hemi_1.75.npz`, p_Lamb 0.307). Joint-fit top-Ag nk, d_ETL 129-138, d_HTL 69-92:
Ag 7/8/10/12 nominal: Ag/Al0.2 0.627/0.624/0.618/0.614, Ag/Liq 0.643/0.641/0.636/0.634. Upper Al2O3 30 nm costs < 0.1 %p.

ZnS >= 100 nm on top Ag under the high-index TFE (`zns_hi2.py` -> `zns_hi2.csv`, 27 wavelengths 560-716 nm;
`zns_hi.py` used 6 wavelengths and overrates resonant optima by up to 8 %p - do not use zns_hi_opt.csv / zns_hi_map.csv).
Ag 8 joint, Ag/Al0.2/Liq, HATCN 10: best ZnS ~100 nm with a thick HTL:
d_ETL<=80: 0.608 | <=100: 0.630 (HTL 252) | <=130: 0.654 (HTL 254) | <=160: 0.663 | <=200: 0.667 (ETL 195, HTL 135).
References at ETL 130: no ZnS/HTL 78 0.622; ZnS 100/HTL 78 0.612; no ZnS/HTL 254 0.597.
Tolerance at (130, 254, ZnS 100): +-10 nm HTL or ZnS <= 0.5 %p, +-10 nm ETL -0.6/+0.4 %p.

ZnS 100-250 nm with thin organics (`zns_thin.py` -> `zns_thin_grid.csv`, 27 wavelengths, 10 nm grid, Ag 8 joint, Ag/Al0.2/Liq):
best (ETL<=130, HTL<=150) / (ETL<=100, HTL<=120) / (ETL<=80, HTL<=100):
ZnS 100: 0.634 (130/150) / 0.614 (100/100) / 0.603 | 125: 0.619 / 0.602 / 0.588 | 150: 0.603 / 0.593 / 0.580 |
175: 0.599 / 0.597 / 0.586 | 200: 0.610 (110/80) / 0.610 (100/70) / 0.601 | 225: 0.621 (110/80) / 0.620 (100/80) / 0.608 |
250: 0.626 (130/110) / 0.618 (100/90) / 0.603.  ZnS phase repeats every ~130 nm (lambda/2n): ZnS ~225-250 restores the
first-order HTL optimum that ZnS 100 pushes to HTL ~250.

Practical TFE (`pary16.py` -> `pary16_grid.csv`, `pary16_best.csv`): ... HATCN 10 / Ag 8 joint / ZnS / Al2O3 30 / ParyC 500 /
Al2O3 30 / ParyC 500 / Al2O3 30 / n1.6 hemisphere MLA, Ag/Al0.2/Liq, ETL 60-140, HTL 40-160, ZnS 100-250, 27 wavelengths.
Best: ZnS 200-225, ETL 80, HTL 60: MLA 0.581 | before TFE (planar, air) 0.32 | TFE no lens 0.29.
Within 1 %p: ETL 70-100, HTL 50-80. Per ZnS best: 100 0.575, 125 0.565, 150 0.561, 175 0.569, 200 0.581, 225 0.581, 250 0.573.
Ag/Liq instead of Ag/Al0.2/Liq at (225, 80, 60): 0.599 / 0.324 / 0.291.

Top electrode IZO instead of HATCN/Ag/ZnS (`izo_top.py` -> `izo_top_grid.csv`, `izo_top_best.csv`; TAPC / HATCN 10 / IZO / same
ParyC TFE / n1.6 MLA, l_IZO). Best within ETL 60-140, HTL 40-160 (HTL always at the 160 bound):
IZO 30: 0.614 (ETL 90) | 50: 0.612 | 70: 0.613 | 100: 0.617 | 130: 0.618 (ETL 100). Before TFE 0.11-0.21, TFE no lens 0.15-0.21.
HTL beyond the bound: IZO 30/ETL 100: HTL 200 0.622, 240 0.630, 280 0.634; IZO 100/ETL 100: HTL 200-240 0.626.
l_IZO_SJ_RTP instead of l_IZO changes the result by < 0.5 %p. Optics only - sputter damage to the hole contact not included.

Lens-free planar TEOLED limit, emission straight into air (`lensfree_limit.py` -> `lensfree_limit.csv`, cross-seeded
`lensfree_refine.py` -> `lensfree_refined.csv`; the per-case DE fell into 2nd-order local optima for Theta 1, refinement
fixes them). Ag 100 | ETL | EML 2 nm sheet | HTL | top Ag | cap | air, PLQY 1 (EQE = eta_out), lossless non-dispersive organics.
Red 620 nm, bulk (McPeak) Ag with free thickness, n 1.8, ZnS cap; rows Theta, columns FWHM (delta/10/20/30/50 nm/Ir 68 nm):
0.67: 0.614 0.610 0.598 0.583 0.537 0.473 | 0.82: 0.713 0.709 0.700 0.686 0.640 0.575 | 0.90: 0.763 0.760 0.753 0.740 0.697 0.636
1.00: 0.827 0.825 0.821 0.811 0.778 0.733.  Optima are 1st order (ETL 57-67, HTL 68-76, top Ag 23-33 nm, ZnS 46-73).
Colours (Theta 1, delta / FWHM 20): 460 nm 0.760/0.726, 530 nm 0.804/0.789, 620 nm 0.827/0.821.
Variants (Theta 1 delta | Theta 0.82 Ir): measured HATCN/Ag 8 nm 0.633 | 0.470; transport n 1.5 0.894 | 0.780 (EML 2 nm sheet,
optimistic); cap n 1.8 0.798 | 0.555; 1st order only = free.  Forward-Lambertian reading of the Theta 1 optimum: 109 %.
Ladder to the real stack (`lensfree_ladder.py`, Theta 0.82, Ir): ideal 0.575 -> Ag 7.87 bulk 0.496 -> joint nk 0.470 ->
real stack (Bebq2/TAPC o,e/HATCN/Al0.2/Liq) 0.486 (ETL 61.5, HTL 51.4, ZnS 87.2) -> x0.96 = 0.466; Octave v2 0.4617.
`lensfree_variants.py`: ZnS >= 100 0.460 (61.1/47.3/100); Al 1 nm 0.447; both 0.438.
