# Directional selectivity of multi-pass OLED extraction

How much of the light can be put into a forward cone when the outcoupling structure is
replaced by **a randomizing layer + an angular filter**, and how this depends on the
round-trip loss A' of the OLED, on the filter (ideal vs. realizable DBR) and on how well
the returned light is randomized.

Inputs follow the round-trip-loss paper (Kim, Lee, Yoo, Gather): CPS dipole model with
birefringent organics and the group's dispersive optical constants (nk library, not
redistributed; `extract_inputs.py` pulls the subset used), trans-scale multi-pass series
(Suppl. Note 1 of that paper), stacks of its Fig. 2 and Ext. Data Table 1, the ZnS/LiF
DBR technology of its Suppl. Note 2, and measured emission spectra.

## Model

`device | randomizer | filter | air`, solved per wavelength (400-700 nm, 2 nm), 0.5 deg
substrate-angle bins, PLQY = 1, Theta = 0.865:

    u  = (I - Rs Df)^-1 Ts v,   BT = diag(Tf) U,   BR = Rs + Ts Df U
    x  = (I - diag(R_LED) BR)^-1 P_sub,           air = BT x

* devices (`sel.STACKS`): Al or Ag cathode on glass (platform of the MLA study, ITO 150 nm),
  ITO cathode + chirped 10-pair ZnS/LiF rear DBR, low-loss Ag/ITO 50 nm stack on n 1.51 and 1.80.
  A' (Lambertian, 530 nm): 0.170, 0.049, 0.050, 0.027, 0.031.
* randomizers: idealized partial mixing m; index-matched Henyey-Greenstein slab (`slab.py`,
  Monte Carlo, g, tau* = S(1-g), albedo); none.
* filters: ideal cone (T = 1 inside, R = 1 outside), bare Fresnel exit, quarter-wave edge
  DBRs (TiO2/SiO2, SiNx/SiO2, ZnS/LiF; 4-16 pairs; edge 420-820 nm), and layer-optimized
  filters (`optfilter.py`; with a full randomizer the series has an exact per-wavelength
  closed form, which makes differential evolution over all thicknesses cheap).
* references: hemispherical MLA (LightTools BSDF, trans-scale series), strong-microcavity
  bottom emitter (20 nm Ag anode, ETL/HTL scanned), planar.
* metrics: EQE into an air cone of half-angle beta; on-axis luminance in Lambertian-equivalent
  EQE units, EQE(<5 deg)/sin^2(5 deg).

Checks: matrix model vs closed form p_b/[p_b+(1-p_b)A'] within 1-3 %; MLA reference
EQE 0.562 (Al) and 0.876 (low-loss n 1.8) as in the papers; fast closed form vs matrix
model for a 34-layer DBR: 0.24386 vs 0.24387.

## Run order

    python3 extract_inputs.py <path>/nk_JH_total.mat
    python3 run_devices.py      # devices + slabs            (~1 min)
    python3 run_sweep.py        # 1224 filters x 5 stacks x 3 randomizers (~30 min, 4 cores)
    python3 run_refs.py         # MLA and strong-cavity references (~7 min)
    python3 run_random.py       # randomization sweep (~2 min)
    python3 run_kito.py         # transparent-electrode loss sweep
    python3 optfilter.py '<json list of cases>'   # layer-optimized filters
    python3 analyze.py; python3 make_figs.py

Outputs: `results/*.csv`, `figures/figSel1-4.{png,pdf,xlsx}`.
