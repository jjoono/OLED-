# Fig. 2(a) mock — ITO thickness, Al versus Ag

Question answered: does the ITO-thickness penalty show with the conventional Al reflector,
or only with Ag?

Stack: reflector (Al, JO n,k / Ag, McPeak n,k; 100 nm) / ETL 200 nm / EML 20 nm / HTL 200 nm /
ITO (n = 1.9 + 0.02i, 30-200 nm) / substrate; 550 nm, isotropic dipole, PLQY = 1, u grid 3000.
`dr4e.m` is `dr4d.m` plus a `TOPMAT` switch for the reflector.

The round-trip loss A' is split into a mirror part, A' recomputed with k_TCO = 0
(`tco3_*_k0.csv`), and a TCO part, the remainder (`tco2_*.csv` minus `tco3_*`).
CSV columns: d_ITO, eta_sub, A', wg, spp, abs, abs_top, abs_bottom, eta_ext, EQE.

| n_sub | reflector | A'_mirror | A'_TCO 50 -> 150 nm | TCO share | eta_ext 50 -> 150 nm |
|---|---|---|---|---|---|
| 1.5 (p 0.38) | Al | 0.130, flat | 0.046 -> 0.129 | 26 % -> 50 % | 0.772 -> 0.703 (-6.9 pp) |
| 1.5 (p 0.38) | Ag | 0.018, flat | 0.054 -> 0.147 | 75 % -> 89 % | 0.894 -> 0.788 (-10.6 pp) |
| 1.8 (p 0.30) | Al | 0.064, flat (dipole abs) | 0.042 -> 0.092 | 40 % -> 58 % | 0.675 -> 0.600 (-7.5 pp) |
| 1.8 (p 0.30) | Ag | 0.007, flat (dipole abs) | 0.049 -> 0.104 | 88 % -> 94 % | 0.819 -> 0.694 (-12.5 pp) |

So the trend is visible with Al, but "the dominant loss path" is only true with Ag: with Al the
TCO reaches parity with the mirror, with Ag it is essentially the only loss. The figure shows
both, which is the point of the section. In EQE the high-index substrate roughly doubles the
penalty (Al -5.3 -> -10.2 pp, Ag -9.9 -> -17.0 pp), which supports the following sentence in
the manuscript.

## Workbook

`fig2a_rawdata.xlsx` holds the data for the three panels, built from the four CSVs here.
Three sheets: README, `nsub_1.5` (the panel as drawn) and `nsub_1.8` (the companion sweep).

Columns per reflector: A′ mirror, A′ TCO, A′ total, eta_ext, eta_sub. The mirror part is A′
recomputed at k_ITO = 0 and the TCO part is the remainder, so the two stack to the total and
can be plotted directly as a stacked area. eta_ext is a live formula driven by the escape
probability on the README sheet (0.38 at n_sub = 1.5, 0.30 at 1.8), so a revised p updates the
column without re-running anything. All 72 formulas were recalculated and every value checked
against the simulation output to 5e-6, with the mirror-plus-TCO split verified to close on the
total.
