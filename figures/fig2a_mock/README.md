# Fig. 2(a) mock — ITO thickness, Al versus Ag

Question answered: does the ITO-thickness penalty show with the conventional Al reflector,
or only with Ag?

Stack: reflector (Al, JO n,k / Ag, McPeak n,k; 100 nm) / ETL 200 nm / EML 20 nm / HTL 200 nm /
ITO (n = 1.8636 + 0.0032285i at 550 nm, Koenig et al. 2014; 30-200 nm) / substrate; 550 nm,
isotropic dipole, PLQY = 1, u grid 3000.
`dr4e.m` is `dr4d.m` plus a `TOPMAT` switch for the reflector.

The round-trip loss A' is split into a mirror part, A' recomputed with k_TCO = 0
(`kn3_*_k0.csv`), and a TCO part, the remainder (`kn2_*.csv` minus `kn3_*`).
CSV columns: d_ITO, eta_sub, A', wg, spp, abs, abs_top, abs_bottom, eta_ext, EQE.

| n_sub | reflector | A'_mirror | A'_TCO 50 -> 150 nm | TCO share | eta_ext 50 -> 150 nm |
|---|---|---|---|---|---|
| 1.5 (p 0.38) | Al | 0.134, flat | 0.008 -> 0.023 | 6 % -> 15 % | 0.811 -> 0.800 (-1.1 pp) |
| 1.5 (p 0.38) | Ag | 0.018, flat | 0.009 -> 0.026 | 34 % -> 60 % | 0.957 -> 0.934 (-2.2 pp) |
| 1.8 (p 0.30) | Al | 0.137, flat | 0.013 -> 0.031 | 8 % -> 18 % | 0.741 -> 0.721 (-2.0 pp) |
| 1.8 (p 0.30) | Ag | 0.017, flat | 0.014 -> 0.034 | 46 % -> 67 % | 0.933 -> 0.894 (-3.9 pp) |

So the reading is now a two-step one. With Al the mirror alone costs 13.4 % per round trip
and the TCO is a rounding error (6 % of A' at 50 nm, 15 % at 150 nm), so ITO thickness barely
moves eta_ext: -1.1 pp over 50 -> 150 nm. Drop the mirror loss to 1.8 % with Ag and the same
ITO film becomes the main loss path (34 % -> 60 %), and eta_ext falls twice as fast, -2.2 pp.
That is the sentence the manuscript now carries: fix the mirror first, and the transparent
electrode is what is left.

A useful check: at 150 nm ITO on glass with an Ag reflector the model gives eta_ext = 0.934
against the 0.916 measured on the green device. The earlier k = 0.02 assumption gave 0.788,
well below the measurement — on its own a good reason to prefer the Koenig constants.

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


## ITO optical constants — settled

The TCO is **n = 1.8636 + 0.0032285i at 550 nm**, interpolated from the Koenig et al. (2014)
tabulation on refractiveindex.info that the author supplied (`Konig.csv`). It coincides with
the `l_ITO` entry of the project library to four decimals, so the library already held this
dataset. The earlier runs used a flat n = 1.9 + 0.02i; those files (`tco2_*`, `tco3_*`) are
kept for comparison, and `plot_fig2a.py` can redraw them with `PRE2=tco2 PRE3=tco3`.

What the project library holds at 550 nm, for reference: `ITO_SNU` 1.998 + 0.0013i,
`l_ITO` 1.864 + 0.0032i, `ITO` 1.809 + 0.0045i, `ITO_test` 1.827 + 0.0038i,
`etri_ITO` 1.922 + 0.0481i, `IZO` 2.062 + 0.0012i, `l_IZO` 2.044 + 0.0053i,
`IZO_NIR` 1.921 + 0.0128i. Koenig sits with the clean films; the ETRI film is an order of
magnitude lossier.

How much the choice matters, from `sim/design_rule4/aprime.py` (a plain TMM for A' alone,
validated against the Octave output to 1e-4). "Parity" is the ITO thickness at which the TCO
absorption equals the mirror's ohmic loss; eta_ext at p = 0.305, n_sub = 1.8, n_TCO = 1.9.

| k_ITO | Al: TCO share at 50 / 150 nm | Al: parity | Ag: TCO share at 50 / 150 nm | Ag: parity |
|---|---|---|---|---|
| 0.0032 (Koenig) | 8 % / 17 % | — | 46 % / 65 % | 64 nm |
| 0.005 | 12 % / 24 % | 482 nm | 57 % / 74 % | 37 nm |
| 0.01 | 21 % / 38 % | 258 nm | 72 % / 85 % | 19 nm |
| 0.02 (former) | 34 % / 54 % | 126 nm | 83 % / 91 % | 10 nm |
| 0.048 (ETRI) | 52 % / 71 % | 45 nm | 91 % / 96 % | 5 nm |
| 0.08 | 62 % / 77 % | 26 nm | 94 % / 97 % | 5 nm |

The Ag statement survives every entry in that table — the TCO is the dominant round-trip loss
for any ITO anyone would deposit. The Al statement does not: "TCO absorption grows to match
the mirror loss" needed k >= about 0.015 and has been dropped from the manuscript.

n_TCO barely touches A' (0.206 / 0.195 at n = 1.8 / 2.0, k = 0.02, 50 nm); its real effect is
on eta_sub, where going 1.8 -> 2.0 costs about 3.5 %p because the TCO becomes its own
waveguide.
