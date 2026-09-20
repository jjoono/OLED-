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


## ITO optical constants — what was used, and how much it matters

The TCO is entered as a single non-dispersive value at 550 nm: **n = 1.9, k = 0.02**
(`NTCO`, `KITO` in `dr4*.m`). n = 1.9 was chosen over 1.8 earlier because it matches the
measured sets in the project library and costs little in A'; k = 0.02 was the value the
design-rule scripts were first set up with and it sits at a quarter of Fig. 3(a)'s sweep
(k = 0 … 0.08).

What the project's own library (`nk_JH_total.mat`) holds at 550 nm:

| dataset | n | k |
|---|---|---|
| `ITO_SNU` | 1.998 | 0.0013 |
| `l_ITO` | 1.864 | 0.0032 |
| `ITO` | 1.809 | 0.0045 |
| `ITO_test` | 1.827 | 0.0038 |
| `etri_ITO` | 1.922 | **0.0481** |
| `IZO` | 2.062 | 0.0012 |
| `l_IZO` | 2.044 | 0.0053 |
| `IZO_NIR` | 1.921 | 0.0128 |

So n = 1.9 is squarely inside the measured range (1.81–2.00), but k = 0.02 is four to six
times the library's clean ITO films and roughly half the ETRI film. It is a "typical
sputtered ITO" choice, not a measured one.

How much the choice matters, from `aprime.py` (a plain TMM for A' alone, validated against
the Octave output to 1e-4). "Parity" is the ITO thickness at which the TCO absorption equals
the mirror's ohmic loss; eta_ext is evaluated at p = 0.305, n_sub = 1.8.

| k_ITO | Al: TCO share at 50 / 150 nm | Al: η_ext at 50 / 150 nm | Al: parity | Ag: TCO share at 50 / 150 nm | Ag: η_ext at 50 / 150 nm | Ag: parity |
|---|---|---|---|---|---|---|
| 0.0032 | 8 % / 17 % | 0.746 / 0.734 | — | 46 % / 65 % | 0.935 / 0.904 | 64 nm |
| 0.005 | 12 % / 24 % | 0.738 / 0.716 | 482 nm | 57 % / 74 % | 0.920 / 0.874 | 37 nm |
| 0.01 | 21 % / 38 % | 0.717 / 0.673 | 258 nm | 72 % / 85 % | 0.883 / 0.803 | 19 nm |
| 0.02 | 34 % / 54 % | 0.680 / 0.606 | 126 nm | 83 % / 91 % | 0.822 / 0.699 | 10 nm |
| 0.03 | 42 % / 62 % | 0.650 / 0.556 | 82 nm | 87 % / 94 % | 0.774 / 0.626 | 7 nm |
| 0.048 | 52 % / 71 % | 0.607 / 0.495 | 45 nm | 91 % / 96 % | 0.707 / 0.540 | 5 nm |
| 0.08 | 62 % / 77 % | 0.552 / 0.430 | 26 nm | 94 % / 97 % | 0.627 / 0.455 | 5 nm |

Two things follow.

* The **Ag statement is robust**: the TCO is the dominant round-trip loss for every k in the
  table, 46 % of A' at 50 nm even with the library's cleanest film. "With a low-loss mirror
  the TCO is essentially the only loss path" holds whatever ITO is assumed.
* The **Al statement is not**: "TCO absorption grows to match the mirror loss" needs
  k >= about 0.015 for parity to fall inside the plotted thickness range. At k = 0.005 parity
  is at 482 nm, i.e. never in practice.

n_TCO barely touches A' (0.206 / 0.195 at n = 1.8 / 2.0, k = 0.02, 50 nm); its real effect is
on eta_sub, where going 1.8 -> 2.0 costs about 3.5 %p because the TCO becomes its own
waveguide.

Open decision for the author: keep k = 0.02 and state it in the caption as a typical
device-grade film (and then quote the parity thickness, 126 nm, so the claim is checkable),
or drop to the library's own measured ITO (k ~ 0.003–0.005), which is the more defensible
provenance but leaves the Al device's round-trip loss almost entirely mirror-driven. The
second reads as a cleaner two-step story — with Al the mirror is the problem, and once the
mirror is fixed the TCO becomes the problem — and it does not weaken any Ag claim. Either
way the value belongs in Methods, since Fig. 3(a) sweeps this exact parameter.
