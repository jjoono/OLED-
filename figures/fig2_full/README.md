# Figure 2, composed — three rows, ten panels

`plot_fig2.py` draws the whole of Fig. 2 at full journal width (180 mm), `fig2_full.png` at
320 dpi and `fig2_full.pdf` as vector. It reads what is already in the repository —
`figures/fig2a_mock/kn*.csv`, `figures/fig2b_mock/fig2b_curves_konig.csv` and
`sim/design_rule4/fig2d.py` — so nothing is recomputed by hand. `plot_fig2_v1.py` is the
earlier two-row version, kept for comparison.

## Layout

| row | panels | subject |
|---|---|---|
| 1 | a, b, c | ITO thickness: the loss split for Al, for Ag, and the extraction efficiency |
| 2 | d, e, f | substrate index: the two inputs of eq. (2), η_ext, then η_sub^(0) and EQE |
| 3 | g, h, i, j | the three reflectors at a fixed 150 nm ITO, in angle and wavelength |

Row 2 was one panel in the earlier draft and is now three, which is what makes the rows
balance. The split follows eq. (2) rather than the drawing: **d** shows the two quantities
that go in — p falls by 2.7× over the range while A′ barely moves, so the whole decline of
η_ext comes from the escape probability, not from the mirror getting worse — **e** shows what
comes out, and **f** the product with the substrate-delivered power.

Colour is the reflector throughout: vermillion Al, blue Ag, purple DBR, with grey reserved for
the TCO part of the loss.

## Draft caption

**Fig. 2 | Parasitic absorption is the binding loss.**
**a–c**, ITO thickness. Round-trip loss A′ = 1 − ⟨R_LED⟩ split into the mirror's ohmic part
(colour; A′ recomputed at k_TCO = 0) and the TCO part (grey), for an Al (**a**) and an Ag
(**b**) reflector — note the three times finer axis in **b** — and the resulting extraction
efficiency (**c**). The mirror part is flat in thickness while the TCO part grows, so the TCO
carries 6 % of A′ at 50 nm with Al and 60 % at 150 nm with Ag. Star: the 0.916 measured on the
green device.
**d–f**, The index of the substrate and of the index-matched outcoupling structure. **d**, the
two quantities eq. (2) takes: the single-pass escape probability p falls by 2.7× across the
range while A′ barely moves. **e**, the extraction efficiency that follows — the same
round-trip loss costs more when p is small, so the two reflectors peel apart by 21 percentage
points. **f**, substrate-delivered power η_sub^(0) (thin) and the product EQE (bold); the
shaded gap is light that reaches the substrate and never escapes. With Al the gain in
η_sub^(0) is spent on the loss in η_ext and EQE turns over at 0.63; with Ag it keeps rising to
0.89. The step at n_sub = n_EML is the waveguide cutoff.
**g–j**, The three reflectors at a fixed 150 nm ITO, resolved in angle and wavelength: Al
(**g**), Ag (**h**), and ten pairs of ZnS/LiF whose thicknesses were optimised over all
substrate angles and over the emission spectrum (**i**). **j**, the same averaged over the
green emission spectrum, with the cos θ sin θ weight shaded; the flux- and spectrum-weighted
losses are 15.6, 4.5 and 5.5 %. The dielectric mirror absorbs less than Ag (3.4 %) but
transmits outside its stopband, and only for angles inside the escape cone to air; designing
the same ten pairs as a plain quarter-wave stack at 550 nm costs 8.8 % (dashed).

All panels: 550 nm family, isotropic dipole, PLQY = 1, generic stack of substrate /
transparent electrode / 420 nm non-absorbing organics (n = 1.8) / 100 nm reflector; McPeak Ag,
Johnson–Christy-type Al, Koenig ITO, measured IZO, ZnS and LiF. **a–c** and **g–j** on glass
(n_sub = 1.5, p = 0.38); **d–f** sweep the substrate index.

## If it has to be smaller

Drop **g** and **h**: panel **j** carries the numbers and the angular story on its own, and the
two metal maps are the least surprising part of the figure. That leaves 3 + 3 + 2 and frees
roughly a fifth of the height.

## A note on "optimised"

The manuscript says the DBR thicknesses came from a genetic algorithm, which is what Methods
records for the orange device. The optimum drawn in this figure was found by an exhaustive
grid search over the two thicknesses (`sim/design_rule4/fig2d.py` plus the search recorded in
`figures/fig2d_mock/README.md`) — with only two free parameters that search is global, and it
lands on the same point the device's genetic algorithm found: ZnS 73 / LiF 110 nm against the
device's 70 / 115 nm, 6.04 % against 6.05 % at n_sub = 1.8. Worth keeping the wording in
Methods aligned with whichever method is described there.

## Raw data for (g)-(j)

`fig2ghij_rawdata.xlsx` (built by `export_ghij.py`): the three colour maps as matrices — rows
wavelength in 5 nm steps, columns the angle in the substrate in 1° steps — plus the
transmitted part of the DBR map on its own, the four spectrum-averaged angle curves of (j)
with the cos·sin weight, and a summary sheet with the weighted numbers at both substrate
indices. `ghij_csv/` holds the same arrays as plain CSV.

`fig2i_DBR_fine.xlsx` (`export_i_fine.py`) is panel (i) on a four times finer grid — every
0.25° from 0 to 89.75° and every 1 nm — with three sheets for the dielectric stack (the
round-trip loss, the transmitted part alone, and absorption alone) and a fourth with Ag's
absorption on the same grid for a direct comparison. That resolution is comfortably
converged: the stack's three angular features are 1.3°, 2.4° and 2.7° wide at half maximum, so
0.25° puts five to ten points across each, and in the wavelength direction the value moves by
at most 0.40 %p per nm. Exact grazing is left out — it carries zero flux weight and the
transfer matrix is singular there.

Spot-checked against the plotted curves: at normal incidence 14.41 / 3.75 / 10.48 / 2.64 %
for Al, Ag, the optimised DBR and the plain quarter-wave stack, matching the figure.

The raw data for (a)-(c) is in `figures/fig2a_mock/fig2a_rawdata.xlsx` and for (d)-(f) in
`figures/fig2b_mock/fig2bc_rawdata.xlsx`.

## Why the DBR map looks worse than it is

Asked whether the dielectric mirror's absorption beyond the escape cone is really as high as
panel (i) suggests. It is not, and the numbers say so.

Beyond arcsin(1/n_sub) = 41.8° the transmitted fraction is identically zero — checked, the
maximum T in the region 45–90°, 470–560 nm is 0.0 % — so there the map already shows
absorption alone. In that region the DBR's median loss is **4.74 %** against Ag's **5.47 %**,
and it exceeds Ag in only 29 % of the region, on the narrow resonance fringes. Angle by angle,
spectrum-averaged: 3.78 vs 4.59 % at 45°, 4.23 vs 4.99 % at 60°, 4.50 vs 5.41 % at 70°. The
DBR is the better mirror everywhere except within a few degrees of grazing.

Weighted over the whole hemisphere: absorption 3.38 % for the DBR against 4.49 % for Ag. The
difference between the two is entirely the metal's ohmic loss (at normal incidence 2.46 % vs
3.75 %, and the 1.3 pp gap is the Ag mirror). The DBR's total only exceeds Ag's because of the
2.09 % it transmits **inside** the escape cone.

The impression comes from the display: Ag's map is smooth while the DBR's carries sharp
guided-mode resonances of the ITO/organic slab, and the eye reads the fringes rather than the
median. The numbers are converged — the hemisphere average moves by 0.005 pp between 451 and
3601 angles — so the fringes are physical, not sampling. Each map now carries its own
flux- and spectrum-weighted number underneath, and panel (i) marks the escape cone with the
two regimes labelled, so the comparison cannot be read off the fringes.
