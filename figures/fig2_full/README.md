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
(**g**), Ag (**h**), and ten pairs of ZnS/LiF designed as a quarter-wave at 550 nm (**i**).
**j**, the same averaged over the green emission spectrum, with the cos θ sin θ weight shaded;
the flux- and spectrum-weighted losses are 15.6, 4.5 and 8.8 %. The dielectric mirror absorbs
less than Ag (3.4 %) but transmits outside its stopband, and only for angles inside the escape
cone to air; re-optimising the two thicknesses for angularly randomised light brings it to
5.5 % (dashed).

All panels: 550 nm family, isotropic dipole, PLQY = 1, generic stack of substrate /
transparent electrode / 420 nm non-absorbing organics (n = 1.8) / 100 nm reflector; McPeak Ag,
Johnson–Christy-type Al, Koenig ITO, measured IZO, ZnS and LiF. **a–c** and **g–j** on glass
(n_sub = 1.5, p = 0.38); **d–f** sweep the substrate index.

## If it has to be smaller

Drop **g** and **h**: panel **j** carries the numbers and the angular story on its own, and the
two metal maps are the least surprising part of the figure. That leaves 3 + 3 + 2 and frees
roughly a fifth of the height.
