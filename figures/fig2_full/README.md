# Figure 2, composed

`plot_fig2.py` draws the whole of Fig. 2 as one figure at full journal width (180 mm,
`fig2_full.png` at 320 dpi and `fig2_full.pdf` as vector). It reads the data already in the
repository — `figures/fig2a_mock/kn*.csv`, `figures/fig2b_mock/fig2b_curves_konig.csv` and
`sim/design_rule4/fig2d.py` — so nothing is recomputed by hand.

## Why this layout

The natural panel counts (3 + 1 + 1 + 4) do not tile. Two changes make them:

* **(a) goes from three panels to two.** The two stacked-area plots were one per reflector
  because the Al and Ag stacks both start at zero and would overlap. Drawn as grouped stacked
  bars at five thicknesses, both fit in one panel, and "the mirror part is flat, the TCO part
  grows" reads better from bars than from an area whose lower edge never moves.
* **(d) goes from four panels to three.** The DBR map is dropped; the structure is still in
  the angle plot, and its stopband behaviour is an SI topic. The two maps that stay are the
  conventional stack and the design rule, which is the comparison the text makes.

That gives two rows of equal height: a | a | b | c on top, and the two maps, the angle plot
and a weighted-average key below. Row 1 is four equal cells; row 2 spends the same width on
two narrow maps, a wide line plot and a text key.

**Colour is the reflector everywhere**: vermillion = Al, blue = Ag, and grey is reserved for
the TCO part of the loss in (a). The earlier mock-ups used orange for the TCO part in (a) and
orange for Al in (b)–(c), which collided as soon as the panels were put side by side.

## Draft caption

**Fig. 2 | Parasitic absorption is the binding loss.**
**a**, Round-trip loss A′ = 1 − ⟨R_LED⟩ of the generic stack against ITO thickness, split into
the mirror's ohmic part (colour; A′ recomputed at k_TCO = 0) and the TCO part (grey; the
remainder), and the resulting extraction efficiency η_ext = p/[p + (1 − p)A′]. With Al the
mirror alone costs 13.4 % per round trip and the TCO is 6 % of A′ at 50 nm; with Ag the mirror
costs 1.8 % and the same film becomes the main loss path, 34 % at 50 nm and 60 % at 150 nm.
Star: the 0.916 measured on the green device.
**b**, η_ext against the index of the substrate and the index-matched outcoupling structure,
with the single-pass escape probability p (dashed). A higher index lowers p, so each photon
makes more round trips and the penalty for a lossy mirror grows from 5 to 21 percentage
points.
**c**, Substrate-delivered power η_sub^(0) (thin) and the product EQE = η_sub^(0) η_ext
(bold); the shaded gap is light that reaches the substrate and never escapes. With Al the gain
in η_sub^(0) is spent on the loss in η_ext and EQE turns over at 0.63 near n_sub = 1.8; with Ag
it keeps rising to 0.89. The step at n_sub = n_EML is the waveguide cutoff.
**d**, Round-trip loss 1 − R resolved in angle and wavelength for the conventional stack and
for the design rule, and averaged over the green emission spectrum for all five candidate
electrode structures. Shading: the cos θ sin θ weight an angularly randomising outcoupling
structure enforces. The dielectric mirror absorbs as little as the best metal but transmits
outside its stopband, and only for angles inside the escape cone to air, so it loses 11.5 % on
glass and 6.4 % at n_sub = 1.8.

All panels: 550 nm family, isotropic dipole, PLQY = 1, generic stack of substrate /
transparent electrode / 420 nm non-absorbing organics (n = 1.8) / 100 nm reflector; McPeak Ag,
Johnson–Christy-type Al, Koenig ITO, measured IZO, ZnS and LiF. (a) and (d) on glass
(n_sub = 1.5, p = 0.38); (b) and (c) sweep the substrate index.

## If a panel has to go

(d)'s two maps are the first thing to cut — the angle plot carries the numbers and the
angular story on its own, and the maps could move to the SI. That would leave a 4 + 2 layout
and free roughly a third of the figure height.
