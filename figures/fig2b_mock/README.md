# Fig. 2(b) mock-up — substrate index, decomposed

The panel the parasitic-absorption paragraph cites, replacing the former Fig. 1(c). It answers
why raising the substrate index stalled the field near 60 %: η_sub^(0) keeps rising with the
index while η_ext keeps falling, and with a lossy mirror the second cancels the first.

Stack: reflector 100 nm (Al, JO n,k / Ag, McPeak n,k) / ETL 200 nm / EML 20 nm / HTL 200 nm /
ITO 50 nm (n = 1.864 + 0.0032i at 550 nm, Koenig et al. 2014) / substrate, index swept 1.30–2.00 and index-matched to the
outcoupling structure. 550 nm, isotropic dipole, PLQY = 1, u grid 3000 points.
`dr4f.m` is `dr4e.m` plus an `nsub` sweep mode. Curves assembled in `fig2b_curves.csv`
(n_sub, p, then η_sub / A′ / η_ext / EQE for each reflector).

η_ext = p/[p + (1−p)A′], with A′ = 1 − ⟨R_LED⟩ from the model and p interpolated from the
single-pass escape probability of the outcoupling structure, which falls roughly as 1/n_sub².

| n_sub | p | Al: η_sub / η_ext / EQE | Ag: η_sub / η_ext / EQE |
|---|---|---|---|
| 1.5 | 0.380 | 0.397 / 0.811 / 0.322 | 0.465 / 0.957 / 0.445 |
| 1.7 | 0.335 | 0.646 / 0.782 / 0.505 | 0.715 / 0.949 / 0.679 |
| 1.8 | 0.305 | 0.842 / 0.746 / 0.628 | 0.956 / 0.934 / 0.893 |
| 2.0 | 0.248 | 0.851 / 0.719 / 0.612 | 0.959 / 0.926 / 0.888 |

The Al device peaks at 0.628 near n_sub = 1.8 and then falls; the Ag device reaches 0.893 and
holds it. The Al ceiling now lands where the field actually stalled, just above 60 %.

## Layout variants

`fig2b_mock.png` (`plot_fig2b.py`) — two panels, one per reflector, four curves each.
These three layout studies were drawn before the Koenig re-run and still carry the
n = 1.9 + 0.02i numbers in their footnotes; they are kept for the layout argument, not for
the values. `fig2bc_squares.png` is the current panel.

`fig2b_single_vs_eqe.png` (`plot_fig2b_single.py`) — the author's intended layout (A: one
panel, p + η_ext(Al) + η_ext(Ag)) beside the same panel with the two EQE curves added (B).
A is self-contained: p falls with n_sub, so each photon makes more round trips and the two
η_ext curves peel apart, from a few %p at 1.3 to 15 %p at 2.0 — the cost of a lossy mirror
grows with the substrate index. What A cannot say is that the high index is worth having at
all, because every curve in it falls; that needs η_sub^(0) or the product.

`fig2b_single_inset.png` (`plot_fig2b_inset.py`) — variant C: A as the hero panel with the
two EQE curves in an inset. Keeps the drafted three-curve reading and still carries the
turnover at 0.54 / the climb to 0.75.

η_ext gap between the reflectors, from `fig2b_curves_konig.csv`: 4.8 %p at n_sub = 1.3,
14.6 %p at 1.5, 18.8 %p at 1.8, 20.7 %p at 2.0.

Still to settle with the author: whether to mark the reported record devices on the EQE curve,
as the earlier draft did with a shaded "previous works" circle. Individual literature points
with reference numbers would be stronger than a shaded region.

## Raw data

`fig2bc_rawdata.xlsx` — the numbers behind `fig2bc_squares.png` (and behind every other
layout variant, since they all plot the same table). Three sheets: `fig2b_2c` with the
plotted curves, where η_ext and EQE are live formulas driven by the editable yellow p
column, and `Al_modes` / `Ag_modes` with the full five-channel budget as the model returns
it. Built by `make_fig2bc_xlsx.py`; the formulas were checked by recalculating the workbook
in LibreOffice and comparing all 60 formula cells against `fig2b_curves.csv` (exact match).

Note that column 9 of `nsub_al.csv` / `nsub_ag.csv` is η_ext at the script's default p = 0.4
and is *not* the plotted curve — the workbook recomputes η_ext from the p column.

## Optical constants

Re-run with the Koenig et al. (2014) ITO (n = 1.8636 + 0.0032285i at 550 nm) after the author
supplied the tabulation: `nsub_{al,ag}_konig.csv`, assembled into `fig2b_curves_konig.csv`,
which is what the panels and the workbook now use. The earlier n = 1.9 + 0.02i run is kept as
`nsub_{al,ag}.csv` / `fig2b_curves.csv`; `DATA=fig2b_curves.csv python3 plot_fig2bc.py`
redraws it.
