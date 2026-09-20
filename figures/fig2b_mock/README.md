# Fig. 2(b) mock-up — substrate index, decomposed

The panel the parasitic-absorption paragraph cites, replacing the former Fig. 1(c). It answers
why raising the substrate index stalled the field near 60 %: η_sub^(0) keeps rising with the
index while η_ext keeps falling, and with a lossy mirror the second cancels the first.

Stack: reflector 100 nm (Al, JO n,k / Ag, McPeak n,k) / ETL 200 nm / EML 20 nm / HTL 200 nm /
ITO 50 nm (n = 1.9 + 0.02i) / substrate, index swept 1.30–2.00 and index-matched to the
outcoupling structure. 550 nm, isotropic dipole, PLQY = 1, u grid 3000 points.
`dr4f.m` is `dr4e.m` plus an `nsub` sweep mode. Curves assembled in `fig2b_curves.csv`
(n_sub, p, then η_sub / A′ / η_ext / EQE for each reflector).

η_ext = p/[p + (1−p)A′], with A′ = 1 − ⟨R_LED⟩ from the model and p interpolated from the
single-pass escape probability of the outcoupling structure, which falls roughly as 1/n_sub².

| n_sub | p | Al: η_sub / η_ext / EQE | Ag: η_sub / η_ext / EQE |
|---|---|---|---|
| 1.5 | 0.380 | 0.384 / 0.772 / 0.297 | 0.455 / 0.894 / 0.407 |
| 1.7 | 0.335 | 0.631 / 0.731 / 0.462 | 0.704 / 0.866 / 0.609 |
| 1.8 | 0.305 | 0.798 / 0.680 / 0.543 | 0.906 / 0.822 / 0.745 |
| 2.0 | 0.248 | 0.827 / 0.637 / 0.527 | 0.931 / 0.786 / 0.732 |

The Al device peaks at 0.54 near n_sub = 1.8 and then falls; the Ag device reaches 0.75 and
holds it.

## Layout variants

`fig2b_mock.png` (`plot_fig2b.py`) — two panels, one per reflector, four curves each.

`fig2b_single_vs_eqe.png` (`plot_fig2b_single.py`) — the author's intended layout (A: one
panel, p + η_ext(Al) + η_ext(Ag)) beside the same panel with the two EQE curves added (B).
A is self-contained: p falls with n_sub, so each photon makes more round trips and the two
η_ext curves peel apart, from a few %p at 1.3 to 15 %p at 2.0 — the cost of a lossy mirror
grows with the substrate index. What A cannot say is that the high index is worth having at
all, because every curve in it falls; that needs η_sub^(0) or the product.

`fig2b_single_inset.png` (`plot_fig2b_inset.py`) — variant C: A as the hero panel with the
two EQE curves in an inset. Keeps the drafted three-curve reading and still carries the
turnover at 0.54 / the climb to 0.75.

η_ext gap between the reflectors, from `fig2b_curves.csv`: 4.4 %p at n_sub = 1.3, 12.2 %p at
1.5, 14.2 %p at 1.8, 14.9 %p at 2.0 (ratio 1.05 → 1.23).

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
