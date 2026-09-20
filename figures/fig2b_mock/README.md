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
holds it. Four curves per panel is close to the limit of what one panel carries — if it reads
as crowded, p can move to the caption, since the text already gives its value and origin.

Still to settle with the author: whether to mark the reported record devices on the EQE curve,
as the earlier draft did with a shaded "previous works" circle. Individual literature points
with reference numbers would be stronger than a shaded region.
