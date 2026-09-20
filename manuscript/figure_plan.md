# Figure plan — unityEQE v12

Citation order in v12, strictly monotonic apart from one deliberate retrospective:

1(a) → 1(b) → 1(c) → 2(a)–(c) → 2(d),(e) → 2(f) → 2(g)–(j) → 3(a) → 3(b) → 3(c) → 4(a) →
4(b) → 4(c) → [3(a) retrospective] → 5(a) → 5(b) → 5(c)

Changes from v9: every Fig. 2 panel is recomputed with the Koenig et al. (2014) ITO
constants (n = 1.864 + 0.0032i at 550 nm) in place of the flat n = 1.9 + 0.02i, at the
author's request. The numbers below are the new ones.

Changes from v8: Fig. 2's index sweep is split across two small square panels — (b) carries p
and η_ext, (c) carries η_sub^(0) and EQE — following the author's revision of the
parasitic-absorption paragraph. The absorption maps move from (c) to (d), and the
Poynting-vector panel moves to the SI.

Changes from v6: Fig. 1 has three panels again. The escape-cone panel stays dropped — the
ray-optics limit is textbook material and now takes one sentence with no figure — but the
decay curve returns as (b) and the master curve moves to (c). The two are not redundant: the
decay curve shows the recycling *process* and belongs at this point in the paper, where the
comparison is still qualitative (large versus small round-trip loss) and no Al/Ag numbers have
been introduced; the master curve is the *conclusion*, with a quantitative R_LED axis.

Changes from v4: the non-existent Fig. 1(c) citation now points at 2(b). Sections 3 and 4 are
swapped, so the transparent-electrode rule follows the parasitic-absorption diagnosis
directly; Fig. 3 panels are renumbered accordingly (old c → a, a → b, b → c).

---

## Fig. 1 — The photon-recycling law

Three panels.

**(a) Schematic, two halves.** Left: a conventional OLED with an external outcoupling
structure. Light delivered to the substrate is partly transmitted at the structure/air
boundary and partly returned, reflects off the OLED stack and is recycled; each round trip
loses A′ to mirror ohmic absorption and TCO absorption. Right: the same picture with the
round-trip loss suppressed by the design rule of this work, so almost all substrate light
eventually escapes. Label P₀, P_sub, B_T, B_R, R_LED and A′ = 1 − R_LED.

**(b) Decay with round-trip number.** Light remaining in the substrate against the number of
round trips, for a stack with a large round-trip loss and one with a small loss. Qualitative at
this point: no Al/Ag labels and no numbers, since the quantitative comparison is Fig. 2's job.
This is the panel that makes the recycling process visible; the master curve alone hides it.

Axes and settings. The vertical axis is the **cumulative fraction of substrate-delivered light
that has escaped** after n round trips, so η_sub drops out of the calculation and only p and A′
remain; the value it saturates at is η_ext, which is the vertical axis of (c). The curve is
p[1 − (1−p)ⁿR_LEDⁿ]/[1 − (1−p)R_LED]. Draw three curves labelled by **A′ = 0.02, 0.10, 0.25**
with no material names — Al and Ag are introduced in Fig. 2 — and state **p** on the panel.
Six round trips on the horizontal axis is enough for all three to flatten at p = 0.4.

The escape probability is not a free parameter. For angularly randomised light inside the
substrate the geometric-optics estimate is T̄/n_sub², where T̄ is the transmittance averaged
over the escape cone (Yablonovitch, *J. Opt. Soc. Am.* **72**, 899 (1982), Sec. 4). Note the
factor of two: Yablonovitch quotes 1/(2n²) because his internal intensity is two-sided, which
suits a slab with a rear mirror; here each round trip involves one encounter with the
extraction structure, so the one-sided T̄/n_sub² applies.

| n_sub | 1/n_sub² | T̄/n_sub² | ray trace (Fig. 2(b)) |
|---|---|---|---|
| 1.3 | 0.592 | 0.556 | 0.665 |
| 1.5 | 0.444 | **0.404** | 0.380 |
| 1.8 | 0.309 | **0.267** | 0.305 |
| 2.0 | 0.250 | 0.210 | 0.248 |

The ray trace sits on the estimate at n_sub = 1.8 and 2.0 and exceeds it at 1.3, where the
curved microlens surface beats the flat escape cone — which is why the text calls T̄/n_sub² an
estimate rather than the value. p = 0.4 for the curve family of (c) is the n_sub = 1.5 entry.

Eq. (2) also assumes angular randomisation at every round trip, which is Yablonovitch's
ergodicity condition: the surface slope must exceed ½ arcsin(1/n_sub), i.e. 20.9° at n_sub =
1.5 and 16.9° at 1.8. A hemispherical microlens spans 0–90° and satisfies it comfortably,
whereas a plane-parallel slab is his explicit counter-example and gets no enhancement.

**(c) Master curve.** η_ext against the round-trip reflectance, with the wide view and a
zoom on the high-reflectance end (the author's existing draft). Settled conventions:

* x axis is **R_LED = 1 − A′** in both the main panel and the zoom, matching the definition
  in the text. Do not use R_LED(1 − A′_para) in one panel and R_LED in the other.
* the value of **p is stated**, and a small family of curves (e.g. p = 0.3, 0.4, 0.6) is drawn
  so the reader sees this is a one-parameter family; the p = 0.3 and 0.4 members tie directly
  to the escape-probability curve of Fig. 2(b).
* the horizontal "w/o outcoupling" and "w/ microcavity or horizontal emitter" bands are
  **removed**. Those are whole-device outcoupling values, not substrate-to-air extraction, so
  they do not belong on an η_ext axis; Fig. 2(b) already makes that comparison on an EQE axis.
  The two boxes that stay are the ones in the zoom: outcoupling with an Al electrode, and the
  proposed low-loss structure.
* the text states the reflective-electrode assumption T_LED = 0, so the caption should too.

## Fig. 2 — Parasitic absorption is the binding loss

Composed as one full-width figure in `figures/fig2_full/` (`plot_fig2.py`, with a draft
caption in that folder's README). Three rows of panels lettered a–j: **a,b,c** the ITO
thickness (loss split for Al, for Ag, then η_ext); **d,e,f** the substrate index (the two
inputs of eq. (2), then η_ext, then η_sub^(0) and EQE); **g,h,i,j** the three reflectors at a
fixed 150 nm ITO. The middle row was one panel and is now three, which is what makes the rows
balance; the split follows eq. (2), with (d) showing that p falls by 2.7× while A′ barely
moves, so the decline of η_ext is the escape probability's doing and not the mirror's.
Colour is the reflector throughout (vermillion Al, blue Ag, purple DBR) with grey reserved for
the TCO part of the loss.

**(a) ITO thickness, Al vs Ag.** Three panels, mocked up in `figures/fig2a_mock/`.
Left and centre: round-trip loss A' against ITO thickness, split into the mirror ohmic part
(A' recomputed at k_TCO = 0) and the TCO part (the remainder), for Al and for Ag. Right: the
resulting η_ext for both. With Al the mirror alone costs 13.4 % per round trip and the TCO is
a rounding error (6 % of A' at 50 nm, 15 % at 150 nm), so η_ext barely moves with ITO
thickness: −1.1 pp over 50 → 150 nm. With Ag the mirror loss drops to 1.8 % and the same ITO
becomes the main loss path (34 % → 60 %), with η_ext falling twice as fast, −2.2 pp. The Ag
panel needs a zoom inset, since its whole stack is a fifth of the Al one. Marking the green
device's measured η_ext = 0.916 at 150 nm is worth it: the model gives 0.934 there.

**(b) The cost: escape probability and extraction efficiency.** Small square panel.
η_ext against the substrate/MLA index for a lossy (Al) and a low-loss (Ag) reflector, with the
single-pass escape probability **p** on the same axes. p falls with n_sub, so each photon makes
more round trips and the two η_ext curves peel apart — 5 %p at n_sub = 1.3, 15 %p at 1.5,
21 %p at 2.0. This is the panel the parasitic-absorption paragraph
cites for "the trend grows stronger when a high-index substrate and an outcoupling structure
are used". p is the substrate angular distribution weighted by the outcoupling structure's
BSDF transmittance, confirmed by the author; label it **p** so the symbol matches eq. (2) and
leave the construction to the caption and Methods. p = 0.30 at n_sub = 1.8 is the value every
extraction-efficiency calculation in `sim/design_rule4/` rests on, and p = 0.30 and 0.40 are
two members of the curve family in Fig. 1(b).

**(c) The result: substrate-delivered power and EQE.** Small square panel, same x axis as (b).
η_sub^(0) (dashed) and EQE = η_sub^(0) η_ext (solid) for both reflectors. η_sub rises
monotonically with the index, but with Al the gain is spent on the loss in η_ext, so EQE turns
over at 0.63 near n_sub = 1.8; with Ag it keeps climbing to 0.89 and holds. The Al ceiling
lands just above 60 %, where the field actually stalled, which is the comparison the section
is making. This is the panel
for the closing sentence of the paragraph — with large parasitic absorption the gain in
substrate-delivered power is outrun by the loss in extraction, which is what limited the field
near 60 %. Mocked up together with (b) in `figures/fig2b_mock/fig2bc_squares.png`.

Open: whether to mark reported record devices on the EQE curve, as the earlier draft did with
a shaded "previous works" circle. Individual literature points with reference numbers would be
stronger than a shaded region.

**(g)–(j) Angle- and wavelength-resolved round-trip loss.** The transparent electrode is
fixed at 150 nm of ITO and only the reflector changes, so the panel compares Al, Ag and a
ten-pair ZnS/LiF quarter-wave stack at 550 nm. Flux- and spectrum-weighted 1 − R at
n_sub = 1.5: 15.6, 4.5 and 8.8 %, the last being 3.4 % absorbed plus 5.4 % leaked. (j) adds
the cos·sin weight as shading, which is how the panel carries "the angular distribution
matters as much as the reflectance", and a dashed curve for the DBR re-optimised for
randomised light (5.5 %). The leak stops at arcsin(1/n_sub) — marked — because beyond the
escape cone the light is trapped by total internal reflection at the back surface.

**Supplementary (was 2(d)): Poynting-vector loss accounting.** Layer-resolved dissipation
through one round trip: flat through the transparent organics, linear in thickness through the
TCO, and a step at the metal that is much larger for Al than for Ag. It is the evidence behind
the mirror/TCO split already drawn in 2(a), so as a main-text panel it repeats (a); the text
now describes it in one sentence and cites the SI.

## Fig. 3 — Design rules

**(a) Electrode quality.** η_sub^(0), η_ext and EQE against the TCO extinction coefficient
(left) and against the real index of a 10 nm Ag electrode (right), with the material bands of
the two electrode families. Mocked up in `figures/design_rule4/`. η_sub falls linearly, η_ext
hyperbolically. This is the panel the design-freedom section refers back to.

**(b) Waveguide and SPP against organic thickness and substrate index.** Mode-resolved power
fractions. Two axes are wanted: raising n_sub to n_EML moves the waveguided power into the
substrate mode, and increasing the organic thickness suppresses the SPP. Shows that with a
uniform organic index roughly 500 nm is needed for sub-10 % SPP loss.

**(c) Isotropic vs anisotropic ETL.** Substrate-delivered power against ETL thickness for an
isotropic ETL (n = 1.8) and low-n_e ETLs. Mark the threshold: the benefit appears only once
n_SPP drops below n_sub, which on Ag with n_o = 1.8 needs n_e near 1.5, moving the threshold
thickness from about 200 nm to about 80 nm. Data in `figures/fig3_Ag_mock/`.

## Fig. 4 — Experiment

**(a) Device structures.** Both stacks side by side: the green device (glass, MLA film,
Ag rear reflector on parylene) and the orange device (high-index MLA substrate, Al cathode,
auxiliary ZnS/LiF DBR).

**(b) Green device.** EQE against luminance with and without the MLA film: 20 % → 55.8 %.
State η_sub^(0) = 0.609 and η_ext = 91.6 %.

**(c) Orange device.** EQE against luminance for the planar reference, the MLA substrate with
Al only, and with the auxiliary DBR: 23 % / 61 % / 77 %. Add the back-leakage measurement
(18.6 % → 0.9 %) as an inset or adjacent panel — it is the only model-free measurement of the
recycling process in the paper and is currently buried in a supplementary table.

## Fig. 5 — Design freedom

**(a) Transport-layer thickness.** EQE against ETL and HTL thickness with and without the
outcoupling structure. Without it the Fabry–Pérot orders are enforced; with it the map
flattens. [[confirm the thickness above which it flattens for this stack — the generic model
puts it near 200 nm, not 100 nm]]

**(b) Dipole orientation.** EQE against horizontal dipole ratio, with and without the
outcoupling structure. Nearly linear without; nearly flat with, so an isotropic emitter still
exceeds 80 %.

**(c) Outcoupling-structure tolerance.** EQE against microlens aspect ratio, and against the
scattering parameters S and g for a scattering layer, for the proposed structure and an Al
reference. Shows the wider tolerance of the low-loss design.
