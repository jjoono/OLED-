# Figure plan — unityEQE v5

Citation order in v5, strictly monotonic apart from one deliberate retrospective:

1(a) → 1(b) → 1(c) → 2(a) → 2(b) → 2(c) → 2(d) → 3(a) → 3(b) → 3(c) → 4(a) → 4(b) → 4(c) →
[3(a) retrospective] → 5(a) → 5(b) → 5(c)

Changes from v4: Fig. 1 has three panels, not four — the old 1(d) becomes 1(c) and the
non-existent 1(c) citation now points at 2(b). Sections 3 and 4 are swapped, so the
transparent-electrode rule follows the parasitic-absorption diagnosis directly; Fig. 3
panels are renumbered accordingly (old c → a, a → b, b → c).

---

## Fig. 1 — The photon-recycling law

**(a) Schematic, two halves.** Left: a conventional OLED with an external outcoupling
structure. Light delivered to the substrate is partly transmitted at the structure/air
boundary and partly returned, reflects off the OLED stack and is recycled; each round trip
loses A′ to mirror ohmic absorption and TCO absorption. Right: the same picture with the
round-trip loss suppressed by the design rule of this work, so almost all substrate light
eventually escapes. Label P₀, P_sub, B_T, B_R, R_LED and A′ = 1 − R_LED.

**(b) Escape cone without an outcoupling structure.** Ray-optics limit: the escape cone in
the substrate and the resulting extraction, scaling as 1/n_sub². Supports the definition of
η_ext and the contrast with the recycling case.

**(c) Decay with round-trip number.** Light remaining in the substrate against the number of
round trips, for a lossy stack (Al, A′ > 0.2) and a low-loss stack. The Al curve decays fast;
the low-loss curve stays high until the light has escaped. Only the decay curve — the master
curve η_ext vs R_LED and the EQE-vs-n_sub curve are not in Fig. 1.

## Fig. 2 — Parasitic absorption is the binding loss

**(a) ITO thickness, Al vs Ag.** Three panels, already mocked up in `figures/fig2a_mock/`.
Left and centre: round-trip loss A′ against ITO thickness, split into the mirror ohmic part
(A′ recomputed at k_TCO = 0) and the TCO part (the remainder), for Al and for Ag. Right: the
resulting η_ext for both. Shows that the mirror loss is flat in ITO thickness while the TCO
loss grows linearly, so with Al the TCO reaches parity with the mirror (26 % → 50 % of A′
over 50–150 nm) and with Ag it is essentially the only loss (75 % → 89 %).

**(b) Substrate index sweep, decomposed.** η_sub^(0), η_ext and their product EQE against the
substrate/MLA index, for a lossy and a low-loss reflector. η_sub rises monotonically while
η_ext falls, so a device with large parasitic absorption saturates around 60 % near n ≈ 1.65,
whereas suppressing the absorption removes the trade-off and EQE rises monotonically. This is
the panel the parasitic-absorption paragraph cites; it replaces the former Fig. 1(c).

**(c) Angle- and wavelength-resolved round-trip absorption.** Maps of the absorption for light
incident from the substrate onto the stack, for the candidate bottom-electrode structures.
Al averages 15–25 %, Ag 5–15 %. Include the angular distribution of substrate-delivered power
alongside, since the weighting matters as much as the reflectance.

**(d) Poynting-vector loss accounting.** Layer-resolved dissipation through one round trip:
flat through the transparent organics, linear in thickness through the TCO, and a step at the
metal that is much larger for Al than for Ag.

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
