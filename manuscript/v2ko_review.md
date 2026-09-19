# unityEQE_v2_ko — storyline review and applied corrections

Revised file: `manuscript/unityEQE_v2_ko_rev.docx` (16 of 46 paragraphs touched; XSD-validated,
paragraph count unchanged, rendered and checked). The author's register and sentence style were
kept; only the statements below were changed.

## Corrected — factual

1. **Mirror reflection loss is ohmic, not plasmonic** (intro, Fig. 1a paragraph, round-trip
   paragraph). Absorption when light reflects off the metal electrode is Joule dissipation in
   the metal. Surface-plasmon coupling is a separate near-field channel at the emitter. The
   manuscript already used "ohmic 손실" correctly in the Fig. 2(d) paragraph, so this was also
   an internal inconsistency.
2. **The low-n_e ETL works through a threshold, not a gradual evanescent cut-off.** The old text
   said large in-plane wavevectors "become evanescent in the transport layer and never reach the
   metal". What actually happens is that a low n_e shifts the plasmon dispersion: once n_SPP
   falls below the substrate index the mode becomes leaky and its power returns as substrate
   light. Simulations in `sim/design_rule4/` show the u>1 bin jumping from 0.013 to 0.047 across
   that threshold, and the threshold moving from n_e = 1.60 to 1.68 when n_sub goes 1.8 to 1.9.
3. **The anisotropic-ETL benefit needs the threshold qualifier** (Fig. 3b paragraph). On Ag with
   n_o = 1.8, dropping n_e to 1.7 or 1.6 leaves the plasmon bound and buys nothing; only near
   n_e = 1.5 does the threshold thickness fall from about 200 nm to about 80 nm. B3PyMPM and
   B4PyMPM sit close to that threshold, so the material criterion is n_SPP against n_sub, not
   n_e alone.
4. **The index ladder needs n_TCO <= n_sub.** Matching the substrate to the EML removes the
   organic waveguide, but a transparent electrode at n = 1.9-2.0 above a substrate at 1.8 guides
   light of its own. Worth 3.5 %p of substrate-delivered power between n_TCO = 1.8 and 2.0.
5. **Thick transport layers are limited by transport and driving voltage, not by injection.**
   Injection happens at the interfaces and does not care about bulk thickness.
6. **1/(2n^2) is the emitter-to-air estimate, not the substrate-to-air one.** For light already
   isotropic inside the substrate the escape-cone fraction scales as 1/n_sub^2. The paragraph
   defines eta_ext immediately before, so as written it was low by a factor of two.
7. **Electrical loss removed from the eta_sub definition.** eta_sub is an optical power budget;
   charge balance and PLQY belong to eta_int.
8. **Absorption scales with the imaginary permittivity.** "2*n*k" now reads as eps_2 = 2nk and
   "흡수" as "소산", which is what makes one figure of merit cover both the TCO (k varies, n does
   not) and the thin Ag (n varies, k does not).

## Corrected — claim strength

9. **"80 % 이상의 초고효율을 갖는 OLED를 구현하였다" in the introduction.** The measurements are
   55.8 % and 77 %; above 80 % is a prediction. Reworded to state the measured 77 % and present
   80 % as the route the design rule opens.
10. **Abstract conflated two devices.** The 90 %-plus extraction is the Ag device on an ordinary
    substrate (91.6 % measured); the 77 % EQE is the high-index-substrate device, whose
    extraction is 87.5 %. As written it could read as one device with both numbers.
11. **The DBR is an auxiliary rear reflector, not a replacement electrode.** It was evaporated
    over the back of the finished device including the area around the Al cathode. CONFIRM this
    against the device record before submission.

## Also changed

- Symbols: `h_sub`, `h_ext`, `h_EQE` were plain "h" with no Symbol font, so they rendered as h.
  Changed to eta.
- The author's inline query "(광자수?power분율?)" resolved to photon-number fraction, which is
  what the simulation code computes.
- Typos: 일반적인유기, 기반의OLED, "는 더 joule", 구현에 가능하다, trade-off사라지고, Yablonobich.

## Flagged, not changed — needs the author's judgement

- **"ETL의 두께를 100 nm 이상으로 형성하는 것만으로도 충분한 고효율 달성이 가능"** (Fig. 5a
  paragraph). In the generic stack, 100 nm transport layers leave 22 % in the plasmon and only
  0.74 of substrate-delivered power; 200 nm is where it saturates at 0.92. If the 100 nm figure
  comes from the specific SpiroAC-TRZ stack, say so; otherwise raise the number.
- **The 500 nm organic thickness for sub-10 % plasmon loss** is consistent with the generic model
  (about 420 nm total), so it was left as is.
- The stray note at the end of the file about the Yablonovitch escape limit is still a to-do.
- Equation (1) is referenced but never written out.

## v3 — full rewrite for flow (`manuscript/unityEQE_v3_ko.docx`)

Requested after the v2 edit still read as repetitive. Built by `rewrite_v3.py` on the
original document's skeleton (same styles, namespaces and section properties; body runs
carry only the east-Asian font hint, as in the original). 38 paragraphs against 50, and
11.2 k characters against 13.8 k in the edited v2, with no claim, figure reference,
supplementary pointer or author placeholder dropped.

Duplications removed:

- "organics are transparent, so parasitic absorption is TCO + metal" was stated three
  times (round-trip, loss-accounting and electrode paragraphs); now once, in the
  loss-accounting paragraph.
- the lambda/4, lambda/2 ITO-thickness rationale appeared twice; now once, in the
  parasitic-absorption paragraph, with a back-reference from the electrode paragraph.
- the TCO absorption-versus-conductivity trade-off appeared twice with a "to be discussed
  later" forward reference; now once, where the 50 nm ITO / 10 nm Ag choice is made.
- "light bounces many times, so a small loss per round trip matters" was in both the
  Fig. 1(a) and the round-trip paragraphs; now once.
- the 20 % Al round-trip loss was quoted twice; the Fig. 1(d) sentence now refers back.
- the two Fig. 5 paragraphs both opened with the "light had to be crammed into the escape
  cone" framing; merged into one opening.
- the Discussion no longer re-lists the results numbers except the two headline ones.

Structure: the Results section is split into six short sub-headed blocks (recycling law,
parasitic absorption, substrate-delivered power, electrode rule, experiment, design
freedom). The sub-headings are plain paragraphs and can be deleted without side effects.
The Fig. 2(a) sentence now carries the two-step wording (comparable to the mirror loss
with Al; the only remaining loss with a low-loss reflector) supported by
`sim/design_rule4/tco2_*.csv`. The stray Yablonovitch note at the end was dropped.

All v2 corrections are carried over. Still flagged, unchanged: the "ETL above 100 nm is
enough" figure in the design-freedom paragraph.
