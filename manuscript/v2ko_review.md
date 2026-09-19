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

## v4 — Methods and back matter added (`manuscript/unityEQE_v4_ko.docx`)

Adds the author's Methods (optical simulation; OLED fabrication and characterisation) in the
same register, with `[[ ]]` placeholders where a Nature-family referee will ask: film PLQY and
orientation factor per emitter and the source of the optical constants; thickness and
deposition of the top Ag reflector of the green device; whether the DBR covers the area
around the Al cathode; angular range and EQE integration method, device count and spread.
Back-matter headings (Acknowledgement, Conflicts of Interest, Supporting Information, Data
Availability Statement, Author contributions) with placeholders.

Body tweaks: the two experiments are now named green (glass + MLA film + Ag) and orange
(high-index MLA substrate + Al + auxiliary DBR) to tie the Results to the Methods; the
Fig. 2(a) sentence now describes the figure as built in `figures/fig2a_mock/` (round-trip
loss split into mirror ohmic and TCO parts, Al versus Ag, with eta_ext).

Suggested Fig. 2(a) legend: "ITO 두께에 따른 왕복 손실 A'와 기판→공기 추출 효율. A'는 기판에서
본 OLED 스택의 반사율 결손으로, 회색은 거울의 ohmic 손실(k_TCO = 0으로 재계산), 주황은 TCO 흡수
(나머지)이다. Al 반사판(왼쪽)에서는 TCO 흡수가 거울 손실과 맞먹는 수준까지 커지고, Ag 반사판
(가운데)에서는 사실상 유일한 손실이 된다. 오른쪽: 그에 따른 η_ext. 유리 기판(n = 1.5) + microlens
film, p = 0.38, 550 nm."

## v5 — figure order and section swap (`manuscript/unityEQE_v5_ko.docx`)

Four changes, all at the author's direction; `manuscript/figure_plan.md` holds the resulting
per-panel specification.

1. Fig. 1 has three panels. The former 1(d) becomes 1(c) and now carries only the decay curve
   against round-trip number; the sentence no longer also promises an extraction-versus-
   reflectance plot. The citation of the non-existent 1(c) in the parasitic-absorption
   paragraph now points at 2(b), which already carries that content.
2. Sections 3 and 4 are swapped: the transparent-electrode rule now follows the
   parasitic-absorption diagnosis directly, and the substrate-delivered-power section follows.
   Fig. 3 panels renumbered old c → a, a → b, b → c. Two bridge sentences were added so the
   sections still hand off cleanly, and the design-freedom retrospective points at 3(a).
3. Fig. 3(b) is introduced once as covering both organic thickness and substrate index, so the
   second citation is no longer a surprise.
4. The Fig. 1(a) paragraph now has a sentence for the proposed half of the schematic.

Also: the thesis sentence about the trade-off disappearing moved to the end of the
parasitic-absorption section, so that section closes on the paper's own claim.

Citation order is now monotonic: 1(a) 1(b) 1(c) 2(a) 2(b) 2(c) 2(d) 3(a) 3(b) 3(c) 4(a) 4(b)
4(c) 5(a) 5(b) 5(c), with one deliberate backward reference to 3(a) in the design-freedom
section.

## v6 — Fig. 1 reduced to two panels (`manuscript/unityEQE_v6_ko.docx`)

The escape-cone panel is dropped: the ray-optics limit is textbook material, so the sentence
stays and the figure citation goes. The master curve (η_ext against R_LED = 1 − A′) takes the
(b) slot, which is what the round-trip paragraph argues in the first place, and the
decay-versus-round-trip-number curve is dropped as redundant with it. The round-trip paragraph
now reads off the master curve: an Al device with A′ > 0.2 sits on the flat part, suppressing
A′ to a few per cent moves it onto the steep part, and a smaller p shifts the whole curve down.

Figure conventions settled with the author: x axis is R_LED = 1 − A′ in both the wide view and
the zoom; p is stated on the panel; the horizontal "w/o outcoupling" and "w/ microcavity"
bands are removed because they are whole-device outcoupling values on an η_ext axis.

Also recorded in `figure_plan.md`: the blue curve of the Fig. 2(b) draft is the single-pass
escape probability p, confirmed by the author, which is the value the η_ext calculations in
`sim/design_rule4/` rest on.

## v7 — Fig. 1 back to three panels (`manuscript/unityEQE_v7_ko.docx`)

The decay curve returns as 1(b) and the master curve moves to 1(c). Dropping it in v6 was a
misjudgement: the two panels do not say the same thing. The decay curve shows the recycling
process, which is the paper's physical premise and is hidden inside the shape of the master
curve; and at this point in the paper the comparison should still be qualitative, since no
Al/Ag numbers have been introduced. The master curve, with a quantitative R_LED axis, is the
conclusion that follows. Section 1 now runs schematic → process → law.

The author's section-1 revisions are carried in: the round-trip paragraph no longer quotes the
20 % figure or cites Supplementary Fig. 2, since that is Fig. 2(c)'s result — it now says only
that the transport layers are transparent, so the loss is set by the mirror's ohmic absorption
and the transparent electrode; the reflective-electrode assumption T_LED = 0 is stated; and
the index-matching and "power remaining" sentences take the author's wording.
