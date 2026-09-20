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

## v8 — a one-line basis for p (`manuscript/unityEQE_v8_ko.docx`)

The Fig. 1(b) paragraph gains one sentence: a round trip is one encounter with the extraction
structure, and the escape probability there is roughly 1/n_sub² for angularly randomised light
inside the substrate, about 0.4 for ordinary glass, citing Yablonovitch (1982). That is all
the text needs — it says why the panel is drawn at p = 0.4 without turning into a derivation.

Two points kept in `figure_plan.md` rather than the text. Yablonovitch quotes 1/(2n²) because
his internal intensity is two-sided, which suits a slab with a rear mirror; each round trip
here involves a single encounter with the extraction structure, so the one-sided 1/n_sub²
applies and glass gives 0.4 rather than 0.22. And eq. (2) rests on his ergodicity condition,
a surface slope above ½ arcsin(1/n_sub), which a hemispherical microlens satisfies — worth
having to hand if a referee asks why every round trip carries the same p.

No figure is cited at this point, so the citation order is undisturbed, and the ray-traced
escape probability of Fig. 2(b) is left to speak for itself when that panel arrives.

## v9 — the author's Fig. 2(a) paragraph, and the index sweep split in two (`manuscript/unityEQE_v9_ko.docx`)

The parasitic-absorption paragraph is replaced with the author's own rewrite. Three things
change in substance. The extraction-efficiency sentence no longer says the drop is steeper
for Ag than for Al *in passing* — it now says η_ext is higher with Ag yet still falls steeply
with TCO absorption, which is what the data show: over 30 → 200 nm of ITO at n_sub = 1.8,
η_ext goes 0.863 → 0.646 with Ag against 0.701 → 0.570 with Al. The low-loss mirror is the
more sensitive one, because dη_ext/dA′ = −p(1−p)/[p + (1−p)A′]² is steepest where A′ is
smallest. Second, the mechanism behind "the trend grows stronger with a high-index structure"
is now stated through eq. (2) rather than through the number of round trips: a smaller p gives
a lower η_ext at the same absorption, because more of the circulating light is lost before it
escapes. Third, the closing sentence — substrate-delivered power rises but extraction falls,
and with large parasitic absorption the second outruns the first — now cites Fig. 2(c).

One wording change beyond the author's text: "기판과 공기 사이의 굴절률 차이" became
"광추출 구조와 공기 사이", since the escape happens at the structure/air boundary and the
structure is index-matched to the substrate. Spacing normalised in two places.

Figure consequences. The substrate-index sweep is now two small square panels sharing an x
axis: (b) p and η_ext, (c) η_sub^(0) and EQE. The absorption maps move from (c) to (d). The
Poynting-vector panel leaves the main text for the SI — it is the evidence behind the
mirror/TCO split that 2(a) already draws, so as a fourth main-text panel it repeated (a); the
text keeps it as one sentence with an SI citation and notes that it agrees with 2(a).
The citation order stays strictly monotonic.

## v10 — Koenig ITO constants, and what they do to Fig. 2 (`manuscript/unityEQE_v10_ko.docx`)

The author supplied the Koenig et al. (2014) ITO tabulation, so the TCO is now
n = 1.8636 + 0.0032285i at 550 nm instead of the flat n = 1.9 + 0.02i the design-rule scripts
had been set up with. It matches the `l_ITO` entry of the project library to four decimals,
so the library already held this dataset. Every Fig. 2 panel was recomputed.

The consequence is not cosmetic. A transparent ITO absorbs six times less per pass, so the
Al stack now loses essentially everything at the mirror: 13.4 % per round trip against a TCO
contribution of 6 % of A' at 50 nm and 15 % at 150 nm. The sentence "통상적인 Al 반사판에서는
TCO 흡수가 거울 손실과 맞먹는 수준까지 커지고" is no longer true — parity would need
k >= 0.015 — and has been replaced by the two-step reading the new numbers give: with Al the
mirror is the whole loss and ITO thickness hardly matters (η_ext −1.1 pp over 50 → 150 nm),
while dropping the mirror loss to 1.8 % with Ag makes the same film the main loss path
(34 % → 60 % of A') and doubles the thickness penalty (−2.2 pp). The Ag half of the claim was
never at risk: the TCO dominates the round trip for every k anyone would deposit, from the
library's cleanest film to the ETRI one.

Two things improve. The model now agrees with the measurement: 150 nm ITO on glass with an Ag
reflector gives η_ext = 0.934 against the 0.916 measured on the green device, where the
k = 0.02 assumption gave 0.788. And the Al ceiling in Fig. 2(c) moves from 0.54 to 0.63,
landing just above the 60 % where the field actually stalled — which is the comparison that
section is making.

Methods now names the optical constants the design-rule calculations rest on (McPeak Ag,
Johnson–Christy-type Al, Koenig ITO, non-absorbing n = 1.8 organics); the remaining bracket
there is only the per-emitter PLQY and orientation factor.

Still to check: Fig. 2(d)'s quoted round-trip absorption ranges (Al 15–25 %, Ag 5–15 %) come
from the earlier slide, not from data in this repository. With the Koenig ITO the generic
stack gives 14–16 % for Al and 2.4–5.3 % for Ag at n_sub = 1.5, so the Ag range in the text
is probably too high and should be recomputed on whatever electrode structures that panel
ends up showing.

## v11 — Fig. 2(d) recomputed (`manuscript/unityEQE_v11_ko.docx`)

The quoted ranges (Al 15–25 %, Ag 5–15 %) came from an old slide and had no data behind them
in this repository. Five electrode structures are now settled and computed with a dispersive
TMM on the generic Fig. 2 stack: Al / ITO 150 nm, Ag / ITO 150 nm, Ag / IZO 50 nm,
Ag / Ag 10 nm and a metal-free ZnS/LiF DBR with an IZO cathode. Flux- and spectrum-weighted
1 − R at n_sub = 1.5: 15.6, 4.5, 3.2, 5.6 and 11.5 %. Al was about right; the Ag range was
roughly twice too high and the text now gives a number per structure.

The DBR sentence gains the new result: its absorption, 3.0 %, matches the best metal, but
outside the stopband it transmits, and that leak occurs only for angles inside the escape cone
to air, so the same stack loses 8.5 % on glass and 3.3 % at n_sub = 1.8.

Not changed, and needing the author: the orange device's rear DBR cannot affect the round trip
through the stack, since 100 nm of Al transmits 1e-7 at 550 nm. The 60 % → 77 % it produced
is light that bypasses the electrode, which is what the back-leakage measurement says; the
text attributes it to joule dissipation at the Al.

## v12 — Fig. 2 at a fixed 150 nm ITO, and relettered a–j (`manuscript/unityEQE_v12_ko.docx`)

Two changes, both the author's.

The electrode panel now fixes the transparent electrode at 150 nm of ITO and varies only the
reflector, so it compares Al, Ag and a ten-pair ZnS/LiF quarter-wave stack at 550 nm: 15.6,
4.5 and 8.8 % of flux- and spectrum-weighted round-trip loss on glass. The thin-Ag and IZO
variants leave the main text. The DBR sentence gains the result that came out of the
recalculation: a quarter-wave stack designed at the emission wavelength is not the right
design for randomised light, because at oblique incidence the stopband moves to shorter
wavelengths and the two materials stop sharing a quarter-wave condition. Re-optimising the two
thicknesses for the angular and spectral average brings the same ten pairs to 5.5 % on glass,
and the orange device's own ZnS 70 / LiF 115 nm turns out to be that optimum to 0.01 %p at
n_sub = 1.8 — which is what the genetic-algorithm design in Methods was for, and is worth
saying, since the λ/4-at-550 design would have lost 12.8 % there.

Fig. 2 is now composed as one full-width figure of ten panels in three rows, so the citations
become ranges: 2(a)–(c) for the ITO thickness, 2(d),(e) for the index sweep, 2(f) for the
product, 2(g)–(j) for the reflectors. The order stays strictly monotonic. The middle row was
one panel and is now three; the extra one shows that the escape probability falls by 2.7×
across the index range while A′ barely moves, which is the cleanest statement of why η_ext
declines with substrate index — it is the escape probability, not a worsening mirror.
