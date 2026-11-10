# DRAFT letters to editors — TE-OLED "EQE approaching 60%" concern
# Status: DRAFT for user review. Send from institutional email.
# Recipients: Angew. Chem. Int. Ed. editorial office (angewandte@wiley-vch.de)
#             Advanced Materials editorial office (advmat@wiley-vch.de)
# Before sending: (1) discuss with advisor, (2) fill in [NAME/AFFILIATION],
# (3) verify both DOIs resolve to the versions cited here.

---

## LETTER 1 — to Angewandte Chemie Int. Ed. (re: anie.2518285)

Subject: Concern regarding the headline efficiency claim in anie.2518285
("Top-Emitting OLED Realizes EQE Approaching 60%...")

Dear Editors,

I am writing to raise a documented methodological concern about the
headline claim of the following article, and to request that it be
evaluated editorially:

  J. Liu et al., "Triphenylene-Involved pi-Extension Combining With
  Phenyl-Blocking Enhances the Stability of MR-TADF Emitter: Top-Emitting
  OLED Realizes EQE Approaching 60% With BT.2020 Green Gamut and Long
  Lifetime", Angew. Chem. Int. Ed. 2026, 65, e2518285.

The concern, in brief: the external quantum efficiency of 59.2% claimed
in the title for the top-emitting (TE) device is not supported by the
measurement described in the paper, and the discrepancy is quantifiable
from the paper's own numbers.

1. Instrument asymmetry, disclosed in the Supporting Information. The
   bottom-emitting devices were measured with an absolute (integrating-
   sphere) EQE system (Hamamatsu C9920-12). The TE device — the one
   carrying the title claim — was instead measured with a forward-
   direction I-V-L tester ("F-Star optical measurement systems,
   FS-1000GA3-IVL"). No angle-resolved electroluminescence, angular
   intensity distribution, or stated angular assumption appears anywhere
   in the main text or Supporting Information (Figs. S1-S39).

2. Why this matters for this specific device. The TE device is a strong
   microcavity (Ag 140 nm mirror anode; Yb/Ag 12 nm semi-transparent
   cathode; EL FWHM narrowed to 19 nm). Microcavity emission is
   forward-peaked, i.e., strongly non-Lambertian. Converting a forward-
   only luminance measurement to EQE requires an angular emission
   profile; assuming a Lambertian profile systematically inflates the
   EQE of forward-peaked devices. This is a long-established point in
   the field's measurement standards (S. R. Forrest, D. D. C. Bradley,
   M. E. Thompson, Adv. Mater. 2003, 15, 1043; quantified for modern
   devices by D. G. Archer et al., Adv. Optical Mater. 2021, 9, 2000838).

3. The paper's own numbers show the Lambertian assumption was used.
   For emission at 520 nm, a Lambertian emitter obeys
   CE/EQE = E_photon x K(lambda)/pi = 3.68 cd/A per % EQE.
   The reported TE values give 218.0 cd/A / 59.2% = 3.682 — exactly the
   Lambertian ratio. For a forward-peaked cavity with independently
   measured CE and EQE, this ratio must exceed the Lambertian value.
   Its exact coincidence demonstrates that EQE and CE are not two
   independent measurements: one is an arithmetic Lambertian derivative
   of the other. For a device whose emission profile is manifestly
   non-Lambertian, the derived member of that pair — here the headline
   EQE — is unsupported.

4. Internal comparison. The same group's integrating-sphere values in
   this and a companion paper (bottom-emitting: 33.8%, 40.1%, 41.9%;
   the latter two in Adv. Mater. 2026, e72684, supported there by
   measured horizontal dipole ratios) are consistent with established
   physics. The only value requiring a switch away from the sphere is
   the 59.2% in this title.

I want to be explicit that I am not alleging data fabrication: the
instruments are disclosed in the SI, and the practice of forward-only
TE characterization is regrettably widespread. My concern is that a
quantitative headline claim — in the article title — rests on a
conversion assumption that is invalid for the device class in question,
and that this is not disclosed to readers.

I would therefore ask the editors to invite the authors to:
  (a) state explicitly how the TE EQE was obtained, including the
      angular assumption used;
  (b) provide angle-resolved EL intensity and spectra for the TE device,
      from which the angle-integrated EQE follows without assumption; and
  (c) correct the reported EQE and the title if the angle-integrated
      value does not support "approaching 60%".

I am happy to provide a detailed version of the arithmetic in point 3
or any further clarification.

Sincerely,
[NAME]
[AFFILIATION, POSITION]
[ORCID / institutional email]

---

## LETTER 2 — to Advanced Materials (re: adma.74082)

Subject: Undisclosed measurement methodology behind the title claim of
adma.74082 ("...Top-Emitting OLEDs with EQE Approaching 60%...")

Dear Editors,

I am writing to raise a documented concern about the following article:

  Z. Xue et al., "Ultra-Narrow Pure-Green MR-TADF Emitter with an FWHM
  of 12 nm Enables Superior-Performance Top-Emitting OLEDs with EQE
  Approaching 60%, Power Efficiency over 300 lm W-1, and CIE Coordinates
  of (0.14, 0.79)", Adv. Mater. 2026, e74082.

The title claims EQE approaching 60% (58.9% reported) for a top-emitting
(TE) microcavity OLED. My concern has two parts.

1. The measurement behind the title claim is not disclosed. The SI's
   "Device Fabrication and Measurement" paragraph describes only an
   integrating-sphere system (Hamamatsu C9920-12), in wording that
   reproduces, essentially verbatim, the measurement paragraph of the
   group's bottom-emitting-only paper (Adv. Mater. 2026, e72684). No
   instrument or procedure for the TE device is named anywhere. For
   comparison, the same laboratory's contemporaneous TE paper in Angew.
   Chem. Int. Ed. (2026, e2518285) explicitly states that its TE device
   was measured with a forward-direction I-V-L tester (F-Star
   FS-1000GA3-IVL), not the sphere. Readers of the present article
   cannot determine which applies here.

2. The reported numbers show that EQE and CE are not independent
   measurements. For emission at 521 nm a Lambertian emitter obeys
   CE/EQE ~ 3.8 cd/A per % EQE. The reported TE values give
   232.5 / 58.9 = 3.95 — on the Lambertian line. A strong microcavity
   device is forward-peaked; if CE (forward) and EQE (angle-integrated)
   had been measured independently, this ratio would necessarily exceed
   the Lambertian value. Its coincidence shows one of the two headline
   numbers is an arithmetic Lambertian derivative of the other. Either
   the EQE was inferred from forward luminance under a Lambertian
   assumption (inflating it for this device class), or the CE and the
   power efficiency "over 300 lm W-1" were back-calculated from a sphere
   EQE under the same invalid assumption. Under either reading, the
   title's combination {EQE ~60%, PE >300 lm W-1} is not supported
   without angle-resolved data — of which the paper contains none
   (the angle-resolved PL of Fig. S25 characterizes film dipole
   orientation, not device electroluminescence).

The measurement standards of this field (Forrest, Bradley, Thompson,
Adv. Mater. 2003, 15, 1043 — published in this journal — and Archer
et al., Adv. Optical Mater. 2021, 9, 2000838) are explicit that
non-Lambertian devices require angle-resolved characterization for a
valid EQE. I am not alleging fabrication; I am stating that a
quantitative title claim rests on an undisclosed measurement and an
assumption invalid for the device class.

I ask the editors to invite the authors to (a) disclose the TE
measurement instrument and conversion assumptions, (b) provide
angle-resolved EL data for the TE device, and (c) correct the EQE, PE,
and title if the angle-integrated values do not support them.

Sincerely,
[NAME]
[AFFILIATION, POSITION]
[ORCID / institutional email]

---

## Sending notes (for user)

- Send each letter separately to its journal. Do not CC the authors;
  the editorial office will contact them.
- Use your institutional address; anonymous concerns carry less weight.
  Discuss with your advisor first — the Yang group is large and this is
  your immediate field.
- Realistic outcomes: author response with angular data (concern
  resolved), an author correction/erratum, or an editorial "no action".
  Keep the correspondence; if no substantive response in ~4-6 weeks, a
  PubPeer post with the same facts is the standard next step.
- Everything asserted above is verifiable from the archived primary
  sources in data/ (both PDFs + both SIs) and EQE60_CLAIM_EVIDENCE.md.
