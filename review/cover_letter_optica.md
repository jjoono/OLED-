# Cover letter (draft for Optica)

[Date]

Dear Editor,

We submit our manuscript, "Reciprocal Angular Transport Explains the Limited Benefit of Mixing in OLED Light Extraction," for consideration as a Research Article in *Optica*.

**The question.** Randomizing the direction of trapped light is the standard route to better light extraction. It is the assumption behind statistical ray optics and behind microlens arrays, scattering films, and textured surfaces. Yet, as we show, strongly mixing structures gain little over simple periodic lenses. We ask how much angular mixing can improve extraction at all, and answer it from optical reciprocity.

**The main result.** We treat repeated encounters between an emitter and a passive extraction structure as a reciprocal Markov chain on polar angle. In a fully mixing, reciprocal return model, recirculating light escapes with probability p − σ²/(1 − p), where σ² is the étendue-weighted variance of the angle-dependent escape probability. It does not escape with p, as the fully randomized model of statistical ray optics assumes, and the efficiency follows in closed form.

**What follows from it.**
- Across 37 scattering, microlens, randomized-lens, and textured structures, complete mixing would add no more than about two percentage points to the spectrally weighted extraction of the devices studied. It adds more (up to 9 points) only at wavelengths where the cavity sends light into angles the structure transmits poorly, and there randomized lenses and textures realize most of that gain. The fully randomized model predicts up to 19 points, and 41 points for a planar interface.
- The result explains why faster mixing gives limited gains: the same angular variance of escape that leaves light behind also limits what mixing can recover.
- We state its scope explicitly. It is a physically constrained reference for a rank-one reciprocal model, not a universal bound; the Supplement gives a reciprocal counterexample.

**Supporting results.**
- The spectrum of the return matrix gives two descriptors: an angular mixing number and a late-pass escape ratio.
- A reduced model built on them reproduces full transport calculations with a mean absolute error of 0.30 percentage points over 342 cases.
- We test it independently. Direct ray tracing in (θ, φ) agrees with the full series to within 0.6 points. Out-of-sample structures are predicted for seven of nine families. The model fails for cone arrays, and we trace that failure to slow angular modes that the escaping light does not occupy.
- A controlled freeform-versus-hemisphere benchmark, with independent searches and a fixed-cavity analysis, shows that lens-shape freedom changes mainly the mixing rate and the angular composition, not the escape probability.

**Why Optica.** The constraint that reciprocity places on angular redistribution is not specific to OLEDs. It applies to any passive structure that recycles light between a source or absorber and an escape interface, including LEDs, light-trapping textures, and scintillators. The closed-form reference and the spectral descriptors give a common, quantitative way to compare such structures. We believe this conceptual result, together with its explicit limits, will interest readers across optical design, photonic devices, and light transport.

**Related work by the authors.** The device model and the LightTools pipeline used here were also used in a related manuscript, which is currently under review at another journal and available as a preprint (ref. 14, arXiv:2610.12065). That work asks why OLED efficiency is limited and how minimizing the round-trip loss of returned light approaches unity quantum efficiency. The present manuscript asks a different question: how much angular redistribution can contribute, and what reciprocity implies for it. [Authors: verify against the submitted manuscript.] None of the results on complete reciprocal mixing, the spectral descriptors, the reduced model, or the validation tests appear in that manuscript. Shared inputs (stack definitions, optical constants, and the round-trip-loss formalism) are cited.

This manuscript has not been published and is not under consideration elsewhere. All authors have approved the submission. Data and code will be made available in a public repository [DOI].

Suggested reviewers: [names, affiliations, e-mail; expertise in OLED outcoupling / light trapping / radiative transport]

Opposed reviewers (if any): [ ]

Sincerely,

[Corresponding author, affiliation, e-mail]
