# English draft: sentences changed from v1 to v2

Source of each change: (P) the author's edits in the Word-exported PDF of v1; (T) Table 1 computed over the full spectrum; (M) Methods addition.

| # | Where | v1 | v2 | Source |
|---|---|---|---|---|
| 1 | Title | …: synergetic photon management | …: Synergetic photon management | P |
| 2 | Front matter | [[corresponding author(s)]] | Correspondence: Seunghyup Yoo (syoo.ee@kaist.edu), Malte C. Gather (malte.gather@uni-koeln.de) | P |
| 3 | Abstract, last sentence | … EQEs in the high 80s with isotropic emitters and above 90% when combined with a lower-absorption transparent electrode … | … and 90% when combined with a lower-absorption transparent electrode … | T |
| 4 | Introduction, last paragraph | … in the high 80s with an isotropic emitter and above 90% with a low-absorption transparent electrode … | … and 90% with a low-absorption transparent electrode … | T |
| 5 | Results, "Quantifying the parasitic absorption", 2nd paragraph | (Supplementary Fig. 2: BSDF_T and BSDF_R versus the refractive index of the outcoupling structure) | (Supplementary Fig. 2: BTDF and BRDF versus …) | P |
| 6 | Results, "Design freedom", 2nd paragraph | The same single-wavelength (550 nm) calculation for a real material stack … gives 88–89% with an isotropic emitter (ETL 200–300 nm; 4–5 percentage points above the generic stack), essentially the same 88–89% with the emitters of Θ ≈ 0.8–0.9 …, and 90–91% when the k of the ITO is lowered from 0.0032 to 0.002 (Table 1) [[to be confirmed by a full-spectrum, multi-wavelength calculation]]. | The same calculation for a real material stack … gives 87–89% with an isotropic emitter (ETL 200–300 nm; 3–5 percentage points above the generic stack), essentially the same 87–89% with the emitters of Θ ≈ 0.8–0.9 …, and 89–91% when the k of the ITO is lowered from 0.0032 to 0.002 (Table 1; values over the full emission spectrum and at 550 nm). | T |
| 7 | Discussion, 2nd sentence | … in the high 80s with an isotropic emitter and above 90% with a lower-absorption transparent electrode … | … and 90% with a lower-absorption transparent electrode … | T |
| 8 | Methods, Fig. 5a,b conditions | The wavelengths span 530–580 nm in 1-nm steps, … | The wavelengths span 400–700 nm in 1-nm steps, … | P |
| 9 | Methods, optical constants | (… n = 1.864 and k = 0.0032 at 550 nm, n = 1.83–1.88 over 530–580 nm) | (… n = 1.864 and k = 0.0032 at 550 nm, n = 1.79–2.08 over 400–700 nm) | P (range recomputed) |
| 10 | Methods, Table 1 calculation | The real-stack calculation of Table 1 obtains P_sub and R_LED at 550 nm from the measured optical constants of each layer (…) and sums the series with the BSDF of the same hemispherical microlens array (n_MLA = 1.8) as in Fig. 5. | The real-stack calculation of Table 1 obtains P_sub and R_LED from the measured, dispersive optical constants of each layer (…) and sums the series with the BSDF …; both the single-wavelength value at 550 nm and the value weighted over the 400–700 nm emission spectrum on a photon-number basis are given. | T |
| 11 | Methods, new sentence after 10 | — | The Python implementation used for Figs 2 and 3, Fig. 5c,d and Table 1 and the MATLAB implementation used for Fig. 4 and Fig. 5a,b,e agree for the stack of Table 1 to within 0.2 percentage points in η_sub and 1 percentage point in η_ext. | M |
| 12 | Figure 4 caption | (maximum EQE 20% → 55.8%) … (reference device 23%, microlens-array substrate + Al 61%, microlens-array substrate + Al + DBR 77%) | (maximum EQE 20% to 55.8%) … (reference device 23%, microlens-array substrate with Al reflector 61%, microlens-array substrate with Al electrode and DBR reflector 77%) | P |
| 13 | Figure 5 caption, last sentence | (c) and (d) … with PLQY = 1 and Θ = 2/3 …; (a) and (b) with the 530–580 nm emission spectrum, PLQY = 1 and Θ = 2/3; (e) with the 400–800 nm emission spectrum, PLQY = 0.98 and Θ = 0.76 (Methods). | (c) and (d) … with PLQY = 1, Θ = 2/3 …; (a) and (b) with the 400–700 nm emission spectrum, PLQY = 1, Θ = 2/3; (e) with the 400–800 nm emission spectrum, PLQY = 0.98, Θ = 0.76 (Methods). | P |
| 14 | Table 1 caption, optical constants | Measured optical constants (B3PyMPM n_o/n_e = 1.821/1.609, TCTA:B3PyMPM 1.833/1.671, TAPC 1.691/1.664); … | Measured optical constants (B3PyMPM: n_o = 1.821, n_e = 1.609; TCTA:B3PyMPM co-host: n_o = 1.833, n_e = 1.671; TAPC: n_o = 1.691, n_e = 1.664); … | P |
| 15 | Table 1 caption, entries | Each entry is the single-wavelength value at 550 nm; the full-spectrum value (530–580 nm, Ir(ppy)₂acac) will be given in the same cell as "full spectrum / 550 nm" [[…]]. | Each entry reads "full emission spectrum (400–700 nm, Ir(ppy)₂acac, photon-number weighting) / single wavelength at 550 nm"; the full-spectrum values are 0.6–0.9 percentage points below the 550-nm values because the ITO and Ag absorb more at shorter wavelengths. A′ is the 550-nm value. | T |
| 16 | Table 1 body | single 550-nm values | "full / 550 nm" pairs in η_sub, SPP loss, η_ext and EQE (A′ unchanged) | T |
| 17 | End of document | "Items still open" list | removed (the PDF export had dropped it) | P |

Not changed on purpose: the Word equation objects and italic maths symbols of the PDF are formatting; the markdown source keeps plain-text symbols (η_sub, R_LED, …).

# Author's PDF (Word export of v1 with edits) versus v2

Content differences only; equation objects and italic symbols in the PDF are formatting and are not listed.

**Changed in v2 after the PDF (Table 1 over the full spectrum, Methods)**

| # | Where | PDF | v2 |
|---|---|---|---|
| A1 | Abstract, last sentence | … and above 90% when combined with a lower-absorption transparent electrode | … and 90% when combined with … |
| A2 | Introduction, last paragraph | … and above 90% with a low-absorption transparent electrode … | … and 90% with … |
| A3 | Design freedom, 2nd paragraph | The same single-wavelength (550 nm) calculation … gives 88–89% … (4–5 percentage points …), … 88–89% …, and 90–91% … (Table 1) [[to be confirmed …]] | The same calculation … gives 87–89% … (3–5 percentage points …), … 87–89% …, and 89–91% … (Table 1; values over the full emission spectrum and at 550 nm) |
| A4 | Discussion, 2nd sentence | … and above 90% with a lower-absorption transparent electrode … | … and 90% with … |
| A5 | Methods, optical constants | n = 1.83–1.88 over 400–700 nm (the PDF changed the range to 400–700 nm but kept the 530–580 nm index span) | n = 1.79–2.08 over 400–700 nm |
| A6 | Methods, Table 1 calculation | … obtains P_sub and R_LED at 550 nm from the measured optical constants … as in Fig. 5. | … from the measured, dispersive optical constants … as in Fig. 5; both the single-wavelength value at 550 nm and the value weighted over the 400–700 nm emission spectrum on a photon-number basis are given. + new sentence: The Python implementation used for Figs 2 and 3, Fig. 5c,d and Table 1 and the MATLAB implementation used for Fig. 4 and Fig. 5a,b,e agree for the stack of Table 1 to within 0.2 percentage points in η_sub and 1 percentage point in η_ext. |
| A7 | Table 1 caption, entries | Each entry is the single-wavelength value at 550 nm; the full-spectrum value (400–700 nm, Ir(ppy)₂acac) will be given … [[…]] | Each entry reads "full emission spectrum (400–700 nm, Ir(ppy)₂acac, photon-number weighting) / single wavelength at 550 nm"; the full-spectrum values are 0.6–0.9 percentage points below the 550-nm values because the ITO and Ag absorb more at shorter wavelengths. A′ is the 550-nm value. |
| A8 | Table 1 body | single 550-nm values | "full / 550 nm" pairs (A′ unchanged) |

**PDF edits carried into v2 with a correction**

| # | Where | PDF | v2 |
|---|---|---|---|
| B1 | Methods, Fig. 5a,b conditions | … and an normalized in-plane wavevector grid of 1,000 points up to u = 3 | … and a normalised in-plane wavevector grid (u = k_∥/(n_EML k₀)) of 1,000 points per unit up to u = 3 (the grid has 1,000 points per unit of u, 3,000 in total) |

**Points in the PDF (Word file) to check on the author's side**

| # | Where | PDF | Note |
|---|---|---|---|
| C1 | Definition of η_ext; eq. (3) paragraph | "this fraction scales as 1/η_sub²", "p is of order 1/η_sub² [15]", "(η_sub = 1.5)" | these three should be n_sub (refractive index), not η_sub; v2 has n_sub |
| C2 | Transparent-electrode paragraph | ε_e = 2nk | v2 has ε₂ = 2nk (imaginary part of the permittivity) |
| C3 | Methods, matrix form | "Writing the power P₀ delivered by the OLED to the substrate and angular distribution …" | "its" was lost; v2 keeps "and its angular distribution" |
