# Evidence file: "TE-OLED EQE approaching 60%" methodology question

Target paper: Liu et al., Angew. Chem. Int. Ed. 2026, 65, e2518285
(BN-Tpl-Ph, TE-OLED EQE_max 59.2%, claim is in the TITLE).
Comparison paper (same group): Xue et al., Adv. Mater. 2026, 38, e72684
(p/m-DBFSi, bottom-emitting only, EQE_max 40.1-41.9%).
Both PDFs + SIs archived in data/ (Liu_AngewChem_*, Xue_AdvMater_*).

## Verified facts (primary sources, read directly)

1. Liu Angew SI, "Device Fabrication and Characterization":
   - BE-OLEDs: Keithley 2400 + "absolute EQE measurement system
     (C9920-12, Hamamatsu Photonics)" -> integrating sphere, absolute.
   - TE-OLEDs: "F-Star optical measurement systems (FS-1000GA3-IVL)"
     -> forward-direction IVL tester; no sphere.
   - The ENTIRE SI (Figs S1-S39, Tables S1-S5) contains zero
     angle-resolved EL spectra, zero angular intensity data, and the
     words "Lambertian/angular/goniometer" never appear.
   - CE/EQE ratio is 3.68 cd/A per %EQE for the TE device vs 3.77 for
     BE -> the same Lambertian conversion was applied to a strongly
     forward-peaked microcavity device (FWHM 19 nm = strong cavity).

2. Xue Adv Mater 2026 (same group, three months earlier):
   - Bottom-emitting ONLY; no TE device anywhere in paper or SI.
   - Measured with the SAME sphere system (C9920-12); SI wording is
     nearly verbatim identical to the Liu SI's BE sentence.
   - They additionally measured horizontal dipole ratios (Fig S31) to
     justify the high (40-42%) sphere-measured EQE. So the group knows
     how to support an unusual efficiency claim with optics data when
     it has that data.

## Reading

Within this group's own practice, BE numbers (33.8%, 40.1%, 41.9%) are
integrating-sphere absolute values. The only number that required a
switch to a forward-only instrument is the 59.2% headline. The 1.75x
TE/BE jump matches what Lambertian conversion of a forward-peaked
cavity produces; genuine sphere-verified TE enhancements are typically
1.2-1.4x. Angle-integrated truth likely mid-40s %.

## Anchors for a formal comment

- Forrest, Bradley, Thompson, Adv. Mater. 2003, 15, 1043: EQE of
  non-Lambertian devices must be measured angle-resolved or in a sphere.
- Archer et al., Adv. Optical Mater. 2021, 9, 2000838: quantifies the
  error of forward-only EQE for non-Lambertian OLEDs.

## Correction log

- An earlier web-search summary claimed Xue Adv Mater 2026 also reports
  a "TE-OLED EQE approaching 60%". WRONG - primary source shows BE-only,
  41.9% max. Corrected 2026-08-14 after reading the actual PDF+SI.

## Status / open items

- Not yet checked: Liu et al. Adv. Mater. 2025 (adma.202411610,
  deuteration green) and other group TE papers - do they also use
  F-Star forward-only for TE? (User to supply SIs; Wiley egress-blocked
  from this environment.)
- Recommended escalation order: (1) email corresponding authors asking
  for angular EL data + conversion assumption; (2) PubPeer with the
  facts above; (3) formal Comment to Angew. Frame as methodology, not
  misconduct: instruments were disclosed, so fabrication is not on the
  table.
