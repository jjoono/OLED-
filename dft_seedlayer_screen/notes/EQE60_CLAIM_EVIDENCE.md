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

## Paper 3 checked (2026-08-14): Xue et al., Adv. Mater. 2026, adma.74082
("Ultra-Narrow ... FWHM of 12 nm ... TE OLEDs with EQE Approaching 60%",
DBFCN emitter; archived as Xue_AdvMater_12nm.pdf / _SI.docx)

- BE device: EQE_max 38.5%; TE device: EQE_max 58.9%, CE_max 232.5 cd/A,
  PE > 300 lm/W. The 60% claim is again in the TITLE.
- The SI "Device Fabrication and Measurement" paragraph mentions ONLY the
  Hamamatsu C9920-12 sphere - and it is verbatim boilerplate copied from
  the group's BE-only TPS paper (adma.72684), with no TE instrument
  named, unlike the Liu Angew SI which explicitly names the F-Star for TE.
  So this paper does not disclose how the TE device was measured at all.
- No angle-resolved EL anywhere. (Angle-resolved PL for dipole
  orientation, Fig. S25, is film PL - not device EL.)

## The Lambertian-ratio arithmetic (applies to both TE papers)

For emission at lambda with luminous efficacy K(lambda), a Lambertian
device obeys CE [cd/A] = EQE x E_photon[eV] x K / pi. At 520-521 nm
(V ~ 0.71, K ~ 485-490 lm/W) that is 3.68-3.8 cd/A per %EQE.

  Liu Angew TE:  218.0 / 59.2 = 3.682   (Lambertian prediction: 3.68)
  Xue 12nm TE:   232.5 / 58.9 = 3.947   (prediction at 521 nm: ~3.8)
  (BE reference: Liu binary BE 110.7 / 29.4 = 3.77 - same family.)

A strong-microcavity TE device is forward-peaked, so its true CE/EQE
ratio must sit ABOVE the Lambertian value if both quantities were
measured independently. Landing exactly ON the Lambertian line proves
that in each TE dataset only ONE of {EQE, CE} is an independent
measurement and the other is an arithmetic Lambertian derivative:
either forward luminance -> EQE (inflated EQE), or sphere EQE ->
CE/PE (fabricated-by-assumption CE and PE for a non-Lambertian
device). Under either reading, the headline pair {EQE ~60%, PE >300
lm/W} is not supported without angular data - which neither paper has.

## Optical simulation of the published stack (scripts/65_te_cavity_eqe.py)

CPS dipole-cavity model of the exact Liu-Angew TE stack (all thicknesses
as published; Ag J&C, Yb 1.1+2.6i, organics n=1.75; gamma = exciton
utilisation = 1, q_PL = 0.97; best dipole plane in the 35 nm EML):

  Horizontal fraction 0.85:  true EQE <= 35.8 %,  Lambertian-read 81.7 %
  Horizontal fraction 1.00:  true EQE <= 42.1 %,  Lambertian-read 96.8 %
  SPP/evanescent loss 12-15 %; metal/ITO absorption ~50 % of dipole power.

The claimed 59.2 % exceeds the perfect-orientation upper bound by 17
points -> not reachable as an angle-integrated EQE in this structure.
A forward-luminance x Lambertian reading of 59.2 % corresponds to a true
EQE in the mid-20s to ~30 % (model inflation factor ~2.3), i.e. BELOW the
group's own sphere-measured bottom-emitting 33.8 %. This is the third,
physics-based leg of the claim (independent of the instrument-disclosure
and CE/EQE-ratio legs), and it survives even if the authors assert the
sphere was used.

## Spacer-optimisation extension (scripts/66_te_cavity_optimize.py)

Freeing HTL/ETL/capping up to 300 nm each (same model, Theta_par 0.85):
optimum lands at 145/55/80 nm - within nanometres of the published
135/45/75 - giving true EQE 44.1 % (eta_out 45 %). The thick-spacer route
never wins: higher cavity orders lose more to per-pass metal absorption
than they gain from SPP suppression. So (a) the published DESIGN is
near-optimal and the claimed NUMBER is the anomaly, and (b) no spacer
configuration with these materials reaches 59.2 % as a true EQE (~50 %
even at perfect dipole orientation). Bonus finding: at the true-EQE
optimum the forward-Lambertian reading collapses to 26.3 % - the
forward-only metric actively mis-ranks TE designs (it rewards the
forward-peaked geometry that reads 81.7 % while delivering 35.8 %).

## Cross-check against the user's own MATLAB framework (scripts/67)

The user's lab dipole model (birefringent CPS, PSO-optimised red TE stack,
no Yb, capping n = 2.3, Theta_par 0.95, eta_rad 1) reports ~60 % true EQE
at ETL/HTL/cap = 229/252/129. Replicated here at 54 % - the two
implementations agree within model detail. Assumption ladder from there to
the Liu conditions: full red spectrum -3.5 pts; adding Yb 3 nm -20 pts
(the single biggest loss); green + Yb with every spacer re-optimised but
keeping the n=2.3 cap and no-ITO advantages: 53.9 % (Tpar 0.95) / 49.3 %
(0.85). The Liu device additionally uses a low-index BPBPA cap (~1.75),
ITO 15 nm inside the cavity, and q = 0.97, landing at the 44.1/35.8 %
bounds of scripts/66. So "~60 % is physically reachable in a red, Yb-free,
high-index-capped TE-OLED" and "59.2 % is unsupported for THIS green
device" are both true - two independent codes agree.

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
