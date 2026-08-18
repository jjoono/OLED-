# Can a rough Ag TOP surface be flattened at ~5 nm? Literature/patent survey (2026-08-18)

Question: is post-growth planarization of the top face of an ultrathin Ag
electrode essentially nonexistent? ANSWER: NO - three real routes exist,
two demonstrated at/near the target thickness. The naive assumption that
"etch-back preserves roughness" is WRONG at low ion energy: ballistic
smoothing beats it.

## Routes that exist

1. ION-BEAM THINNING-BACK (closest match to the question)
   Deposit Ag thick (10-15 nm, self-smoothed and closed), then thin back
   with low-energy ions to ~4.5 nm with an ULTRASMOOTH top surface.
   "Pushing the thinness limit of silver films ... via ion-beam
   thinning-back process" (PMC10933474). Companion result on Au: ~100 eV
   Ar+ ballistic smoothing + grain-boundary restructuring keeps 4 nm Au
   compact (PMC11935187). Equipment: ordinary ion mill run at <=200 eV.
   Synergy with our physics: A(d) is flat in thickness (scripts/78), so
   depositing thick costs nothing in absorption; thinning back recovers T.

2. GAS CLUSTER ION BEAM (GCIB) SMOOTHING - patented, demonstrated post-fab
   US6613240 "Method and apparatus for smoothing thin conductive films by
   gas cluster ion beam". Applied on 12 nm Ag plasmonic films POST-
   fabrication: grains widened ~3x, up to 4-fold optical/electrical
   improvement (Appl. Phys. A 117, 2014; s00339-014-8728-1). Per-atom
   energy ~eV scale despite 15-30 keV clusters -> damage depth ~1 nm; the
   Ag itself shields the organics underneath. Equipment rare (Exogenesis/
   ANAB class).

3. TEMPLATE-STRIP + TRANSFER (avoid growing the top face at all)
   Al2O3/Ag stripped from Si: exposed Ag RMS < 0.21 nm, Rs 7.4 Ohm/sq,
   T 93.9 % (Sci. Rep. 7, 44576). Transfer onto ORGANICS is established by
   the cold-welding patent family: US6294398 (stamp pressed onto
   unpatterned OLEDs), US6895667 (transfer of patterned metal by
   cold-welding), US8222072 / US8637345 (LOW-pressure cold welding).
   Our E_b screen doubles as the adhesion-layer ranking for the receiving
   side (HATCN 1.03 > TPBi 0.89 > phen 0.87 eV).

## Routes that do NOT work at 5 nm

- Laser smoothing: ns-laser melting of 2-9.5 nm Ag produces spontaneous
  bicontinuous DEWETTING (arXiv:1002.2133) - the laser is a dewetting
  tool at this thickness, not a smoother.
- Plain thermal annealing (uncapped): dewetting, established.
- CMP / electropolish: solution/slurry, incompatible and unscaled to 5 nm.
- Mechanical flat-stamp pressing of a continuous film: no credible
  demonstration found at this scale (cold welding transfers, it does not
  planarize a grown film).

## Compatibility on TOP of an OLED

Both ion routes are vacuum, room-temperature; 100 eV Ar+ range in Ag is
<1 nm and GCIB damage ~1 nm, so the metal shields the stack. Heat load and
charging manageable. Direct on-OLED demonstrations: cold-welding stamping
was done on finished OLEDs (US6294398); ion-thinning-back demonstrated on
flexible optoelectronic substrates - on-OLED validation would be novel.

## Practical ranking for this lab

1. Ion-beam thinning-back - needs only an ion mill at low energy; test on
   glass/HATCN/Ag(12 nm) -> thin to 5-6 nm; verify with Rs (rho = Rs*d
   should approach the Fuchs-Sondheimer floor as the surface smooths).
2. Capped mild anneal (80-120 C) - zero new equipment, small gain.
3. Template-strip transfer - highest ceiling (sub-0.3 nm RMS), highest
   integration risk; our seed chemistry transfers to the adhesion layer.
4. GCIB - only if a tool can be accessed (rare).

Coverage caveat: web-snippet search, not exhaustive Scopus/patent-family
due diligence; run 42_lit_harvest locally with these seeds before quoting
"first" claims in a paper.
