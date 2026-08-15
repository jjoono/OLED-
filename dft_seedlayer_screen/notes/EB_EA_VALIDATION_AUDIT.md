# E_b / EA literature-validation audit (2026-08-15)

Question: are all seed-candidate E_b and EA values literature-validated?
Answer: NO - and half of them cannot be, in principle. Three tiers.

## A. Experimentally anchored (method-level validation)
- Ag2 dimer BE: 1.86 eV calc (PBE-D3/def2-SVP+CP) vs 1.65 eV exp (+13%).
  The only experimental anchor for the E_b machinery. REPORT.md sec.1.
- Ag bulk cohesion 2.95 eV (exp) as upper reference.
- EA(F4TCNQ)-EA(TCNQ) = 0.48 eV calc vs ~0.55 eV lit (scripts/45).
- F4TCNQ solid EA 5.24 eV (UPS/IPES) referenced with the gas->solid
  polarisation shift carried as an explicit separate term.
- xTB deliberately invalidated for energies (Ag2 5.04 vs 1.65) and
  restricted to geometry generation.

## B. Qualitative literature consistency only (not numerical)
- phen chelate 0.87 eV <-> BCP:Ag coordination (Org.Electron.2014),
  p-bPPhenB used as wetting inducer in sandwich electrodes.
- ZnS 1.60 eV <-> Kim AFM 2015 ZnS/Ag priors.
- MoO3 weak binding <-> MoO3/Ag wetting reports.

## C. Not validated - mostly impossible to validate
1. ALL per-molecule absolute E_b (HATCN 1.03, TPBi 0.89, F4TCNQ 0.97...):
   single-adatom adsorption energies are not experimentally measurable
   quantities; no literature counterpart exists. Paper must quote
   ordering + -0.2..0.3 eV model error, not absolutes.
2. Cluster vs slab gap (HATCN 1.03 vs 1.346): documented, but no external
   referee; slab side has internal lattice cross-check only.
3. New-candidate EAs (TCPM 0.87, SBF2CN 0.54 @SVP): novel molecules, no
   literature values exist. Internal ranking only (correct usage).
   ** TCPSi EA is MISSING from runs/candidate_EA.json - incomplete. **
4. HATCN EA 3.38 (gas calc): no gas-phase experiment exists; solid IPES
   4.4-5.0 not directly comparable (polarisation).
5. TCPM slab E_b: number itself unconfirmed (SCF state-hopping; local OT
   diagnostic queued).
6. F6TCNNQ E_b: not computable at PBE (SCF oscillation); left blank
   honestly (scripts/52 header).
7. EA->E_b proxy (r=0.962, slope 0.055): internal 4-fragment correlation,
   not literature - calibration tool only.

## Action items
- FREE WIN: add TCNQ absolute gas EA benchmark row - modern photoelectron
  experiment gives 3.383 eV vs our 3.54 (wB97X/TZVP, +0.15 eV / 4%).
  Calculation already cached (runs/ea_cache.json); just record the
  comparison in the benchmark table.
- Complete TCPSi EA (local queue, scripts/54).
- Resolve TCPM slab state consistency (local queue, LOCAL_HANDOFF sec.5).
- 42_lit_harvest.py (local) may surface prior Ag-on-organic DFT numbers
  (e.g. Ag/PTCDA-class) for additional cross-checks.
