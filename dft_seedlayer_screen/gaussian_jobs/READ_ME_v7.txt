v7 -- PBE0 route (as v6) + three geometry/reference fixes
=========================================================

v6 fixed the SCF route. The 54 jobs that finished under v5 turned out to have
three further faults that had nothing to do with SCF settings. All three are
fixed here, so DELETE the v6 folder as well and use this one.

FAULT 1  The first job of each folder landed on the wrong SCF solution.
  It is the only job with no Guess=Read, so it converges from a Harris guess.
  Three of seven folders collapsed:
      benzene_t0_z0   +3.41 eV   S**2 = 0.7567
      DMABN_t0_z0     +2.20 eV   S**2 = 0.7534
      Bphen_t0_z0     +0.15 eV   S**2 = 0.7551
  and that point is the barrier's zero, so the whole folder was poisoned.
  FIX: each folder now ends with <name>_t0_z0_ref.gjf -- the same geometry
  re-run from the orbitals the chain finished on. One cheap extra job puts the
  reference on the same electronic state as the path.

FAULT 2  Cu4I4's "1.5 eV barrier" was Pauli repulsion, not diffusion.
  Ag sat 3.10 A from iodine at t=0 and 2.24 A at t=2, because the code used one
  2.2 A contact distance for every element pair. Ag-I needs 2.84 A.
  FIX: contact distance is now per element pair from covalent radii.
  Cu4I4 now runs at 2.85-3.10 A across the whole path instead of 2.24-3.10.
  Every folder's geometry changed, so v5/v6 outputs cannot be mixed with these.

FAULT 3  Three path points cannot show a barrier.
  With t = 0, 0.5, 1 there is no interior point, so the maximum is whichever
  endpoint is higher. Cs2CO3 (0.023, 0.000, 0.135 eV) and Bphen (0.003, 0.093,
  0.103 eV) were both still rising at the last point -- lower bounds, not
  barriers.
  FIX: 5 path points minimum (7 where the hop is long).

COST
  246 jobs in 22 folders, up from 168. The v5 set of 54 took about 3 hours,
  so budget roughly 12-15 hours, plus whatever the hard ones cost.

STILL GOOD NEWS FROM v5
  All 54 jobs terminated normally. The 90-degree saddle that hung HATCN did
  not appear anywhere else.

WHAT TO REPORT BACK
  - HATCN_t0_z0: "SCF Done" energy and iteration count
  - HATCN_t0_z1: iteration count (tells us Guess=Read is working)
