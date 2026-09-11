v7 -- PBE0 route (as v6) + four fixes the v5 results exposed
============================================================

DELETE the v5 and v6 folders. Every geometry in here is different from both,
so their outputs cannot be mixed with these.

FAULT 1  The path never left its starting atom.
  destination() put the adatom above the destination atom along the STARTING
  atom's outward normal. On anything three-dimensional that direction points
  back into the molecule, so place() pushed the point straight out again and it
  landed next to where it began. Al4O6 v5, nearest atom at every point:

      t0 O4:2.20   t1 O4:2.23   t2 O4:2.21   t3 O4:2.24
      (the destination oxygen never came closer than 3.00 A)

  There was no hop in that path at all -- the 0.230 -> 0.047 eV "downhill
  barrier" was the adatom sliding around one oxygen.
  FIX: each site gets its own outward direction, chosen as the one that puts
  the adatom furthest from every substrate atom at its own contact distance.
  The path now interpolates the SITE and the NORMAL separately instead of
  drawing a chord between two adatom positions, so the adatom rides over the
  surface. Al4O6 v7:

      t0 O4:2.20   t1 O4:2.23   t2 Al0:2.72   t3 O6:2.53   t4 O6:2.20

  A site2site path is refused if it still ends over the atom it started on.

  BONUS: the two ends now sit at the same height over equivalent atoms
  (HATCN 2.30 -> 2.30, F4TCNQ 2.30 -> 2.30, Cs2CO3 2.40 -> 2.40,
  Cu4I4 3.10 -> 3.10, LiF32 2.50 -> 2.50). By symmetry their energies must be
  equal, so end-minus-start is a free error bar on the barrier. The reader
  prints it next to E_d and flags any candidate whose ends disagree.

FAULT 2  The first job of each folder landed on the wrong SCF solution.
  It is the only job with no Guess=Read. Three of seven v5 folders collapsed:
      benzene_t0_z0  +3.41 eV  S**2 = 0.7567
      DMABN_t0_z0    +2.20 eV  S**2 = 0.7534
      Bphen_t0_z0    +0.15 eV  S**2 = 0.7551
  and that point is the barrier's zero.
  FIX: each folder ends with <name>_t0_z0_ref.gjf, the same geometry re-run
  from the orbitals the chain finished on. The reader also rejects any point
  with S**2 > 0.7525 -- on the v5 data that picks out exactly those three.

FAULT 3  Cu4I4's 1.5 eV barrier was Pauli repulsion.
  Ag sat 3.10 A from iodine at t=0 and 2.24 A at t=2, because one 2.2 A contact
  distance served every element pair. Ag-I needs 2.84 A.
  FIX: contact distance per element pair from covalent radii. Cu4I4 now runs
  at 2.87-3.10 A across the whole path.

FAULT 4  Three path points cannot show a barrier.
  No interior point means the maximum is whichever endpoint is higher. Cs2CO3
  (0.023, 0.000, 0.135 eV) and Bphen (0.003, 0.093, 0.103 eV) were both still
  rising at the last point.
  FIX: 5 points minimum, sized by how far the site moves.

COST
  246 jobs in 22 folders, up from 168. The v5 set of 54 took about 3 hours, so
  budget roughly 12-15 hours plus the hard ones.

STILL GOOD NEWS FROM v5
  All 54 jobs terminated normally. The 90-degree saddle that hung HATCN did not
  appear anywhere else.

WHAT TO REPORT BACK
  - HATCN_t0_z0: "SCF Done" energy and iteration count
  - HATCN_t0_z1: iteration count (tells us Guess=Read is working)
