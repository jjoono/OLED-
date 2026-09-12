v10 -- every path is a real hop, and the three big dimers are in
================================================================

WHAT TO RUN
  RUN_ALL.bat   (or)   python run_campaign.py
  Stop and restart freely; finished jobs are skipped.
  402 jobs in 21 folders. Expect about a day -- see SCHEDULING below.

WHAT WAS WRONG IN v8
  A diffusion barrier is the climb between two EQUIVALENT binding sites. Ten of
  the twenty-two v8 paths were not that. They dragged Ag from its site onto the
  ring centre, and those profiles climb monotonically to the end:

      Bphen     0.101  0.189  0.471  0.971  0.790   (eV)
      pyridine  0.010  0.042  0.228  0.661  0.750

  That is the binding-energy difference between two INEQUIVALENT sites. It put
  Bphen at the top of the ranking -- the molecule this project's own model has
  closing last, at 28 nm.

FOUR FIXES

1. Equivalence is tested, not assumed.
   "The nearest other atom of the same element" is not an equivalent site; it
   paired the two ends of the Al4O6 and Cu4I4 paths across different chemical
   environments. Atoms are labelled by iterating each label over its bonded
   neighbours (Morgan relabelling, 3 rounds), and two sites count as equivalent
   only when those labels match. This also caught Bphen and benzene, which had
   equivalent sites all along and were being sent down the ring-centre path by
   a hand-set rule.

2. A molecule with a unique binding site gets a neighbour molecule.
   Eleven of them do: Cs2CO3, pyridine, Me3P=O, Ph3P=O, thiophene, benzene,
   PhCN, DMABN, TPA, Mo3O8, Mo3O9, plus Bphen by fix 3. The neighbour is the
   copy produced by turning the molecule 180 degrees about an axis along the
   surface normal through the midpoint of the two sites, pushed in until every
   cross-molecule atom pair is at van der Waals contact.

   That rotation maps the whole complex at t=0 onto the whole complex at t=1 --
   checked to 1e-15 A on all of them -- so the two ends are the SAME STATE by
   construction and any energy difference between them is purely numerical.
   The error bar stops being an estimate.

   PhCN now runs
       N0 2.50 -> 2.74 -> 3.35 -> 4.18 -> 3.35 -> 2.74 -> 2.50 N13
   a hop with a desorption maximum in the middle, symmetric on its face.

3. Two atoms 2.7 A apart are one site, not two.
   Bphen's equivalent nitrogens are the two halves of one bidentate pocket and
   Ag bridges both at once; its v8 "hop" never left that pocket. An equivalent
   atom counts as a second site only if Ag can sit on one without sitting on
   the other. Bphen therefore gets a dimer too.

4. The adatom stays on the face it started on.
   "Furthest from every atom" on a small flat molecule is off the edge past a
   C-H, not above the ring, and benzene's path left the ring and crossed over a
   hydrogen. Directions more than 60 degrees off the molecular surface normal
   are no longer offered. Paths now end where they should:
       HATCN   N0 2.50 -> ... -> 2.50 N29
       F4TCNQ  N0 2.50 -> ... -> 2.50 N4
       Al4O6   O4 2.40 -> ... -> 2.40 O6

CARRIED OVER FROM v9
  Three adatom heights (dz = -0.20, +0.20, +0.60 A) instead of two, so the
  reader fits a parabola and quotes the energy at the RELAXED height. Two
  heights cannot bracket a minimum, and on the v8 data that was worth up to
  0.456 eV -- more than any barrier in the table.

SCHEDULING -- this is why it is a day and not six hours
  Bphen's dimer is 85 atoms, Ph3P=O's 71, TPA's 69, against 9 to 43 for
  everything else. Two things keep that from doubling the wall time:

  - Those three get %NProcShared=16 and %Mem=96GB; the rest keep 8 and 48GB.
    A bigger Fock build has more work per thread, so it is the large jobs that
    can actually use them.
  - The runner now claims folders against a THREAD budget (32 physical cores)
    rather than counting folders, and starts the longest folder first. Taken
    alphabetically, Bphen would grind on one worker for most of a day after
    everything else had finished.

  Bphen sets the floor at roughly 18 hours; the rest fits underneath it.
  If you want it shorter, drop Bphen and run it separately later.

WHAT IS STILL OUT
  PhCz -- its input complex has Ag resting on a hydrogen at 2.51 A, so it was
  never on a binding site. That needs re-optimising, not re-running.
  p-bPPhenB, B3PyMPM, Liq -- same problem, worse (Ag inside a bond length).

STILL NOT FIXED
  The substrate is one or two molecules, not a film -- no extended packing, no
  bulk polarisation. The adatom cannot move sideways off the path, so the
  barrier stays an upper bound. F4TCNQ is spin contaminated throughout
  (S**2 0.776-0.846). The output is a screening ORDER, not a table of barriers,
  and it is the correlation against measured closure thickness that licenses it.
