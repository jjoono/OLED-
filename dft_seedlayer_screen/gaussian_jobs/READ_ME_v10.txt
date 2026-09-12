v10 -- every path is now a real hop between two equivalent sites
================================================================

WHAT TO RUN
  RUN_ALL.bat   (or)   python run_campaign.py
  Stop and restart freely; finished jobs are skipped.
  336 jobs in 18 folders.

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

1. Equivalence is now tested, not assumed.
   "The nearest other atom of the same element" is not an equivalent site. It
   paired the two ends of the Al4O6 and Cu4I4 paths across different chemical
   environments. Atoms are now labelled by iterating each one's label over its
   bonded neighbours (Morgan relabelling, 3 rounds), and two sites count as
   equivalent only when those labels match.

2. A molecule with a unique binding site gets a neighbour molecule.
   Cs2CO3, pyridine, Me3P=O, thiophene, benzene, PhCN, DMABN, Mo3O8 and Mo3O9
   have nowhere to hop to on one molecule, so the neighbour is built: the copy
   produced by turning the molecule 180 degrees about an axis along the surface
   normal, through the midpoint of the two sites, pushed in until every
   cross-molecule atom pair is at van der Waals contact.

   That rotation is not cosmetic. It maps the whole complex at t=0 onto the
   whole complex at t=1 -- verified to 1e-15 A on all nine -- so the two ends
   are the SAME STATE by construction, and any energy difference between them
   is purely numerical. The error bar stops being an estimate.

   PhCN, for instance, now runs
       N0 2.50 -> 2.74 -> 3.35 -> 4.18 -> 3.35 -> 2.74 -> 2.50 N13
   which is a hop with a desorption maximum in the middle, symmetric on its
   face.

3. Two atoms 2.7 A apart are one site, not two.
   Bphen's equivalent nitrogens are the two halves of one bidentate pocket and
   Ag bridges both at once; its v8 "hop" never left that pocket. An equivalent
   atom now counts as a second site only if Ag can sit on one without sitting
   on the other.

4. The adatom has to stay on the face it started on.
   "Furthest from every atom" on a small flat molecule is off the edge past a
   C-H, not above the ring -- benzene's path left the ring and crossed over a
   hydrogen. Directions more than 60 degrees off the molecular surface normal
   are no longer offered.

   Paths now come out symmetric on both ends, which they should be:
       HATCN   N0 2.50 -> ... -> 2.50 N29
       F4TCNQ  N0 2.50 -> ... -> 2.50 N4
       Al4O6   O4 2.40 -> ... -> 2.40 O6

ALSO CARRIED OVER FROM v9
  Three adatom heights (dz = -0.20, +0.20, +0.60 A) instead of two, so the
  reader can fit a parabola and quote the energy at the RELAXED height. Two
  heights cannot bracket a minimum, and on the v8 data that was worth up to
  0.456 eV -- more than any barrier in the table.

WHAT IS DEFERRED, AND WHY
  Bphen (85 atoms as a dimer), Ph3P=O (71), TPA (69) all need the neighbour
  molecule and are over the 50-atom budget. Raise it with GAUSS_DIMER_MAX if
  you want them, but expect each to cost more than the rest of the campaign.
  PhCz's input complex has Ag resting on a hydrogen, so it was never on a
  binding site; it needs re-optimising, not re-running.

STILL NOT FIXED
  The substrate is one or two molecules, not a film -- no extended packing, no
  bulk polarisation. The adatom still cannot move sideways off the path, so the
  barrier stays an upper bound. F4TCNQ is spin contaminated throughout
  (S**2 0.776-0.846). The output is a screening ORDER, not a table of barriers,
  and it is the correlation against measured closure thickness that will
  license it.
