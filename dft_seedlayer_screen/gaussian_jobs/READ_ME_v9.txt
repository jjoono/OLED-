v9 -- make the barrier bigger than its own error bar
====================================================

WHAT TO RUN
  RUN_ALL.bat   (or)   python run_campaign.py
  Same as v8. Stop and restart freely; finished jobs are skipped.

WHY A v9
  The v8 numbers are not usable as barriers, and the reason is not SCF, basis
  or functional. It is that the adatom's height was never relaxed.

  Two heights cannot bracket a minimum. Whichever of the two is lower is an
  endpoint of the scan, so nothing was relaxed -- the height was guessed
  between two arbitrary values. Measured on the v8 output, that guess was worth

      DMABN t3   0.456 eV        Me3PO t4   0.424 eV
      F4TCNQ     0.350 eV        Al4O6 t3   0.176 eV

  which is more than any barrier in the v8 table. Worse, it varies along the
  path (DMABN: 0.008 eV at t0, 0.456 eV at t3), so it distorts the shape of the
  curve rather than shifting it, and which height wins flips mid-path:

      BTD      z0  z0  z1  z1  z1
      Cs2CO3   z0  z0  z0  z1  z0

  Part of every v8 "barrier" is just which of two arbitrary heights won.

WHAT CHANGED

1. Three heights instead of two: dz = -0.20, +0.20, +0.60 A from contact.
   The contact height won 26 of 42 v8 points, so the minimum sits near it on
   both sides and the scan has to reach below it too. A point slightly inside
   the repulsive wall is what lets the reader fit a parabola and report the
   energy at the relaxed height instead of the lower of two. A height that
   would put Ag inside 85% of its contact distance is dropped rather than
   fitted, so the fit is never pulled by a wall.

   Cost: 16 jobs per candidate instead of 11. 358 jobs, about 1.45x of v8.

2. The reader now quotes each barrier with its own error bar.
   A hop between equivalent sites is symmetric about its midpoint, so E(t) and
   E(1-t) must agree. Half the largest disagreement is an error bar measured on
   exactly the same footing as the barrier. A candidate whose barrier is not at
   least 3x that is marked "barrier below its own error" -- on the v8 data that
   is six of eight.

3. It also reports the barrier relative to a common reference molecule.
   Screening needs the ordering, not the absolute value, and the absolute value
   is what this protocol gets least right: a rigid substrate, a small basis and
   a single molecule instead of a film all push the same way on every
   candidate. Differences against a reference keep what screening uses.

WHAT IS STILL NOT FIXED (and cannot be, at this cost)
  - the substrate is one molecule, not a film: no neighbouring molecules, no
    bulk polarisation. A slab would change these numbers, possibly by 2x.
  - the adatom still cannot move sideways off the constructed path, so the
    barrier remains an upper bound on the true minimum-energy-path value.
  - F4TCNQ is spin contaminated at every point (S**2 0.776-0.846), a real
    open-shell Ag-to-acceptor charge transfer that a single-reference doublet
    describes poorly. It is the top-ranked candidate, so this matters.

  These are why the output is a screening ORDER, not a table of barriers.
