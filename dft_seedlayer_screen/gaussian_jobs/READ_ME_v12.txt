v12 -- the two errors v11 isolated, both fixed. 402 jobs, 21 folders.
=====================================================================

WHAT TO RUN
  RUN_ALL.bat   (or)   python run_campaign.py
  This is a FULL run -- the geometries changed, so v10/v11 outputs cannot be
  mixed in. Use a fresh folder. About 14-18 hours (Bphen is back in).

WHAT v11 ESTABLISHED
  Two things, one good and one that pointed straight at the remaining error.

  The SCF is not the problem any more. A dimer path's two ends are the same
  geometry computed from opposite ends of a 22-job orbital chain, and they agree
  to at most 0.43 meV against barriers of 50-500 meV. That closes out everything
  v4 through v8 worried about.

  Relaxing the adatom height did NOT shrink the monomer asymmetry -- HATCN stayed
  at 0.179 eV, Cu4I4 went 0.176 -> 0.203. So the monomer error is not the height
  and not the SCF. By elimination it is the construction, and the dimers show
  exactly what the difference is: their ends are related by an exact rotation,
  the monomers' ends were built by two independent outward() calls.

FIX 1 -- the far end is now the near end's symmetry image
  Every monomer candidate turns out to have a real symmetry operation carrying
  its starting site onto its destination site:

      HATCN, F4TCNQ, BTD, pDCNB, oDCNB, LiF32   180 degrees
      Al4O6, Cu4I4                              120 degrees
      triazine                                  177.8 (its structure is slightly
                                                 distorted)

  The adatom's end position is that operation applied to its start position, and
  the far surface normal is the near one carried over by it. For the two-fold
  cases the whole path is equivariant, not just the endpoints: the point at 1-t
  is the image of the point at t. Checked on all 21 folders -- every path now
  starts and ends at an identical distance from an identical atom, to the last
  digit printed.

  Expect HATCN, Cu4I4, BTD and pDCNB, the four that failed the error test, to
  come back with dimer-like error bars.

FIX 2 -- the height scan is centred on the binding distance
  place() lifts until EVERY atom clears its contact distance, so one bulky atom
  set the height for the whole complex: on Cs2CO3 the caesium, whose contact
  distance is 3.89 A, held Ag well off the carbonate oxygen it is bound to. The
  scan window then sat above the binding distance by more than its own width,
  and no number of extra heights could reach it -- Cs2CO3, Mo3O8 and F4TCNQ had
  not one bracketed point out of nineteen even after the v11 top-up.

  Each path point is now placed at the same tightness -- distance over contact
  distance -- that the adatom has at its own binding site in the input complex.
  Along a path that number is now flat and symmetric where it should be:

      HATCN   1.136  1.128  1.088  1.128  1.136
      Cs2CO3  1.129  1.134  1.132  1.122  1.132  1.134  1.129

  And where the input already has Ag tighter than its contact distance (Mo3O8,
  0.936) the scan is lifted as a whole rather than having its lowest height
  dropped, so every point gets three heights instead of two.

BPHEN IS BACK IN
  As a dimer, 85 atoms, at 16 threads and 96 GB. It is the long pole at roughly
  18 hours; the runner starts it first so the rest fits underneath.

STILL NOT FIXED
  One or two molecules, not a film. No lateral relaxation, so the barrier stays
  an upper bound. F4TCNQ is spin contaminated throughout (S**2 ~ 0.81), and it
  is near the top of the ranking. PhCz, p-bPPhenB, B3PyMPM and Liq still have
  Ag parked somewhere that is not a binding site in their input structures and
  need re-optimising, not re-running.
