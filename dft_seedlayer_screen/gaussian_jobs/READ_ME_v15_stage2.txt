v15 / stage 2 -- paths rebuilt from the relaxed sites. 454 jobs, 25 folders.
============================================================================

HOW TO RUN
  FRESH folder. Every geometry changed, so nothing from v10-v14 can be reused.
  Extract, RUN_ALL.bat. Expect roughly a day: Bphen (85 atoms as a dimer) and
  p-bPPhenB (73) set the floor, and the runner starts them first.

WHAT STAGE 1 FOUND
  Every input complex was relaxed Ag-only at PBE0-D3/def2-SVP, the level the
  barriers are computed at. How far Ag moved from where the input file had it:

      physisorbed         Mo3O9 0.02  Al4O6 0.10  Cu4I4 0.18  pyridine 0.20
                          Mo3O8 0.29  LiF32 0.33  benzene 0.51
      N / O donors        Bphen 0.82  F4TCNQ 1.37  BTD 1.37  B3PyMPM 1.48
                          HATCN 1.59  DMABN 1.59  p-bPPhenB 1.58  PhCN 1.87
                          pDCNB 1.91  oDCNB 1.92  Ph3PO 2.10  Cs2CO3 2.21
                          thiophene 2.36  Me3PO 2.41  Liq 3.08

  The physisorbed sites were right to within a tenth of an angstrom and their
  earlier barriers stand. The nitrogen and oxygen donors -- the molecules this
  project is about -- were off by more than a bond length, and every path built
  on them started in the wrong place and ended at the symmetry image of the
  wrong place. On HATCN that cost 0.165 eV at the site and inverted which point
  on the path is the well. HATCN reproduces the independent v14 probe to six
  decimals (-1487.228289 Ha); so does Cu4I4.

WHAT IS DIFFERENT NOW
  - The path's t=0 IS the relaxed geometry. Checked: at dz=+0.2 every endpoint
    sits exactly 0.2 A above its relaxed Ag-X distance (HATCN N 2.37, F4TCNQ N
    2.29, BTD N 2.47, Cu4I4 I 3.48). The one exception is Cs2CO3, 0.2 A high,
    because the neighbouring molecule's caesium is a real steric constraint in
    the dimer that the monomer relaxation did not see; the dz=-0.2 height
    covers it.
  - The far end is still the exact symmetry image of the near end, so the two
    ends stay one state by construction.
  - Re-centring can no longer press mid-path points into the gap: the straight
    line between the two relaxed endpoints is now a floor.
  - Four candidates rejoin: p-bPPhenB, B3PyMPM, Liq, PhCz. Their relaxations
    found real sites (N 2.35, N 2.66, O 2.42, ring face 3.12 A). TPA stays in
    too: its relaxed Ag sits at the periphery of a phenyl, nearest a hydrogen
    at 3.07 A, which after a relaxation is a result rather than a placement
    error. 25 candidates, up from 21.

WHAT THIS RUN SETTLES
  This is the last geometry fix available at this level of model. After it,
  the numbers are what PBE0-D3/def2-SVP says for an adatom on one or two
  molecules with the substrate frozen. What remains is the model itself -- no
  film, no lateral relaxation along the path (the barrier is an upper bound),
  F4TCNQ's spin contamination -- and none of that is fixed by another run.

WHEN IT IS BACK
  Send the whole folder. If any point still misses its bracket the reader will
  say so and one top-up will finish it, as before.
