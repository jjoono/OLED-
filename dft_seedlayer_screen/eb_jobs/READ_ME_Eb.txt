PACKAGE C -- E_b at the same level and geometry as E_d. 26 single points.
=========================================================================

HOW TO RUN
  NEW folder, independent of everything else. Extract, RUN_ALL.bat.
  About 1-2 hours. Can run alongside package B if you want -- it is small, but
  they will share the machine, so B finishes later by roughly that much.

WHAT IT IS
      E_b = E(molecule) + E(Ag) - E(complex)

  The complex energies already exist: stage 1 relaxed Ag on all 25 candidates
  at PBE0-D3/def2-SVP. This is the two pieces they are measured against -- 25
  bare molecules and one Ag atom.

WHY IT IS WORTH AN HOUR
  The E_b table this project has been carrying came from the PBE/DIIS route the
  audit disqualified -- the route where one geometry reached SCF solutions
  0.51 eV apart depending on how it converged. So the two halves of the
  screening metric are currently on different, and one of them on discredited,
  footing.

  It also matters for what the barriers have been saying. HATCN's Ag sits
  chemisorbed at 2.17 A from a nitrile nitrogen, and its intramolecular hop is
  0.065 eV -- it moves easily BETWEEN nitrogens while being held tightly ON
  them. That is a binding-energy story, and E_d alone cannot tell it. Venables
  puts nucleation density under both terms, not one.

WHAT WAS DONE TO KEEP IT CLEAN
  - The molecule is taken at the geometry it has IN the complex. The substrate
    was frozen through the stage-1 relaxation, so that is the input structure
    unchanged and no relaxation energy of the molecule leaks into E_b. Checked
    against the relaxed outputs: 25 of 25 match to under 0.01 A.
  - Dispersion is kept on all three terms. It is a large part of the binding
    for the physisorbed candidates, and taking it off one term only would be
    worse than taking it off all.
  - Every molecule is run as a closed-shell singlet, Ag as a doublet.

  Not corrected for basis set superposition. On def2-SVP that overbinds by
  something like 0.1-0.3 eV, roughly in proportion to contact area, so it
  inflates the strong binders more than the weak ones. Treat the column as
  ordering rather than absolute, the same as E_d.

WHEN IT IS BACK
  Send the eb_jobs folder. Harvest reads it against the stage-1 outputs and
  writes runs/binding_energies_pbe0.json.
