v6 -- PBE0 hybrid + level shift
===============================

WHAT CHANGED FROM v5
  route:  #P PBE1PBE/Def2SVP EmpiricalDispersion=GD3BJ SP
          SCF=(XQC,MaxCycle=128,VShift=300) NoSymm

  v5 used PBEPBE with SCF=(QC,MaxCycle=64,Conver=6).

WHY
  The v5 HATCN job did not fail on a threshold. It failed on the shape of
  the SCF surface. Gaussian printed, 33 times:

      Angle between quadratic step and gradient=  90.01 degrees.
      Incorrect curvature in search direction -- initial direction reversed.

  Exactly 90 degrees means the Newton step is orthogonal to the gradient:
  the residual gradient lies in a direction of negative curvature. That is a
  saddle point in orbital-rotation space -- an unstable wavefunction. No step
  length lowers the energy along it, so the gradient sat at 5.96e-3 from
  iteration 13 to 43 while the energy was flat to 1e-9 Hartree. Conver=6
  could not help: QC tests the gradient, and the gradient was stuck for a
  structural reason.

  The cause is the 0.22 eV alpha HOMO-LUMO gap, which comes from PBE's
  delocalisation error smearing the Ag 5s electron onto the acceptor. The
  same error is what put two SCF solutions 0.08 eV apart (-1487.257852 under
  DIIS vs -1487.254892 here) -- 28% of the 0.29 eV barrier being measured.

  PBE0 carries 25% exact exchange, which opens the gap and removes both the
  saddle and the two-solution ambiguity. VShift=300 (0.3 Hartree) separates
  the frontier orbitals during the DIIS phase; at the fixed point a level
  shift moves only virtual eigenvalues, so the converged density is a genuine
  SCF solution.

COST
  Exact exchange costs more per Fock build, but a QC iteration here ran ~60
  linear-equation micro-iterations. Trading 43 stalled QC iterations for
  ~20 DIIS iterations should be faster in wall time, not slower.

HOW TO RUN
  1. stop the current run
  2. delete the old gaussian_jobs folder
  3. extract this one, run RUN_ALL.bat

  The first job in each folder is the expensive one; the rest read the
  previous point's orbitals from the folder's .chk.

WHAT TO REPORT BACK
  - HATCN_t0_z0: the "SCF Done" energy and the iteration count
  - HATCN_t0_z1: the iteration count (tells us if Guess=Read is working)

  168 jobs, 22 folders. TPBi excluded.
  p-bPPhenB / B3PyMPM / Liq still refused -- their complexes need
  re-optimising before a barrier means anything.
