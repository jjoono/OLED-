v8 -- sized for the 32-core / 768 GB workstation
================================================

WHAT TO RUN
  Extract this folder, then double-click            RUN_ALL.bat
  or, from the Miniforge Prompt inside the folder:  python run_campaign.py

  Options:  python run_campaign.py --workers 4 --g16 C:\G16W
  Stop it any time (Ctrl+C or close the window) and run it again -- finished
  jobs are skipped, so it continues where it left off.

  When every folder says "folder complete", zip the whole gaussian_jobs folder
  and send it back.


WHY THE CPU SAT AT 30%
  Nothing was wrong. %NProcShared=16 on a 32-core machine with SMT is 16 of 64
  logical processors = 25%, plus OS overhead. Gaussian was using nearly all of
  what it had been given; the other three quarters of the machine was idle.

  Do NOT read anything into "Job cpu time = Elapsed time" in the .out file.
  G16W reports cpu with one-second granularity -- every value in the file is a
  whole number -- so it cannot tell you how many threads ran.

  Raising %NProcShared would not have helped either. A 475-basis-function
  single point stops scaling well before 16 threads: the Fock build
  parallelises, the diagonalisation and DIIS bookkeeping do not. The way to use
  the machine is several jobs at once.


WHAT CHANGED

1. Eight threads per job, several folders at a time.
   Each candidate folder is already an independent orbital chain with its own
   .chk, so folders are the natural unit of parallelism and nothing has to be
   split. run_campaign.py runs (physical cores / 8) folders at once -- four on
   this machine -- each with its OWN Gaussian scratch directory. One shared
   scratch is what stops several G16W processes from coexisting; that was the
   real blocker, not the thread count.

   Expected: about 4x the throughput, with the CPU near 90% instead of 30%.

2. %Mem=48GB per job (4 x 48 = 192 GB of 768).

3. The level shift now applies only where it earns its cost.
   Bphen_t3_z1 needed 38 SCF cycles under VShift=300 while READING the previous
   point's orbitals -- the same point took 5 cycles in v5. A level shift damps
   the update by construction: that is what you want from a Harris guess and
   pure overhead once you are already sitting on the answer. So:
       first job of a folder   SCF=(XQC,MaxCycle=128,VShift=300)
       every chained job       SCF=(XQC,MaxCycle=128)
   This is an inference from one point, not a measurement of both halves --
   PBE0 and VShift changed together in v6. If the chained jobs still take tens
   of cycles, say so and the shift is not the cause.

4. The job order is written down, not inferred.
   Sorted alphabetically, "<name>_t0_z0_ref" lands immediately after
   "<name>_t0_z0", so the old batch loop would have run the reference re-run
   SECOND -- reading only its own orbitals, which is exactly the guess it
   exists to escape. Every folder now carries ORDER.txt and the runner follows
   it. This was a real bug in the v7 package.

5. A folder stops at its first failure.
   Continuing down a chain whose orbitals were never written means every later
   point reads a stale .chk, and the barrier ends up made of two different
   electronic states. Other folders carry on.


ABOUT IN-CORE INTEGRALS
  Two-electron integrals for 475 basis functions are roughly 51 GB, which does
  fit. Holding them in memory instead of recomputing them each cycle could help
  -- but PBE0's cost here is split between the exchange integrals and the XC
  quadrature, and only the first would benefit. It is not assumed. If you want
  to know, run one point both ways:

      python run_campaign.py --workers 1            (as shipped)
      then add SCF=(...,InCore) to one .gjf by hand and time the same point

  The per-job elapsed time is printed for every job, so the comparison is
  there to read. Job parallelism is the 4x; this would be a few percent on top
  and is not worth risking a memory failure across four concurrent jobs.


CONTENTS
  246 jobs in 22 folders. Route unchanged from v6/v7 apart from item 3:
  PBE1PBE/Def2SVP with GD3BJ dispersion, NoSymm, orbitals chained along each
  path. Geometry as v7 (per-element contact distances, per-site outward
  normals, 5-point paths).


WHAT TO REPORT BACK
  - the first few lines of run_campaign.py's output (how many workers it chose)
  - HATCN_t0_z0: "SCF Done" energy and cycle count
  - HATCN_t0_z1: cycle count -- if it is still ~38, item 3 was the wrong call
