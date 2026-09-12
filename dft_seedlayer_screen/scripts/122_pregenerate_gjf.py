"""Write every Gaussian input the barrier campaign needs, up front.

The workstation has Gaussian and nothing else. Installing Python there to drive
it is a dependency, a PATH problem and a support burden, and none of it is
necessary: the inputs do not depend on any result, so they can all be written
here and shipped as files. The workstation then runs one batch file that loops
g16 over them, which needs nothing that is not already installed.

One folder per candidate, one .gjf per point on its path, plus a RUN_ALL.bat
that walks them in order and skips anything already finished, so the job
survives being interrupted. Two details that matter on Windows: G16W writes
.out rather than .log, so the batch names the output explicitly and checks for
.out; and "already finished" means the output ends in Normal termination, so a
crashed job is retried rather than skipped forever. Reading the answers back is
a separate script that runs here.

Geometry comes from pathgeom, unchanged, so these inputs describe exactly the
same paths the psi4 driver would take. Candidates whose stored complex is not a
bonding geometry are refused here rather than silently producing a confident
zero.
"""
import os, sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathgeom import (CANDIDATES, NPATH_MAX, NPATH_MIN, SPACING, STRUCT,
                      ZSCAN, contact, destination, dimer_frames, equivalents,
                      frames, geometry, read_xyz, sanity)

# A dimer doubles the atoms, and the cost of a hybrid single point runs far
# faster than that. Above this size the neighbour-molecule hop is deferred
# rather than generated, so one candidate cannot swallow the campaign.
DIMER_MAX_ATOMS = int(os.environ.get("GAUSS_DIMER_MAX", "90"))

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                   "gaussian_jobs")
# The functional is a lever on convergence, not only on accuracy. Under PBE the
# Ag/HATCN complex has an alpha HOMO-LUMO gap of 0.22 eV, because GGA
# delocalisation error smears the Ag 5s electron across the acceptor and leaves
# the frontier space nearly degenerate. Quadratic convergence then stalls: the
# Newton step comes out 90.0 degrees from the gradient on every iteration (33 of
# them in the v5 HATCN run), so no step length lowers the energy and the
# gradient sits at 5.96e-3 forever while the energy is flat to 1e-9 Hartree.
# That is a saddle in orbital-rotation space -- an unstable wavefunction -- and
# no convergence threshold fixes it. A hybrid opens the gap, which removes the
# degeneracy that creates the saddle, and removes the delocalisation error that
# put the two SCF solutions 0.08 eV apart in the first place.
FUNC = os.environ.get("GAUSS_FUNC", "PBE1PBE")
# VShift separates the near-degenerate orbitals during the DIIS phase; at the
# fixed point a level shift moves virtual eigenvalues only, so the converged
# density is a genuine SCF solution. It changes which solution is reached, not
# what that solution is. XQC keeps quadratic convergence as a fallback rather
# than the first resort -- a QC iteration here runs ~60 linear-equation
# micro-iterations, so it costs far more than the hybrid's exact exchange.
# The shift is worth its cost only on the job that has to *find* the solution.
# Bphen_t3_z1 took 38 cycles under VShift=300 reading the previous point's
# orbitals; the same point took 5 under v5. A level shift damps the update by
# construction, which is what you want from a Harris guess and pure overhead
# once you are already sitting on the answer. So the first job of each folder
# keeps the shift and every chained job drops it. Both are overridable.
SCF_FIRST = os.environ.get("GAUSS_SCF_FIRST", "XQC,MaxCycle=128,VShift=300")
SCF_CHAIN = os.environ.get("GAUSS_SCF_CHAIN", "XQC,MaxCycle=128")


def route(read):
    scf = SCF_CHAIN if read else SCF_FIRST
    return (f"#P {FUNC}/Def2SVP EmpiricalDispersion=GD3BJ SP "
            f"SCF=({scf}) NoSymm" + (" Guess=Read" if read else ""))


RUNNER = r'''r"""Run the campaign: several folders at a time, each folder in path order.

    python run_campaign.py [--workers N] [--g16 C:\G16W]

Each candidate folder is one orbital chain with its own .chk, so folders are
the unit of parallelism and no chain has to be split. One worker owns a folder
from start to finish and runs its jobs in the order ORDER.txt gives, which is
not the alphabetical order: the reference re-run <name>_t0_z0_ref must go last,
after the chain it is meant to read.

Each worker gets its own Gaussian scratch directory. Sharing one is what stops
several G16W processes from running at the same time.

A job whose .out already ends in "Normal termination" is skipped, so the script
can be stopped and restarted and will pick up where it left off.
"""
import argparse, os, shutil, subprocess, sys, threading, time

HERE = os.path.dirname(os.path.abspath(__file__))
LOCK = threading.Lock()
STATE = {"done": 0, "skip": 0, "bad": 0}


def finished(out):
    """True if this job already ran to completion."""
    try:
        with open(out, "rb") as f:
            f.seek(max(0, os.path.getsize(out) - 4096))
            return b"Normal termination" in f.read()
    except OSError:
        return False


def say(msg):
    with LOCK:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def jobs_of(folder):
    """Job stems in the order they must run."""
    order = os.path.join(folder, "ORDER.txt")
    if os.path.exists(order):
        names = [l.strip() for l in open(order) if l.strip()]
        return [n for n in names
                if os.path.exists(os.path.join(folder, n + ".gjf"))]
    # No ORDER.txt: fall back to alphabetical but push the reference re-runs to
    # the end, which is the one thing alphabetical order gets wrong.
    stems = sorted(f[:-4] for f in os.listdir(folder) if f.endswith(".gjf"))
    return ([s for s in stems if not s.endswith("_ref")]
            + [s for s in stems if s.endswith("_ref")])


def _threads_of(folder):
    """%NProcShared this folder's jobs ask for."""
    for fn in sorted(os.listdir(folder)):
        if fn.endswith(".gjf"):
            for line in open(os.path.join(folder, fn)):
                if line.lower().startswith("%nprocshared="):
                    try:
                        return int(line.split("=")[1])
                    except ValueError:
                        break
            break
    return 8


def _cost(folder):
    """Rough work in a folder: atoms cubed times the number of jobs."""
    gjfs = [f for f in os.listdir(folder) if f.endswith(".gjf")]
    if not gjfs:
        return 0.0
    atoms = sum(1 for line in open(os.path.join(folder, gjfs[0]))
                if len(line.split()) == 4 and line.split()[0].isalpha())
    return (atoms ** 3) * len(gjfs)


def run_folder(folder, g16, scr):
    os.makedirs(scr, exist_ok=True)
    env = dict(os.environ)
    env["GAUSS_EXEDIR"] = g16
    env["GAUSS_SCRDIR"] = scr
    env["PATH"] = g16 + os.pathsep + env.get("PATH", "")
    name = os.path.basename(folder)
    exe = os.path.join(g16, "g16.exe")

    for stem in jobs_of(folder):
        out = os.path.join(folder, stem + ".out")
        if finished(out):
            with LOCK:
                STATE["skip"] += 1
            continue
        say(f"{name}: {stem}")
        t0 = time.time()
        try:
            subprocess.run([exe, stem + ".gjf", stem + ".out"],
                           cwd=folder, env=env, check=False)
        except OSError as e:
            say(f"{name}: cannot start Gaussian -- {e}")
            with LOCK:
                STATE["bad"] += 1
            return
        dt = time.time() - t0
        ok = finished(out)
        with LOCK:
            STATE["done" if ok else "bad"] += 1
        if not ok:
            # Do not go on down a chain whose orbitals were never written: every
            # later point would read a stale .chk and the barrier would be made
            # of two different electronic states.
            say(f"{name}: {stem} did NOT finish normally after {dt/60:.1f} min "
                f"-- stopping this folder")
            return
        say(f"{name}: {stem} done in {dt/60:.1f} min")
        for junk in os.listdir(folder):
            if junk.startswith("Gau-") or junk == "fort.7":
                try:
                    os.remove(os.path.join(folder, junk))
                except OSError:
                    pass
    say(f"{name}: folder complete")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=0,
                    help="folders to run at once (default: cores // 8)")
    ap.add_argument("--g16", default=r"C:\G16W")
    ap.add_argument("--threads", type=int, default=0,
                    help="physical cores to keep busy (default: cores // 2)")
    ap.add_argument("--scratch", default="")
    a = ap.parse_args()

    if not os.path.exists(os.path.join(a.g16, "g16.exe")):
        print(f"Gaussian not found at {os.path.join(a.g16, 'g16.exe')}")
        print("Pass the right one:  python run_campaign.py --g16 C:\\path\\to\\G16W")
        return 1

    cores = os.cpu_count() or 8
    # os.cpu_count() reports logical processors, and SMT threads do not help a
    # Fock build -- they contend for the same floating-point units. Budget on
    # physical cores, assuming SMT is on, and give each job the 8 threads its
    # %NProcShared asks for.
    workers = a.workers or max(1, (cores // 2) // 8)
    budget = a.threads or max(8, cores // 2)
    scratch = a.scratch or os.path.join(
        os.environ.get("TEMP", HERE), "gauscr")

    folders = [os.path.join(HERE, d) for d in os.listdir(HERE)
               if os.path.isdir(os.path.join(HERE, d))
               and any(f.endswith(".gjf") for f in os.listdir(os.path.join(HERE, d)))]
    # Longest first. Folders taken alphabetically leave the 85-atom Bphen dimer
    # grinding on one worker for a day after everything else has finished;
    # starting the long poles immediately overlaps them with the short ones,
    # which is the whole of the difference in when the campaign ends.
    folders.sort(key=lambda f: -_cost(f))
    if not folders:
        print("no candidate folders next to this script")
        return 1

    print(f"{len(folders)} folders, {budget} physical cores to fill, "
          f"{cores} logical processors visible")
    print("  order: " + ", ".join(os.path.basename(f) for f in folders[:5])
          + ", ...")
    print(f"scratch under {scratch}\n")

    # Folders are claimed against a thread budget, not a worker count, because
    # they no longer all ask for the same number of threads: an 85-atom dimer
    # gets 16 where a 20-atom complex gets 8. Counting folders would either
    # oversubscribe the machine or leave half of it idle.
    queue = list(folders)
    qlock = threading.Condition()
    free = [budget]

    def worker(k):
        scr = os.path.join(scratch, f"w{k}")
        while True:
            with qlock:
                while True:
                    if not queue:
                        return
                    i = next((j for j, f in enumerate(queue)
                              if _threads_of(f) <= free[0]), None)
                    if i is not None:
                        break
                    qlock.wait(30)
                folder = queue.pop(i)
                want = _threads_of(folder)
                free[0] -= want
            try:
                run_folder(folder, a.g16, scr)
            finally:
                with qlock:
                    free[0] += want
                    qlock.notify_all()

    t0 = time.time()
    nthread = max(workers, budget // 8)
    threads = [threading.Thread(target=worker, args=(k,), daemon=True)
               for k in range(nthread)]
    for t in threads:
        t.start()
    try:
        for t in threads:
            t.join()
    except KeyboardInterrupt:
        print("\nstopped -- re-run this script to continue where it left off")
        return 1

    print(f"\n{STATE['done']} ran, {STATE['skip']} already done, "
          f"{STATE['bad']} failed, {(time.time()-t0)/3600:.1f} h wall")
    if STATE["bad"]:
        print("Re-run this script to retry the folders that stopped.")
    shutil.rmtree(scratch, ignore_errors=True)
    print("When everything is done, send the whole gaussian_jobs folder back.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
'''

LAUNCHER = """@echo off
rem  Start the campaign. Everything real happens in run_campaign.py -- this file
rem  only finds a Python and hands over, because nested batch loops cannot run
rem  several Gaussian jobs at once or keep a job order that is not alphabetical.
cd /d "%~dp0"

rem  A real interpreter is looked for BEFORE whatever "python" is on PATH.
rem  Windows ships a stub at
rem  %LOCALAPPDATA%\\Microsoft\\WindowsApps\\python.exe that exists, is found by
rem  "where python", prints the word Python and does nothing. Testing that the
rem  candidate can actually run code is the only check that tells the two apart.
setlocal enabledelayedexpansion
set "PY="
for %%P in (
  "%USERPROFILE%\\miniforge3\\python.exe"
  "%LOCALAPPDATA%\\miniforge3\\python.exe"
  "%ProgramData%\\miniforge3\\python.exe"
  "C:\\Miniforge3\\python.exe"
  "%USERPROFILE%\\anaconda3\\python.exe"
  "%USERPROFILE%\\miniconda3\\python.exe"
) do if not defined PY if exist %%P set "PY=%%~P"

if not defined PY (
  for /f "delims=" %%P in ('where python 2^>nul') do (
    if not defined PY (
      echo %%P | find /i "WindowsApps" >nul || set "PY=%%P"
    )
  )
)

if defined PY (
  "!PY!" -c "import sys" >nul 2>&1
  if errorlevel 1 (
    echo   Found "!PY!" but it cannot run Python code.
    set "PY="
  )
)

if defined PY (
  echo   using !PY!
  echo.
  "!PY!" run_campaign.py %*
  goto :done
)
echo   No working Python found.
echo   The "python" on PATH is the Microsoft Store stub, which only prints the
echo   word Python. Open the Miniforge Prompt and run there:
echo       cd /d "%~dp0"
echo       python run_campaign.py

:done
echo.
pause
"""


def gjf(syms, xyz, title, nproc, mem, mult=2, chk=None, read=False):
    """One Gaussian input. `chk` names a checkpoint shared across the folder and
    `read` starts from the orbitals already in it.

    Chaining the orbitals along the path is not an optimisation, it is what
    makes the barrier mean anything. These Ag complexes have more than one SCF
    solution -- the same Bphen geometry converged to two states 0.51 eV apart
    under two different convergence routes -- and the barrier being measured is
    0.29 eV. Points that land on different solutions produce a number made of
    the gap between states rather than of the path. Reading the previous point's
    orbitals keeps every point on the same electronic state, and as a side
    effect converges in a few cycles instead of a hundred.
    """
    out = []
    if chk:
        out.append(f"%Chk={chk}")
    out += [f"%NProcShared={nproc}", f"%Mem={mem}GB",
            route(read), "", title, "",
            f"0 {mult}"]
    for a, c in zip(syms, xyz):
        out.append(f" {a:<2s} {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
    return "\n".join(out) + "\n\n"


def main():
    # Eight threads per job, not all of them. A 475-basis-function SP stops
    # scaling well before 16 threads -- the Fock build parallelises, the
    # diagonalisation and DIIS bookkeeping do not -- so the way to use a
    # 32-core machine is several jobs at once, not one wide one. Each folder
    # is an independent orbital chain with its own .chk, so folders are the
    # natural unit of parallelism and nothing has to be split.
    nproc = int(os.environ.get("GAUSS_NPROC", "8"))
    mem = int(os.environ.get("GAUSS_MEM_GB", "48"))
    # Eight threads is right for a 20-atom complex and wasteful for an 85-atom
    # one: a bigger Fock build has more work per thread, so it keeps more of
    # them busy. The three dimers that dominate the campaign get double.
    big_n = int(os.environ.get("GAUSS_BIG_ATOMS", "60"))
    os.makedirs(OUT, exist_ok=True)
    made, refused, index = 0, [], []

    for tag, fn, rule, mult in CANDIDATES:
        p = os.path.join(STRUCT, fn)
        if not os.path.exists(p):
            refused.append((tag, f"missing {fn}"))
            continue
        syms, xyz = read_xyz(p)
        why = sanity(syms, xyz, tag)
        if why:
            refused.append((tag, why))
            continue

        sub_s, sub_x, ag, anchor, nrm = geometry(syms, xyz)

        # A hop needs somewhere to hop to. If the binding site has no
        # equivalent partner on this molecule, the monomer cannot define a
        # barrier at all -- the v8 run silently measured the climb onto the ring
        # centre instead, which is the binding difference between two
        # inequivalent sites and is what put Bphen (0.101 -> 0.971 eV, uphill
        # all the way) at the top of a table where this project's own prediction
        # has it closing last. The hop that happens in a film is onto the
        # neighbouring molecule, so that is what gets built.
        # An equivalent atom is only a separate site if the adatom can sit on
        # one without sitting on the other. Bphen's two nitrogens are 2.7 A
        # apart and Ag bridges both at once -- that is one bidentate pocket, not
        # two sites, and the "hop" between them stays inside the same pocket.
        near_ok = []
        for j in equivalents(sub_s, sub_x, anchor):
            if np.linalg.norm(ag - sub_x[j]) > 1.5 * contact(sub_s[j]):
                near_ok.append(j)
        if near_ok:
            # The rule column in CANDIDATES was set by hand and it is not
            # allowed to veto a real hop: "face" on Bphen, benzene and PhCz sent
            # three molecules that do have equivalent sites down the ring-centre
            # path, which does not measure a barrier.
            near, dest, cls = destination(sub_s, sub_x, ag, anchor, "auto")
            # Size the path by how far the *site* moves, not by how far the
            # adatom's start and end points are apart. On a cage those differ by
            # arc: Al4O6 hops between oxygens 2.83 A apart, but the two adatom
            # positions sit 5.7 A from each other because each is lifted along
            # its own outward normal. Spacing by the chord would sample the hop
            # twice as finely as it needs and double the cost of the campaign.
            span = float(np.linalg.norm(sub_x[near] - sub_x[anchor]))
            npath = int(np.clip(round(span / SPACING) + 1,
                                NPATH_MIN, NPATH_MAX))
            pts, cls, _ = frames(sub_s, sub_x, ag, anchor, "auto", npath)

            # A hop that never leaves its starting atom is not a hop. On Al4O6
            # every v5 point stayed 2.20-2.24 A from the same oxygen while the
            # destination oxygen never came closer than 3.00 A, so the "barrier"
            # was only the adatom sliding downhill around one site.
            end = pts[-1][0]
            got = int(np.linalg.norm(sub_x - end, axis=1).argmin())
            if got == anchor:
                refused.append((tag, f"path ends over {sub_s[anchor]}{anchor} "
                                     f"again -- it never reaches "
                                     f"{sub_s[near]}{near}"))
                continue
        else:
            n_dimer = 2 * len(sub_s) + 1
            if n_dimer > DIMER_MAX_ATOMS:
                refused.append((tag, f"binding site is unique, so the hop needs "
                                     f"a neighbour molecule -- {n_dimer} atoms, "
                                     f"over the {DIMER_MAX_ATOMS}-atom budget"))
                continue
            _, _, _, span = dimer_frames(sub_s, sub_x, ag, anchor, nrm,
                                         NPATH_MIN)
            npath = int(np.clip(round(span / SPACING) + 1,
                                NPATH_MIN, NPATH_MAX))
            sub_s, sub_x, pts, span = dimer_frames(sub_s, sub_x, ag, anchor,
                                                   nrm, npath)
            cls = "dimer"

        safe = tag.replace("=", "").replace("-", "")
        d = os.path.join(OUT, safe)
        os.makedirs(d, exist_ok=True)
        n = 0
        first = None
        order = []
        np_ = nproc * 2 if len(sub_s) + 1 >= big_n else nproc
        mem_ = mem * 2 if len(sub_s) + 1 >= big_n else mem
        lim = np.array([contact(x) for x in sub_s])
        for i, (pos, nrm) in enumerate(pts):
            for j, dz in enumerate(ZSCAN):
                agx = pos + dz * nrm
                # The lowest height deliberately sits inside the contact
                # distance -- a point on the repulsive wall is what brackets the
                # minimum from below. Inside 85% of it the energy is no longer
                # parabolic and the fit would be pulled by a wall, so that
                # height is dropped for this point rather than fitted.
                if float((np.linalg.norm(sub_x - agx, axis=1) / lim).min()) < 0.85:
                    continue
                name = f"{safe}_t{i}_z{j}"
                body = (sub_s + ["Ag"], np.vstack([sub_x, agx]),
                        f"{tag} t={i/(len(pts)-1):.3f} dz={dz:+.2f} class={cls}")
                with open(os.path.join(d, name + ".gjf"), "w") as f:
                    f.write(gjf(*body, np_, mem_, mult,
                                chk=f"{safe}.chk", read=(n > 0)))
                order.append(name)
                if n == 0:
                    first = body
                n += 1

        # Re-run the path's first point last, from the orbitals the chain ended
        # on. Only the first job has no previous orbitals to read, so it is the
        # one point converged from a Harris guess -- and in the v5 run that guess
        # collapsed onto a far higher solution three times out of seven: benzene
        # by 3.41 eV, DMABN by 2.20 eV, Bphen by 0.15 eV, each flagged by a spin
        # contamination the rest of its folder did not have. That point is the
        # barrier's zero, so a bad one poisons every energy difference in the
        # folder. The re-run costs one cheap job and puts the reference on the
        # same electronic state as the path it is subtracted from.
        with open(os.path.join(d, f"{safe}_t0_z0_ref.gjf"), "w") as f:
            f.write(gjf(*first, np_, mem_, mult,
                        chk=f"{safe}.chk", read=True))
        order.append(f"{safe}_t0_z0_ref")
        n += 1
        # The runner must not infer the order from the filenames. Sorted
        # alphabetically "<name>_t0_z0_ref" lands immediately after
        # "<name>_t0_z0", which would run the reference re-run second -- reading
        # only its own orbitals, which is exactly the guess it exists to escape.
        # The order is written down instead.
        with open(os.path.join(d, "ORDER.txt"), "w") as f:
            f.write("\n".join(order) + "\n")
        index.append((tag, safe, cls, len(sub_s) + 1, npath, len(ZSCAN), n))
        made += n
        print(f"  {tag:<12} {cls:<10} {len(sub_s)+1:>3} atoms  "
              f"{npath} x {len(ZSCAN)} = {n:>2} jobs")

    with open(os.path.join(OUT, "INDEX.csv"), "w") as f:
        f.write("candidate,folder,path_class,atoms,path_points,heights,jobs\n")
        for row in index:
            f.write(",".join(str(x) for x in row) + "\n")

    with open(os.path.join(OUT, "run_campaign.py"), "w") as f:
        f.write(RUNNER)
    with open(os.path.join(OUT, "RUN_ALL.bat"), "w", newline="\r\n") as f:
        f.write(LAUNCHER)

    print(f"\n{made} jobs in {len(index)} folders -> {os.path.relpath(OUT)}")
    if refused:
        print("\nrefused:")
        for tag, why in refused:
            print(f"  {tag:<12} {why}")
        print("  (these need their complex re-optimised before a barrier means anything)")


if __name__ == "__main__":
    main()
