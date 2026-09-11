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
                      ZSCAN, destination, frames, geometry, read_xyz, sanity)

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


RUNNER = r'''"""Run the campaign: several folders at a time, each folder in path order.

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
    scratch = a.scratch or os.path.join(
        os.environ.get("TEMP", HERE), "gauscr")

    folders = sorted(
        os.path.join(HERE, d) for d in os.listdir(HERE)
        if os.path.isdir(os.path.join(HERE, d))
        and any(f.endswith(".gjf") for f in os.listdir(os.path.join(HERE, d))))
    if not folders:
        print("no candidate folders next to this script")
        return 1

    print(f"{len(folders)} folders, {workers} at a time, "
          f"{cores} logical processors visible")
    print(f"scratch under {scratch}\n")

    queue = list(folders)
    qlock = threading.Lock()

    def worker(k):
        scr = os.path.join(scratch, f"w{k}")
        while True:
            with qlock:
                if not queue:
                    return
                folder = queue.pop(0)
            run_folder(folder, a.g16, scr)

    t0 = time.time()
    threads = [threading.Thread(target=worker, args=(k,), daemon=True)
               for k in range(workers)]
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

where python >nul 2>&1
if not errorlevel 1 (
  python run_campaign.py %*
  goto :done
)
if exist "%USERPROFILE%\\miniforge3\\python.exe" (
  "%USERPROFILE%\\miniforge3\\python.exe" run_campaign.py %*
  goto :done
)
echo   No Python found on PATH.
echo   Open the Miniforge Prompt and run:  python run_campaign.py

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
        near, dest, cls = destination(sub_s, sub_x, ag, anchor, rule)
        # Size the path by how far the *site* moves, not by how far the adatom's
        # start and end points are apart. On a cage those two differ by metres
        # of arc: Al4O6 hops between oxygens 2.83 A apart, but the two adatom
        # positions sit 5.7 A from each other because each is lifted along its
        # own outward normal. Spacing the path by the chord would sample the hop
        # twice as finely as it needs and double the cost of the campaign.
        heavy = np.array([x for sym, x in zip(sub_s, sub_x) if sym != "H"])
        end_site = sub_x[near] if near is not None else heavy.mean(axis=0)
        span = float(np.linalg.norm(end_site - sub_x[anchor]))
        npath = int(np.clip(round(span / SPACING) + 1, NPATH_MIN, NPATH_MAX))
        pts, cls, _ = frames(sub_s, sub_x, ag, anchor, rule, npath)

        # A hop that never leaves its starting atom is not a hop. On Al4O6
        # every v5 point stayed 2.20-2.24 A from the same oxygen while the
        # supposed destination oxygen never came closer than 3.00 A, so the
        # "barrier" was only the adatom sliding downhill around one site. Test
        # the thing that actually failed -- which atom the adatom ends up over.
        # Only site2site paths have a destination atom; a toface path is meant
        # to end over the ring centre, a short move that this test would
        # misread as a stalled hop.
        if cls == "site2site" and near is not None:
            end = pts[-1][0]
            got = int(np.linalg.norm(sub_x - end, axis=1).argmin())
            if got == anchor:
                refused.append((tag, f"path ends over {sub_s[anchor]}{anchor} "
                                     f"again -- it never reaches "
                                     f"{sub_s[near]}{near}"))
                continue

        safe = tag.replace("=", "").replace("-", "")
        d = os.path.join(OUT, safe)
        os.makedirs(d, exist_ok=True)
        n = 0
        first = None
        order = []
        for i, (pos, nrm) in enumerate(pts):
            for j, dz in enumerate(ZSCAN):
                name = f"{safe}_t{i}_z{j}"
                body = (sub_s + ["Ag"], np.vstack([sub_x, pos + dz * nrm]),
                        f"{tag} t={i/(len(pts)-1):.3f} dz={dz:+.2f} class={cls}")
                with open(os.path.join(d, name + ".gjf"), "w") as f:
                    f.write(gjf(*body, nproc, mem, mult,
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
            f.write(gjf(*first, nproc, mem, mult,
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
