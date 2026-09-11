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
SCF = os.environ.get("GAUSS_SCF", "XQC,MaxCycle=128,VShift=300")
ROUTE = (f"#P {FUNC}/Def2SVP EmpiricalDispersion=GD3BJ SP "
         f"SCF=({SCF}) NoSymm")


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
            ROUTE + (" Guess=Read" if read else ""), "", title, "",
            f"0 {mult}"]
    for a, c in zip(syms, xyz):
        out.append(f" {a:<2s} {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
    return "\n".join(out) + "\n\n"


def main():
    nproc = int(os.environ.get("GAUSS_NPROC", "16"))
    mem = int(os.environ.get("GAUSS_MEM_GB", "32"))
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
        for i, (pos, nrm) in enumerate(pts):
            for j, dz in enumerate(ZSCAN):
                name = f"{safe}_t{i}_z{j}"
                body = (sub_s + ["Ag"], np.vstack([sub_x, pos + dz * nrm]),
                        f"{tag} t={i/(len(pts)-1):.3f} dz={dz:+.2f} class={cls}")
                with open(os.path.join(d, name + ".gjf"), "w") as f:
                    f.write(gjf(*body, nproc, mem, mult,
                                chk=f"{safe}.chk", read=(n > 0)))
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
        n += 1
        index.append((tag, safe, cls, len(sub_s) + 1, npath, len(ZSCAN), n))
        made += n
        print(f"  {tag:<12} {cls:<10} {len(sub_s)+1:>3} atoms  "
              f"{npath} x {len(ZSCAN)} = {n:>2} jobs")

    with open(os.path.join(OUT, "INDEX.csv"), "w") as f:
        f.write("candidate,folder,path_class,atoms,path_points,heights,jobs\n")
        for row in index:
            f.write(",".join(str(x) for x in row) + "\n")

    bat = ["@echo off",
           "rem  Run every Gaussian job in this folder tree, in order.",
           "rem",
           "rem  The output file is named explicitly rather than left to Gaussian.",
           "rem  G16W on Windows writes .out, not .log, and a skip check written",
           "rem  against .log never fires -- every job would be re-run on restart.",
           "rem",
           "rem  A job counts as done only if its output ends in Normal termination,",
           "rem  so a crash is retried on the next pass instead of being skipped",
           "rem  forever.",
           "setlocal enabledelayedexpansion", "",
           'set "G16DIR=C:\\G16W"',
           'set "PATH=%SystemRoot%\\System32;%SystemRoot%;%G16DIR%;%PATH%"',
           'set "GAUSS_EXEDIR=%G16DIR%"',
           'set "GAUSS_SCRDIR=%TEMP%\\gauscr"',
           'if not exist "%GAUSS_SCRDIR%" mkdir "%GAUSS_SCRDIR%"', "",
           'if not exist "%G16DIR%\\g16.exe" (',
           '  echo   Gaussian not found at %G16DIR%\\g16.exe',
           '  echo   Edit G16DIR at the top of this file.',
           '  pause & exit /b 1', ')', "",
           'cd /d "%~dp0"',
           "set /a done=0", "set /a skip=0", "set /a bad=0", "",
           'for /r %%F in (*.gjf) do (',
           '  set "OUT=%%~dpnF.out"',
           '  set "SKIPME="',
           '  if exist "!OUT!" (',
           '    findstr /c:"Normal termination" "!OUT!" >nul 2>&1',
           '    if not errorlevel 1 set "SKIPME=1"',
           '  )',
           '  if defined SKIPME (',
           '    set /a skip+=1',
           '  ) else (',
           '    echo [!time!] %%~nF',
           '    pushd "%%~dpF"',
           '    "%G16DIR%\\g16.exe" "%%~nxF" "%%~nF.out"',
           '    popd',
           '    findstr /c:"Normal termination" "!OUT!" >nul 2>&1',
           '    if errorlevel 1 (',
           '      set /a bad+=1',
           '      echo        ... did NOT finish normally',
           '    ) else (',
           '      set /a done+=1',
           '    )',
           '  )',
           '  del /q "%%~dpF\\Gau-*.*" 2>nul',
           '  del /q "%%~dpF\\fort.7" 2>nul',
           ')', "",
           "echo.",
           "echo   ran !done! successfully, !skip! already done, !bad! failed",
           "echo.",
           "if !bad! gtr 0 echo   Re-run this file to retry the failed ones.",
           "echo   When all are done, send the whole gaussian_jobs folder back.",
           "pause"]
    with open(os.path.join(OUT, "RUN_ALL.bat"), "w", newline="\r\n") as f:
        f.write("\n".join(bat) + "\n")

    print(f"\n{made} jobs in {len(index)} folders -> {os.path.relpath(OUT)}")
    if refused:
        print("\nrefused:")
        for tag, why in refused:
            print(f"  {tag:<12} {why}")
        print("  (these need their complex re-optimised before a barrier means anything)")


if __name__ == "__main__":
    main()
