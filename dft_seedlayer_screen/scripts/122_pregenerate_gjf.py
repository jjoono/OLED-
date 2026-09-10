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
from pathgeom import (CANDIDATES, NPATH_MAX, NPATH_MIN, SPACING, STRUCT, ZSCAN,
                      destination, geometry, place, read_xyz, sanity)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                   "gaussian_jobs")
ROUTE = ("#P PBEPBE/Def2SVP EmpiricalDispersion=GD3BJ SP "
         "SCF=(XQC,MaxCycle=128) NoSymm")


def gjf(syms, xyz, title, nproc, mem, mult=2):
    out = [f"%NProcShared={nproc}", f"%Mem={mem}GB", ROUTE, "", title, "",
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
        dest, cls = destination(sub_s, sub_x, ag, anchor, rule)
        span = float(np.linalg.norm(dest - ag))
        npath = int(np.clip(round(span / SPACING) + 1, NPATH_MIN, NPATH_MAX))

        safe = tag.replace("=", "").replace("-", "")
        d = os.path.join(OUT, safe)
        os.makedirs(d, exist_ok=True)
        n = 0
        for i, t in enumerate(np.linspace(0.0, 1.0, npath)):
            pos = place(sub_x, (1 - t) * ag + t * dest, nrm)
            for j, dz in enumerate(ZSCAN):
                name = f"{safe}_t{i}_z{j}"
                with open(os.path.join(d, name + ".gjf"), "w") as f:
                    f.write(gjf(sub_s + ["Ag"],
                                np.vstack([sub_x, pos + dz * nrm]),
                                f"{tag} t={t:.3f} dz={dz:+.2f} class={cls}",
                                nproc, mem, mult))
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
