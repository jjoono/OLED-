"""Let Gaussian relax the adatom at a few path points, substrate frozen.

    python scripts/125_relax_probe.py <returned gaussian_jobs folder> HATCN:0 F4TCNQ:0 Cu4I4:2

The height scan approximates a relaxation with three points and a parabola.
On HATCN and F4TCNQ the energy at the nitrile site is still falling steeply at
the lowest height the scan reaches, so the scan is not measuring the site at
all -- it is measuring a point on the way down to it. Nothing in a scan can say
how far down the bottom is. An optimisation can.

Each probe is the same complex with every substrate atom frozen (freeze code
-1) and Ag free (0), started from the scanned point's orbitals. Gaussian moves
Ag in all three coordinates, so the probe also relaxes the lateral position the
scan never does. The result is the true site energy at this level of theory,
and the difference from the fitted scan value is the size of the error the scan
was carrying at that point.

Cu4I4:2 is the control. Its minimum IS bracketed by the scan, so the fitted
value is what the probe should reproduce; agreement there says the parabola
fit is trustworthy where it does bracket, and disagreement says it is not.
"""
import os, re, sys, glob

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
gen = importlib.import_module("122_pregenerate_gjf")
topup = importlib.import_module("124_topup_heights")
from pathgeom import CANDIDATES

E_RE = re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)")
DZ_RE = re.compile(r"dz=([+-]?\d+\.\d+)")


def lowest_height(folder, i):
    """(dz, E) of the lowest scanned height at path point i."""
    best = None
    for f in glob.glob(os.path.join(folder, f"*_t{i}_z*.out")):
        if f.endswith("_ref.out"):
            continue
        txt = open(f, errors="replace").read()
        if "Normal termination" not in txt:
            continue
        e, dz = E_RE.findall(txt), DZ_RE.findall(txt)
        if e and dz and (best is None or float(e[-1]) < best[1]):
            best = (float(dz[0]), float(e[-1]))
    return best


def main():
    root = sys.argv[1]
    want = [(a.split(":")[0], int(a.split(":")[1])) for a in sys.argv[2:]]
    nproc = int(os.environ.get("GAUSS_NPROC", "8"))
    mem = int(os.environ.get("GAUSS_MEM_GB", "48"))
    big = int(os.environ.get("GAUSS_BIG_ATOMS", "60"))

    for tag, fn, rule, mult in CANDIDATES:
        safe = tag.replace("=", "").replace("-", "")
        pts_wanted = [i for t, i in want if t == safe]
        if not pts_wanted:
            continue
        d = os.path.join(root, safe)
        built = topup.rebuild(tag, fn, rule, mult)
        if built is None:
            print(f"  {tag}: cannot rebuild the path")
            continue
        sub_s, sub_x, pts, mult = built
        np_ = nproc * 2 if len(sub_s) + 1 >= big else nproc
        mem_ = mem * 2 if len(sub_s) + 1 >= big else mem
        for i in pts_wanted:
            pos, nrm = pts[i][0], pts[i][1]
            start = lowest_height(d, i)
            dz0 = start[0] if start else 0.0
            agx = pos + dz0 * nrm
            name = f"{safe}_t{i}_relax"
            lines = [f"%Chk={safe}.chk", f"%NProcShared={np_}", f"%Mem={mem_}GB",
                     f"#P {gen.FUNC}/Def2SVP EmpiricalDispersion=GD3BJ "
                     f"Opt=(MaxCycles=60) SCF=({gen.SCF_CHAIN}) NoSymm "
                     f"Guess=Read", "",
                     f"{tag} t={i/(len(pts)-1):.3f} relax Ag from dz={dz0:+.2f}; "
                     f"substrate frozen", "", f"0 {mult}"]
            for a, c in zip(sub_s, sub_x):
                lines.append(f" {a:<2s} -1 {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
            lines.append(f" Ag  0 {agx[0]:14.8f} {agx[1]:14.8f} {agx[2]:14.8f}")
            with open(os.path.join(d, name + ".gjf"), "w") as f:
                f.write("\n".join(lines) + "\n\n")
            print(f"  {name:<22} from dz={dz0:+.2f}"
                  + (f"  (scan E={start[1]:.6f})" if start else ""))


if __name__ == "__main__":
    main()
