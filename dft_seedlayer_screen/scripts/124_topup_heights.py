"""Add the one height each path point still needs, reusing the run already done.

    python scripts/124_topup_heights.py <returned gaussian_jobs folder>

The v10 scan used dz = -0.20, +0.20, +0.60 A. On 75 of its 120 path points the
lowest of those three is an END of the scan, so the minimum is outside it and
the reader cannot fit anything -- it falls back to the lowest raw point, which
is what the three heights were meant to stop. The misses are lopsided: 53 want
to go lower and 22 higher, and the energy is still falling by a median 60 meV
(worst 389 meV) at the edge.

Widening the window for everything would cost another two thirds of a campaign
to re-measure heights that are already bracketed. Instead each point gets ONE
more height, on the side its own data points to. That is about 125 jobs rather
than 670, and it brackets every point that a fourth height can bracket.

The extra jobs go into the folders they belong to and read the same .chk, so
they continue the orbital chain rather than starting a new one. Keep the .chk
files from the finished run in place.
"""
import os, re, sys, glob, collections

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
gen = importlib.import_module("122_pregenerate_gjf")
from pathgeom import (CANDIDATES, STRUCT, ZSCAN, contact, destination,
                      dimer_frames, equivalents, frames, geometry, read_xyz,
                      sanity)

STEP = 0.4                      # same spacing as the existing scan
E_RE = re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)")
DZ_RE = re.compile(r"dz=([+-]?\d+\.\d+)")


def done_heights(folder):
    """{path point: {dz: energy}} from the outputs already in this folder."""
    pts = collections.defaultdict(dict)
    for f in glob.glob(os.path.join(folder, "*.out")):
        stem = os.path.basename(f)[:-4]
        if stem.endswith("_ref"):
            continue
        txt = open(f, errors="replace").read()
        if "Normal termination" not in txt:
            continue
        e, dz = E_RE.findall(txt), DZ_RE.findall(txt)
        if not e or not dz:
            continue
        i = int(stem.split("_t")[1].split("_z")[0])
        pts[i][round(float(dz[0]), 2)] = float(e[-1])
    return pts


def rebuild(tag, fn, rule, mult):
    """The same geometry the campaign used. The generator is deterministic."""
    syms, xyz = read_xyz(os.path.join(STRUCT, fn))
    if sanity(syms, xyz, tag):
        return None
    sub_s, sub_x, ag, anchor, nrm = geometry(syms, xyz)
    near_ok = [j for j in equivalents(sub_s, sub_x, anchor)
               if np.linalg.norm(ag - sub_x[j]) > 1.5 * contact(sub_s[j])]
    if near_ok:
        near, _, _ = destination(sub_s, sub_x, ag, anchor, "auto")
        span = float(np.linalg.norm(sub_x[near] - sub_x[anchor]))
        npath = int(np.clip(round(span / gen.SPACING) + 1,
                            gen.NPATH_MIN, gen.NPATH_MAX))
        pts, cls, _ = frames(sub_s, sub_x, ag, anchor, "auto", npath)
    else:
        if 2 * len(sub_s) + 1 > gen.DIMER_MAX_ATOMS:
            return None
        _, _, _, span = dimer_frames(sub_s, sub_x, ag, anchor, nrm,
                                     gen.NPATH_MIN)
        npath = int(np.clip(round(span / gen.SPACING) + 1,
                            gen.NPATH_MIN, gen.NPATH_MAX))
        sub_s, sub_x, pts, span = dimer_frames(sub_s, sub_x, ag, anchor,
                                               nrm, npath)
    return sub_s, sub_x, pts, mult


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "gaussian_jobs")
    nproc = int(os.environ.get("GAUSS_NPROC", "8"))
    mem = int(os.environ.get("GAUSS_MEM_GB", "48"))
    big = int(os.environ.get("GAUSS_BIG_ATOMS", "60"))

    made, skipped, unfixable = 0, 0, []
    for tag, fn, rule, mult in CANDIDATES:
        safe = tag.replace("=", "").replace("-", "")
        d = os.path.join(root, safe)
        if not os.path.isdir(d):
            continue
        built = rebuild(tag, fn, rule, mult)
        if built is None:
            continue
        sub_s, sub_x, pts, mult = built
        have = done_heights(d)
        if len(have) != len(pts):
            print(f"  {tag}: {len(have)} points measured but {len(pts)} in the "
                  f"rebuilt path -- skipped, the geometry does not match")
            continue

        np_ = nproc * 2 if len(sub_s) + 1 >= big else nproc
        mem_ = mem * 2 if len(sub_s) + 1 >= big else mem
        lim = np.array([contact(x) for x in sub_s])
        order, wall = [], []
        for i, (pos, nrm) in enumerate(pts):
            b = have[i]
            lo = min(b, key=b.get)
            if min(b) < lo < max(b):
                skipped += 1
                continue                      # already bracketed
            dz = round(lo - STEP if lo == min(b) else lo + STEP, 2)
            agx = pos + dz * nrm
            if float((np.linalg.norm(sub_x - agx, axis=1) / lim).min()) < 0.85:
                # Going lower would press Ag into the repulsive wall, where the
                # curve is not a parabola. Nothing a fourth height can do.
                wall.append(i)
                continue
            name = f"{safe}_t{i}_z{'m' if dz < 0 else 'p'}{abs(dz):.1f}".replace(".", "")
            with open(os.path.join(d, name + ".gjf"), "w") as f:
                f.write(gen.gjf(list(sub_s) + ["Ag"],
                                np.vstack([sub_x, agx]),
                                f"{tag} t={i/(len(pts)-1):.3f} dz={dz:+.2f} topup",
                                np_, mem_, mult, chk=f"{safe}.chk", read=True))
            order.append(name)
            made += 1
        if wall:
            unfixable.append((tag, wall))
        if order:
            with open(os.path.join(d, "ORDER.txt"), "a") as f:
                f.write("\n".join(order) + "\n")
            print(f"  {safe:<10} +{len(order):>2} heights")

    print(f"\n{made} top-up jobs; {skipped} points were already bracketed")
    if unfixable:
        print("\nat the repulsive wall -- a fourth height cannot bracket these:")
        for tag, w in unfixable:
            print(f"  {tag:<10} points {w}")


if __name__ == "__main__":
    main()
