"""Ag adatom surface diffusion barrier E_d across the whole candidate set.

E_b has been available for twenty-odd candidates for a while; E_d existed for
four. A seed layer needs both -- strong binding alone still lets the adatom walk
until it meets another and nucleates an island -- so the screening cannot be
read without this axis.

PROTOCOL, identical for every candidate.

  start   the relaxed Ag complex, structures/<tag>.xyz
  place   a straight line between two sites cuts THROUGH the molecule at its
          midpoint, which drives Ag to within a bond length of the substrate and
          returns barriers of hundreds of eV. So at every step the adatom is
          pushed back out along the surface normal until nothing is closer than
          D_MIN, and the height scan starts from there.
  path    Ag is dragged laterally from its bound position to a destination on
          the same molecule -- the NEAREST equivalent site, not any equivalent
          site, since diffusion proceeds by the cheapest hop available. At every
          lateral position its height is scanned along the local surface normal
          and the minimum kept, which is the relaxed-Ag, frozen-substrate
          approximation used in script 24. Points are placed about 1 A apart, so
          a long path gets more of them rather than a coarser description.
  grid    deliberately coarse -- three to five lateral points, two heights. One
          SCF on a 31-atom Ag complex runs about ten minutes on the four cores
          available, so the fine grid script 24 used for a single system is out
          of reach for twenty-six. The barrier this returns is a lower bound:
          a coarse path can miss the true saddle but cannot invent one.
  E_d     max(E along path) - E(start), with the START GIVEN THE SAME TREATMENT
          as every other point. The structures are DFT minima, so a method that
          did not produce them does not have its own minimum there; leaving the
          start unrelaxed while every other point gets a height scan makes the
          reference the highest point on the path and returns a barrier of
          exactly zero.

The destination depends on what the molecule offers, and the class is reported
alongside the number because the three are not interchangeable:

  site2site   another equivalent anchor exists on the molecule (HATCN's six
              nitriles, F4TCNQ's four). This is a true intramolecular hop and
              the only class where E_d means what Venables means by it.
  toface      one anchor only, so the adatom's escape route is across the
              molecular pi face. The barrier is real but the destination is a
              shallow or repulsive site, and what it measures is closer to
              partial detachment than to a hop.
  chelate     the adatom sits in a bidentate pocket; the path runs out of the
              pocket toward the nearest ring centroid.

Comparing a site2site number against a toface number is the same error as
putting a monodentate E_b next to a bidentate one, which is why the class
travels with the value everywhere it is used.

SCF is the practical obstacle, not the arithmetic. An open-shell Ag on a strong
acceptor has near-degenerate states, and plain DIIS oscillates indefinitely
rather than converging. Every point therefore walks a ladder of settings and
takes the first rung that converges -- level-shifted and damped, then harder
damped, then second-order SCF, then SOSCF from a GWH guess. maxiter is low on
the early rungs so a hopeless attempt fails in a minute instead of grinding
through four hundred iterations, and the rungs actually used are recorded with
the result so a barrier obtained only by the desperate settings can be spotted.
"""
import json, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import psi4

# The path geometry lives in pathgeom so this driver and the Gaussian one
# (script 102) cannot drift apart.
from pathgeom import (CANDIDATES, D_MIN, H2EV, NPATH_MAX, NPATH_MIN, RUNS,
                      SPACING, STRUCT, ZSCAN, destination, geometry, place,
                      read_xyz, sanity)

OUT = os.path.join(RUNS, "diffusion_barriers_all.json")




def n_cores():
    """Cores actually available. sched_getaffinity respects a scheduler's or a
    container's mask and is what SLURM sets, but it does not exist on Windows,
    where the workstation route lands -- so fall back rather than crash on
    import."""
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 4


# Both read from the machine rather than being hard-coded, so the same file runs
# on a Windows workstation, in a container, and under a scheduler. PSI4_MEM_GB is
# set by scripts/run_Ed.slurm from the SLURM request.
psi4.set_memory(f"{os.environ.get('PSI4_MEM_GB', '10')} GB")
psi4.set_num_threads(n_cores())
psi4.core.set_output_file(os.path.join(RUNS, "psi4_diff_all.out"), False)

BASE_SCF = {"basis": "def2-svp", "scf_type": "df", "reference": "uks"}
LADDER = [
    dict(BASE_SCF, maxiter=150, guess="sad", level_shift=1.0,
         level_shift_cutoff=1e-3, damping_percentage=15.0),
    dict(BASE_SCF, maxiter=200, guess="sad", level_shift=2.0,
         level_shift_cutoff=1e-2, damping_percentage=40.0),
    dict(BASE_SCF, maxiter=200, guess="sad", soscf=True,
         soscf_start_convergence=1e-2),
    dict(BASE_SCF, maxiter=250, guess="gwh", soscf=True,
         soscf_start_convergence=1e-2, damping_percentage=25.0),
]

def gstr(syms, xyz, mult):
    st = f"0 {mult}\n"
    for a, c in zip(syms, xyz):
        st += f"{a} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n"
    return st + "symmetry c1\nno_reorient\nno_com\n"


def sp(syms, xyz, mult):
    """First converging rung of the SCF ladder; (None, -1) if every rung fails."""
    g = gstr(syms, xyz, mult)
    for rung, opts in enumerate(LADDER):
        try:
            psi4.set_options(opts)
            e = psi4.energy("pbe-d3bj", molecule=psi4.geometry(g))
            psi4.core.clean()
            return e, rung
        except Exception:
            psi4.core.clean()
            continue
    return None, -1


def barrier(tag, fn, rule, mult):
    syms, xyz = read_xyz(os.path.join(STRUCT, fn))
    why = sanity(syms, xyz, tag)
    if why:
        print(f"[{tag}] REFUSED: {why}", flush=True)
        print(f"[{tag}] re-optimise the complex before asking for a barrier",
              flush=True)
        return None
    sub_s, sub_x, ag, anchor, nrm = geometry(syms, xyz)
    dest, cls = destination(sub_s, sub_x, ag, anchor, rule)
    print(f"[{tag}] {len(sub_s)} substrate atoms, anchor {sub_s[anchor]}{anchor}, "
          f"path {cls}, span {np.linalg.norm(dest-ag):.2f} A", flush=True)

    span = float(np.linalg.norm(dest - ag))
    npath = int(np.clip(round(span / SPACING) + 1, NPATH_MIN, NPATH_MAX))
    E, rungs = [], []
    for t in np.linspace(0.0, 1.0, npath):
        pos = place(sub_x, (1 - t) * ag + t * dest, nrm)
        best = None
        for dz in ZSCAN:
            trial = pos + dz * nrm
            e, rung = sp(sub_s + ["Ag"], np.vstack([sub_x, trial]), mult)
            if e is None:
                print(f"    t={t:.2f} dz={dz:+.2f} SCF failed on every rung", flush=True)
                continue
            rungs.append(rung)
            if best is None or e < best:
                best = e
        E.append(best)
        print(f"    t={t:.2f}  E={'FAIL' if best is None else f'{best:.6f}'}", flush=True)

    ok = [(i, e) for i, e in enumerate(E) if e is not None]
    if E[0] is None or len(ok) < 3:
        print(f"[{tag}] insufficient converged points", flush=True)
        return None
    ed = (max(e for _, e in ok) - E[0]) * H2EV
    print(f"[{tag}] E_d = {ed:.3f} eV  ({cls})  SCF rungs {sorted(set(rungs))}",
          flush=True)
    return {"E_d_eV": round(ed, 4), "class": cls, "anchor": sub_s[anchor],
            "scf_rungs": sorted(set(rungs)),
            "n_atoms": len(sub_s), "path_A": round(float(np.linalg.norm(dest - ag)), 3),
            "points": [None if e is None else round(e, 8) for e in E]}


def main():
    only = sys.argv[1:] or None
    res = {}
    if os.path.exists(OUT):
        res = json.load(open(OUT))
    for tag, fn, rule, mult in CANDIDATES:
        if only and tag not in only:
            continue
        if tag in res and res[tag] is not None:
            print(f"[{tag}] already done, skipping", flush=True)
            continue
        if not os.path.exists(os.path.join(STRUCT, fn)):
            print(f"[{tag}] missing {fn}", flush=True)
            continue
        t0 = time.time()
        try:
            res[tag] = barrier(tag, fn, rule, mult)
        except Exception as exc:
            print(f"[{tag}] aborted: {type(exc).__name__}: {exc}", flush=True)
            res[tag] = None
        print(f"[{tag}] {time.time()-t0:.0f} s\n", flush=True)
        json.dump(res, open(OUT, "w"), indent=1)      # checkpoint after every system
    print(json.dumps({k: (v or {}).get("E_d_eV") for k, v in res.items()}, indent=1))


if __name__ == "__main__":
    main()
