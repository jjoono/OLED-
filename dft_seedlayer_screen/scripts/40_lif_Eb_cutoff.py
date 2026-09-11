"""Does E_b converge with CUTOFF even though E_total does not?

The LiF CUTOFF scan (scripts/38 output, slabs/cutoff_scan) showed the total
energy drifting by ~5 mHa per 100 Ry all the way to 700 Ry, without settling.
That is not a basis-set error: with NGRIDS 5 and progression factor 3, changing
CUTOFF moves the whole ladder of grid levels (500 -> 500/167/56/19/6 Ry,
700 -> 700/233/78/26/9), so Gaussians get reassigned between levels and the
total energy picks up a discretisation offset.

The claim to test is that this offset is common to all systems sharing a cell and
a grid setting, so it cancels in

    E_b = E(slab) + E(Ag) - E(slab+Ag)

This script measures that directly on LiF(001), which is cheap enough to run at
several cutoffs: if E_b is stable to ~0.01 eV across cutoffs where E_total moves
by 5 mHa (0.14 eV), the cancellation is real and the HATCN production runs can
use a modest cutoff.

It also gives a second useful number: the slab E_b of Ag on LiF, comparable to
the cluster value of 0.25 eV.
"""
import os, sys, subprocess, json
import numpy as np
from ase.io import read as ase_read
from ase import Atoms

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "b", os.path.join(os.path.dirname(os.path.abspath(__file__)), "38_build_slabs.py"))
B = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(B)

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(BASE, "slabs", "lif_Eb_cutoff")
HA = 27.211386
CUTOFFS = (300, 400, 500, 600)


def energy(path):
    if not os.path.exists(path):
        return None
    txt = open(path, errors="ignore").read().splitlines()
    for line in reversed(txt):
        if "extrapolated to T->0" in line:
            return float(line.split()[-1])
    for line in reversed(txt):
        if "Total FORCE_EVAL" in line:
            return float(line.split()[-1])
    return None


def build():
    os.makedirs(OUT, exist_ok=True)
    slab = B.lif_slab()
    slab_ag = B.ag_on_lif(slab)
    # isolated Ag in the SAME cell as the slab -- this is what makes the grid
    # offset cancel; a different cell would reintroduce it.
    ag = Atoms("Ag", positions=[[0.0, 0.0, slab.get_cell()[2, 2] / 2]],
               cell=slab.get_cell(), pbc=[True, True, True])
    for c in CUTOFFS:
        B.write_cp2k(slab, f"slab_c{c}", OUT, cutoff=c,
                     title=f"LiF(001) clean, CUTOFF {c}")
        B.write_cp2k(slab_ag, f"slabag_c{c}", OUT, cutoff=c, uks=True,
                     title=f"LiF(001) + Ag, CUTOFF {c}")
        B.write_cp2k(ag, f"ag_c{c}", OUT, cutoff=c, uks=True,
                     title=f"isolated Ag in the LiF cell, CUTOFF {c}")
    print(f"built {3*len(CUTOFFS)} inputs in {OUT}")


def run():
    env = dict(os.environ,
               CP2K_DATA_DIR="/root/miniforge3/envs/cp2k/share/cp2k/data",
               OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "4"))
    cp2k = "/root/miniforge3/envs/cp2k/bin/cp2k.psmp"
    for c in CUTOFFS:
        for tag in ("slab", "ag", "slabag"):
            n = f"{tag}_c{c}"
            if energy(os.path.join(OUT, f"{n}.out")) is not None:
                print(f"[skip] {n}", flush=True)
                continue
            print(f"[run ] {n}", flush=True)
            subprocess.run([cp2k, "-i", f"{n}.inp", "-o", f"{n}.out"],
                           cwd=OUT, env=env,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print(f"       E = {energy(os.path.join(OUT, f'{n}.out'))}", flush=True)


def report():
    print(f"\n{'CUTOFF':>8}{'E_slab (Ha)':>18}{'E_slab+Ag (Ha)':>18}"
          f"{'E_Ag (Ha)':>16}{'E_b (eV)':>11}")
    rows = []
    for c in CUTOFFS:
        e_s = energy(os.path.join(OUT, f"slab_c{c}.out"))
        e_sa = energy(os.path.join(OUT, f"slabag_c{c}.out"))
        e_a = energy(os.path.join(OUT, f"ag_c{c}.out"))
        if None in (e_s, e_sa, e_a):
            print(f"{c:>8}{'incomplete':>18}")
            continue
        eb = (e_s + e_a - e_sa) * HA
        rows.append((c, e_s, eb))
        print(f"{c:>8}{e_s:>18.6f}{e_sa:>18.6f}{e_a:>16.6f}{eb:>11.3f}")

    if len(rows) >= 2:
        de_tot = (max(r[1] for r in rows) - min(r[1] for r in rows)) * HA
        de_b = max(r[2] for r in rows) - min(r[2] for r in rows)
        print(f"\n  spread of E_total across cutoffs : {de_tot:.3f} eV")
        print(f"  spread of E_b     across cutoffs : {de_b:.3f} eV")
        if de_b < 0.1 * max(de_tot, 1e-9):
            print("  -> the grid offset cancels in E_b, as expected. A modest")
            print("     cutoff is fine for the production runs provided every")
            print("     system in the difference uses the same cell and grid.")
        else:
            print("  -> E_b is NOT protected by cancellation; the cutoff has to")
            print("     be converged on E_b itself before trusting any number.")
        print(f"\n  cluster reference for Ag on LiF: 0.25 eV")
        json.dump({"cutoffs": [r[0] for r in rows],
                   "E_slab_Ha": [r[1] for r in rows],
                   "E_b_eV": [r[2] for r in rows],
                   "spread_Etot_eV": de_tot, "spread_Eb_eV": de_b},
                  open(os.path.join(OUT, "lif_Eb_cutoff.json"), "w"), indent=2)


if __name__ == "__main__":
    build()
    if "--build-only" not in sys.argv:
        run()
    report()
