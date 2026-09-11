"""Where does the slab's extra 0.3 eV come from -- neighbouring molecules or BSSE?

The Ag-N scan on the HATCN monolayer (a = 15.0 A) gave E_b = 1.35 eV at 2.2 A,
against 1.03 eV from the counterpoise-corrected molecular cluster. Two candidate
explanations were raised:

  (1) in the periodic monolayer the adatom sits in a pocket and also contacts a
      NEIGHBOURING molecule -- geometry check at the minimum puts the nearest
      neighbour-molecule N at 3.06 A, inside van der Waals contact (~3.3 A), so
      this is possible;
  (2) basis-set superposition error, since the slab carries no counterpoise
      correction.

They are separable: stretch the 2D lattice so the molecules no longer touch and
recompute the SAME Ag-N distance. Neighbour contact scales away, BSSE does not.
Runs at a = 15 (reference) and a = 19 A with a reduced cutoff and vacuum so the
larger cell still fits the memory budget; the cutoff cancels in E_b (scripts/40),
and both lattice constants use identical settings so the comparison is clean.
"""
import os, sys, subprocess
import numpy as np
from ase import Atoms
from ase.io import write as ase_write
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location("b", os.path.join(HERE, "38_build_slabs.py"))
B = importlib.util.module_from_spec(_s)
_s.loader.exec_module(B)

OUT = os.path.join(B.BASE, "slabs", "lattice_test")
HA = 27.211386
R_AG_N = 2.2          # the scan minimum
CUTOFF = 300
VACUUM = 16.0
LATTICES = (15.0, 19.0)


def nearest_neighbour_distance(at):
    """Closest approach between the adatom and any atom of a periodic image."""
    p, cell = at.get_positions(), at.get_cell()
    ag, rest = p[-1], p[:-1]
    best = np.inf
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            if i == 0 and j == 0:
                continue
            d = np.linalg.norm(rest + i * np.asarray(cell[0]) + j * np.asarray(cell[1]) - ag, axis=1)
            best = min(best, d.min())
    return float(best)


def build():
    os.makedirs(OUT, exist_ok=True)
    info = {}
    for a in LATTICES:
        ml = B.hatcn_monolayer(a=a, vacuum=VACUUM)
        base = B.ag_on_hatcn(ml, "nitrile")
        # put Ag at exactly R_AG_N from its nitrile N, along the same axis
        p = base.get_positions()
        ns = [k for k, s in enumerate(base.get_chemical_symbols()[:-1]) if s == "N"]
        cen = p[:-1, :2].mean(axis=0)
        i_n = max(ns, key=lambda k: np.linalg.norm(p[k, :2] - cen))
        axis = p[-1] - p[i_n]
        axis /= np.linalg.norm(axis)
        p[-1] = p[i_n] + R_AG_N * axis
        base.set_positions(p)

        tag = f"a{a:.0f}"
        info[a] = nearest_neighbour_distance(base)
        ase_write(os.path.join(OUT, f"cx_{tag}.xyz"), base)
        B.write_cp2k(base, f"cx_{tag}", OUT, cutoff=CUTOFF, uks=True,
                     title=f"Ag on HATCN monolayer, a = {a} A, Ag-N = {R_AG_N} A")
        B.write_cp2k(ml, f"ml_{tag}", OUT, cutoff=CUTOFF,
                     title=f"clean HATCN monolayer, a = {a} A")
        ag = Atoms("Ag", positions=[[0.0, 0.0, ml.get_cell()[2, 2] / 2]],
                   cell=ml.get_cell(), pbc=[True, True, True])
        B.write_cp2k(ag, f"ag_{tag}", OUT, cutoff=CUTOFF, uks=True,
                     title=f"isolated Ag in the a = {a} A cell")
    print(f"{'a (A)':>7}{'nearest neighbour-molecule atom (A)':>38}")
    for a, d in info.items():
        note = "contact" if d < 3.6 else "no contact"
        print(f"{a:>7.0f}{d:>30.2f}   {note}")
    return info


def energy(p):
    if not os.path.exists(p):
        return None
    txt = open(p, errors="ignore").read().splitlines()
    for l in reversed(txt):
        if "extrapolated to T->0" in l:
            return float(l.split()[-1])
    for l in reversed(txt):
        if "Total FORCE_EVAL" in l:
            return float(l.split()[-1])
    return None


def run():
    env = dict(os.environ,
               CP2K_DATA_DIR="/root/miniforge3/envs/cp2k/share/cp2k/data",
               OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "4"))
    cp2k = "/root/miniforge3/envs/cp2k/bin/cp2k.psmp"
    for a in LATTICES:
        for pre in ("ml", "ag", "cx"):
            n = f"{pre}_a{a:.0f}"
            if energy(os.path.join(OUT, f"{n}.out")) is not None:
                print(f"[skip] {n}", flush=True); continue
            print(f"[run ] {n}", flush=True)
            subprocess.run([cp2k, "-i", f"{n}.inp", "-o", f"{n}.out"], cwd=OUT,
                           env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print(f"       {energy(os.path.join(OUT, f'{n}.out'))}", flush=True)


def report(info):
    print(f"\n{'a (A)':>7}{'nn dist':>10}{'E_b (eV)':>11}")
    ebs = {}
    for a in LATTICES:
        e_ml = energy(os.path.join(OUT, f"ml_a{a:.0f}.out"))
        e_ag = energy(os.path.join(OUT, f"ag_a{a:.0f}.out"))
        e_cx = energy(os.path.join(OUT, f"cx_a{a:.0f}.out"))
        if None in (e_ml, e_ag, e_cx):
            print(f"{a:>7.0f}{info[a]:>10.2f}{'incomplete':>11}"); continue
        eb = (e_ml + e_ag - e_cx) * HA
        ebs[a] = eb
        print(f"{a:>7.0f}{info[a]:>10.2f}{eb:>11.3f}")
    if len(ebs) == 2:
        d = ebs[15.0] - ebs[19.0]
        print(f"\n  neighbour contribution = {d:+.3f} eV")
        print(f"  isolated-molecule slab E_b = {ebs[19.0]:.3f} eV")
        print(f"  counterpoise-corrected cluster reference = 1.03 eV")
        print(f"  residual (basis / BSSE / method) = {ebs[19.0]-1.03:+.3f} eV")
        if d > 0.15:
            print("\n  -> A real part of the monolayer's extra binding comes from the")
            print("     adatom contacting more than one molecule. That is physics, not")
            print("     an artefact: on a real HATCN film the adatom sits in an")
            print("     inter-molecular pocket and is anchored by several nitriles.")
        else:
            print("\n  -> Neighbour contact is not the explanation; the gap to the")
            print("     cluster value is basis-set/BSSE, and the slab number should")
            print("     be treated as an upper bound.")


if __name__ == "__main__":
    info = build()
    if "--build-only" not in sys.argv:
        run()
        report(info)
