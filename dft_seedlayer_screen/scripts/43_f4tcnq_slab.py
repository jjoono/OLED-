"""F4TCNQ monolayer: does it overtake HATCN once the monolayer effect is included?

WHY THIS ONE. The cluster screen put F4TCNQ at 0.97 eV against HATCN's 1.03 eV --
a 0.06 eV gap, well inside the method error. The slab test (scripts/41) then found
that a periodic HATCN monolayer adds +0.53 eV over the isolated molecule, because
the adatom sits in an inter-molecular pocket and coordinates to nitriles of more
than one molecule. F4TCNQ has the same ingredients: planar, terminal nitriles on
the molecular rim, dense packing. So the correction should be comparable and the
ranking inside the nitrile family is not decided by the cluster numbers.

FAIR COMPARISON. The size of the pocket effect depends on how close the
neighbouring molecule is, so the two systems must be compared at the SAME
inter-molecular contact distance, not at the same lattice constant. HATCN
(a = 15.0 A, N-to-N molecular extent 12.03 A) put the nearest neighbour-molecule
atom 3.06 A from the adatom. F4TCNQ has a 9.53 A N-to-N extent, so the cell is
chosen to reproduce that 3.06 A and the achieved value is measured and reported.

  a: along the long molecular axis, the direction the scanned nitrile points
  b: along the short axis

Runs the same protocol as scripts/39+41: Ag-N rigid scan for the binding curve,
plus a stretched-lattice control that switches the neighbour contact off.
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

OUT = os.path.join(B.BASE, "slabs", "f4tcnq")
HA = 27.211386
VACUUM = 18.0
CUTOFF = 350
SCAN_R = (3.5, 3.0, 2.6, 2.4, 2.2, 2.0)
TARGET_GAP = 3.06        # A, matched to the HATCN monolayer (scripts/41)


def _oriented():
    """F4TCNQ with the long axis on x, short on y, molecular plane = xy."""
    sym, x = B.read_xyz(os.path.join(B.STR, "F4TCNQ.xyz"))
    x = x - x.mean(axis=0)
    u, s, vt = np.linalg.svd(x)
    x = x @ vt.T
    x[:, 2] -= x[:, 2].mean()
    return sym, x


def _cell_at(scale, vacuum=VACUUM):
    """Rectangular cell whose in-plane size is `scale` x the molecular bbox."""
    sym, x = _oriented()
    ext_a = x[:, 0].max() - x[:, 0].min()
    ext_b = x[:, 1].max() - x[:, 1].min()
    a, b = ext_a * scale, ext_b * scale
    cell = np.array([[a, 0, 0], [0, b, 0], [0, 0, vacuum]])
    xx = x.copy(); xx[:, 2] += vacuum / 2
    return Atoms(symbols=sym, positions=xx, cell=cell, pbc=[True, True, True]), a, b


def monolayer(target_nn=TARGET_GAP, r_probe=2.2, vacuum=VACUUM):
    """Cell sized so the adatom-to-neighbouring-molecule distance equals
    `target_nn` at Ag-N = `r_probe`.

    Setting the lattice constant from a geometric estimate does not work for
    F4TCNQ: its nitriles splay away from the long axis (the two most distant N
    are diagonal, 9.53 A, while the long-axis extent is only 8.49 A), so an
    a = extent + gap rule put the adatom 1.5-2.5 A from the next molecule --
    inside it. Solve for the cell that reproduces the measured HATCN contact
    instead, which is also what makes the two systems comparable.
    """
    lo, hi = 1.0, 4.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        at, a, b = _cell_at(mid, vacuum)
        d = nn_dist(place_ag(at, r_probe))
        if d < target_nn:
            lo = mid
        else:
            hi = mid
    at, a, b = _cell_at(0.5 * (lo + hi), vacuum)
    return at, a, b


def place_ag(at, r):
    """Ag on the outermost nitrile N, along the C-N axis, at distance r."""
    at = at.copy()
    sym, p = at.get_chemical_symbols(), at.get_positions()
    ns = [i for i, t in enumerate(sym) if t == "N"]
    cen = p[:, :2].mean(axis=0)
    i_n = max(ns, key=lambda k: np.linalg.norm(p[k, :2] - cen))
    # C-N axis of that nitrile
    d = np.linalg.norm(p - p[i_n], axis=1); d[i_n] = 9e9
    i_c = int(np.argmin(d))
    ax = p[i_n] - p[i_c]; ax /= np.linalg.norm(ax)
    at.append(Atoms("Ag", positions=[p[i_n] + r * ax])[0])
    return at


def nn_dist(at):
    p, cell = at.get_positions(), at.get_cell()
    ag, rest = p[-1], p[:-1]
    best = np.inf
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            if i == 0 and j == 0:
                continue
            sh = i * np.asarray(cell[0]) + j * np.asarray(cell[1])
            best = min(best, np.linalg.norm(rest + sh - ag, axis=1).min())
    return float(best)


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


def build():
    os.makedirs(OUT, exist_ok=True)
    ml, a, b = monolayer()
    print(f"monolayer cell: a = {a:.2f} A, b = {b:.2f} A, vacuum = {VACUUM} A")

    # references
    B.write_cp2k(ml, "ml", OUT, cutoff=CUTOFF, title="clean F4TCNQ monolayer")
    ag = Atoms("Ag", positions=[[0.0, 0.0, ml.get_cell()[2, 2] / 2]],
               cell=ml.get_cell(), pbc=[True, True, True])
    B.write_cp2k(ag, "ag", OUT, cutoff=CUTOFF, uks=True,
                 title="isolated Ag in the F4TCNQ cell")

    # binding curve, far -> near with wavefunction carry-over
    print(f"\n{'r (A)':>7}{'nearest neighbour-molecule atom (A)':>38}")
    for r in SCAN_R:
        cx = place_ag(ml, r)
        nm = f"cx_r{r:.1f}".replace(".", "p")
        ase_write(os.path.join(OUT, f"{nm}.xyz"), cx)
        B.write_cp2k(cx, nm, OUT, cutoff=CUTOFF, uks=True,
                     guess="ATOMIC" if r == SCAN_R[0] else "RESTART",
                     title=f"Ag on F4TCNQ nitrile, Ag-N = {r} A")
        print(f"{r:>7.1f}{nn_dist(cx):>30.2f}")

    # stretched-lattice control: neighbour contact switched off
    big, ab, bb = monolayer(target_nn=7.0)
    print(f"\ncontrol cell: a = {ab:.2f} A, b = {bb:.2f} A")
    B.write_cp2k(big, "ml_big", OUT, cutoff=300, title="clean F4TCNQ, stretched lattice")
    agb = Atoms("Ag", positions=[[0.0, 0.0, big.get_cell()[2, 2] / 2]],
                cell=big.get_cell(), pbc=[True, True, True])
    B.write_cp2k(agb, "ag_big", OUT, cutoff=300, uks=True,
                 title="isolated Ag in the stretched cell")
    cxb = place_ag(big, 2.2)
    ase_write(os.path.join(OUT, "cx_big.xyz"), cxb)
    B.write_cp2k(cxb, "cx_big", OUT, cutoff=300, uks=True,
                 title="Ag on F4TCNQ, stretched lattice, Ag-N = 2.2 A")
    print(f"control nearest neighbour-molecule atom: {nn_dist(cxb):.2f} A")
    return ml


def run():
    env = dict(os.environ,
               CP2K_DATA_DIR="/root/miniforge3/envs/cp2k/share/cp2k/data",
               OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "4"))
    cp2k = "/root/miniforge3/envs/cp2k/bin/cp2k.psmp"
    order = ["ml", "ag"] + [f"cx_r{r:.1f}".replace(".", "p") for r in SCAN_R] \
            + ["ml_big", "ag_big", "cx_big"]
    prev = None
    for n in order:
        if energy(os.path.join(OUT, f"{n}.out")) is not None:
            print(f"[skip] {n}", flush=True)
            prev = n if n.startswith("cx_r") else prev
            continue
        if n.startswith("cx_r") and prev:
            src = os.path.join(OUT, f"{prev}-RESTART.wfn")
            if os.path.exists(src):
                import shutil
                shutil.copy(src, os.path.join(OUT, f"{n}-RESTART.wfn"))
        print(f"[run ] {n}", flush=True)
        subprocess.run([cp2k, "-i", f"{n}.inp", "-o", f"{n}.out"], cwd=OUT, env=env,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"       {energy(os.path.join(OUT, f'{n}.out'))}", flush=True)
        if n.startswith("cx_r"):
            prev = n


def report():
    e_ml, e_ag = energy(f"{OUT}/ml.out"), energy(f"{OUT}/ag.out")
    print(f"\nE(F4TCNQ monolayer) = {e_ml}\nE(isolated Ag)      = {e_ag}")
    if None in (e_ml, e_ag):
        print("references missing"); return
    print(f"\n{'Ag-N (A)':>10}{'E_b (eV)':>11}")
    rows = []
    for r in sorted(SCAN_R):
        e = energy(f"{OUT}/cx_r{r:.1f}".replace(".", "p") + ".out")
        if e is None:
            print(f"{r:>10.1f}{'--':>11}"); continue
        eb = (e_ml + e_ag - e) * HA
        rows.append((r, eb))
        print(f"{r:>10.1f}{eb:>11.3f}")
    if not rows:
        return
    rb, eb_max = max(rows, key=lambda t: t[1])
    print(f"\nF4TCNQ monolayer E_b = {eb_max:.3f} eV at Ag-N = {rb:.1f} A")
    print(f"HATCN  monolayer E_b = 1.346 eV at Ag-N = 2.2 A   (scripts/39)")
    print(f"cluster values:  F4TCNQ 0.97 eV,  HATCN 1.03 eV")

    eb_big = None
    e3 = [energy(f"{OUT}/{k}.out") for k in ("ml_big", "ag_big", "cx_big")]
    if None not in e3:
        eb_big = (e3[0] + e3[1] - e3[2]) * HA
        print(f"\nstretched-lattice (isolated molecule) E_b = {eb_big:.3f} eV")
        print(f"neighbour contribution = {eb_max - eb_big:+.3f} eV")
        print(f"  HATCN for comparison  = +0.534 eV")

    print("\nVERDICT")
    if eb_max > 1.346:
        print(f"  F4TCNQ OVERTAKES HATCN by {eb_max-1.346:.3f} eV in the monolayer.")
        print("  The cluster ranking (HATCN > F4TCNQ) does not survive, and the")
        print("  paper's lead material has to be reconsidered.")
    else:
        print(f"  HATCN stays ahead by {1.346-eb_max:.3f} eV in the monolayer, so the")
        print("  cluster ranking survives inside the nitrile family.")
    print("  NOTE: absorption decides this too -- F4TCNQ's radical anion is")
    print("  strongly coloured (see scripts/44). Anchoring is not the only axis.")


if __name__ == "__main__":
    build()
    if "--build-only" not in sys.argv:
        run()
        report()
