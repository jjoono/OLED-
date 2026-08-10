"""Does a tetrahedral seed form the inter-molecular pocket? TCPM monolayer.

THE QUESTION. HATCN's monolayer binds Ag at 1.346 eV, of which +0.53 eV comes from
the adatom sitting in an inter-molecular POCKET and touching nitriles of more than
one molecule (scripts/41: stretching the lattice from 15 to 19 A removes exactly
that much). The isolated molecule is worth only 0.81 eV.

TCPM was picked because an sp3 centre keeps four nitrile anchors while halving the
crystallisation risk (92 -> 35, scripts/53) at a cost of only ~0.10 eV in core
acceptor strength (EA 0.872 vs HATCN 2.711, scripts/54, via a proxy slope of
0.055 eV/eV). But the pocket is the part that proxy cannot see, and it is the part
most at risk: HATCN is a flat disc whose six nitriles all lie in the packing plane,
so a landing adatom is surrounded. TCPM's four nitriles point in four directions at
96-117 degrees, so three of them point AWAY from any adatom on the fourth. The
whole case for this molecule turns on whether a pocket survives that geometry.

Two prior results say it might not. An intramolecular o-dinitrile chelate binds
Ag WORSE than a single nitrile (0.130 vs 0.292 eV, scripts/51), so nitriles do not
help each other just by being nearby. And the pocket is a packing effect, which a
bulky 3D molecule is exactly the kind of thing to prevent.

HOW THE CELL IS CHOSEN, and why not the way scripts/43 chose it. For F4TCNQ the
lattice was tuned to reproduce HATCN's measured 3.06 A adatom-to-neighbour contact,
which was right there: both molecules are planar and the question was which binds
better in the SAME pocket. Here the question is whether a pocket forms at all, so
tuning the cell to produce one would assume the answer. Instead the monolayer is
packed to physical van der Waals contact -- the lattice constant is bisected until
the closest inter-molecular atom pair sits at 3.5 A -- and whatever adatom-to-
neighbour distance results is MEASURED and reported.

Protocol otherwise identical to scripts/39/41/43 so the number is comparable:
CP2K GPW, DZVP-MOLOPT/GTH-PBE, Ag-N rigid scan far to near with wavefunction
carry-over, plus a stretched-lattice control that switches the neighbour contact
off. Every completed CP2K run is checkpointed and pushed, because this container
is snapshot-restored and a slab scan is hours of work.
"""
import os, sys, json, subprocess, shutil
import numpy as np
from ase import Atoms
from ase.io import write as ase_write
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _ckpt import Checkpoint

_s = importlib.util.spec_from_file_location("b", os.path.join(HERE, "38_build_slabs.py"))
B = importlib.util.module_from_spec(_s)
_s.loader.exec_module(B)

OUT = os.path.join(B.BASE, "slabs", "tcpm")
RUNS = os.path.join(B.BASE, "runs")
HA = 27.211386
VACUUM = 16.0
CUTOFF = 350
SCAN_R = (4.0, 3.2, 2.8, 2.5, 2.3, 2.2)
VDW_CONTACT = 3.5        # A, C...C van der Waals contact for the packed monolayer
CONTROL_GAP = 8.0        # A, lattice stretched until the neighbour is out of reach
# The control cell comes out at a = 18.7 A, whose grid needs ~17.7 GB at 350 Ry --
# past the 9 GB budget. Dropping the control triple to 220 Ry brings it under
# (grid memory goes as cutoff^1.5). That is safe here for the same reason
# scripts/40 established: the grid offset cancels in E_b provided every system in
# the difference shares a cell and a grid, and ml_big/ag_big/cx_big do. The
# measured E_b spread across 300-600 Ry there was 0.014 eV against a 0.376 eV
# spread in E_total, so the control's E_b is not meaningfully cutoff-dependent.
CONTROL_CUTOFF = 220

REF = {"HATCN_monolayer": 1.346, "HATCN_isolated": 0.812,
       "HATCN_pocket": 0.534, "F4TCNQ_monolayer": 1.556}


def oriented():
    """TCPM with one nitrile pointing along +z -- the arm an adatom would meet."""
    sym, x = B.read_xyz(os.path.join(RUNS, "TCPM", "xtbopt.xyz"))
    x = x - x.mean(axis=0)
    ns = [i for i, s in enumerate(sym) if s == "N"]
    cen = x.mean(axis=0)
    i_n = max(ns, key=lambda k: np.linalg.norm(x[k] - cen))
    d = np.linalg.norm(x - x[i_n], axis=1); d[i_n] = 9e9
    i_c = int(np.argmin(d))
    axis = x[i_n] - x[i_c]; axis /= np.linalg.norm(axis)
    # rotate `axis` onto +z
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(axis, z); c = float(np.dot(axis, z))
    if np.linalg.norm(v) < 1e-8:
        R = np.eye(3) if c > 0 else -np.eye(3)
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = np.eye(3) + vx + vx @ vx * (1 / (1 + c))
    x = x @ R.T
    return sym, x, ns


def cell_at(a, sym, x, vacuum=VACUUM):
    ext_z = x[:, 2].max() - x[:, 2].min()
    xx = x.copy(); xx[:, 2] -= xx[:, 2].min() - vacuum / 2
    cell = np.array([[a, 0, 0], [0, a, 0], [0, 0, ext_z + vacuum]])
    return Atoms(symbols=sym, positions=xx, cell=cell, pbc=[True, True, True])


def min_intermolecular(at):
    """Closest atom pair between the molecule and any periodic image."""
    p, cell = at.get_positions(), np.asarray(at.get_cell())
    best = np.inf
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            if i == 0 and j == 0:
                continue
            sh = i * cell[0] + j * cell[1]
            dm = np.linalg.norm(p[:, None, :] + sh - p[None, :, :], axis=-1)
            best = min(best, float(dm.min()))
    return best


def pack(sym, x, target=VDW_CONTACT):
    """Bisect the lattice constant until neighbours touch at `target`."""
    lo, hi = 6.0, 30.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if min_intermolecular(cell_at(mid, sym, x)) < target:
            lo = mid
        else:
            hi = mid
    a = 0.5 * (lo + hi)
    return a, cell_at(a, sym, x)


def place_ag(at, i_n, r):
    """Ag on the +z nitrile N, along the C-N axis."""
    at = at.copy()
    p = at.get_positions()
    d = np.linalg.norm(p - p[i_n], axis=1); d[i_n] = 9e9
    i_c = int(np.argmin(d))
    ax = p[i_n] - p[i_c]; ax /= np.linalg.norm(ax)
    at.append(Atoms("Ag", positions=[p[i_n] + r * ax])[0])
    return at


def ag_to_neighbour(at):
    """Distance from the adatom to the nearest atom of a NEIGHBOURING molecule."""
    p, cell = at.get_positions(), np.asarray(at.get_cell())
    ag, rest = p[-1], p[:-1]
    best = np.inf
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            if i == 0 and j == 0:
                continue
            sh = i * cell[0] + j * cell[1]
            best = min(best, float(np.linalg.norm(rest + sh - ag, axis=1).min()))
    return best


def energy(path):
    if not os.path.exists(path):
        return None
    txt = open(path, errors="ignore").read().splitlines()
    for l in reversed(txt):
        if "extrapolated to T->0" in l:
            return float(l.split()[-1])
    for l in reversed(txt):
        if "Total FORCE_EVAL" in l:
            return float(l.split()[-1])
    return None


def build():
    os.makedirs(OUT, exist_ok=True)
    sym, x, ns = oriented()
    i_n = int(np.argmax(x[:, 2] * np.array([1 if s == "N" else -1e9
                                            for s in sym])))
    a, ml = pack(sym, x)
    contact = min_intermolecular(ml)
    print(f"packed monolayer: a = {a:.2f} A, closest inter-molecular pair "
          f"{contact:.2f} A")

    B.write_cp2k(ml, "ml", OUT, cutoff=CUTOFF, title="TCPM monolayer, packed")
    ag = Atoms("Ag", positions=[[0.0, 0.0, ml.get_cell()[2, 2] / 2]],
               cell=ml.get_cell(), pbc=[True, True, True])
    B.write_cp2k(ag, "ag", OUT, cutoff=CUTOFF, uks=True,
                 title="isolated Ag in the TCPM cell")

    print(f"\n{'Ag-N (A)':>9}{'Ag to neighbour molecule (A)':>32}")
    geom = {}
    for r in SCAN_R:
        cx = place_ag(ml, i_n, r)
        nm = f"cx_r{r:.1f}".replace(".", "p")
        ase_write(os.path.join(OUT, f"{nm}.xyz"), cx)
        B.write_cp2k(cx, nm, OUT, cutoff=CUTOFF, uks=True,
                     guess="ATOMIC" if r == SCAN_R[0] else "RESTART",
                     title=f"Ag on TCPM monolayer, Ag-N = {r} A")
        geom[nm] = ag_to_neighbour(cx)
        print(f"{r:>9.1f}{geom[nm]:>26.2f}")

    # stretched control
    ab, big = pack(sym, x, target=CONTROL_GAP)
    print(f"\ncontrol lattice: a = {ab:.2f} A")
    B.write_cp2k(big, "ml_big", OUT, cutoff=CONTROL_CUTOFF,
                 title="TCPM, stretched lattice")
    agb = Atoms("Ag", positions=[[0.0, 0.0, big.get_cell()[2, 2] / 2]],
                cell=big.get_cell(), pbc=[True, True, True])
    B.write_cp2k(agb, "ag_big", OUT, cutoff=CONTROL_CUTOFF, uks=True,
                 title="isolated Ag in the stretched cell")
    cxb = place_ag(big, i_n, 2.2)
    ase_write(os.path.join(OUT, "cx_big.xyz"), cxb)
    B.write_cp2k(cxb, "cx_big", OUT, cutoff=CONTROL_CUTOFF, uks=True,
                 title="Ag on TCPM, stretched lattice, Ag-N = 2.2 A")
    print(f"control Ag to neighbour: {ag_to_neighbour(cxb):.2f} A")

    json.dump({"a_packed": a, "contact": contact, "a_control": ab,
               "ag_to_neighbour": geom,
               "ag_to_neighbour_control": ag_to_neighbour(cxb)},
              open(os.path.join(OUT, "geometry.json"), "w"), indent=2)
    return a, contact, geom


def run():
    env = dict(os.environ,
               CP2K_DATA_DIR="/root/miniforge3/envs/cp2k/share/cp2k/data",
               OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "4"))
    cp2k = "/root/miniforge3/envs/cp2k/bin/cp2k.psmp"
    ck = Checkpoint("dft_seedlayer_screen/runs/tcpm_slab.json", label="TCPM slab")
    order = (["ml", "ag"] + [f"cx_r{r:.1f}".replace(".", "p") for r in SCAN_R]
             + ["ml_big", "ag_big", "cx_big"])
    prev = None
    for n in order:
        if ck.has(n):
            print(f"[skip] {n} = {ck.get(n)}", flush=True)
            prev = n if n.startswith("cx_r") else prev
            continue
        if n.startswith("cx_r") and prev:
            src = os.path.join(OUT, f"{prev}-RESTART.wfn")
            if os.path.exists(src):
                shutil.copy(src, os.path.join(OUT, f"{n}-RESTART.wfn"))
        print(f"[run ] {n}", flush=True)
        subprocess.run([cp2k, "-i", f"{n}.inp", "-o", f"{n}.out"], cwd=OUT, env=env,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        e = energy(os.path.join(OUT, f"{n}.out"))
        if e is None:
            print(f"       FAILED", flush=True)
            continue
        ck.put(n, e)
        print(f"       {e}  [checkpointed]", flush=True)
        if n.startswith("cx_r"):
            prev = n


def report():
    ck = Checkpoint("dft_seedlayer_screen/runs/tcpm_slab.json", push=False)
    geo = {}
    gp = os.path.join(OUT, "geometry.json")
    if os.path.exists(gp):
        geo = json.load(open(gp))
    e_ml, e_ag = ck.get("ml"), ck.get("ag")
    if e_ml is None or e_ag is None:
        print("references not done yet"); return
    print(f"\n{'Ag-N (A)':>9}{'E_b (eV)':>11}{'Ag-neighbour (A)':>19}")
    rows = []
    for r in sorted(SCAN_R):
        nm = f"cx_r{r:.1f}".replace(".", "p")
        e = ck.get(nm)
        if e is None:
            print(f"{r:>9.1f}{'--':>11}"); continue
        eb = (e_ml + e_ag - e) * HA
        rows.append((r, eb))
        d = (geo.get("ag_to_neighbour") or {}).get(nm)
        print(f"{r:>9.1f}{eb:>11.3f}{(f'{d:.2f}' if d else '--'):>19}")
    if not rows:
        return
    rb, eb_max = max(rows, key=lambda t: t[1])

    print(f"\nTCPM monolayer E_b = {eb_max:.3f} eV at Ag-N = {rb:.1f} A")
    print(f"HATCN monolayer    = {REF['HATCN_monolayer']:.3f} eV")
    print(f"F4TCNQ monolayer   = {REF['F4TCNQ_monolayer']:.3f} eV")

    e3 = [ck.get(k) for k in ("ml_big", "ag_big", "cx_big")]
    if None not in e3:
        eb_big = (e3[0] + e3[1] - e3[2]) * HA
        pocket = eb_max - eb_big
        print(f"\nstretched lattice (isolated molecule) E_b = {eb_big:.3f} eV")
        print(f"POCKET CONTRIBUTION = {pocket:+.3f} eV")
        print(f"  HATCN for comparison = {REF['HATCN_pocket']:+.3f} eV")
        print("\nVERDICT")
        if pocket > 0.35:
            print("  A tetrahedral molecule DOES form the pocket. The sp3 core buys")
            print("  glass formation without giving up the packing effect, and TCPM")
            print("  is a real candidate for a crack-free seed.")
        elif pocket > 0.15:
            print("  Partial pocket: weaker than HATCN's but not absent. TCPM trades")
            print("  some packing binding for glass formation; whether that is worth")
            print("  it depends on how much binding the kMC says is actually needed.")
        else:
            print("  NO pocket. The bulky 3D shape that suppresses crystallisation")
            print("  also stops the adatom reaching neighbouring molecules, so TCPM")
            print("  keeps only its single-nitrile binding. The crystallisation and")
            print("  anchoring axes are then NOT independent after all, and the")
            print("  crack-free design has to come from somewhere other than an sp3")
            print("  core -- this is the result that would kill the approach.")
    print("\n  Reminder: E_b here carries no counterpoise correction, as in")
    print("  scripts/39/41/43, so it is an upper bound. The POCKET number is a")
    print("  difference between two runs in the same cell, so it is the more")
    print("  trustworthy of the two.")


if __name__ == "__main__":
    build()
    if "--build-only" not in sys.argv:
        run()
        report()
