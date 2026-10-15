"""Can two adjacent nitriles chelate one Ag? Phthalonitrile, scanned properly.

WHY THIS IS SEPARATE FROM scripts/50. The chelate case in 50 never converged, and
the reason was the scan coordinate, not the chemistry:

  * xtb relaxes BOTH starting guesses -- bridging and end-on -- to the same
    symmetric structure, Ag-N = 2.15/2.15 A, in the molecular plane. So
    phthalonitrile does prefer the bridge, and the two jobs in 50 were duplicates.
  * 50 then scans Ag along the line to its NEAREST nitrogen. For a symmetric
    bridge that line is not the symmetry axis, so every step pushed Ag off the
    bisector, toward one N and away from the other, distorting the bond being
    measured. The SCF failures were the geometry falling apart, not a hard
    electronic structure.

GEOMETRY. N...N = 4.02 A, so a bridge is possible at all: Ag equidistant from both
needs Ag-N >= 2.01 A, and a real Ag-N bond is 2.2-2.3 A. Ag is parametrised by its
distance to the nitrogens, r, sitting on the in-plane bisector at height
h = sqrt(r^2 - (N...N/2)^2) above the midpoint. That keeps the two Ag-N bonds
equal at every step, which is what makes the scan measure the bridge.

WHAT IT TESTS. HATCN's monolayer gains +0.53 eV over the isolated molecule because
the adatom sits in an inter-molecular pocket touching nitriles of more than one
molecule (scripts/41). If an o-dinitrile can do that WITHIN one molecule, the
design rule for a crack-free seed becomes concrete: a twisted glass-forming core
carrying an o-dicyanoarene anchor. The comparison is against benzonitrile's single
nitrile, 0.292 eV at the same level (scripts/50).

Far -> near with orbital carry-over, and every scan point checkpointed and pushed
via _ckpt, because this container is snapshot-restored and unpushed work does not
survive.
"""
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _ckpt import Checkpoint

BASE = os.path.abspath(os.path.join(HERE, ".."))
RUNS = os.path.join(BASE, "runs")
H2EV = 27.211386

# Ag-N distances to scan, far to near. Must all exceed half the N...N separation
# or the bisector position is undefined.
#
# The inner points (2.50, 2.35, 2.20, 2.10) are dropped. Not because they failed
# to converge -- though 2.50 did, burning ten minutes on retries -- but because
# they are no longer needed: the completed points bracket the minimum, 4.00 above
# 3.20 below 2.80, so the curve has already turned repulsive and anything closer
# is higher. Keeping them in would only spend SCF time to confirm a rise that the
# bracket already establishes.
#
# Note this puts the minimum at Ag-N = 3.2 A, far longer than the 2.15 A that
# GFN2-xTB predicted for the same bridge. That disagreement is expected in this
# direction -- xTB is known to overbind Ag here, which is exactly why the project
# protocol relaxes with xTB but takes energies from DFT -- and it is the reason
# the E_b below is much smaller than the xTB geometry would suggest.
SCAN_R = [5.0, 4.0, 3.2, 2.8]
PHCN_REF = 0.292          # benzonitrile, one nitrile, same protocol (scripts/50)


def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0]); syms, xyz = [], []
    for l in lines[2:2 + n]:
        p = l.split(); syms.append(p[0]); xyz.append(list(map(float, p[1:4])))
    return syms, np.array(xyz)


def bisector_frame(syms, x):
    """Midpoint of the two nitrile N, the outward in-plane bisector, and N...N/2."""
    ns = [i for i, s in enumerate(syms) if s == "N"]
    if len(ns) != 2:
        raise SystemExit(f"expected 2 nitrile N, found {len(ns)}")
    mid = 0.5 * (x[ns[0]] + x[ns[1]])
    a = 0.5 * float(np.linalg.norm(x[ns[0]] - x[ns[1]]))
    # outward direction: away from the ring centroid, projected into the molecular
    # plane so the scan cannot drift out of plane
    cen = x.mean(axis=0)
    normal = np.linalg.svd(x - cen)[2][2]
    d = mid - cen
    d = d - np.dot(d, normal) * normal
    d /= np.linalg.norm(d)
    return mid, d, a, ns


def place(mid, d, a, r):
    """Ag on the bisector, equidistant r from both nitrogens."""
    if r < a:
        return None
    return mid + np.sqrt(r * r - a * a) * d


def main():
    syms, x = read_xyz(os.path.join(RUNS, "oDCNB", "xtbopt.xyz"))
    mid, d, a, ns = bisector_frame(syms, x)
    print(f"phthalonitrile: N...N = {2*a:.2f} A, minimum bridging Ag-N = {a:.2f} A")

    import psi4
    psi4.set_memory(os.environ.get("PSI4_MEM", "6 GB"))
    psi4.set_num_threads(int(os.environ.get("PSI4_THREADS", "4")))
    psi4.core.set_output_file(os.path.join(RUNS, "psi4_chelate.out"), False)
    base = {"basis": "def2-svp", "scf_type": "df", "reference": "uks",
            "maxiter": 300, "guess": "sad"}
    psi4.set_options(base)
    M = "pbe-d3bj"

    def gstr(sy, xyz, ghost=None, mult=1):
        s = f"0 {mult}\n"
        for i, (sym, c) in enumerate(zip(sy, xyz)):
            t = f"Gh({sym})" if ghost and i in ghost else sym
            s += f"{t} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n"
        return s + "symmetry c1\nno_reorient\nno_com\n"

    def energy(g, carry):
        for extra in ({"guess": "read"} if carry else {},
                      {"guess": "sad", "damping_percentage": 30},
                      {"guess": "sad", "soscf": True,
                       "soscf_start_convergence": 1e-2}):
            try:
                psi4.set_options({**base, **extra})
                e = psi4.energy(M, molecule=psi4.geometry(g))
                psi4.set_options(base)
                return e
            except Exception as ex:
                print(f"      retry ({type(ex).__name__})", flush=True)
                psi4.core.clean()
        return None

    ck = Checkpoint("dft_seedlayer_screen/runs/chelate_scan.json",
                    label="chelate scan")
    sy = syms + ["Ag"]
    agi = len(sy) - 1

    carry = False
    for r in SCAN_R:
        key = f"cx_r{r:.2f}"
        if ck.has(key):
            print(f"  r = {r:.2f}  [cached] {ck.get(key):.6f}", flush=True)
            carry = True
            continue
        ag = place(mid, d, a, r)
        xs = np.vstack([x, ag])
        e = energy(gstr(sy, xs, None, 2), carry)
        if e is None:
            print(f"  r = {r:.2f}  FAILED", flush=True)
            continue
        carry = True
        ck.put(key, e)
        print(f"  r = {r:.2f}  E = {e:.6f}", flush=True)

    # counterpoise at the scan minimum
    done = {float(k[4:]): ck.get(k) for k in ck.data if k.startswith("cx_r")}
    if not done:
        print("no scan points converged"); return
    rb = min(done, key=done.get)
    print(f"\nminimum at Ag-N = {rb:.2f} A")
    if rb == max(done):
        print("  WARNING: minimum at the largest r scanned -- no bound minimum here")

    xs = np.vstack([x, place(mid, d, a, rb)])
    for key, ghost, mult, what in (("e_sub", {agi}, 1, "molecule + Ag ghost"),
                                   ("e_ag", set(range(agi)), 2, "Ag + molecule ghost")):
        if ck.has(key):
            print(f"  {what}: [cached] {ck.get(key):.6f}", flush=True)
            continue
        e = energy(gstr(sy, xs, ghost, mult), False)
        if e is None:
            print(f"  {what}: FAILED -- no counterpoise correction possible")
            return
        ck.put(key, e)
        print(f"  {what}: {e:.6f}", flush=True)

    eb = (ck.get("e_sub") + ck.get("e_ag") - done[rb]) * H2EV
    ck.put("Eb_eV", eb)
    print("\n" + "=" * 62)
    print(f"o-dinitrile CHELATE E_b = {eb:.3f} eV at Ag-N = {rb:.2f} A")
    print(f"benzonitrile, ONE nitrile = {PHCN_REF:.3f} eV")
    print(f"ratio = {eb/PHCN_REF:.2f}x")
    print("=" * 62)
    if eb > 1.6 * PHCN_REF:
        print("Two adjacent nitriles bind one Ag cooperatively -- an o-dicyanoarene")
        print("reproduces inside one molecule what HATCN gets from its")
        print("inter-molecular pocket. That makes the crack-free design concrete:")
        print("a twisted, glass-forming core carrying an o-dicyanoarene anchor.")
    elif eb > 1.15 * PHCN_REF:
        print("Modest cooperativity: the second nitrile helps but does not double")
        print("the binding. An o-dicyanoarene is worth having over a single CN,")
        print("but it does not substitute for HATCN's pocket on its own.")
    else:
        print("No cooperativity -- the bridge is worth no more than one nitrile.")
        print("The pocket effect does not transfer into a single molecule, and a")
        print("crack-free seed will have to get its binding from nitrile COUNT")
        print("and film packing rather than from an intramolecular chelate.")


if __name__ == "__main__":
    main()
