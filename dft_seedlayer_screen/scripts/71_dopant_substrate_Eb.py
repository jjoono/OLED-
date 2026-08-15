"""Dopant-substrate binding: do co-deposited Mg / Yb atoms bind the organic
(nitrile) substrate much more strongly than Ag does?  Mechanism #2 of the
alloy-smoothing story: if yes, dopant atoms act as in-situ nucleation seeds.

Systems (same level as the whole screen: PBE-D3BJ/def2-SVP, CP-corrected,
rigid distance scan along the anchor axis, UKS where open shell):
  - Mg2 dimer              (comparison row vs Ag-Ag 1.86 / Ag-Mg 0.96 eV)
  - Mg  on PhCN nitrile N  (PhCN = cheap single-nitrile HATCN proxy; Ag: 0.29 eV)
  - Ag-Yb dimer            (scan 2.6-3.8 A)
  - Yb  on PhCN nitrile N  (Yb: def2-SVP + auto def2-ECP; [Xe]4f14 6s2 closed
                            shell -> RKS first, UKS 1 then 3 as fallback)

Checkpointed per SCF point via _ckpt.Checkpoint -> runs/dopant_substrate_Eb.json.
Ordered cheapest/most-useful-first so a timeout leaves the best subset.
"""
import os, sys, json, glob, shutil
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _ckpt import Checkpoint

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")
H2EV = 27.211386
M = "pbe-d3bj"

import psi4
psi4.set_memory("20 GB")
psi4.set_num_threads(12)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_dopant_substrate.out"), False)
base_opts = {"basis": "def2-svp", "scf_type": "df", "maxiter": 80, "guess": "sad"}


def clear_scratch():
    psi4.core.clean()
    for p in glob.glob("/tmp/psi.*"):
        try:
            os.remove(p) if os.path.isfile(p) else shutil.rmtree(p, ignore_errors=True)
        except OSError:
            pass


def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0]); syms, xyz = [], []
    for l in lines[2:2 + n]:
        p = l.split(); syms.append(p[0]); xyz.append(list(map(float, p[1:4])))
    return syms, np.array(xyz)


def gstr(syms, xyz, ghost=None, mult=1):
    s = f"0 {mult}\n"
    for i, (sym, c) in enumerate(zip(syms, xyz)):
        t = f"Gh({sym})" if ghost and i in ghost else sym
        s += f"{t} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n"
    return s + "symmetry c1\nno_reorient\nno_com\n"


def robust_energy(geom_str, mult):
    """PBE-D3BJ energy with escalating SCF rescue. ~4 attempts, then None.
    Closed shell: RKS first (project rule for Yb), then UKS."""
    attempts = []
    if mult == 1:
        attempts.append(("rks", {}))
    else:
        attempts.append(("uks", {}))
    attempts += [("uks", {"level_shift": 2.0, "level_shift_cutoff": 0.005,
                          "damping_percentage": 30}),
                 ("uks", {"soscf": True, "soscf_max_iter": 30,
                          "damping_percentage": 20})]
    for i, (ref, extra) in enumerate(attempts):
        try:
            psi4.set_options({**base_opts, "reference": ref, **extra})
            e = psi4.energy(M, molecule=psi4.geometry(geom_str))
            psi4.set_options({**base_opts, "reference": "uks"})
            return e
        except Exception as ex:
            print(f"    scf attempt {i} ({ref},{list(extra) or 'plain'}) failed: "
                  f"{type(ex).__name__}", flush=True)
            psi4.core.clean()
    return None


# ---------------- geometry builders ----------------
def phcn_complex(metal, r):
    """Metal along the C#N axis, r Angstrom beyond N. PhCN xtb geometry:
    atom0 = N(nitrile), atom1 = C(nitrile)."""
    syms, x = read_xyz(os.path.join(RUNS, "PhCN", "xtbopt.xyz"))
    assert syms[0] == "N" and syms[1] == "C"
    ax = x[0] - x[1]; ax /= np.linalg.norm(ax)
    return syms + [metal], np.vstack([x, x[0] + r * ax])


def dimer(a, b, r):
    return [a, b], np.array([[0.0, 0.0, 0.0], [0.0, 0.0, r]])


# ---------------- system table (cheapest-first) ----------------
# name: (builder(r)->syms,xyz, complex_mult(s), fragA_indices, multA, multB, r-list)
# For PhCN systems fragA = substrate (all but last atom), fragB = metal atom.
SYSTEMS = [
    ("Mg2",     lambda r: dimer("Mg", "Mg", r), [1], 1, 1,
     [3.0, 3.2, 3.4, 3.6, 3.9, 4.2]),
    ("Mg_PhCN", lambda r: phcn_complex("Mg", r), [1], 1, 1,
     [2.0, 2.2, 2.4, 2.7, 3.0, 3.4, 3.8]),
    ("AgYb",    lambda r: dimer("Ag", "Yb", r), [2], 2, 1,
     [2.6, 2.8, 3.0, 3.2, 3.5, 3.8]),
    ("Yb_PhCN", lambda r: phcn_complex("Yb", r), [1, 3], 1, 1,
     [2.3, 2.5, 2.7, 3.0, 3.3, 3.7]),
]
# multA = fragment-A multiplicity (substrate or first atom); multB = metal atom.

ck = Checkpoint(os.path.join(RUNS, "dopant_substrate_Eb.json"),
                label="dopant-substrate Eb")


def run_system(name, build, cmults, multA, multB, rs):
    if ck.has(name):
        print(f"[skip] {name} done: {ck.get(name)}", flush=True)
        return
    scan = {}
    failed_pts = 0
    for r in rs:
        key = f"{name}|r{r:.2f}"
        if ck.has(key):
            scan[r] = ck.get(key); continue
        syms, xyz = build(r)
        entry = None
        for cm in cmults:
            e = robust_energy(gstr(syms, xyz, None, cm), cm)
            if e is not None and (entry is None or e < entry["E_h"]):
                entry = {"E_h": e, "mult": cm}
            print(f"  [{name}] r={r:.2f} mult={cm} "
                  f"E={'FAILED' if e is None else f'{e:.6f}'}", flush=True)
        if entry is None:
            failed_pts += 1
            ck.put(key, {"error": "scf_unconverged"})
            continue
        ck.put(key, entry)
        scan[r] = entry
    ok = {r: v for r, v in scan.items() if isinstance(v, dict) and "E_h" in v}
    if not ok:
        ck.put(name, {"error": "all_scan_points_unconverged"})
        print(f"[RESULT] {name}: FAILED (no converged scan point)", flush=True)
        return
    rbest = min(ok, key=lambda r: ok[r]["E_h"])
    e_cx, cmult = ok[rbest]["E_h"], ok[rbest]["mult"]
    edge = (rbest == min(rs)) or (rbest == max(rs))
    # CP fragments at the scan minimum geometry
    syms, xyz = build(rbest)
    nA = len(syms) - 1  # fragment A = all but last atom (metal B)
    kA, kB = f"{name}|cpA", f"{name}|cpB"
    if ck.has(kA):
        e_A = ck.get(kA)["E_h"]
    else:
        e_A = robust_energy(gstr(syms, xyz, {nA}, multA), multA)
        if e_A is not None: ck.put(kA, {"E_h": e_A})
    if ck.has(kB):
        e_B = ck.get(kB)["E_h"]
    else:
        e_B = robust_energy(gstr(syms, xyz, set(range(nA)), multB), multB)
        if e_B is not None: ck.put(kB, {"E_h": e_B})
    if e_A is None or e_B is None:
        ck.put(name, {"error": "cp_fragment_unconverged", "r_A": rbest,
                      "scan_eV_shape": {f"{r:.2f}": v["E_h"] for r, v in ok.items()}})
        print(f"[RESULT] {name}: CP fragment failed", flush=True)
        return
    eb = (e_A + e_B - e_cx) * H2EV
    res = {"Eb_eV": round(eb, 4), "r_A": rbest, "complex_mult": cmult,
           "min_at_scan_edge": edge, "n_failed_points": failed_pts,
           "scan_E_h": {f"{r:.2f}": ok[r]["E_h"] for r in sorted(ok)},
           "level": "PBE-D3BJ/def2-SVP(+def2-ECP for Yb), CP-corrected rigid scan"}
    ck.put(name, res)
    print(f"[RESULT] {name}: E_b(CP) = {eb:.3f} eV at r = {rbest:.2f} A "
          f"(mult {cmult}{', EDGE MINIMUM' if edge else ''})", flush=True)


if __name__ == "__main__":
    for name, build, cmults, multA, multB, rs in SYSTEMS:
        run_system(name, build, cmults, multA, multB, rs)
        clear_scratch()
    print("ALL DONE", flush=True)
    print(json.dumps({k: v for k, v in ck.data.items() if "|" not in k}, indent=2))
