"""Ag adatom surface diffusion barrier E_d across the whole candidate set.

E_b has been available for twenty-odd candidates for a while; E_d existed for
four. A seed layer needs both -- strong binding alone still lets the adatom walk
until it meets another and nucleates an island -- so the screening cannot be
read without this axis.

PROTOCOL, identical for every candidate.

  start   the relaxed Ag complex, structures/<tag>.xyz
  path    Ag is dragged laterally from its bound position to a destination on
          the same molecule -- the NEAREST equivalent site, not any equivalent
          site, since diffusion proceeds by the cheapest hop available. At every
          lateral position its height is scanned along the local surface normal
          and the minimum kept, which is the relaxed-Ag, frozen-substrate
          approximation used in script 24. Points are placed about 1 A apart, so
          a long path gets more of them rather than a coarser description.
  E_d     max(E along path) - E(start)

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

SCF needs the level-shifted, damped settings -- with plain DIIS these open-shell
Ag complexes oscillate indefinitely rather than converge.
"""
import json, os, sys, time
import numpy as np
import psi4

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STRUCT = os.path.join(BASE, "structures")
RUNS = os.path.join(BASE, "runs")
OUT = os.path.join(RUNS, "diffusion_barriers_all.json")
H2EV = 27.211386

psi4.set_memory("10 GB")
psi4.set_num_threads(len(os.sched_getaffinity(0)))
psi4.core.set_output_file(os.path.join(RUNS, "psi4_diff_all.out"), False)

ROB = {"basis": "def2-svp", "scf_type": "df", "reference": "uks", "maxiter": 400,
       "guess": "sad", "level_shift": 1.0, "level_shift_cutoff": 1e-3,
       "damping_percentage": 15.0}

SPACING = 1.0                   # A between lateral points
NPATH_MIN, NPATH_MAX = 5, 11
ZSCAN = (-0.4, -0.2, 0.0, 0.25, 0.5)

# tag -> (structure file, destination rule, multiplicity)
#   "auto"  hop to the nearest other atom of the same element as the anchor
#   "face"  ring/molecular centroid of the heavy atoms
CANDIDATES = [
    ("HATCN",      "HATCN_Ag_CN.xyz",          "auto", 2),
    ("F4TCNQ",     "F4TCNQ_Ag.xyz",            "auto", 2),
    ("TPBi",       "TPBi_Ag_N.xyz",            "auto", 2),
    ("p-bPPhenB",  "pbPPhenB_Ag_chelate.xyz",  "face", 2),
    ("Cs2CO3",     "Cs2CO3_Ag.xyz",            "auto", 2),
    ("B3PyMPM",    "B3PyMPM_Ag.xyz",           "auto", 2),
    ("Bphen",      "Bphen_Ag.xyz",             "face", 2),
    ("pyridine",   "pyridine_Ag.xyz",          "face", 2),
    ("triazine",   "triazine_Ag_refined.xyz",  "auto", 2),
    ("BTD",        "BTD_Ag_refined.xyz",       "auto", 2),
    ("Me3P=O",     "Me3PO_Ag_refined.xyz",     "face", 2),
    ("Ph3P=O",     "Ph3PO_Ag.xyz",             "face", 2),
    ("thiophene",  "thiophene_Ag_refined.xyz", "face", 2),
    ("benzene",    "benzene_Ag.xyz",           "face", 2),
    ("PhCN",       "PhCN_Ag.xyz",              "face", 2),
    ("DMABN",      "DMABN_Ag.xyz",             "face", 2),
    ("pDCNB",      "pDCNB_Ag.xyz",             "auto", 2),
    ("oDCNB",      "oDCNB_Ag.xyz",             "auto", 2),
    ("PhCz",       "PhCz_Ag.xyz",              "face", 2),
    ("TPA",        "TPA_Ag.xyz",               "face", 2),
    ("Liq",        "Liq_Ag.xyz",               "face", 2),
    ("Al4O6",      "Al4O6_Ag.xyz",             "auto", 2),
    ("Mo3O9",      "Mo3O9_Ag.xyz",             "auto", 2),
    ("Mo3O8",      "Mo3O8_Ag.xyz",             "auto", 2),
    ("Cu4I4",      "Cu4I4_Ag_refined.xyz",     "auto", 2),
    ("LiF32",      "LiF32_Ag.xyz",             "auto", 2),
]


def read_xyz(p):
    L = open(p).read().strip().splitlines()
    n = int(L[0])
    s, x = [], []
    for l in L[2:2 + n]:
        q = l.split()
        s.append(q[0])
        x.append(list(map(float, q[1:4])))
    return s, np.array(x)


def gstr(syms, xyz, mult):
    st = f"0 {mult}\n"
    for a, c in zip(syms, xyz):
        st += f"{a} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n"
    return st + "symmetry c1\nno_reorient\nno_com\n"


def sp(syms, xyz, mult):
    psi4.set_options(ROB)
    e = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, mult)))
    psi4.core.clean()
    return e


def geometry(syms, xyz):
    """Split off Ag; return substrate, Ag position, anchor index, surface normal."""
    i_ag = [i for i, s in enumerate(syms) if s == "Ag"]
    if len(i_ag) != 1:
        raise ValueError(f"expected exactly one Ag, found {len(i_ag)}")
    i = i_ag[0]
    sub_s = [s for j, s in enumerate(syms) if j != i]
    sub_x = np.delete(xyz, i, axis=0)
    ag = xyz[i]
    d = np.linalg.norm(sub_x - ag, axis=1)
    anchor = int(np.argmin(d))
    n = ag - sub_x[anchor]
    n /= np.linalg.norm(n)
    return sub_s, sub_x, ag, anchor, n


def destination(sub_s, sub_x, ag, anchor, rule):
    """Where the adatom is dragged to, and the class of that path."""
    if rule == "auto":
        el = sub_s[anchor]
        same = [j for j, s in enumerate(sub_s) if s == el and j != anchor]
        if same:
            # nearest equivalent site: diffusion takes the cheapest hop, and a
            # drag across the whole molecule is a different process entirely
            near = min(same, key=lambda j: np.linalg.norm(sub_x[j] - sub_x[anchor]))
            # sit the same distance above the destination atom as above the anchor
            h = np.linalg.norm(ag - sub_x[anchor])
            nrm = ag - sub_x[anchor]
            return sub_x[near] + h * nrm / np.linalg.norm(nrm), "site2site"
    heavy = np.array([x for s, x in zip(sub_s, sub_x) if s != "H"])
    cen = heavy.mean(axis=0)
    h = np.linalg.norm(ag - sub_x[anchor])
    nrm = ag - sub_x[anchor]
    return cen + h * nrm / np.linalg.norm(nrm), ("chelate" if rule == "chelate" else "toface")


def barrier(tag, fn, rule, mult):
    syms, xyz = read_xyz(os.path.join(STRUCT, fn))
    sub_s, sub_x, ag, anchor, nrm = geometry(syms, xyz)
    dest, cls = destination(sub_s, sub_x, ag, anchor, rule)
    print(f"[{tag}] {len(sub_s)} substrate atoms, anchor {sub_s[anchor]}{anchor}, "
          f"path {cls}, span {np.linalg.norm(dest-ag):.2f} A", flush=True)

    span = float(np.linalg.norm(dest - ag))
    npath = int(np.clip(round(span / SPACING) + 1, NPATH_MIN, NPATH_MAX))
    E = []
    for t in np.linspace(0.0, 1.0, npath):
        pos = (1 - t) * ag + t * dest
        best = None
        for dz in (ZSCAN if t > 0 else (0.0,)):     # start point is already relaxed
            trial = pos + dz * nrm
            try:
                e = sp(sub_s + ["Ag"], np.vstack([sub_x, trial]), mult)
            except Exception as exc:
                print(f"    t={t:.2f} dz={dz:+.2f} SCF failed: {type(exc).__name__}",
                      flush=True)
                continue
            if best is None or e < best:
                best = e
        E.append(best)
        print(f"    t={t:.2f}  E={'FAIL' if best is None else f'{best:.6f}'}", flush=True)

    ok = [(i, e) for i, e in enumerate(E) if e is not None]
    if E[0] is None or len(ok) < 3:
        print(f"[{tag}] insufficient converged points", flush=True)
        return None
    ed = (max(e for _, e in ok) - E[0]) * H2EV
    print(f"[{tag}] E_d = {ed:.3f} eV  ({cls})", flush=True)
    return {"E_d_eV": round(ed, 4), "class": cls, "anchor": sub_s[anchor],
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
