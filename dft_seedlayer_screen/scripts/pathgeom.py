"""Path construction shared by every E_d driver.

The geometry half of the barrier calculation -- where the adatom starts, where
it is dragged to, and how it is kept off the substrate -- is identical whatever
program computes the energies. It lives here so the psi4 driver (script 97) and
the Gaussian driver (script 102) cannot drift apart, and so a machine with only
one of the two can still import it.

Nothing in this file imports a quantum chemistry package.
"""
import os

import numpy as np

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STRUCT = os.path.join(BASE, "structures")
RUNS = os.path.join(BASE, "runs")
H2EV = 27.211386

SPACING = 1.0                   # A between lateral points
NPATH_MIN, NPATH_MAX = 3, 5
ZSCAN = (0.0, 0.4)           # measured up from the D_MIN contact height
D_MIN = 2.2                  # A, closest approach allowed to any substrate atom

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


BAD_GEOM = {"p-bPPhenB", "B3PyMPM", "Liq"}
"""Structures whose Ag sits closer than a bond length to the substrate.

These files are initial placements that were never replaced by the relaxed
geometry -- Ag at 1.39 A from nitrogen on p-bPPhenB, 1.85 A from hydrogen on
B3PyMPM, 0.91 A from lithium on Liq. A barrier computed from one of them is
meaningless: the path starts inside a repulsive wall, so the "site" energy is
enormous and every step downhill from it. Drivers refuse them rather than
returning a number.
"""


def sanity(syms, xyz, tag=""):
    """Reject a complex whose Ag is not at a plausible bonding distance."""
    import numpy as _np
    i = [k for k, s in enumerate(syms) if s == "Ag"]
    if len(i) != 1:
        return f"expected one Ag, found {len(i)}"
    i = i[0]
    d = _np.linalg.norm(_np.delete(xyz, i, 0) - xyz[i], axis=1)
    if d.min() < 2.0:
        near = [s for k, s in enumerate(syms) if k != i][int(_np.argmin(d))]
        return (f"Ag is {d.min():.2f} A from {near} -- an unrelaxed placement, "
                f"not a minimum")
    return None


def read_xyz(p):
    L = open(p).read().strip().splitlines()
    n = int(L[0])
    s, x = [], []
    for l in L[2:2 + n]:
        q = l.split()
        s.append(q[0])
        x.append(list(map(float, q[1:4])))
    return s, np.array(x)


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


def place(sub_x, pos, nrm):
    """Lift `pos` along `nrm` until no substrate atom is closer than D_MIN.

    Interpolating the adatom's absolute position between two binding sites sends
    it straight through the molecule. Sliding it back out along the normal keeps
    the path on the surface, which is what a diffusing adatom actually does.
    """
    p = np.array(pos, float)
    for _ in range(60):
        dmin = float(np.min(np.linalg.norm(sub_x - p, axis=1)))
        if dmin >= D_MIN:
            return p
        p = p + 0.1 * nrm
    return p
