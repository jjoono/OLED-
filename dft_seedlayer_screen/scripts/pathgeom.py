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
# Three points cannot show a barrier. A path sampled at t=0, 0.5, 1 has no
# interior point, so the reported maximum is whichever endpoint is higher and the
# number is a lower bound, not a barrier. The v5 run made this concrete: Cs2CO3
# (0.023, 0.000, 0.135 eV) and Bphen (0.003, 0.093, 0.103 eV) were both still
# rising at their last point.
NPATH_MIN, NPATH_MAX = 5, 7
# Three heights, straddling the contact distance in both directions.
#
# Two heights cannot bracket a minimum: whichever is lower is an endpoint of the
# scan, so the height is never relaxed, only guessed between two arbitrary
# values. In the v8 run that guess was worth up to 0.456 eV (DMABN t3) -- more
# than any barrier in the table -- and it varied along the path, so it distorted
# the shape of the curve rather than shifting it. Which of the two won even
# flipped mid-path (BTD z0,z0,z1,z1,z1), which is what put spurious structure
# into the barrier.
#
# The contact height won 26 of 42 points, so the minimum sits near it on both
# sides; the scan has to reach below the contact distance as well as above it.
# A point slightly inside the repulsive wall is what makes the parabola fit
# possible, so the reader gets a real minimum instead of a lower of two.
ZSCAN = (-0.2, 0.2, 0.6)     # measured from the contact height, along the normal
D_MIN = 2.2                  # A, floor on the closest approach to any substrate atom

# One contact distance for every element pair is wrong by an angstrom on the
# heavy ones. In the v5 run Cu4I4 sat 3.10 A from iodine at t=0 and 2.24 A at
# t=2, because 2.2 A is the floor for Ag-I as much as for Ag-C. The resulting
# 1.5 eV "barrier" was Pauli repulsion from pushing Ag into the iodine core, and
# it collapsed to 0.29 eV half an angstrom further out -- which no chemisorption
# barrier does. The contact distance is set per pair instead, from the covalent
# radii (Cordero 2008), so Ag-I gets 2.84 A while Ag-C keeps 2.21 A.
R_COV = {"H": 0.31, "Li": 1.28, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
         "Al": 1.21, "P": 1.07, "S": 1.05, "Cl": 1.02, "Cu": 1.32, "Br": 1.20,
         "Mo": 1.54, "Ag": 1.45, "I": 1.39, "Cs": 2.44, "Ba": 2.15}
R_AG = R_COV["Ag"]


def contact(sym):
    """Closest approach allowed between Ag and one substrate element."""
    return max(D_MIN, R_AG + R_COV.get(sym, 1.0))


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


def _sphere(n=1200):
    """n roughly uniform unit vectors (Fibonacci sphere)."""
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * k / n)
    th = np.pi * (1 + 5 ** 0.5) * k
    return np.column_stack([np.cos(th) * np.sin(phi),
                            np.sin(th) * np.sin(phi), np.cos(phi)])


def outward(sub_s, sub_x, site, h, hint=None):
    """Direction from `site` that puts an adatom at height h furthest from the
    substrate.

    Every site needs its own outward direction. The v5 run reused the starting
    site's normal at the destination, which on a three-dimensional molecule
    points back into it: on Al4O6 the destination came out 2.0 A from the atom
    the path started on, so place() pushed it back out to where it began and the
    adatom orbited its starting oxygen for the whole path. That is why the
    energy fell monotonically -- there was no hop in the path at all.

    Clearance is measured against each element's own contact distance, so the
    direction chosen is the one that is genuinely open rather than the one that
    happens to point away from the heaviest neighbour.
    """
    u = _sphere()
    if hint is not None:
        h_ = np.asarray(hint, float)
        nh = np.linalg.norm(h_)
        if nh > 1e-9:
            u = u[u @ (h_ / nh) > 0.0]
    lim = np.array([contact(x) for x in sub_s])
    d = np.linalg.norm((site + h * u[:, None, :]) - sub_x[None, :, :], axis=2)
    return u[int((d - lim).min(axis=1).argmax())]


def destination(sub_s, sub_x, ag, anchor, rule):
    """Where the adatom is dragged to, and the class of that path.

    Returns the destination *site atom* index (None for a face target), the
    adatom position above it, and the path class.
    """
    heavy = np.array([x for s, x in zip(sub_s, sub_x) if s != "H"])
    cen = heavy.mean(axis=0)
    h = float(np.linalg.norm(ag - sub_x[anchor]))

    if rule == "auto":
        el = sub_s[anchor]
        same = [j for j, s in enumerate(sub_s) if s == el and j != anchor]
        if same:
            # nearest equivalent site: diffusion takes the cheapest hop, and a
            # drag across the whole molecule is a different process entirely
            near = min(same,
                       key=lambda j: np.linalg.norm(sub_x[j] - sub_x[anchor]))
            u = outward(sub_s, sub_x, sub_x[near], h, hint=sub_x[near] - cen)
            return near, sub_x[near] + h * u, "site2site"

    u = outward(sub_s, sub_x, cen, h, hint=ag - cen)
    return None, cen + h * u, ("chelate" if rule == "chelate" else "toface")


def frames(sub_s, sub_x, ag, anchor, rule, npath):
    """Adatom positions and surface normals along the hop.

    The site and the normal are interpolated separately. Interpolating the
    adatom's absolute position instead draws a chord, which on anything that is
    not flat cuts through the molecule and has to be pushed back out along one
    fixed normal -- and pushing it out along the *starting* site's normal is what
    collapsed the Al4O6 path onto its own start. Walking the site across the
    surface while the normal rotates with it keeps the adatom above the surface
    the whole way, which is the path a diffusing atom actually takes.
    """
    near, dest, cls = destination(sub_s, sub_x, ag, anchor, rule)
    h = float(np.linalg.norm(ag - sub_x[anchor]))
    a_site = sub_x[anchor]
    d_site = sub_x[near] if near is not None else dest - (dest - a_site) * 0.0
    if near is None:
        heavy = np.array([x for s, x in zip(sub_s, sub_x) if s != "H"])
        d_site = heavy.mean(axis=0)
    n_a = (ag - a_site) / np.linalg.norm(ag - a_site)
    n_d = (dest - d_site) / np.linalg.norm(dest - d_site)

    out = []
    for t in np.linspace(0.0, 1.0, npath):
        site = (1 - t) * a_site + t * d_site
        u = (1 - t) * n_a + t * n_d
        nu = np.linalg.norm(u)
        u = n_a if nu < 1e-6 else u / nu
        out.append((place(sub_x, site + h * u, u, sub_s), u))
    return out, cls, float(np.linalg.norm(dest - ag))


def place(sub_x, pos, nrm, sub_s=None):
    """Lift `pos` along `nrm` until no substrate atom is closer to Ag than its
    own contact distance.

    Interpolating the adatom's absolute position between two binding sites sends
    it straight through the molecule. Sliding it back out along the normal keeps
    the path on the surface, which is what a diffusing adatom actually does.

    `sub_s` gives the substrate element symbols, so each atom gets its own
    contact distance. Without it every atom falls back to the flat D_MIN, which
    is the behaviour that produced Cu4I4's spurious 1.5 eV barrier.
    """
    p = np.array(pos, float)
    lim = (np.array([contact(x) for x in sub_s]) if sub_s is not None
           else np.full(len(sub_x), D_MIN))
    for _ in range(80):
        if float((np.linalg.norm(sub_x - p, axis=1) - lim).min()) >= 0.0:
            return p
        p = p + 0.1 * nrm
    return p
