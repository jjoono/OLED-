"""A composite descriptor for a seed layer, and what it must not be asked to do.

E_b and E_d describe one adsorption site. Whether a film closes at 5 nm or 8
depends on how many silver islands nucleate per unit area, and that is a
different quantity: a molecule that binds twice as hard but presents half as
many sites per square nanometre is a wash, and one that binds hard but sits at
an angle so only one of its three sites faces the surface is worse than the
per-site number suggests. TPBi is exactly that case.

THE CHAIN, from what is computed to what is measured.

  n_s      accessible sites per nm^2, counted the way a substrate would count
           them. The molecule is laid against a flat surface in every
           orientation on a spherical grid; for each, the surface is the tangent
           plane to the molecule's van der Waals envelope, and a site counts if
           it reaches within CONTACT of that plane. The orientation that puts
           the most sites in contact is the one the molecule will adopt, and its
           count is what is used. A propeller loses the sites that point away
           without needing a fudge factor, and a flat molecule keeps all of
           them -- but neither outcome is assumed.
  N_sat    Venables saturation island density. Adatoms land at rate R, hop with
           barrier E_d, and stick where they meet. N_sat ~ n_s exp(chi E_d/kT),
           with chi = 0.284 as fitted against the kinetic Monte Carlo runs.
  d_c      islands touch when their spacing is N_sat^-1/2, and the thickness
           needed to reach that scales the same way, so d_c ~ N_sat^-1/2.

One free constant, fixed on HATCN, where the closure thickness is measured at
5 nm by transport and by optics independently. Everything else is a prediction.

WHAT THIS NUMBER CANNOT DECIDE. It ranks nucleation, nothing else. Four hard
filters sit outside it and have vetoed candidates before:

  electronic role   a seed under a top anode is in the injection path. HATCN is
                    a hole-injection material, so it is free there; TPBi is an
                    electron transporter and would block it.
  refractive index  the seed is also an optical layer. p-bPPhenB was chosen over
                    better binders in the literature because its n = 2.21
                    matches MoO3.
  film formation    HATCN is small, rigid and flat -- the worst case for glass
                    formation, and it crystallises and cracks above ~30 nm. Fine
                    at 5 nm, not a general-purpose layer.
  deposition temp   MoO3 sublimes near 556 C against 345 C for p-bPPhenB, which
                    is what rules it out on plastic.

A candidate has to pass all four before its nucleation rank means anything.
"""
import os, sys

import numpy as np
from scipy.spatial import ConvexHull

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathgeom import STRUCT, read_xyz

KT = 0.025852          # eV at 300 K
CHI = 0.284            # Venables exponent, fitted against the kMC runs
CONTACT = 1.2          # A; how close a site must come to the substrate plane
N_ORIENT = 4000        # orientations tried when seating the molecule
VDW = 1.7              # A, carbon van der Waals radius
D_HATCN = 5.0          # nm, measured closure thickness -- the single calibration
ED_FROM_EB = 0.28      # monodentate correlation, used only where E_d is unknown

# name, file, E_b per site, E_d if computed, mode, site element, site kind
#
# "kind" picks out which atoms of that element actually bind, by counting heavy
# neighbours. A nitrile nitrogen has one; an aromatic ring nitrogen has two.
# HATCN carries both -- six nitriles that bind at 1.029 eV and six core aza
# nitrogens on a face that is REPULSIVE at -0.024 eV -- so counting every
# nitrogen would double its site density on the strength of sites that push
# silver away.
SEEDS = [
    ("HATCN",     "HATCN_Ag_CN.xyz",          1.029, 0.286, "mono", "N", "nitrile"),
    ("F4TCNQ",    "F4TCNQ_Ag.xyz",            0.966, None,  "mono", "N", "nitrile"),
    ("TPBi",      "TPBi_Ag_N.xyz",            0.889, None,  "bi",   "N", "ring"),
    ("p-bPPhenB", "pbPPhenB_Ag_chelate.xyz",  0.870, None,  "bi",   "N", "ring"),
    ("B3PyMPM",   "B3PyMPM_Ag.xyz",           0.626, None,  "mono", "N", "ring"),
    ("Bphen",     "Bphen_Ag.xyz",             0.490, None,  "bi",   "N", "ring"),
]


def sphere(n):
    """n roughly uniform directions, Fibonacci lattice."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5**0.5) * i
    return np.c_[np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi),
                 np.cos(phi)]


def geometry_terms(fn, elem, kind):
    """Sites reaching the substrate in the molecule's best orientation, and area."""
    s, x = read_xyz(os.path.join(STRUCT, fn))
    heavy_i = [i for i, q in enumerate(s) if q not in ("H", "Ag")]
    heavy = x[heavy_i]
    c0 = heavy.mean(0)

    # which atoms of `elem` are binding sites, by heavy-neighbour count
    want = 1 if kind == "nitrile" else 2
    sites = []
    for i, q in enumerate(s):
        if q != elem:
            continue
        d = np.linalg.norm(x - x[i], axis=1)
        d[i] = 9
        nb = sum(1 for j, dd in enumerate(d) if dd < 1.75 and s[j] not in ("H", "Ag"))
        if nb == want:
            sites.append(i)
    if not sites:
        return 0, 0, 0.0, 0.0

    # seat the molecule: for each orientation the substrate is the tangent plane
    # to the heavy-atom envelope, and a site counts if it reaches CONTACT of it
    best = (0, None)
    for n in sphere(N_ORIENT):
        floor = np.min((heavy - c0) @ n)
        reach = (x[sites] - c0) @ n
        k = int(np.sum(reach - floor <= CONTACT))
        if k > best[0]:
            best = (k, n)
    facing, normal = best
    if normal is None:
        normal = np.array([0.0, 0.0, 1.0])

    # footprint in the plane of that same orientation
    e1 = np.cross(normal, [1.0, 0.0, 0.0])
    if np.linalg.norm(e1) < 1e-6:
        e1 = np.cross(normal, [0.0, 1.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(normal, e1)
    proj = np.c_[(heavy - c0) @ e1, (heavy - c0) @ e2]
    h = ConvexHull(proj)
    per = sum(np.linalg.norm(proj[a] - proj[b])
              for a, b in zip(h.vertices, np.roll(h.vertices, -1)))
    area = (h.volume + per * VDW + np.pi * VDW**2) / 100.0

    cc = heavy - c0
    _, _, vt = np.linalg.svd(cc, full_matrices=False)
    rms = float(np.sqrt(((cc @ vt[-1])**2).mean()))
    return len(sites), facing, area, rms


def main():
    rows = []
    for name, fn, eb, ed, mode, elem, kind in SEEDS:
        free, facing, area, rms = geometry_terms(fn, elem, kind)
        # a bidentate pocket is ONE site that uses two atoms, not two sites
        sites = facing / 2.0 if mode == "bi" else float(facing)
        n_s = sites / area if area > 0 else 0.0
        ed_use = ed if ed is not None else ED_FROM_EB * eb
        rows.append(dict(name=name, eb=eb, ed=ed_use, known=ed is not None,
                         mode=mode, free=free, facing=facing, sites=sites,
                         area=area, rms=rms, n_s=n_s,
                         nsat=n_s * np.exp(CHI * ed_use / KT)))

    ref = [r for r in rows if r["name"] == "HATCN"][0]
    for r in rows:
        r["dc"] = (D_HATCN * np.sqrt(ref["nsat"] / r["nsat"])
                   if r["nsat"] > 0 else float("inf"))

    print(f"{'seed':<11} {'E_b':>6} {'E_d':>7} {'mode':<5} {'sites':>7} "
          f"{'contact':>8} {'area':>7} {'n_s':>7} {'d_c pred':>9}")
    print(f"{'':<11} {'eV':>6} {'eV':>7} {'':<5} {'total':>7} {'in reach':>8} "
          f"{'nm2':>7} {'/nm2':>7} {'nm':>9}")
    print("-" * 78)
    for r in sorted(rows, key=lambda r: r["dc"]):
        ed = f"{r['ed']:.3f}" + ("" if r["known"] else "*")
        print(f"{r['name']:<11} {r['eb']:>6.3f} {ed:>7} {r['mode']:<5} "
              f"{r['free']:>8} {r['facing']:>7} {r['area']:>7.2f} "
              f"{r['n_s']:>7.2f} {r['dc']:>9.1f}")
    print("  * E_d estimated from E_b by the monodentate correlation, not computed")

    print("\nWhat the composite changes, against ranking on E_b alone:")
    by_eb = [r["name"] for r in sorted(rows, key=lambda r: -r["eb"])]
    by_dc = [r["name"] for r in sorted(rows, key=lambda r: r["dc"])]
    print(f"  by E_b : {' > '.join(by_eb)}")
    print(f"  by d_c : {' > '.join(by_dc)}")

    t = [r for r in rows if r["name"] == "TPBi"][0]
    h = ref
    print(f"\nTPBi is the case that makes the point. Its E_b is {t['eb']/h['eb']*100:.0f} % "
          f"of HATCN's,")
    print(f"but it is a propeller -- heavy atoms {t['rms']:.2f} A rms off their own")
    print(f"plane against {h['rms']:.2f} for HATCN -- so only {t['facing']} of its "
          f"{t['free']} nitrogens face")
    print(f"the surface, and being bidentate those {t['facing']} make {t['sites']:.0f} site. "
          f"Per nm2 it offers")
    print(f"{t['n_s']:.2f} against {h['n_s']:.2f}, and the predicted closure thickness is "
          f"{t['dc']:.1f} nm against {h['dc']:.1f}.")

    print("\nUSE IT AS A RANKING, NOT A NUMBER. One constant is fitted, on one")
    print("measured film. The exponent chi comes from kMC. Half the E_d column is")
    print("an estimate. What it is good for is ordering candidates and saying")
    print("which measurement would discriminate them -- and the discriminating")
    print("measurement is closure thickness, which is cheap: evaporate, probe,")
    print("look for the knee.")


if __name__ == "__main__":
    main()
