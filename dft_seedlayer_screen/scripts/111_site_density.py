"""Accessible site density in a real film, and a correction to how it was counted.

Script 110 counted binding sites by seating the molecule against a flat plane
and asking which sites reach it. That answers the wrong question twice over.

The sites do not point at the substrate; they point outward. On HATCN the
silver binds END-ON to a nitrile at 1.029 eV and is REPELLED by the molecular
face at -0.024, so a molecule lying flat presents its nitriles edgewise, not
face up. And a seed layer is an amorphous film, not one molecule on a plane --
molecules sit in every orientation, and what the arriving silver meets is a cut
through that film.

So the right quantity is orientation-free: how many binding sites there are per
unit area of an arbitrary cut through the bulk. For molecules of volume V packed
at random, a cut crosses one molecule per V^(2/3) of area, giving

    n_s = sites per molecule / V^(2/3)

with V the van der Waals volume, integrated here by Monte Carlo over the union
of atomic spheres rather than assumed from a density.

This changes the numbers from script 110 but not its conclusion, which is the
useful thing: TPBi still lands far below HATCN, for the same reason -- it has
half the sites in a larger molecule.
"""
import os, sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathgeom import STRUCT, read_xyz

# Bondi van der Waals radii, A
RVDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "S": 1.80,
        "P": 1.80, "Li": 1.82, "Al": 1.84, "Mo": 2.10, "Cs": 3.43, "Cu": 1.40,
        "I": 1.98, "Ag": 1.72, "Zn": 1.39}
NSAMP = 400000
KT, CHI = 0.025852, 0.284
D_HATCN = 5.0
ED_FROM_EB = 0.28

# name, file, sites/molecule, E_b per site, E_d (None = estimate), mode
SEEDS = [
    ("HATCN",     "HATCN_Ag_CN.xyz",         6, 1.029, 0.286, "mono"),
    ("NDP-9",     "NDP9.xyz",                6, None,  None,  "mono"),
    ("F4TCNQ",    "F4TCNQ_Ag.xyz",           4, 0.966, None,  "mono"),
    ("TPBi",      "TPBi_Ag_N.xyz",           3, 0.889, None,  "bi"),
    ("p-bPPhenB", "pbPPhenB_Ag_chelate.xyz", 4, 0.870, None,  "bi"),
    ("B3PyMPM",   "B3PyMPM_Ag.xyz",          6, 0.626, None,  "mono"),
    ("Bphen",     "Bphen_Ag.xyz",            2, 0.490, None,  "bi"),
]


def vdw_volume(fn, seed=0):
    s, x = read_xyz(os.path.join(STRUCT, fn))
    keep = [i for i, q in enumerate(s) if q != "Ag"]
    x = x[keep]
    r = np.array([RVDW.get(s[i], 1.70) for i in keep])
    lo, hi = (x - r[:, None]).min(0), (x + r[:, None]).max(0)
    rng = np.random.default_rng(seed)
    p = rng.uniform(lo, hi, size=(NSAMP, 3))
    inside = np.zeros(NSAMP, bool)
    for c, rr in zip(x, r):
        inside |= np.sum((p - c)**2, axis=1) <= rr * rr
    return float(inside.mean() * np.prod(hi - lo))


def main():
    rows = []
    for name, fn, ns, eb, ed, mode in SEEDS:
        if not os.path.exists(os.path.join(STRUCT, fn)):
            continue
        v = vdw_volume(fn)
        area = (v**(2.0 / 3.0)) / 100.0                      # A^2 -> nm^2
        sites = ns / 2.0 if mode == "bi" else float(ns)
        rows.append(dict(name=name, eb=eb, ed=ed, mode=mode, ns=ns,
                         sites=sites, v=v, area=area, dens=sites / area))

    print(f"{'seed':<11} {'sites':>6} {'mode':<5} {'V vdW':>9} {'cut area':>9} "
          f"{'n_s':>7}")
    print(f"{'':<11} {'/mol':>6} {'':<5} {'A^3':>9} {'nm2':>9} {'/nm2':>7}")
    print("-" * 56)
    for r in rows:
        print(f"{r['name']:<11} {r['ns']:>6} {r['mode']:<5} {r['v']:>9.0f} "
              f"{r['area']:>9.3f} {r['dens']:>7.1f}")

    print("\nwith the binding energy folded in (Venables, calibrated on HATCN):\n")
    print(f"{'seed':<11} {'E_b':>7} {'E_d':>8} {'n_s':>7} {'d_c pred':>9}")
    print("-" * 46)
    ref = None
    for r in rows:
        if r["eb"] is None:
            continue
        r["edu"] = r["ed"] if r["ed"] is not None else ED_FROM_EB * r["eb"]
        r["nsat"] = r["dens"] * np.exp(CHI * r["edu"] / KT)
        if r["name"] == "HATCN":
            ref = r
    for r in sorted((r for r in rows if r["eb"] is not None),
                    key=lambda r: -r["nsat"]):
        dc = D_HATCN * np.sqrt(ref["nsat"] / r["nsat"])
        star = "" if r["ed"] is not None else "*"
        print(f"{r['name']:<11} {r['eb']:>7.3f} {r['edu']:>7.3f}{star} "
              f"{r['dens']:>7.1f} {dc:>8.1f} nm")
    print("  * E_d estimated from E_b, not computed")

    n9 = [r for r in rows if r["name"] == "NDP-9"][0]
    h = ref
    print(f"\nNDP-9 carries the same six nitriles as HATCN in a molecule "
          f"{n9['v']/h['v']:.2f} times")
    print(f"the volume, so its site density is {n9['dens']:.1f} against "
          f"{h['dens']:.1f} per nm2 -- a factor")
    print(f"{h['dens']/n9['dens']:.2f} behind before binding energy is considered at "
          f"all. To come out")
    print(f"level it would need an E_d larger by "
          f"{KT/CHI*np.log(h['dens']/n9['dens'])*1000:.0f} meV, i.e. "
          f"{h['edu'] + KT/CHI*np.log(h['dens']/n9['dens']):.3f} eV against HATCN's "
          f"{h['edu']:.3f}.")
    print("Its E_b is being computed; that is what decides it.")


if __name__ == "__main__":
    main()
