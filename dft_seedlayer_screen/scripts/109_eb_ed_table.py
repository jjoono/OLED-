"""The screening as it actually stands: E_b and E_d, side by side.

Assembled from runs/*.json rather than transcribed, and the binding mode is read
off the optimised geometry rather than off the name of the file -- Ag is counted
as coordinated to every substrate atom within COORD_MAX of it. That is the check
that caught TPBi being reported as a single-nitrogen contact when it relaxes
onto two, and it is applied here to every candidate uniformly.

Two axes, and they are not equally complete:

  E_b   counterpoise-corrected PBE-D3BJ/def2-SVP, twenty-odd candidates
  E_d   the same level, four candidates

E_d is the gap. The campaign is written (scripts 97 and 102) and running on a
workstation with Gaussian; a GFN2-xTB screen was tried and rejected, since it
returns 0.000 eV on both systems where a DFT barrier exists.

READ ACROSS A ROW, NOT DOWN A COLUMN. A monodentate E_b and a bidentate one are
not the same quantity, and neither are a site-to-site E_d and a to-face one; the
class is printed with every number for that reason. Entries with no real anchor
are listed separately at the end rather than ranked.
"""
import json, os, sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathgeom import RUNS, STRUCT, read_xyz

COORD_MAX = 2.75          # A; beyond this Ag is not bonded, it is nearby

# candidate -> (E_b eV, structure file with the optimised complex, note)
EB = [
    ("HATCN (CN)",  1.029, "HATCN_Ag_CN.xyz",         ""),
    ("F4TCNQ",      0.966, "F4TCNQ_Ag.xyz",           ""),
    ("Cs2CO3",      0.900, "Cs2CO3_Ag.xyz",           ""),
    ("TPBi",        0.889, "TPBi_Ag_N.xyz",           ""),
    ("p-bPPhenB",   0.870, "pbPPhenB_Ag_chelate.xyz", ""),
    ("B3PyMPM",     0.626, "B3PyMPM_Ag.xyz",          ""),
    ("Bphen",       0.490, "Bphen_Ag.xyz",            ""),
    ("Mo3O8",       0.445, "Mo3O8_Ag.xyz",            "oxygen-deficient MoOx"),
    ("Al4O6",       0.422, "Al4O6_Ag.xyz",            "stoichiometric Al2O3"),
    ("pyridine",    0.342, "pyridine_Ag.xyz",         "reference fragment"),
    ("triazine",    0.291, "triazine_Ag_refined.xyz", ""),
    ("Mo3O9",       0.284, "Mo3O9_Ag.xyz",            "stoichiometric MoO3"),
    ("BTD",         0.277, "BTD_Ag_refined.xyz",      ""),
    ("Me3P=O",      0.253, "Me3PO_Ag_refined.xyz",    ""),
    ("LiF32",       0.249, "LiF32_Ag.xyz",            ""),
    ("thiophene",   0.169, "thiophene_Ag_refined.xyz", ""),
    ("Liq",         0.167, "Liq_Ag.xyz",              ""),
    ("benzene",     0.167, "benzene_Ag.xyz",          "dispersion floor"),
    ("HATCN (face)", -0.024, "HATCN_Ag_face.xyz",     "same molecule, wrong site"),
]

# runs/diffusion_barrier_eV.json, keyed differently from the E_b table
ED_KEY = {"HATCN (CN)": "HATCN", "Al4O6": "AlOx"}

# ranked separately: the relaxed geometry shows no coordinating atom at all
NO_ANCHOR = [
    ("PhCz",  0.279, "Ag rests on the pi face, nearest C 2.74 A"),
    ("TPA",   0.254, "Ag sits over a hydrogen, 2.89 A"),
    ("Ag-Ag", 1.858, "cohesion of the metal itself, the number to beat"),
]

RETRACTED = [
    ("Al2O3 0.910", "an Al-rich Al10O10 cluster; stoichiometric Al4O6 gives 0.422"),
    ("MoO3 -0.718", "unrefined placement; the refined value is +0.284"),
    ("ZnS 1.60",    "a literature prior, never computed here"),
    ("PhCN/ZnS 1.218", "ZnS-to-nitrile binding, nothing to do with Ag"),
    ("CuI, CuSCN",  "Ag lands on the copper -- a metallic contact, not a seed site"),
]


def mode_of(fn):
    """Coordination read from the relaxed geometry: what is within COORD_MAX."""
    p = os.path.join(STRUCT, fn)
    if not os.path.exists(p):
        return "?", "", 0
    syms, xyz = read_xyz(p)
    i = [k for k, s in enumerate(syms) if s == "Ag"]
    if len(i) != 1:
        return "?", "", 0
    i = i[0]
    d = np.linalg.norm(np.delete(xyz, i, 0) - xyz[i], axis=1)
    s = [x for k, x in enumerate(syms) if k != i]
    near = sorted([(dd, ss) for dd, ss in zip(d, s) if dd <= COORD_MAX])
    if not near:
        return "none", f"nearest {min(d):.2f} A", 0
    n = len(near)
    name = {1: "mono", 2: "bi"}.get(n, f"{n}-fold")
    return name, " / ".join(f"{ss} {dd:.2f}" for dd, ss in near[:3]), n


def main():
    ed = {}
    p = os.path.join(RUNS, "diffusion_barrier_eV.json")
    if os.path.exists(p):
        ed = json.load(open(p))
    for extra in ("diffusion_barriers_all.json", "diffusion_barriers_gaussian.json",
                  "Ed_merged.json"):
        q = os.path.join(RUNS, extra)
        if os.path.exists(q):
            for k, v in json.load(open(q)).items():
                if isinstance(v, dict) and v.get("E_d_eV") is not None:
                    ed[k] = v["E_d_eV"]

    print(f"{'candidate':<14} {'E_b':>7} {'mode':<7} {'E_d':>7}  contacts within "
          f"{COORD_MAX} A")
    print("-" * 78)
    for name, eb, fn, note in EB:
        m, contacts, _ = mode_of(fn)
        key = ED_KEY.get(name, name.split()[0])
        v = ed.get(key)
        vs = f"{v:7.3f}" if v is not None else "      -"
        tail = contacts + (f"   ({note})" if note else "")
        print(f"{name:<14} {eb:>7.3f} {m:<7} {vs}  {tail}")

    print(f"\n{'':<14} {'E_b':>7}  no coordinating atom in the relaxed geometry")
    print("-" * 78)
    for name, eb, why in NO_ANCHOR:
        print(f"{name:<14} {eb:>7.3f}  {why}")

    have = sum(1 for n, _, _, _ in EB if ed.get(ED_KEY.get(n, n.split()[0])) is not None)
    print(f"\nE_b: {len(EB)} candidates.   E_d: {have}.")
    missing = [n for n, _, _, _ in EB
               if ed.get(ED_KEY.get(n, n.split()[0])) is None]
    print(f"awaiting a barrier: {', '.join(missing)}")

    print("\nretracted, and why:")
    for what, why in RETRACTED:
        print(f"  {what:<18} {why}")

    print("\nTwo cautions on reading this table.")
    print("  A bidentate E_b is not comparable with a monodentate one -- TPBi's")
    print("  0.889 comes from two nitrogens at 2.33 and 2.34 A, while HATCN's")
    print("  1.029 comes from a single nitrile nitrogen. Per site the gap is")
    print("  larger than the column suggests, not smaller.")
    print("  And E_b is per site; what sets nucleation density is sites per unit")
    print("  area. HATCN carries six nitriles on a flat molecule, all exposed.")


if __name__ == "__main__":
    main()
