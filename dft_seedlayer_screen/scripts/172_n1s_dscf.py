"""N 1s core-level shifts of HATCN under Ag: Delta-SCF, for the XPS test.

    QC_THREADS=4 python scripts/172_n1s_dscf.py   (resumable) -> runs/n1s_dscf.json

The charge picture predicts what the substrate side of an XPS experiment sees
(script 139 turns charges into Mo5+ / HATCN-anion fractions). What it does not
give is the N 1s spectrum itself: where the peaks of an Ag-contacted HATCN
appear relative to the pristine film, which is what an ex-situ step-deposition
XPS (no in-situ tool in the lab) can resolve. Literature anchor: HAT-CN 5 nm
N 1s 399.6 eV on ITO, 400.2 eV on Ag2O/Ag (notes/lit_survey, [18]).

Delta-SCF: BE = E(N 1s core-ionised, MOM) - E(ground state), PBE0, UKS.
The target N gets a core-valence basis (cc-pCVTZ), every other atom def2-SVP
(Ag with its ECP), which also localises the core hole on the target among the
symmetry-equivalent N. Systems, all at the package-F PBE0 geometry of
HATCN-Ag1 (structures/HATCN_Ag1_deep.xyz):
  neutral HATCN (Ag removed)     aza N, nitrile N
  HATCN-Ag1 (doublet)            aza N bonded to Ag, aza N far from Ag, nitrile N
  HATCN anion (Ag removed)       aza N, nitrile N
Absolute BEs carry ~0.5 eV method error; the shifts between them are the
prediction (~0.1-0.2 eV).
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "runs", "n1s_dscf.json")
EV = 27.211386


def load():
    return json.load(open(OUT)) if os.path.exists(OUT) else {}


def put(k, v):
    d = load()
    d[k] = v
    json.dump(d, open(OUT, "w"), indent=1)


def build(S, X, charge, spin, target):
    from pyscf import gto
    m = gto.Mole()
    atoms, basis = [], {}
    for k, (s, x) in enumerate(zip(S, X)):
        lab = f"{s}{k}" if k == target else s
        atoms.append((lab, tuple(x)))
        basis[lab] = "cc-pcvtz" if k == target else "def2-svp"
    m.atom = atoms
    m.unit = "Angstrom"
    m.basis = basis
    m.ecp = {"Ag": "def2-svp"}
    m.charge = charge
    m.spin = spin
    m.verbose = 0
    m.max_memory = int(os.environ.get("QC_MEM_MB", "9000"))
    m.build()
    return m


def ground(m):
    from pyscf import dft
    mf = dft.UKS(m).density_fit()
    mf.xc = "PBE0"
    mf.grids.level = 4
    mf.conv_tol = 1e-9
    mf.max_cycle = 60
    mf.level_shift = 0.2
    mf.kernel()
    if not mf.converged:
        mf = mf.newton()
        mf.kernel(mf.make_rdm1())
    return mf


def core_hole(mf0, target):
    """Remove the alpha electron from the occupied orbital with the largest
    weight on the target N's 1s functions, then MOM-relax."""
    from pyscf import dft, scf
    m = mf0.mol
    ao = [i for i, lab in enumerate(m.ao_labels()) if lab.split()[0] == str(target) and "1s" in lab]
    mo_a, occ_a = mf0.mo_coeff[0], mf0.mo_occ[0].copy()
    S = m.intor("int1e_ovlp")
    w = [(np.dot(mo_a[ao, i], (S @ mo_a[:, i])[ao]), i) for i in np.where(occ_a > 0)[0]]
    i_core = max(w)[1]
    occ = [occ_a.copy(), mf0.mo_occ[1].copy()]
    occ[0][i_core] = 0.0
    mf = dft.UKS(m).density_fit()
    mf.xc = "PBE0"
    mf.grids.level = 4
    mf.conv_tol = 1e-8
    mf.max_cycle = 150
    mf = scf.addons.mom_occ(mf, mf0.mo_coeff, occ)
    dm = mf.make_rdm1(mf0.mo_coeff, occ)
    mf.kernel(dm)
    return mf, float(max(w)[0])


def main():
    from pyscf import lib
    from pathgeom import read_xyz
    lib.num_threads(int(os.environ.get("QC_THREADS", "4")))
    S, X = read_xyz(os.path.join(ROOT, "structures", "HATCN_Ag1_deep.xyz"))
    X = np.asarray(X, float)
    kag = S.index("Ag")
    N = [i for i, s in enumerate(S) if s == "N"]
    dN = {i: float(np.linalg.norm(X[i] - X[kag])) for i in N}

    def is_nitrile(i):
        return any(S[j] == "C" and np.linalg.norm(X[j] - X[i]) < 1.25 for j in range(len(S)))

    aza = sorted([i for i in N if not is_nitrile(i)], key=lambda i: dN[i])
    nit = sorted([i for i in N if is_nitrile(i)], key=lambda i: dN[i])
    Sm = [s for i, s in enumerate(S) if i != kag]
    Xm = np.delete(X, kag, axis=0)
    rm = lambda i: i - (1 if i > kag else 0)  # noqa: E731  index after Ag removal
    cases = [
        ("neutral_aza", Sm, Xm, 0, 0, rm(aza[-1])),
        ("neutral_nitrile", Sm, Xm, 0, 0, rm(nit[-1])),
        ("Ag_aza_bonded", S, X, 0, 1, aza[0]),
        ("Ag_aza_far", S, X, 0, 1, aza[-1]),
        ("Ag_nitrile", S, X, 0, 1, nit[0]),
        ("anion_aza", Sm, Xm, -1, 1, rm(aza[-1])),
        ("anion_nitrile", Sm, Xm, -1, 1, rm(nit[-1])),
    ]
    put("_meta", {"aza_bonded_d_Ag": round(dN[aza[0]], 3), "aza_far_d_Ag": round(dN[aza[-1]], 3),
                  "nitrile_d_Ag": round(dN[nit[0]], 3)})
    for tag, Sx, Xx, ch, sp, tgt in cases:
        if tag in load() and "BE" in load()[tag]:
            continue
        t0 = time.time()
        m = build(Sx, Xx, ch, sp, tgt)
        g = ground(m)
        c, wcore = core_hole(g, tgt)
        be = (c.e_tot - g.e_tot) * EV
        put(tag, {"BE": round(float(be), 3), "gs_conv": bool(g.converged), "ch_conv": bool(c.converged),
                  "core_weight": round(wcore, 3), "sec": round(time.time() - t0)})
        print(tag, load()[tag], flush=True)
    d = load()
    ref_a, ref_n = d["neutral_aza"]["BE"], d["neutral_nitrile"]["BE"]
    put("shifts_vs_neutral", {k: round(d[k]["BE"] - (ref_n if "nitrile" in k else ref_a), 3)
                              for k, *_ in cases if k in d and "BE" in d[k]})
    print(json.dumps(load()["shifts_vs_neutral"], indent=1))


if __name__ == "__main__":
    main()
