"""CCSD(T) benchmark of the Ag-adatom binding the screen is built on.

    QC_THREADS=4 python scripts/171_ccsdt_benchmark.py [--test]   (resumable)
    -> runs/ccsdt_benchmark.json

Literature gap 9 (notes/lit_survey): there is no CCSD(T) reference for a
neutral Ag atom on an N-heterocycle or a nitrile, and the only metal-aromatic
one (Ag-benzene, Granatier et al. JCTC 2011) says DFT-D3 overbinds pi
contact 1.7x. The ranking's ETL-class values are dispersion-heavy, so the
question is how much PBE0-D3BJ/def2-SVP (the screen's level, no CP)
overbinds per site TYPE:

  pyridine   N lone pair (the ETL/phen anchor type)
  MeCN       nitrile N (F4TCNQ / HATCN-nitrile type)
  pyrazine   aza N of a diazine (HAT-core type)
  benzene    pi hollow (the host/HTL type; literature check)

Protocol per system:
  1. Ag position: 1-D scan along the binding direction at PBE0-D3BJ/def2-SVP
     (molecule at its GFN2 geometry, rigid), parabola through the minimum.
  2. At that geometry: PBE0-D3BJ/def2-SVP (screen level, no CP), PBE0-D3BJ/
     def2-TZVP (CP), and UCCSD(T)/def2-TZVP (CP, frozen core incl. Ag 4s4p),
     plus UCCSD(T)/def2-SVP at d-0.15/d/d+0.15 A to see whether the CCSD(T)
     minimum sits at a different distance.
CP = Boys-Bernardi counterpoise, ghost basis on the partner.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "runs", "ccsdt_benchmark.json")
EV = 27.211386

SYSTEMS = {  # smiles, binding atom selector
    "MeCN": "CC#N",
    "pyridine": "c1ccncc1",
    "pyrazine": "c1cnccn1",
    "benzene": "c1ccccc1",
}


def load():
    return json.load(open(OUT)) if os.path.exists(OUT) else {}


def put(name, key, val):
    d = load()
    d.setdefault(name, {})[key] = val
    json.dump(d, open(OUT, "w"), indent=1)


def geometry(smi):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from ase import Atoms
    from ase.optimize import BFGS
    from tblite.ase import TBLite
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(m, randomSeed=3)
    AllChem.MMFFOptimizeMolecule(m)
    a = Atoms([x.GetSymbol() for x in m.GetAtoms()], positions=m.GetConformer().GetPositions())
    a.calc = TBLite(method="GFN2-xTB", verbosity=0)
    BFGS(a, logfile=None).run(fmax=0.01, steps=300)
    return list(a.get_chemical_symbols()), a.positions.copy()


def direction(S, X, name):
    """Anchor point and unit vector along which Ag approaches."""
    if name == "benzene":
        c = X[[i for i, s in enumerate(S) if s == "C"]].mean(0)
        u = np.linalg.svd(X[[i for i, s in enumerate(S) if s == "C"]] - c)[2][-1]
        return c, u
    i = [k for k, s in enumerate(S) if s == "N"][0]
    nb = [j for j in range(len(S)) if j != i and np.linalg.norm(X[j] - X[i]) < 1.6]
    u = -sum((X[j] - X[i]) / np.linalg.norm(X[j] - X[i]) for j in nb)
    return X[i], u / np.linalg.norm(u)


def mol(S, X, spin, basis, ghost=()):
    from pyscf import gto
    m = gto.Mole()
    m.atom = [((f"ghost-{s}" if k in ghost else s), tuple(x)) for k, (s, x) in enumerate(zip(S, X))]
    m.unit = "Angstrom"
    m.basis = basis
    m.ecp = {"Ag": basis}
    m.spin = spin
    m.verbose = 0
    m.max_memory = int(os.environ.get("QC_MEM_MB", "9000"))
    m.build()
    return m


def dft_e(m):
    from pyscf import dft
    mf = (dft.UKS(m) if m.spin else dft.RKS(m))
    mf.xc = "PBE0"
    mf.grids.level = 4
    mf.conv_tol = 1e-10
    mf.max_cycle = 200
    mf.level_shift = 0.2
    e = mf.kernel()
    if not mf.converged:
        mf = mf.newton()
        e = mf.kernel(mf.make_rdm1())
    return e


def cc_e(m):
    """UCCSD(T) (RCCSD(T) for closed shell) on an HF reference, frozen core."""
    from pyscf import scf, cc
    mf = (scf.UHF(m) if m.spin else scf.RHF(m))
    mf.conv_tol = 1e-10
    mf.max_cycle = 200
    mf.level_shift = 0.2
    mf.kernel()
    if not mf.converged:
        mf = mf.newton()
        mf.kernel()
    # frozen core: all-electron 1s on C/N; Ag (ECP28) 4s4p = 4 orbitals
    nfz = sum({"C": 1, "N": 1, "Ag": 4}.get(m.atom_pure_symbol(k), 0)
              for k in range(m.natm) if m.atom_charge(k) > 0)
    c = cc.CCSD(mf, frozen=nfz)
    c.direct = True
    c.max_cycle = 100
    c.conv_tol = 1e-8
    c.kernel()
    et = c.ccsd_t()
    return float(mf.e_tot + c.e_corr + et), bool(c.converged)


def d3(S, X):
    import qclocal as q
    return q.d3bj(S, X)


def run(name, smi, test=False):
    d = load().get(name, {})
    if "geom" not in d:
        S, X = geometry(smi)
        put(name, "geom", {"S": S, "X": X.round(5).tolist()})
    g = load()[name]["geom"]
    S, X = g["S"], np.array(g["X"])
    anchor, u = direction(S, X, name)
    SA = S + ["Ag"]
    n = len(S)
    pos = lambda dist: np.vstack([X, anchor + dist * u])  # noqa: E731
    # 1. PBE0-D3/SVP scan
    if "scan" not in load()[name]:
        e_ag = dft_e(mol(["Ag"], [[0, 0, 0]], 1, "def2-svp"))
        e_m = dft_e(mol(S, X, 0, "def2-svp"))
        rows = []
        grid = (2.0, 2.2, 2.4, 2.6, 2.8, 3.1) if name != "benzene" else (2.4, 2.7, 3.0, 3.3, 3.6, 4.0)
        for dist in grid:
            P = pos(dist)
            e = dft_e(mol(SA, P, 1, "def2-svp")) + d3(SA, P)
            rows.append([dist, (e_m + d3(S, X) + e_ag - e) * EV])
            print(name, "scan", dist, round(rows[-1][1], 4), flush=True)
        put(name, "scan", rows)
    rows = np.array(load()[name]["scan"])
    k = int(rows[:, 1].argmax())
    k = min(max(k, 1), len(rows) - 2)
    c = np.polyfit(rows[k - 1:k + 2, 0], rows[k - 1:k + 2, 1], 2)
    dmin = float(-c[1] / (2 * c[0])) if c[0] < 0 else float(rows[k, 0])
    put(name, "d_min", round(dmin, 3))
    P = pos(dmin)

    def eb(method, basis, cp, P=P):
        fn = dft_e if method == "pbe0" else (lambda m: cc_e(m)[0])
        e_cx = fn(mol(SA, P, 1, basis))
        if cp:
            e_m = fn(mol(SA, P, 0, basis, ghost=(n,)))
            e_a = fn(mol(SA, P, 1, basis, ghost=tuple(range(n))))
        else:
            e_m = fn(mol(S, X, 0, basis))
            e_a = fn(mol(["Ag"], [[0, 0, 0]], 1, basis))
        disp = (d3(SA, P) - d3(S, X)) * EV if method == "pbe0" else 0.0
        return round((e_m + e_a - e_cx) * EV - disp, 4)

    tasks = [("pbe0_svp_noCP", "pbe0", "def2-svp", False),
             ("pbe0_tzvp_CP", "pbe0", "def2-tzvp", True),
             ("ccsdt_svp_CP", "cc", "def2-svp", True)]
    if not test:
        tasks.append(("ccsdt_tzvp_CP", "cc", "def2-tzvp", True))
    for key, meth, bas, cp in tasks:
        if key in load()[name]:
            continue
        t0 = time.time()
        put(name, key, eb(meth, bas, cp))
        print(name, key, load()[name][key], f"{time.time()-t0:.0f}s", flush=True)
    if not test and "ccsdt_svp_CP_curve" not in load()[name]:
        curve = []
        for dd in (-0.15, 0.15):
            curve.append([round(dmin + dd, 3), eb("cc", "def2-svp", True, P=pos(dmin + dd))])
        curve.insert(1, [round(dmin, 3), load()[name]["ccsdt_svp_CP"]])
        put(name, "ccsdt_svp_CP_curve", curve)
        print(name, "curve", curve, flush=True)


def main():
    from pyscf import lib
    lib.num_threads(int(os.environ.get("QC_THREADS", "4")))
    test = "--test" in sys.argv
    names = [a for a in sys.argv[1:] if not a.startswith("--")] or list(SYSTEMS)
    for nm in names:
        run(nm, SYSTEMS[nm], test)
    d = load()
    summ = {}
    for nm in SYSTEMS:
        r = d.get(nm, {})
        if "ccsdt_tzvp_CP" in r:
            summ[nm] = {"d_min_pbe0": r["d_min"], "screen_level": r["pbe0_svp_noCP"],
                        "pbe0_tzvp_CP": r["pbe0_tzvp_CP"], "ccsdt_tzvp_CP": r["ccsdt_tzvp_CP"],
                        "screen_over_ccsdt": round(r["pbe0_svp_noCP"] / r["ccsdt_tzvp_CP"], 3)
                        if r["ccsdt_tzvp_CP"] > 0 else None}
    if summ:
        d["summary"] = summ
        json.dump(d, open(OUT, "w"), indent=1)
        print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
