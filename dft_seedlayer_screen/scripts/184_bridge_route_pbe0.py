"""PBE0-D3(BJ)/def2-SVP + ddCOSMO(acetonitrile) single points on the 183 geometries.

    QC_THREADS=4 python scripts/184_bridge_route_pbe0.py [species...]
    -> runs/design/bridge_route_pbe0.json  (energies in Eh, resumable)
    python scripts/184_bridge_route_pbe0.py --report

Why: GFN2-ALPB reaction energies that trade Cl-/CN- for a large delocalised
anion are not trustworthy (Cl -> CN substitution came out at -3.3 eV). The
isomer energies (neutral, same atoms) are the main question:
GFN2 puts the ketenimine 0.62 eV below Hz-C(CN)2-Hz and the carbodiimide
0.64 eV below Hz-N(CN)-Hz.
"""
import importlib.util
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, HERE)
GEO = os.path.join(ROOT, "structures", "bridge_route")
LANE = os.environ.get("LANE", "")        # parallel runs write separate files; --report merges them
OUT = os.path.join(ROOT, "runs", "design", f"bridge_route_pbe0{LANE}.json")
HA = 27.211386
_s = importlib.util.spec_from_file_location("b183", os.path.join(HERE, "183_bridge_route_thermo.py"))
b183 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(b183)


def read(name):
    L = open(os.path.join(GEO, f"{name}.xyz")).read().splitlines()
    ch = int(L[1].split("charge")[1])
    S = [l.split()[0] for l in L[2:2 + int(L[0])]]
    X = np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + int(L[0])]])
    return S, X, ch


def energy(S, X, ch):
    from pyscf import gto, dft
    m = gto.Mole()
    m.atom = [(s, tuple(x)) for s, x in zip(S, X)]
    m.unit = "Angstrom"; m.basis = "def2-svp"; m.charge = ch; m.spin = 0; m.verbose = 0
    m.max_memory = int(os.environ.get("QC_MEM_MB", "3500")); m.build()
    mf = dft.RKS(m).density_fit().ddCOSMO()
    mf.with_solvent.eps = 35.69          # acetonitrile
    mf.xc = "pbe0"; mf.grids.level = 3; mf.conv_tol = 1e-8; mf.max_cycle = 150
    mf.kernel()
    if not mf.converged:
        mf = mf.newton(); mf.kernel(mf.make_rdm1())
    from dftd3.interface import RationalDampingParam, DispersionModel
    d3 = DispersionModel(m.atom_charges(), m.atom_coords()).get_dispersion(
        RationalDampingParam(method="pbe0"), grad=False)["energy"]   # D3(BJ), Bohr in, Eh out
    e = float(mf.e_tot) + float(d3)
    conv = bool(mf.converged)
    try:
        mf.with_df.reset(); mf.with_df._cderi_to_save = None
    except AttributeError:
        pass
    return e, conv


def report():
    import glob
    E = {}
    for f in glob.glob(os.path.join(ROOT, "runs", "design", "bridge_route_pbe0*.json")):
        E.update({k: v["E"] for k, v in json.load(open(f)).items()})
    rows = {}
    for k, (r, p) in b183.RX.items():
        if all(s in E for s in list(r) + list(p)):
            rows[k] = round((sum(n * E[s] for s, n in p.items()) - sum(n * E[s] for s, n in r.items())) * HA, 2)
    g = json.load(open(os.path.join(ROOT, "runs", "design", "bridge_route_thermo.json")))["reactions_eV"]
    for k, v in rows.items():
        print(f"{k:45s} PBE0 {v:+.2f}   GFN2 {g[k]:+.2f}")
    return rows


def main(names):
    from pyscf import lib
    lib.num_threads(int(os.environ.get("QC_THREADS", "4")))
    for name in names:
        d = json.load(open(OUT)) if os.path.exists(OUT) else {}
        if name in d:
            continue
        S, X, ch = read(name)
        e, conv = energy(S, X, ch)
        d = json.load(open(OUT)) if os.path.exists(OUT) else {}
        d[name] = {"E": e, "conv": conv, "charge": ch}
        json.dump(d, open(OUT, "w"), indent=1)
        print(name, "pbe0", round(e, 6), conv, flush=True)


if __name__ == "__main__":
    if sys.argv[1:] == ["--report"]:
        report()
    else:
        main(sys.argv[1:])
