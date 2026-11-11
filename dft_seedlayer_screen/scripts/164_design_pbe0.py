"""PBE0 anchors for the design library: the active-learning step.

    QC_THREADS=4 FRAG_RCUT=7.0 python scripts/164_design_pbe0.py NAME [NAME ...]
        -> runs/design/pbe0_anchor.json (checkpointed per molecule and stage)

The Delta-ML map (script 163) has three training points above 1 eV (HATCN,
F4TCNQ, F6-TCNNQ), and the library's best HAT derivatives sit exactly there.
So the map is not trusted in that region until it has been fed PBE0 numbers
from the region: the molecules picked here are the ones the GP is most
uncertain about among those that would matter (high predicted E_b, passing
the hard constraints), plus HATCN at its GFN2 design geometry as the
consistency check against the production HATCN value (1.631 eV).

Level and recipe are script 149's: PBE0-D3(BJ)/def2-SVP, Ag on the GFN2 site
(deepest site of script 162), fragment within FRAG_RCUT of Ag with H caps,
E_b = elec - [D3(mol+Ag) - D3(mol)] with D3 on the whole molecule. SCF
starts from a converged PBE density (script 157's fix for acceptor-Ag
complexes that would not converge from the default guess).
"""
import importlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
D = os.path.join(ROOT, "runs", "design")
OUT = os.path.join(D, "pbe0_anchor.json")
EV = 27.211386


def pbe_guess(m):
    from pyscf import dft
    mf = (dft.UKS(m) if m.spin else dft.RKS(m)).density_fit()
    mf.xc = "PBE"
    mf.grids.level = 2
    mf.conv_tol = 1e-6
    mf.max_cycle = 150
    mf.level_shift = 0.3
    mf.damp = 0.3
    mf.kernel()
    return mf.make_rdm1()


CK = os.path.join(D, "pbe0_ckpt")


def scf_ck(m, tag):
    """q.scf with a PySCF chkfile, so a container restart resumes from the
    last SCF iterate instead of from scratch (and a PBE guess only once)."""
    from pyscf import dft, scf as pscf
    os.makedirs(CK, exist_ok=True)
    ck = os.path.join(CK, tag + ".chk")
    mf = (dft.UKS(m) if m.spin else dft.RKS(m)).density_fit()
    mf.xc = "PBE0"
    mf.grids.level = 4
    mf.conv_tol = 1e-9
    mf.max_cycle = 200
    mf.level_shift = 0.2
    if os.path.exists(ck):
        dm0 = mf.from_chk(ck)
    else:
        dm0 = pbe_guess(m)
    mf.chkfile = ck
    e = mf.kernel(dm0)
    if not mf.converged:
        mf = mf.newton()
        mf.chkfile = ck
        e = mf.kernel(mf.make_rdm1())
    return mf, e


def load():
    return json.load(open(OUT)) if os.path.exists(OUT) else {}


def save(name, key, val):
    d = load()
    d.setdefault(name, {})[key] = val
    json.dump(d, open(OUT, "w"), indent=1)


def main():
    from pyscf import lib
    lib.num_threads(int(os.environ.get("QC_THREADS", "4")))
    q = importlib.import_module("qclocal")
    s149 = importlib.import_module("149_fragment_pbe0")
    from pathgeom import read_xyz
    sites = {json.loads(l)["name"]: json.loads(l)
             for l in open(os.path.join(D, "sites.jsonl")) if l.strip()}
    d = load()
    if "_E_ag" not in d:
        mf, e = q.scf(q.mol(["Ag"], [[0, 0, 0]], 2))
        save("_E_ag", "E", float(e))
    e_ag = load()["_E_ag"]["E"]
    for name in sys.argv[1:]:
        st = sites[name]
        S, X = read_xyz(os.path.join(ROOT, "structures", "design", f"{name}.xyz"))
        X = np.asarray(X, float)
        ag = np.array(st["top3"][st["deepest"]]["pos"])
        FS, FX, ncap = s149.fragment(S, X, ag)
        rec = load().get(name, {})
        t0 = time.time()
        if "E_frag" not in rec:
            m = q.mol(FS, FX, 1)
            mf, e = scf_ck(m, name + "_frag")
            if not mf.converged:
                raise RuntimeError(f"{name}: fragment SCF not converged")
            save(name, "E_frag", float(e))
            print(f"{name} fragment {len(FS)} atoms done {time.time()-t0:.0f}s", flush=True)
        if "E_cx" not in load()[name]:
            m = q.mol(FS + ["Ag"], np.vstack([FX, ag]), 2)
            mf, e = scf_ck(m, name + "_cx")
            if not mf.converged:
                raise RuntimeError(f"{name}: complex SCF not converged")
            save(name, "E_cx", float(e))
            save(name, "q_Ag_mulliken", round(float(q.mulliken(mf)[-1]), 4))
            print(f"{name} complex done {time.time()-t0:.0f}s", flush=True)
        rec = load()[name]
        elec = (rec["E_frag"] + e_ag - rec["E_cx"]) * EV
        d3 = (q.d3bj(S + ["Ag"], np.vstack([X, ag])) - q.d3bj(S, X)) * EV
        save(name, "result", {"E_b": round(elec - d3, 4), "elec": round(elec, 4),
                              "disp": round(-d3, 4), "frag_atoms": len(FS), "caps": ncap,
                              "mol_atoms": len(S), "site": st["deepest"],
                              "E_b_gfn2": st["E_b_gfn2"], "rcut": s149.R_CUT})
        print(name, load()[name]["result"], flush=True)


if __name__ == "__main__":
    main()
