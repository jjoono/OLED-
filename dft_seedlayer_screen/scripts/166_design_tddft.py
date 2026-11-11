"""Transparency check for design finalists: neutral AND anion, TDA-wB97X/def2-SVP.

    QC_THREADS=4 python scripts/166_design_tddft.py NAME [NAME ...]
        -> runs/design/tddft.json (checkpointed per molecule and charge state)

Spec axis 4b (notes/DESIGN_SPEC_CGL_SEED.md). The first seed monolayer is
reduced by Ag, so the anion spectrum matters as much as the neutral one: the
F4TCNQ anion's 2.4 eV band is what disqualifies it optically, while the
HATCN anion has only a weak 1.96 eV pair (f 0.15 each, script 44).

Level as script 44 (range-separated hybrid for the CT-like anion states),
GFN2 geometry of the neutral for both (vertical), 12 roots. Reported:
  neutral  first state with f > 0.02, and any such state below 3.1 eV
  anion    sum of f between 1.65 and 3.1 eV (the visible), first bright state
HATCN is computed here at the same level and geometry protocol, so the
comparison is internal.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "runs", "design", "tddft.json")
HA = 27.211386
VIS = (1.65, 3.10)


def run(S, X, charge, nroots=12):
    from pyscf import gto, dft, tdscf
    m = gto.Mole()
    m.atom = [(s, tuple(x)) for s, x in zip(S, X)]
    m.unit = "Angstrom"
    m.basis = "def2-svp"
    m.charge = charge
    m.spin = 1 if charge else 0
    m.verbose = 0
    m.max_memory = int(os.environ.get("QC_MEM_MB", "8000"))
    m.build()
    mf = (dft.UKS(m) if m.spin else dft.RKS(m)).density_fit()
    mf.xc = "wb97x"
    mf.grids.level = 3
    mf.conv_tol = 1e-8
    mf.max_cycle = 60
    mf.kernel()
    if not mf.converged:
        mf = mf.newton()
        mf.kernel(mf.make_rdm1())
    td = tdscf.TDA(mf)
    td.nstates = nroots
    td.kernel()
    e = np.array(td.e) * HA
    f = np.array(td.oscillator_strength())
    return [{"eV": round(float(a), 3), "nm": round(1239.84 / float(a), 1), "f": round(float(b), 4)}
            for a, b in zip(e, f)]


def summarise(states, anion):
    bright = [s for s in states if s["f"] > 0.02]
    out = {"first_bright_eV": bright[0]["eV"] if bright else None,
           "first_bright_f": bright[0]["f"] if bright else None}
    if anion:
        out["vis_sum_f"] = round(sum(s["f"] for s in states if VIS[0] <= s["eV"] <= VIS[1]), 3)
    else:
        out["bright_below_3p1"] = [s for s in bright if s["eV"] < VIS[1]]
    return out


def main():
    from pyscf import lib
    from pathgeom import read_xyz
    lib.num_threads(int(os.environ.get("QC_THREADS", "4")))
    for name in sys.argv[1:]:
        S, X = read_xyz(os.path.join(ROOT, "structures", "design", f"{name}.xyz"))
        X = np.asarray(X, float)
        for tag, ch in (("neutral", 0), ("anion", -1)):
            d = json.load(open(OUT)) if os.path.exists(OUT) else {}
            if tag in d.get(name, {}):
                continue
            t0 = time.time()
            st = run(S, X, ch)
            d = json.load(open(OUT)) if os.path.exists(OUT) else {}
            d.setdefault(name, {})[tag] = {"states": st, **summarise(st, ch != 0)}
            json.dump(d, open(OUT, "w"), indent=1)
            print(name, tag, {k: v for k, v in d[name][tag].items() if k != "states"},
                  f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
