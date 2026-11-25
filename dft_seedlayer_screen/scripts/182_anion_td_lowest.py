"""Anion TD-DFT on the lowest SCF solution.

    QC_THREADS=4 python scripts/182_anion_td_lowest.py name [name...]
    -> runs/design/tddft.json, key "anion_n12_lowest"

Why: script 166 converged the anion UKS reference from the default guess at the
neutral design geometry. For HATCN that guess lands on a solution 0.24 eV above
the lowest one at the GFN2 neutral geometry (see 181). The anion absorption
used as the transparency reference (vis f sum 0.30) may therefore belong to
the wrong state.

Protocol: at the same geometry as 166 (structures/design/{name}.xyz), anion UKS
from (a) the default guess with level shift and (b) the density of the
relaxed anion (structures/electrical/{name}_T30_an.xyz). Keep the lower
solution, then TDA-PBE0/def2-SVP with 12 roots, same settings as 166. The
energies of both solutions are stored with the result.
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, HERE)
OUT = os.path.join(ROOT, "runs", "design", "tddft.json")
HA = 27.211386
_s = importlib.util.spec_from_file_location("t166", os.path.join(HERE, "166_design_tddft.py"))
t166 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(t166)


def uks(S, X, dm0=None):
    from pyscf import gto, dft
    m = gto.Mole()
    m.atom = [(s, tuple(x)) for s, x in zip(S, X)]
    m.unit = "Angstrom"; m.basis = "def2-svp"; m.charge = -1; m.spin = 1; m.verbose = 0
    m.max_memory = int(os.environ.get("QC_MEM_MB", "3500")); m.build()
    mf = dft.UKS(m).density_fit()
    mf.xc = "pbe0"; mf.grids.level = 3; mf.conv_tol = 1e-8; mf.max_cycle = 100
    if dm0 is None:
        mf.level_shift = 0.1
    mf.kernel(dm0)
    if not mf.converged:
        mf = mf.newton(); mf.kernel(mf.make_rdm1())
    mf.level_shift = 0
    return mf


def main(names):
    from pyscf import lib, tdscf
    from pathgeom import read_xyz
    lib.num_threads(int(os.environ.get("QC_THREADS", "4")))
    for name in names:
        d = json.load(open(OUT))
        if "anion_n12_lowest" in d.get(name, {}):
            continue
        t0 = time.time()
        S, X = read_xyz(os.path.join(ROOT, "structures", "design", f"{name}.xyz"))
        X = np.asarray(X, float)
        S2, XI = read_xyz(os.path.join(ROOT, "structures", "electrical", f"{name}_T30_an.xyz"))
        assert list(S) == list(S2), "atom order differs"
        a = uks(S, X)
        e_default = a.e_tot
        guess = uks(S, np.asarray(XI, float)).make_rdm1()
        b = uks(S, X, guess)
        e_cross = b.e_tot
        mf = a if a.e_tot <= b.e_tot else b
        del a, b
        td = tdscf.TDA(mf)
        td.nstates = 12; td.conv_tol = 1e-4; td.max_cycle = 4
        td.kernel()
        e = np.array(td.e) * HA
        f = np.array(td.oscillator_strength())
        st = [{"eV": round(float(x), 3), "nm": round(1239.84 / float(x), 1), "f": round(float(y), 4)}
              for x, y in zip(e, f)]
        row = {"states": st, **t166.summarise(st, True),
               "E_default": e_default, "E_cross": e_cross,
               "shift_default_minus_cross_eV": round((e_default - e_cross) * HA, 3)}
        d = json.load(open(OUT)); d.setdefault(name, {})["anion_n12_lowest"] = row
        json.dump(d, open(OUT, "w"), indent=1)
        print(name, "anion_n12_lowest", {k: v for k, v in row.items() if k not in ("states", "E_default", "E_cross")},
              f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
