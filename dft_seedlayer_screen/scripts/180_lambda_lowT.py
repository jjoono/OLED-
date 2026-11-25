"""Reorganisation energy with integer-like occupations in the GFN2 geometries.

    QC_THREADS=4 python scripts/180_lambda_lowT.py [names...]  -> runs/design/electrical_T30.json

Why: script 179 optimised the anions with tblite's default electronic
temperature (300 K). For HATCN (D3h, doubly degenerate LUMO) the Fermi
smearing splits the extra electron over both LUMOs, so the anion stays
symmetric and the Jahn-Teller relaxation is lost; a random kick does not help
(179 with kick: lambda 0.034 eV, unchanged). At 30 K the occupation is close
to integer and the anion distorts. GFN2 lambda (three random starts each):
    HATCN 0.141   Hz-NCN-Hz 0.081   Hz-O-Hz 0.079   (300 K: HATCN 0.03-0.08)
Same protocol for every molecule: GFN2 30 K, fmax 0.01 eV/A, three kicked
anion starts, lowest anion kept; then PBE0/def2-SVP four-point lambda and
EA as in 179.
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "runs", "design", "electrical_T30.json")
GEO = os.path.join(ROOT, "structures", "electrical")
HA = 27.211386
TEL = float(os.environ.get("XTB_TEL", "30"))

_s = importlib.util.spec_from_file_location("e179", os.path.join(HERE, "179_electrical.py"))
e179 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(e179)


def _calc(ch):
    from tblite.ase import TBLite
    return TBLite(method="GFN2-xTB", charge=ch, multiplicity=1 + abs(ch) % 2, verbosity=0,
                  electronic_temperature=TEL, max_iterations=500)


def opt(S, X, ch):
    from ase import Atoms
    from ase.optimize import BFGS
    a = Atoms(S, positions=X)
    a.calc = _calc(ch)
    BFGS(a, logfile=None, maxstep=0.1).run(fmax=0.01, steps=3000)
    return a.positions.copy(), float(a.get_potential_energy())


def sp(S, X, ch):
    from ase import Atoms
    a = Atoms(S, positions=X)
    a.calc = _calc(ch)
    return float(a.get_potential_energy())


def load():
    return json.load(open(OUT)) if os.path.exists(OUT) else {}


def save(name, upd):
    out = load()
    out.setdefault(name, {}).update(upd)
    json.dump(out, open(OUT, "w"), indent=1)


def main(names):
    from pyscf import lib
    lib.num_threads(int(os.environ.get("QC_THREADS", "4")))
    for name in names:
        row = load().get(name, {})
        if "lambda_eV" in row:
            continue
        ion = +1 if name == "NPB" else -1
        gN = os.path.join(GEO, f"{name}_T30_0.xyz")
        gI = os.path.join(GEO, f"{name}_T30_{'cat' if ion > 0 else 'an'}.xyz")
        t0 = time.time()
        if "gfn2_lambda_eV" not in row:
            src = os.path.join(GEO, f"{name}_0.xyz")
            if not os.path.exists(src):   # not run through 179: start from the design structure
                src = os.path.join(ROOT, "structures", "design", f"{name}.xyz")
            S, X0 = e179.read(src)
            XN, EN = opt(S, X0, 0)
            runs = []
            for seed in (1, 2, 3):
                XI, EI = opt(S, XN + np.random.default_rng(seed).normal(0, 0.05, XN.shape), ion)
                runs.append((EI, XI, (sp(S, XN, ion) - EI) + (sp(S, XI, 0) - EN)))
            EI, XI, lam = min(runs, key=lambda r: r[0])
            e179.write(gN, S, XN, f"{name} neutral GFN2 Tel={TEL:g}K")
            e179.write(gI, S, XI, f"{name} ion GFN2 Tel={TEL:g}K lowest of 3 starts")
            row = {"gfn2_lambda_eV": round(lam, 3), "gfn2_lambda_seeds": [round(r[2], 3) for r in runs]}
            save(name, row)
            print(name, "GFN2", row, flush=True)
        S, XN = e179.read(gN)
        XI = e179.read(gI)[1]
        e = {t: row[f"E_{t}"] for t in ("N@N", "I@N", "I@I", "N@I") if f"E_{t}" in row}
        for tag, X, ch in (("N@N", XN, 0), ("I@N", XN, ion), ("I@I", XI, ion), ("N@I", XI, 0)):
            if tag in e:
                continue
            e[tag], conv, pr = e179.pbe0(S, X, ch)
            upd = {f"E_{tag}": e[tag], f"conv_{tag}": conv}
            if pr is not None and tag == "I@I":
                upd["spin_PR_heavy_atoms"] = round(pr, 1)
            save(name, upd)
        r = {}
        if ion < 0:
            r["EA_vert_eV"] = round((e["N@N"] - e["I@N"]) * HA, 3)
            r["EA_ad_eV"] = round((e["N@N"] - e["I@I"]) * HA, 3)
        else:
            r["IP_vert_eV"] = round((e["I@N"] - e["N@N"]) * HA, 3)
            r["IP_ad_eV"] = round((e["I@I"] - e["N@N"]) * HA, 3)
        r["lambda_eV"] = round(((e["I@N"] - e["I@I"]) + (e["N@I"] - e["N@N"])) * HA, 3)
        r["seconds"] = round(time.time() - t0)
        save(name, r)
        print(name, "PBE0", r, flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or ["HATCN", "Hz-NCN-Hz", "Hz-O-Hz"])
