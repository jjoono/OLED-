"""Cheap properties for the design library (synthesizable members only).

    python scripts/161_design_props.py [--workers 4]  ->  runs/design/props.jsonl (resumable)

Per molecule, one CPU each:
  geometry   RDKit ETKDG + MMFF, then GFN2-xTB (tblite) to fmax 0.1 eV/A
  IP, EA     IPEA1-xTB vertical, at the GFN2 geometry. IPEA1 carries a
             constant offset against experiment; only differences to HATCN
             are used, so it cancels.
  gap        GFN2 HOMO-LUMO gap: a transparency proxy, read relative to HATCN
             (a molecule with a smaller gap than HATCN starts absorbing in the
             visible first). Not an optical gap.
  shape      3D-ness and thickness from the inertia tensor; crystallisation
             risk as in script 53 (flat + rigid + symmetric + small).
Mixed-route HAT (synth_class 3, no one-pot route) is skipped.
"""
import json
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
LIB = os.path.join(ROOT, "runs", "design", "library.json")
OUT = os.path.join(ROOT, "runs", "design", "props.jsonl")
XYZ = os.path.join(ROOT, "structures", "design")
BOHR = 0.52917721092
HA = 27.211386
SHIFT = 4.8455
MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Cl": 35.45}


def shape(S, X):
    m = np.array([MASS.get(s, 12.0) for s in S])
    x = X - (X * m[:, None]).sum(0) / m.sum()
    I = np.zeros((3, 3))
    for mi, xi in zip(m, x):
        I += mi * (np.dot(xi, xi) * np.eye(3) - np.outer(xi, xi))
    ev = np.sort(np.linalg.eigvalsh(I))
    thick = np.linalg.svd(x - x.mean(0))[1][2] / np.sqrt(len(S))
    return ev[0] / ev[2] + ev[1] / ev[2], thick


def risk(threeD, nrot, sym_frac, mw):
    """Script 53's crystallisation-risk score (0-100; HATCN 92)."""
    flat = np.clip((1.20 - threeD) / 0.20, 0, 1)
    rigid = np.clip((6 - nrot) / 6, 0, 1)
    sym = np.clip((0.55 - sym_frac) / 0.45, 0, 1)
    small = np.clip((750 - mw) / 500, 0, 1)
    return float(100 * (0.40 * flat + 0.20 * rigid + 0.25 * sym + 0.15 * small))


def sp(method, Z, X, charge, uhf):
    from tblite.interface import Calculator
    c = Calculator(method, Z, X / BOHR, charge=charge, uhf=uhf)
    c.set("verbosity", 0)
    try:
        c.set("max-iter", 300)
    except Exception:
        pass
    return c.singlepoint()


def work(item):
    name, d = item
    t0 = time.time()
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from ase import Atoms
    from ase.optimize import BFGS
    from tblite.ase import TBLite
    m = Chem.AddHs(Chem.MolFromSmiles(d["smiles"]))
    if AllChem.EmbedMolecule(m, randomSeed=7) != 0:
        AllChem.EmbedMolecule(m, randomSeed=7, useRandomCoords=True)
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    S = [a.GetSymbol() for a in m.GetAtoms()]
    X = m.GetConformer().GetPositions()
    a = Atoms(S, positions=X)
    a.calc = TBLite(method="GFN2-xTB", verbosity=0)
    try:
        BFGS(a, logfile=None, maxstep=0.2).run(fmax=0.1, steps=400)
    except Exception as e:
        return {"name": name, "error": f"opt: {e!r}"[:200]}
    X = a.positions.copy()
    Z = np.array(a.numbers)
    r0 = sp("GFN2-xTB", Z, X, 0, 0)
    eps = np.array(r0.get("orbital-energies")) * HA
    occ = np.array(r0.get("orbital-occupations"))
    homo = eps[occ > 0.5].max()
    lumo = eps[occ < 0.5].min()
    e0 = sp("IPEA1-xTB", Z, X, 0, 0).get("energy")
    ep = sp("IPEA1-xTB", Z, X, 1, 1).get("energy")
    em = sp("IPEA1-xTB", Z, X, -1, 1).get("energy")
    td, thick = shape(S, X)
    os.makedirs(XYZ, exist_ok=True)
    with open(os.path.join(XYZ, f"{name}.xyz"), "w") as f:
        f.write(f"{len(S)}\n{name} GFN2\n")
        for s, x in zip(S, X):
            f.write(f"{s} {x[0]:.5f} {x[1]:.5f} {x[2]:.5f}\n")
    # xtb's empirical IPEA shift (4.8455 eV), as in `xtb --vipea`
    return {"name": name, "IP": round(float((ep - e0) * HA - SHIFT), 3),
            "EA": round(float((e0 - em) * HA - SHIFT), 3),
            "gap_gfn2": round(float(lumo - homo), 3), "threeD": round(float(td), 3),
            "thick_A": round(float(thick), 3),
            "cryst_risk": round(risk(td, d["nrot"], d["sym_frac"], d["MW"]), 1),
            "natoms": len(S), "sec": round(time.time() - t0)}


def main():
    from multiprocessing import Pool
    workers = 4
    if "--workers" in sys.argv:
        workers = int(sys.argv[sys.argv.index("--workers") + 1])
    lib = json.load(open(LIB))
    done = set()
    if os.path.exists(OUT):
        done = {json.loads(l)["name"] for l in open(OUT) if l.strip()}
    todo = [(k, v) for k, v in lib.items()
            if k not in done and v.get("synth_class", 1) != 3]
    # small first: early results calibrate the references
    todo.sort(key=lambda kv: (kv[1]["family"] != "reference", kv[1]["heavy"]))
    print(f"{len(todo)} to do, {len(done)} done", flush=True)
    with Pool(workers) as p, open(OUT, "a") as f:
        for r in p.imap_unordered(work, todo):
            f.write(json.dumps(r) + "\n")
            f.flush()
            print(r["name"], r.get("EA"), r.get("gap_gfn2"), r.get("cryst_risk"), r.get("sec"),
                  r.get("error", ""), flush=True)


if __name__ == "__main__":
    main()
