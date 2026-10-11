"""Does the charge-capacity veto separate the oxides the literature has already
ranked?

    python scripts/146_oxide_capacity_xtb.py ZnO        # one oxide
    python scripts/146_oxide_capacity_xtb.py            # all

The literature synthesis (notes/LIT_SEED_SYNTHESIS_20261011.md) left one
descriptor standing as a veto: a seed fails when it keeps oxidising silver as
clusters grow -- MoO3 percolates at 7 nm but stays porous, in the literature and
in this lab. If that is the mechanism, oxides the literature finds GOOD should
hand charge back, and oxides as deep and reducible as MoO3 should not:

    ZnO (10-10)   industry-standard Ag seed, weak acceptor (EA ~4.3 eV)
    TiO2 rutile (110)   reducible, weak acceptor (EA ~4.2 eV); Ag on
                        amorphous TiO2 is discontinuous where ZnO is continuous
    Cu2O (111)    copper oxide seed; oxidised-Cu seeds percolate at ~1-2 nm
    V2O5 (001)    layered, deep acceptor like MoO3 (EA ~6.7 eV); no Ag-seed
                  report found -- this one is a prediction

CuO itself, the ~1 nm seed of the ICFO work, is not computed: it is an
antiferromagnetic Cu(II) Mott insulator, and GFN2-xTB places its empty states
wrongly, so the silver would have the wrong place to put its electron. Cu2O,
closed-shell, is the copper-oxide stand-in, and the substitution is stated
wherever the number is used.

Protocol identical to script 137 (MoO3): periodic slab from a COD structure,
substrate frozen, silver relaxed, 1000 K smearing, Ag1 from five starts over the
surface atoms nearest the cell centre, Ag2..Ag4 grown from the deepest previous
cluster. Charges are GFN2's; the ordering, not the values, is the result.
"""
import os, sys, json, importlib

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.build import surface

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
s137 = importlib.import_module("137_moo3_slab_xtb")

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
ST = os.path.join(ROOT, "structures")
OUT = os.path.join(ROOT, "slabs", "oxides_xtb")
RES = os.path.join(ROOT, "runs", "oxide_capacity_xtb.json")
D_AG = 2.60


def v2o5_layer():
    b = read(os.path.join(ST, "V2O5_COD9012221.cif"))
    f = b.get_scaled_positions()
    f[:, 2] = (f[:, 2] + 0.5) % 1.0 - 0.5          # one layer, centred on z = 0
    a, bb, c = b.cell.lengths()
    pos = f * [a, bb, c]
    pos[:, 2] -= pos[:, 2].min()
    cell = np.diag([a, bb, pos[:, 2].max() + 20.0])
    return Atoms(b.get_chemical_symbols(), positions=pos, cell=cell, pbc=True)


def build(tag):
    if tag == "ZnO":
        s = surface(read(os.path.join(ST, "ZnO_wurtzite_COD9004178.cif")), (1, 0, 0), 4,
                    vacuum=10).repeat((3, 2, 1))
    elif tag == "TiO2":
        s = surface(read(os.path.join(ST, "TiO2_rutile_COD9004141.cif")), (1, 1, 0), 3,
                    vacuum=10).repeat((2, 3, 1))
    elif tag == "Cu2O":
        s = surface(read(os.path.join(ST, "Cu2O_cuprite_COD9007497.cif")), (1, 1, 1), 3,
                    vacuum=10).repeat((2, 2, 1))
    elif tag == "V2O5":
        s = v2o5_layer().repeat((1, 3, 1))
    s.pbc = True
    s.positions[:, 2] += 2.0 - s.positions[:, 2].min()
    return s


def relax(atoms, n_sub, n_ag):
    return s137.relax(atoms, n_sub, n_ag)


def grow(s, tag):
    n_sub = len(s)
    z = s.positions[:, 2]
    ztop = z.max()
    ctr = (s.cell[0][:2] + s.cell[1][:2]) / 2
    top = [i for i in range(n_sub) if z[i] > ztop - 1.2]
    near = sorted(top, key=lambda i: np.linalg.norm(s.positions[i, :2] - ctr))[:4]
    starts = [s.positions[i] + [0, 0, 2.3] for i in near]
    hol = s.positions[near].mean(axis=0)
    starts.append([hol[0], hol[1], ztop + 2.3])
    best = None
    for p in starts:
        r = relax(s + Atoms("Ag", positions=[p]), n_sub, 1)
        if r is not None and (best is None or r[1] < best[1]):
            best = r
    if best is None:
        print(f"  {tag}: every Ag1 start failed", flush=True)
        return {}
    res, cur = {}, best
    for n in range(1, 5):
        if n > 1:
            prev = cur[0]
            ag = prev.positions[n_sub:]
            h = D_AG * np.sqrt(3) / 2
            if n == 2:
                cand_p = [ag[0] + D_AG * np.array([np.cos(t), np.sin(t), 0])
                          for t in (0.0, np.pi / 2, np.pi / 4)]
            else:
                v = ag[1] - ag[0]
                v[2] = 0
                v /= np.linalg.norm(v)
                w = np.cross([0, 0, 1.0], v)
                base = ag[0] if n == 3 else ag[1]
                cand_p = [base + D_AG / 2 * v + sg * h * w for sg in (1, -1)]
            cand = None
            for p in cand_p:
                p = p.copy()
                p[2] = ag[:, 2].mean()
                r = relax(prev + Atoms("Ag", positions=[p]), n_sub, n)
                if r is not None and (cand is None or r[1] < cand[1]):
                    cand = r
            if cand is None:
                print(f"  {tag}: every Ag{n} start failed; stopping", flush=True)
                break
            cur = cand
        a, e, q = cur
        qa = q[n_sub:]
        ag = a.positions[n_sub:]
        dmin = min(float(np.min(np.linalg.norm(a.positions[:n_sub] - x, axis=1))) for x in ag)
        res[n] = {"q_total": round(float(qa.sum()), 4),
                  "q_per_atom": round(float(qa.mean()), 4),
                  "Ag_sub_min_A": round(dmin, 3)}
        write(os.path.join(OUT, f"{tag}_Ag{n}.xyz"), a)
        print(f"  {tag:<5} Ag{n}: q/atom {qa.mean():+.3f}  total {qa.sum():+.3f}  "
              f"Ag-sub {dmin:.2f} A", flush=True)
        d = json.load(open(RES)) if os.path.exists(RES) else {}
        d[tag] = res
        json.dump(d, open(RES, "w"), indent=1)
    return res


def main():
    os.makedirs(OUT, exist_ok=True)
    tags = sys.argv[1:] or ["ZnO", "TiO2", "Cu2O", "V2O5"]
    for t in tags:
        s = build(t)
        write(os.path.join(OUT, f"{t}_substrate.xyz"), s)
        print(f"{t}: {len(s)} atoms {s.get_chemical_formula()}, cell "
              + " x ".join(f"{x:.2f}" for x in s.cell.lengths()), flush=True)
        grow(s, t)


if __name__ == "__main__":
    main()
