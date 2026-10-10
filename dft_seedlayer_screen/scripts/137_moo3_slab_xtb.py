"""Is the oxide's ionised silver an artefact of a 12-atom cluster?

    python scripts/137_moo3_slab_xtb.py

The single most predictable referee objection to package G: Mo3O9 is a ring with
nothing but edges, every oxygen undercoordinated, and package E already showed
the edges bind harder than anything else on it. Perhaps an extended MoO3 surface
lets silver keep its electron.

This puts Ag_n (n = 1..4) on a periodic alpha-MoO3 (010) slab -- one bilayer cut
from Kihlborg's 1963 structure (COD 9014282, the cell the project has used all
along), 3x3 surface cell, 72 substrate atoms, 20 A of vacuum -- and on the same
slab with one terminal oxygen removed, for the substoichiometric MoO(3-x) a
thermally evaporated film actually is. Substrate frozen, silver relaxed, as in
every other package.

The level is GFN2-xTB, not PBE0, because a periodic hybrid is out of reach here.
That is defensible only because it was checked first: on all sixteen cluster
geometries from package G, GFN2 charges reproduce the PBE0 ordering benzene <
HATCN ~ F4TCNQ < Mo3O9 at every n, including HATCN's hand-back (+0.31 -> +0.06)
and Mo3O9's values to within 0.03-0.09 e. The ordering is what this asks about;
the absolute numbers are GFN2's own.

Ag1 is started atop a terminal oxygen, at the two bridges and at the hollow --
the symmetry-distinct sites of the (010) face -- and the deepest kept; each
larger cluster grows from the previous one's relaxed geometry, flat at 2.6 A
in two or three orientations, deepest kept.
"""
import os, sys, json

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from tblite.ase import TBLite

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "slabs", "moo3_xtb")
CIF = os.path.join(ROOT, "structures", "MoO3_alpha_COD9014282.cif")
VAC = 20.0
D_AG = 2.60
FMAX = float(os.environ.get("SLAB_FMAX", "0.08"))


def slab(rep=3, vacancy=False):
    """One bilayer of alpha-MoO3, (010) as the surface, b along z."""
    bulk = read(CIF)
    a, b, c = bulk.cell.lengths()
    f = bulk.get_scaled_positions(wrap=True)
    f[:, 1] = (f[:, 1] + 0.5) % 1.0 - 0.5          # centre a bilayer on y = 0
    keep = np.abs(f[:, 1]) < 0.25
    sym = np.array(bulk.get_chemical_symbols())[keep]
    f = f[keep]
    # new axes: x = a, y = c, z = b (surface normal)
    pos = np.stack([f[:, 0] * a, f[:, 2] * c, f[:, 1] * b], axis=1)
    pos[:, 2] -= pos[:, 2].min()
    thick = pos[:, 2].max()
    cell = np.diag([a, c, thick + VAC])
    unit = Atoms(sym, positions=pos, cell=cell, pbc=True)
    s = unit.repeat((rep, rep, 1))
    s.positions[:, 2] += 2.0
    if vacancy:
        # remove the top terminal oxygen nearest the cell centre
        z = s.positions[:, 2]
        top = [i for i, x in enumerate(s.get_chemical_symbols())
               if x == "O" and z[i] > z.max() - 0.3]
        ctr = s.cell.diagonal()[:2] / 2
        k = min(top, key=lambda i: np.linalg.norm(s.positions[i, :2] - ctr))
        del s[k]
    return s


def calc(n_ag, etemp=1000.0):
    # An adatom over a reducible oxide sits close to a level crossing: at 300 K
    # the SCF oscillates partway through a relaxation and kills the run. Every
    # slab calculation therefore uses the same mild smearing, 1000 K, which
    # moved q(Ag) by 0.02 e at the one geometry where 300 K and 1500 K were
    # compared (0.834 -> 0.812). Same smearing for the bare slab, so energies
    # stay consistent.
    return TBLite(method="GFN2-xTB", verbosity=0, multiplicity=1 + (n_ag % 2),
                  max_iterations=500, mixer_damping=0.05,
                  electronic_temperature=etemp)


def relax(atoms, n_sub, n_ag):
    for etemp in (1000.0, 3000.0):
        a = atoms.copy()
        a.set_constraint(FixAtoms(indices=list(range(n_sub))))
        a.calc = calc(n_ag, etemp)
        try:
            opt = BFGS(a, logfile=None, maxstep=0.2)
            opt.run(fmax=FMAX, steps=250)
            return a, a.get_potential_energy(), a.get_charges()
        except Exception:
            continue
    return None


def e_ag_atom():
    a = Atoms("Ag", positions=[[0, 0, 0]], cell=np.eye(3) * 30, pbc=False)
    a.calc = TBLite(method="GFN2-xTB", verbosity=0, multiplicity=2)
    return a.get_potential_energy()


def substrate_energy(s):
    a = s.copy()
    a.calc = calc(0)
    return a.get_potential_energy()


def anchor(s):
    """The top terminal oxygen nearest the cell centre -- for the vacancy slab,
    the site the removed oxygen occupied."""
    st = slab(vacancy=False)
    z = st.positions[:, 2]
    top = [i for i, x in enumerate(st.get_chemical_symbols())
           if x == "O" and z[i] > z.max() - 0.3]
    ctr = st.cell.diagonal()[:2] / 2
    k = min(top, key=lambda i: np.linalg.norm(st.positions[i, :2] - ctr))
    return st.positions[k].copy()


def grow(s, tag):
    n_sub = len(s)
    ztop = s.positions[:, 2].max()
    e_sub = substrate_energy(s)
    e_ag = e_ag_atom()
    a0, b0 = s.cell.lengths()[0] / 3, s.cell.lengths()[1] / 3
    o = anchor(s)
    # atop the terminal O, the two bridges between neighbouring terminal Os,
    # and the hollow of four: the symmetry-distinct sites of the (010) face
    best = None
    for fx, fy in ((0, 0), (0.5, 0), (0, 0.5), (0.5, 0.5)):
        p = [o[0] + fx * a0, o[1] + fy * b0, ztop + 2.2]
        r = relax(s + Atoms("Ag", positions=[p]), n_sub, 1)
        if r is None:
            print(f"  {tag}: Ag1 start {fx},{fy} failed SCF; skipped", flush=True)
            continue
        if best is None or r[1] < best[1]:
            best = r
    res = {}
    cur = best
    for n in range(1, 5):
        if n > 1:
            prev = cur[0]
            ag = prev.positions[n_sub:]
            h = D_AG * np.sqrt(3) / 2
            starts = []
            if n == 2:
                for th in (0.0, np.pi / 2, np.pi / 4):
                    starts.append(ag[0] + D_AG * np.array([np.cos(th), np.sin(th), 0]))
            else:
                v = ag[1] - ag[0]
                v[2] = 0
                v /= np.linalg.norm(v)
                w = np.cross([0, 0, 1.0], v)
                base = ag[0] if n == 3 else ag[1]
                for sgn in (1, -1):
                    starts.append(base + D_AG / 2 * v + sgn * h * w)
            cand = None
            for p in starts:
                p = p.copy()
                p[2] = ag[:, 2].mean()
                r = relax(prev + Atoms("Ag", positions=[p]), n_sub, n)
                if r is None:
                    continue
                if cand is None or r[1] < cand[1]:
                    cand = r
            if cand is None:
                print(f"  {tag}: every Ag{n} start failed SCF; stopping", flush=True)
                break
            cur = cand
        a, e, q = cur
        qa = q[n_sub:]
        ag = a.positions[n_sub:]
        d_sub = [float(np.min(np.linalg.norm(a.positions[:n_sub] - x, axis=1)))
                 for x in ag]
        dd = sorted(float(np.linalg.norm(ag[i] - ag[j]))
                    for i in range(n) for j in range(i + 1, n))
        res[n] = {"E_bind_total_eV": round(e_sub + n * e_ag - e, 4),
                  "q_total": round(float(qa.sum()), 4),
                  "q_per_atom": round(float(qa.mean()), 4),
                  "Ag_sub_min_A": [round(x, 3) for x in d_sub],
                  "Ag_Ag_A": [round(x, 3) for x in dd[:max(n - 1, 0)]]}
        write(os.path.join(OUT, f"{tag}_Ag{n}.xyz"), a)
        print(f"  {tag:<10} Ag{n}: q/atom {qa.mean():+.3f}  total {qa.sum():+.3f}  "
              f"E_bind {res[n]['E_bind_total_eV']:.3f} eV  "
              f"Ag-sub {min(d_sub):.2f} A", flush=True)
        _save(tag, res)
    return res


def _save(tag, res):
    p = os.path.join(ROOT, "runs", "moo3_slab_xtb.json")
    d = json.load(open(p)) if os.path.exists(p) else {}
    d[tag] = res
    json.dump(d, open(p, "w"), indent=1)


def main():
    os.makedirs(OUT, exist_ok=True)
    out = {}
    for tag, vac in (("MoO3_010", False), ("MoO3-x_010", True)):
        s = slab(vacancy=vac)
        write(os.path.join(OUT, f"{tag}_substrate.xyz"), s)
        print(f"{tag}: {len(s)} atoms, cell "
              + " x ".join(f"{x:.2f}" for x in s.cell.lengths()), flush=True)
        out[tag] = grow(s, tag)
    print("wrote runs/moo3_slab_xtb.json")


if __name__ == "__main__":
    main()
