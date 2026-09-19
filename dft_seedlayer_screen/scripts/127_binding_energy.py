"""E_b at the same level and geometry as E_d.

    python scripts/127_binding_energy.py                # write eb_jobs/
    python scripts/127_binding_energy.py --harvest DIR  # read it back

    E_b = E(molecule) + E(Ag) - E(complex)

The complex energies already exist: stage 1 relaxed Ag on all 25 candidates at
PBE0-D3/def2-SVP, which is the level the barriers are computed at. What is
missing is the two pieces it is measured against. The molecule is taken at the
geometry it has IN the complex -- the substrate was frozen through the
relaxation, so that geometry is the input structure unchanged, and no relaxation
energy of the molecule leaks into E_b.

That matters because the E_b table in the audit came from the PBE/DIIS route
this project disqualified: the same geometry reached SCF solutions 0.51 eV
apart under different convergence paths. These 26 single points replace it with
numbers on the same footing as E_d, which is what pairing them in a Venables
estimate requires.

Ag is a doublet; every molecule here is a closed-shell singlet, so the molecule
jobs are RKS and cheap. The dispersion correction is kept on all three terms --
it is a large part of the binding for the physisorbed candidates, and dropping
it from one term only would be worse than dropping it from all.
"""
import os, re, sys, glob

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
gen = importlib.import_module("122_pregenerate_gjf")
from pathgeom import CANDIDATES, STRUCT, read_xyz, relaxed_file

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "eb_jobs")
SKIP = {"TPBi"}
H2EV = 27.211386
E_RE = re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)")


def sp(name, syms, xyz, mult, nproc, mem, chk=None):
    route = (f"#P {gen.FUNC}/Def2SVP EmpiricalDispersion=GD3BJ SP "
             f"SCF=({gen.SCF_FIRST}) NoSymm")
    out = ([f"%Chk={chk}"] if chk else []) + \
          [f"%NProcShared={nproc}", f"%Mem={mem}GB", route, "", name, "",
           f"0 {mult}"]
    for a, c in zip(syms, xyz):
        out.append(f" {a:<2s} {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
    return "\n".join(out) + "\n\n"


def write():
    nproc = int(os.environ.get("GAUSS_NPROC", "8"))
    mem = int(os.environ.get("GAUSS_MEM_GB", "48"))
    big = int(os.environ.get("GAUSS_BIG_ATOMS", "60"))
    os.makedirs(OUT, exist_ok=True)

    d = os.path.join(OUT, "Ag_atom")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "Ag_atom.gjf"), "w") as f:
        f.write(sp("Ag atom, doublet", ["Ag"], [[0.0, 0.0, 0.0]], 2, nproc, mem))
    n = 1

    for tag, fn, rule, mult in CANDIDATES:
        if tag in SKIP:
            continue
        fn = relaxed_file(fn)
        p = os.path.join(STRUCT, fn)
        if not os.path.exists(p):
            continue
        syms, xyz = read_xyz(p)
        if "Ag" not in syms:
            continue
        i = syms.index("Ag")
        # The molecule exactly as it sits in the complex: the substrate never
        # moved, so this is the complex minus one atom and nothing else.
        m_s = [s for k, s in enumerate(syms) if k != i]
        m_x = np.delete(xyz, i, axis=0)
        safe = tag.replace("=", "").replace("-", "")
        dd = os.path.join(OUT, safe)
        os.makedirs(dd, exist_ok=True)
        np_ = nproc * 2 if len(m_s) >= big else nproc
        mem_ = mem * 2 if len(m_s) >= big else mem
        # A molecule that binds Ag as a doublet complex is a singlet on its own;
        # anything else here would be a different electronic state, not a
        # reference.
        with open(os.path.join(dd, f"{safe}_mol.gjf"), "w") as f:
            f.write(sp(f"{tag} bare molecule at its geometry in the complex",
                       m_s, m_x, 1, np_, mem_))
        n += 1
        print(f"  {tag:<12} {len(m_s):>3} atoms  {np_:>2} threads")

    with open(os.path.join(OUT, "run_campaign.py"), "w") as f:
        f.write(gen.RUNNER)
    with open(os.path.join(OUT, "RUN_ALL.bat"), "w", newline="\r\n") as f:
        f.write(gen.LAUNCHER)
    print(f"\n{n} single points -> {os.path.relpath(OUT)}")


def energy(p):
    """Last converged SCF energy, or None.

    A relaxation that ran out of optimisation steps still terminates with an
    error, and DMABN's did -- after sitting perfectly still for 74 steps with
    the force on Ag 14% over a threshold it was never going to cross. Its
    geometry and energy are converged; only Gaussian's bookkeeping is not. The
    same stalled-is-converged test the stage-1 harvest uses applies here, or
    DMABN silently drops out of the E_b column.
    """
    txt = open(p, errors="replace").read()
    hits = E_RE.findall(txt)
    if not hits:
        return None
    if "Normal termination" in txt or "Stationary point found" in txt:
        return float(hits[-1])
    tail = [float(x) for x in hits[-6:]]
    steps = re.findall(r"Maximum Displacement\s+([\d.]+)", txt)
    if len(hits) >= 8 and max(tail) - min(tail) < 2e-6 \
            and steps and all(float(x) < 1e-4 for x in steps[-5:]):
        return float(hits[-1])
    return None


def harvest(root, relax_root):
    ag = None
    for c in (os.path.join(root, "Ag_atom", "Ag_atom.out"),
              os.path.join(root, "Ag_atom.out")):
        if os.path.exists(c):
            ag = energy(c)
    if ag is None:
        print("no Ag atom energy; cannot compute anything")
        return
    print(f"E(Ag) = {ag:.6f} Ha\n")
    print(f"{'candidate':<12}{'E_b (eV)':>10}{'complex':>16}{'molecule':>16}")
    rows = []
    for tag, fn, rule, mult in CANDIDATES:
        safe = tag.replace("=", "").replace("-", "")
        mol = None
        for c in (os.path.join(root, safe, f"{safe}_mol.out"),
                  os.path.join(root, f"{safe}_mol.out")):
            if os.path.exists(c):
                mol = energy(c)
        cx = None
        for c in (os.path.join(relax_root, safe, f"{safe}_site.out"),
                  os.path.join(relax_root, f"{safe}_site.out")):
            if os.path.exists(c):
                cx = energy(c)
        if mol is None or cx is None:
            continue
        eb = (mol + ag - cx) * H2EV
        rows.append((tag, eb, cx, mol))
    for tag, eb, cx, mol in sorted(rows, key=lambda r: -r[1]):
        print(f"{tag:<12}{eb:>10.3f}{cx:>16.6f}{mol:>16.6f}")
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                       "runs", "binding_energies_pbe0.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    import json
    json.dump({t: round(e, 4) for t, e, _, _ in rows}, open(out, "w"), indent=1)
    print(f"\nwrote {os.path.relpath(out)}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--harvest":
        harvest(sys.argv[2], sys.argv[3])
    else:
        write()
