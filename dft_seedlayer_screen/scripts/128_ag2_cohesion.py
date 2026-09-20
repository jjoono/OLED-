"""Does the seed layer compete with Ag-Ag bonding, or only pin single atoms?

    python scripts/128_ag2_cohesion.py                 # write ag2_jobs/
    python scripts/128_ag2_cohesion.py --harvest D S    # D = ag2_jobs, S = relax_sites

E_b alone cannot explain the experiment. MoOx binds a silver atom 1.2 eV harder
than HATCN does -- exp(2dE/3kT) is 3e13 in Venables -- while HATCN's site
density is only twice as high. Something outside that picture decides it.

The seed layer's job is not to hold one atom. It is to stop atoms finding each
other: a film grows smooth when adatoms stay dispersed long enough to nucleate
densely, and Volmer-Weber islands form because Ag-Ag cohesion (2.95 eV/atom in
the bulk, 1.66 eV in the free dimer) beats anything an organic offers. So the
quantity that matters is how much the substrate WEAKENS the Ag-Ag bond:

    dE_dimer(sub) = E(sub+Ag2) + E(sub) - 2 E(sub+Ag)
    quality       = dE_dimer(sub) - dE_dimer(vacuum)

dE_dimer is the energy of bringing two separately adsorbed atoms together. In
vacuum it is -1.66 eV and clustering always wins. A substrate that holds each
atom in its own pocket has to pay some of that back, and the amount it claws
back is what this measures. A substrate can have a huge E_b and still claw back
nothing -- if it pins Ag as an ion at an isolated site, the pinned atom is out
of the metal film altogether, which is one way an oxide can bind hard and still
leave voids.

Five candidates: HATCN and F4TCNQ against Mo3O9 and Mo3O8, with benzene as the
floor. Substrate frozen, both Ag free, started from the stage-1 site with the
second atom at the free-dimer bond length along the surface.
"""
import os, re, sys, glob

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
gen = importlib.import_module("122_pregenerate_gjf")
from pathgeom import CANDIDATES, STRUCT, contact, read_xyz, relaxed_file

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ag2_jobs")
WANT = ["HATCN", "F4TCNQ", "Mo3O9", "Mo3O8", "benzene"]
D_AG2 = 2.53           # A, PBE0-D3 free Ag2 bond length, refined by the vacuum job
E_RE = re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)")
H2EV = 27.211386


def second_site(syms, xyz, i_ag):
    """Where to put the second atom: beside the first, along the surface, as far
    from the substrate as the first one is."""
    ag = xyz[i_ag]
    sub = np.delete(xyz, i_ag, axis=0)
    s_sub = [s for k, s in enumerate(syms) if k != i_ag]
    d = np.linalg.norm(sub - ag, axis=1)
    n = ag - sub[int(d.argmin())]
    n = n / np.linalg.norm(n)
    # in-plane directions, pick the one that keeps the new atom furthest from
    # the substrate -- beside the first atom, not on top of it or under it
    lim = np.array([contact(x) for x in s_sub])
    best, bestclear = None, -1e9
    for th in np.linspace(0, 2 * np.pi, 72, endpoint=False):
        e1 = np.cross(n, [1.0, 0.0, 0.0])
        if np.linalg.norm(e1) < 1e-6:
            e1 = np.cross(n, [0.0, 1.0, 0.0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(n, e1)
        u = np.cos(th) * e1 + np.sin(th) * e2
        p = ag + D_AG2 * u
        clear = float((np.linalg.norm(sub - p, axis=1) / lim).min())
        if clear > bestclear:
            best, bestclear = p, clear
    return best


def write():
    nproc = int(os.environ.get("GAUSS_NPROC", "8"))
    mem = int(os.environ.get("GAUSS_MEM_GB", "48"))
    os.makedirs(OUT, exist_ok=True)

    d = os.path.join(OUT, "Ag2_vacuum")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "Ag2_vacuum.gjf"), "w") as f:
        f.write(f"%Chk=Ag2_vacuum.chk\n%NProcShared={nproc}\n%Mem={mem}GB\n"
                f"#P {gen.FUNC}/Def2SVP EmpiricalDispersion=GD3BJ "
                f"Opt SCF=({gen.SCF_FIRST}) NoSymm\n\n"
                f"Ag2 free dimer, singlet\n\n0 1\n"
                f" Ag  0.000000  0.000000  0.000000\n"
                f" Ag  0.000000  0.000000  {D_AG2:.6f}\n\n")
    n = 1

    for tag, fn, rule, mult in CANDIDATES:
        if tag not in WANT:
            continue
        fn = relaxed_file(fn)
        p = os.path.join(STRUCT, fn)
        if not os.path.exists(p):
            continue
        syms, xyz = read_xyz(p)
        i = syms.index("Ag")
        p2 = second_site(syms, xyz, i)
        safe = tag.replace("=", "").replace("-", "")
        dd = os.path.join(OUT, safe)
        os.makedirs(dd, exist_ok=True)
        # Two silvers is an even electron count on a closed-shell molecule, so
        # the complex is a singlet -- the atom was a doublet, the dimer is not.
        lines = [f"%Chk={safe}_ag2.chk", f"%NProcShared={nproc}", f"%Mem={mem}GB",
                 f"#P {gen.FUNC}/Def2SVP EmpiricalDispersion=GD3BJ "
                 f"Opt=(MaxCycles=80) SCF=({gen.SCF_FIRST}) NoSymm", "",
                 f"{tag} + Ag2, substrate frozen, both Ag free", "", "0 1"]
        for a, c in zip(syms, xyz):
            lines.append(f" {a:<2s} {0 if a == 'Ag' else -1:>2d} "
                         f"{c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
        lines.append(f" Ag  0 {p2[0]:14.8f} {p2[1]:14.8f} {p2[2]:14.8f}")
        with open(os.path.join(dd, f"{safe}_ag2.gjf"), "w") as f:
            f.write("\n".join(lines) + "\n\n")
        n += 1
        print(f"  {tag:<10} {len(syms)+1:>3} atoms   2nd Ag at "
              f"{np.linalg.norm(p2-xyz[i]):.2f} A from the first")

    with open(os.path.join(OUT, "run_campaign.py"), "w") as f:
        f.write(gen.RUNNER)
    with open(os.path.join(OUT, "RUN_ALL.bat"), "w", newline="\r\n") as f:
        f.write(gen.LAUNCHER)
    print(f"\n{n} jobs -> {os.path.relpath(OUT)}")


def energy(p):
    t = open(p, errors="replace").read()
    h = E_RE.findall(t)
    if not h:
        return None
    if "Normal termination" in t or "Stationary point found" in t:
        return float(h[-1])
    st = re.findall(r"Maximum Displacement\s+([\d.]+)", t)
    if len(h) >= 8 and max(h[-6:], key=float) and st \
            and all(float(x) < 1e-4 for x in st[-5:]):
        return float(h[-1])
    return None


def harvest(root, relax):
    import json
    ag1 = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "..", "runs", "binding_energies_pbe0.json")))
    e_ag = None
    for c in glob.glob(os.path.join(relax, "..", "**", "Ag_atom.out"),
                       recursive=True) + glob.glob("eb_jobs/Ag_atom/Ag_atom.out"):
        e_ag = energy(c)
    d2 = energy(os.path.join(root, "Ag2_vacuum", "Ag2_vacuum.out"))
    if d2 is None or e_ag is None:
        print("need both the Ag atom and the free Ag2")
        return
    vac = (d2 - 2 * e_ag) * H2EV
    print(f"free Ag2 bond = {vac:+.3f} eV\n")
    print(f"{'candidate':<10}{'E_b(1 Ag)':>11}{'dE_dimer':>11}{'clawed back':>13}")
    for tag in WANT:
        safe = tag.replace("=", "").replace("-", "")
        e2 = energy(os.path.join(root, safe, f"{safe}_ag2.out"))
        e1 = energy(os.path.join(relax, safe, f"{safe}_site.out"))
        mol = energy(os.path.join("eb_jobs", safe, f"{safe}_mol.out"))
        if None in (e2, e1, mol):
            print(f"{tag:<10} incomplete")
            continue
        dd = (e2 + mol - 2 * e1) * H2EV
        print(f"{tag:<10}{ag1.get(tag, float('nan')):>11.3f}{dd:>+11.3f}"
              f"{dd - vac:>+13.3f}")


if __name__ == "__main__":
    if len(sys.argv) > 3 and sys.argv[1] == "--harvest":
        harvest(sys.argv[2], sys.argv[3])
    else:
        write()
