"""Does the next silver atom join the cluster, or take a fresh site?

    python scripts/133_cluster_growth.py                  # write cluster_growth/
    python scripts/133_cluster_growth.py --harvest DIR    # read it back

Everything measured so far is a one-atom quantity, and one atom cannot tell a
nucleus from a trapped ion. The kMC on the measured landscapes (script 132)
closed that route: a featureless surface at the oxide's own terrain barrier
reproduces the oxide's closure thickness to the digit, so the barrier landscape
is not what separates HATCN from MoOx. What is left is the charge state -- Ag(0)
at q=+0.07 on most of HATCN against Ag(I) at q=+0.70 on every site of the oxide
-- and the question that turns it into a number is whether a silver atom that
has given up an electron can still start a metal island.

Ag_n for n = 1..4 on the deepest site of each substrate, substrate frozen, all
silver free. Two numbers come out of it per step:

    E_add(n)  = E(Ag_(n-1)/sub) + E(Ag) - E(Ag_n/sub)
                what the nth atom gains by joining the cluster
    E_b(1)    what that same atom would gain by taking an empty deep site

E_add(n) > E_b(1) means growth wins and the deposit coarsens into islands.
E_add(n) < E_b(1) means an atom is better off alone on the substrate, which is a
deposit that stays dispersed -- many nuclei that never become a film. This is
the Volmer-Weber criterion written with the two energies that are actually
measurable here, and it is the first quantity in the project that can rank
substrates by whether they nucleate metal rather than by how hard they pull.

Beside it, the charge on the cluster as it grows. If Ag_4 on the oxide is still
carrying most of a positive charge per atom, it is not a metal nucleus, and the
oxide's 2.16 eV of binding is buying an ionic adlayer instead of a seed.

Clusters start flat on the surface -- dimer, triangle, rhombus, all at the free
Ag2 bond length, laid in the plane parallel to the substrate. Whether they stay
flat is a result: package D found the organics stand a dimer up on one atom,
which is a three-dimensional island seed, while the oxides lay it down with both
atoms bound and the bond stretched.
"""
import os, re, sys, json, glob, shutil

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
gen = importlib.import_module("122_pregenerate_gjf")
r126 = importlib.import_module("126_relax_sites")
from pathgeom import STRUCT, contact, read_xyz

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cluster_growth")
RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "runs")
H2EV = 27.211386
E_RE = re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)")
S2_RE = re.compile(r"S\*\*2 before annihilation\s+(\d+\.\d+)")
D_AG = 2.60                # A, on-surface Ag-Ag from package D (2.57-2.61)
NMAX = 4
OPT = os.environ.get("GAUSS_OPT", "Loose,MaxCycles=80")
# (tag, the deepest site's E_b in eV) -- from packages E and F
TARGETS = [("HATCN", 1.631), ("Mo3O9", 2.164), ("F4TCNQ", 1.208),
           ("benzene", 0.205)]


def basis(sub_s, sub_x, ag):
    """Outward normal at the adatom, and two lateral directions."""
    d = np.linalg.norm(sub_x - ag, axis=1)
    n = ag - sub_x[int(d.argmin())]
    n = n / np.linalg.norm(n)
    u = np.cross(n, [0.0, 0.0, 1.0])
    if np.linalg.norm(u) < 1e-6:
        u = np.cross(n, [1.0, 0.0, 0.0])
    u /= np.linalg.norm(u)
    return n, u, np.cross(n, u)


def lift(sub_s, sub_x, p, n):
    """Back a silver atom off along n until nothing is inside contact."""
    lim = np.array([contact(x) for x in sub_s])
    p = np.array(p, float)
    for _ in range(60):
        if float((np.linalg.norm(sub_x - p, axis=1) - lim).min()) >= 0.0:
            return p
        p = p + 0.1 * n
    return p


def cluster(sub_s, sub_x, ag, k):
    """k silver atoms: the first on its site, the rest flat beside it.

    Dimer, then an equilateral triangle, then a rhombus -- the gas-phase shapes
    of Ag2, Ag3 and Ag4 -- laid in the plane parallel to the surface so that the
    relaxation, not the input, decides whether the cluster stands up.
    """
    n, u, v = basis(sub_s, sub_x, ag)
    h = D_AG * np.sqrt(3) / 2
    offs = [np.zeros(3), D_AG * u, D_AG / 2 * u + h * v, D_AG * 1.5 * u + h * v]
    # the lateral directions are arbitrary in sign; take the pair that keeps the
    # cluster furthest from the substrate
    best, bestclear = None, -1e9
    for su in (1, -1):
        for sv in (1, -1):
            offs = [su * D_AG * u,
                    su * D_AG / 2 * u + sv * h * v,
                    su * D_AG * 1.5 * u + sv * h * v][:k - 1]
            pos = [ag] + [lift(sub_s, sub_x, ag + o, n) for o in offs]
            clear = min(float((np.linalg.norm(sub_x - q, axis=1)
                               - [contact(x) for x in sub_s]).min())
                        for q in pos)
            if clear > bestclear:
                best, bestclear = pos, clear
    return best


def write():
    nproc = int(os.environ.get("GAUSS_NPROC", "3"))
    mem = int(os.environ.get("GAUSS_MEM_GB", "24"))
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT, exist_ok=True)
    # E_add needs the free Ag energy on its own, and no earlier package's
    # outputs are kept in the repository, so it ships with this one. It is a
    # one-atom single point: seconds.
    d = os.path.join(OUT, "Ag_atom")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "Ag_atom.gjf"), "w") as f:
        f.write("\n".join([f"%NProcShared={nproc}", f"%Mem={mem}GB",
                           f"#P {gen.FUNC}/Def2SVP EmpiricalDispersion=GD3BJ SP "
                           f"SCF=({gen.SCF_FIRST}) NoSymm", "",
                           "Ag atom, doublet", "", "0 2",
                           " Ag   0.00000000   0.00000000   0.00000000"]) + "\n\n")
    with open(os.path.join(d, "ORDER.txt"), "w") as f:
        f.write("Ag_atom\n")
    n_job = 1
    for tag, eb1 in TARGETS:
        p = os.path.join(STRUCT, f"{tag}_Ag1_deep.xyz")
        if not os.path.exists(p):
            print(f"  {tag:<9} no {os.path.basename(p)}; skipped")
            continue
        syms, X = read_xyz(p)
        i = syms.index("Ag")
        sub_s = [s for k, s in enumerate(syms) if k != i]
        sub_x = np.delete(X, i, axis=0)
        for k in range(2, NMAX + 1):
            pos = cluster(sub_s, sub_x, X[i], k)
            name = f"{tag}_Ag{k}"
            d = os.path.join(OUT, name)
            os.makedirs(d, exist_ok=True)
            # n silver atoms, n electrons of spin choice: odd n is a doublet,
            # even n a singlet, which is what the free clusters are.
            mult = 1 if k % 2 == 0 else 2
            lines = [f"%Chk={name}.chk", f"%NProcShared={nproc}", f"%Mem={mem}GB",
                     f"#P {gen.FUNC}/Def2SVP EmpiricalDispersion=GD3BJ "
                     f"Opt=({OPT}) SCF=({gen.SCF_FIRST}) NoSymm", "",
                     f"{tag} Ag{k} on the deepest site, substrate frozen",
                     "", f"0 {mult}"]
            for a, c in zip(sub_s, sub_x):
                lines.append(f" {a:<2s} {-1:>2d} {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
            for c in pos:
                lines.append(f" {'Ag':<2s} {0:>2d} {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
            with open(os.path.join(d, f"{name}.gjf"), "w") as f:
                f.write("\n".join(lines) + "\n\n")
            with open(os.path.join(d, "ORDER.txt"), "w") as f:
                f.write(name + "\n")
            n_job += 1
            dd = min(float(np.linalg.norm(sub_x - c, axis=1).min()) for c in pos)
            print(f"  {tag:<9} Ag{k}  mult {mult}  {len(sub_s) + k:>3} atoms  "
                  f"closest Ag-substrate {dd:.2f} A")
    with open(os.path.join(OUT, "run_campaign.py"), "w") as f:
        f.write(gen.RUNNER)
    with open(os.path.join(OUT, "RUN_ALL.bat"), "w", newline="\r\n") as f:
        f.write(gen.LAUNCHER)
    print(f"\n{n_job} relaxations -> {os.path.relpath(OUT)}")


def energy(txt):
    hits = E_RE.findall(txt)
    if not hits:
        return None
    if "Stationary point found" in txt or "Normal termination" in txt:
        return float(hits[-1])
    tail = [float(x) for x in hits[-6:]]
    steps = re.findall(r"Maximum Displacement\s+([\d.]+)", txt)
    if len(hits) >= 8 and max(tail) - min(tail) < 2e-6 and steps \
            and all(float(x) < 1e-4 for x in steps[-5:]):
        return float(hits[-1])
    return None


def q_ag(txt, n_ag):
    """Summed Mulliken charge on the silver, from the last population block."""
    b = txt.split("Mulliken charges")
    if len(b) < 2:
        return float("nan")
    q = [float(x[1]) for x in re.findall(r"\s+(\d+)\s+Ag\s+(-?\d+\.\d+)", b[-1])]
    return sum(q[-n_ag:]) if len(q) >= n_ag else float("nan")


def harvest(root):
    e_ag = None
    for c in glob.glob(os.path.join(root, "Ag_atom", "*.out")) \
            + glob.glob(os.path.join(root, "Ag_atom.out")):
        e_ag = energy(open(c, errors="replace").read())
    res = {}
    for tag, eb1 in TARGETS:
        p = os.path.join(STRUCT, f"{tag}_Ag1_deep.xyz")
        if not os.path.exists(p):
            continue
        head = open(p).readlines()[1]
        e1 = float(re.search(r"E=(-?\d+\.\d+)", head).group(1))
        # E(Ag) is recovered from this substrate's own E_b(1) and E(Ag1/sub),
        # together with E(substrate) -- the same trick the earlier packages use,
        # except here only the SUM E(sub)+E(Ag) is needed and it cancels out of
        # every incremental step but the first.
        ens = {1: e1}
        qs = {}
        for k in range(2, NMAX + 1):
            f = os.path.join(root, f"{tag}_Ag{k}", f"{tag}_Ag{k}.out")
            if not os.path.exists(f):
                continue
            txt = open(f, errors="replace").read()
            e = energy(txt)
            if e is None:
                continue
            ens[k] = e
            qs[k] = q_ag(txt, k)
        if len(ens) < 2 or e_ag is None:
            print(f"{tag}: need the Ag atom energy and at least Ag2")
            continue
        print(f"\n{tag}   deepest site E_b(1) = {eb1:.3f} eV")
        print(f"{'n':>3}{'E_add(n) eV':>13}{'vs E_b(1)':>11}{'q(Ag_n)':>9}"
              f"{'q/atom':>8}{'Ag-Ag (A)':>11}")
        print(f"{1:>3}{eb1:>13.3f}{'--':>11}")
        row = {"E_b1": eb1, "steps": {}}
        for k in range(2, NMAX + 1):
            if k not in ens or (k - 1) not in ens:
                continue
            add = (ens[k - 1] + e_ag - ens[k]) * H2EV
            f = os.path.join(root, f"{tag}_Ag{k}", f"{tag}_Ag{k}.out")
            S, X = r126.final_geometry(open(f, errors="replace").read())
            ag = X[[j for j, s in enumerate(S) if s == "Ag"]]
            dd = [np.linalg.norm(ag[a] - ag[b])
                  for a in range(len(ag)) for b in range(a + 1, len(ag))]
            near = sorted(dd)[:k - 1]
            verdict = "cluster" if add > eb1 else "fresh site"
            print(f"{k:>3}{add:>13.3f}{verdict:>11}{qs[k]:>9.3f}"
                  f"{qs[k] / k:>8.3f}"
                  + ("  " + ", ".join(f"{x:.2f}" for x in near)).rjust(11))
            row["steps"][k] = {"E_add_eV": round(add, 4),
                               "q_total": round(float(qs[k]), 4),
                               "q_per_atom": round(float(qs[k]) / k, 4),
                               "prefers": verdict,
                               "Ag_Ag": [round(float(x), 3) for x in near]}
        res[tag] = row
    if res:
        out = os.path.join(RUNS, "cluster_growth.json")
        json.dump(res, open(out, "w"), indent=1)
        print(f"\nwrote {os.path.relpath(out)}")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--harvest":
        harvest(sys.argv[2])
    else:
        write()
