"""Does the Mo3O9 cluster's binding site represent an MoOx surface?

    python scripts/130_oxide_site_survey.py                  # write oxide_sites/
    python scripts/130_oxide_site_survey.py --harvest DIR    # read it back

Mo3O9 binds Ag at 1.806 eV, more than twice HATCN's 0.604, and that is the
single number that puts the oxide above the organic in the screening table --
against the only experiment we have, where HATCN/Ag closes at 7 nm and MoOx/Ag
is still percolated with voids at 8. Three attempts to explain it away have
failed: charge transfer, site density, and Ag-Ag cohesion competition all
favour the oxide too.

That leaves the cluster itself. Every oxygen in an isolated Mo3O9 ring is
undercoordinated relative to an amorphous MoOx network: six of the nine are
terminal Mo=O with no second metal neighbour, and Ag sits in the hollow over
three of them at once. A real surface exposes mostly two-coordinate bridging O
and saturated Mo. If Ag binds far more weakly at the bridging sites, then 1.806
eV is the cluster's worst case and not the film's, and the oxide's place in the
table is an artefact of the model rather than a fact about the material.

So: the same cluster, the same level, the same frozen substrate, Ag relaxed
from several starting sites -- the hollow it already occupies, atop a terminal
O, atop a bridging O both above the ring plane and in it, and atop Mo. Mo3O8
(one oxygen removed, a reduced Mo exposed) gets the same treatment, because a
real film is substoichiometric and that is the other half of the question.

E_b needs no new reference jobs. The molecule and the Ag atom are unchanged
between sites -- the substrate is frozen at the same geometry -- so their sum
is recovered from the E_b already measured at the hollow and that complex's
energy, and every new site is referenced to it.
"""
import os, re, sys, json, glob

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
gen = importlib.import_module("122_pregenerate_gjf")
r126 = importlib.import_module("126_relax_sites")
from pathgeom import STRUCT, contact, read_xyz

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "oxide_sites")
RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "runs")
H2EV = 27.211386
E_RE = re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)")
CLUSTERS = ["Mo3O9", "Mo3O8"]


def coordination(syms, X, i, cut=2.3):
    """How many Mo the oxygen at i is bonded to: 1 = terminal, 2 = bridging."""
    d = np.linalg.norm(X - X[i], axis=1)
    return sum(1 for k, s in enumerate(syms)
               if k != i and s == "Mo" and d[k] < cut)


def place(syms, X, j, u, scale=1.0):
    """Ag at contact distance from atom j along u, backed off until no other
    atom is inside its own contact distance."""
    u = np.asarray(u, float)
    u = u / np.linalg.norm(u)
    for step in np.arange(0.0, 2.01, 0.1):
        p = X[j] + (contact(syms[j]) * scale + step) * u
        d = np.linalg.norm(X - p, axis=1)
        if all(d[k] >= contact(syms[k]) - 0.15 for k in range(len(syms))):
            return p
    return X[j] + (contact(syms[j]) * scale + 2.0) * u


def sites(tag):
    """(label, Ag position, one line saying what the site is) for one cluster."""
    syms, X = read_xyz(os.path.join(STRUCT, f"{tag}.xyz"))
    cen = X.mean(axis=0)
    # The ring is planar; its normal is the least-spread principal direction.
    nrm = np.linalg.svd(X - cen)[2][-1]
    out = []

    mos = [k for k, s in enumerate(syms) if s == "Mo"]
    oxy = [k for k, s in enumerate(syms) if s == "O"]
    term = [k for k in oxy if coordination(syms, X, k) == 1]
    brid = [k for k in oxy if coordination(syms, X, k) >= 2]

    def radial(k):
        v = X[k] - cen
        v = v - nrm * float(v @ nrm)          # keep it in the ring plane
        return v / np.linalg.norm(v)

    if brid:
        b = brid[0]
        out.append(("bridgeO_perp", place(syms, X, b, nrm),
                    f"atop bridging O{b}, above the ring plane"))
        out.append(("bridgeO_out", place(syms, X, b, radial(b)),
                    f"atop bridging O{b}, outward in the ring plane"))
    if term:
        # The terminal O furthest from the ring plane: the exposed Mo=O tip.
        t = max(term, key=lambda k: abs(float((X[k] - cen) @ nrm)))
        m = min(mos, key=lambda k: np.linalg.norm(X[k] - X[t]))
        u = X[t] - X[m]
        out.append(("termO_top", place(syms, X, t, u),
                    f"atop terminal O{t}, along the Mo={syms[t]} axis"))
    if mos:
        m = mos[0]
        out.append(("Mo_top", place(syms, X, m, radial(m)),
                    f"atop Mo{m}, outward in the ring plane"))
    return syms, X, out


def write():
    nproc = int(os.environ.get("GAUSS_NPROC", "8"))
    mem = int(os.environ.get("GAUSS_MEM_GB", "48"))
    os.makedirs(OUT, exist_ok=True)
    n = 0
    for tag in CLUSTERS:
        syms, X, ss = sites(tag)
        d = os.path.join(OUT, tag)
        os.makedirs(d, exist_ok=True)
        order = []
        for label, pos, what in ss:
            name = f"{tag}_{label}"
            dmin = float(np.linalg.norm(X - pos, axis=1).min())
            lines = [f"%Chk={name}.chk", f"%NProcShared={nproc}",
                     f"%Mem={mem}GB",
                     f"#P {gen.FUNC}/Def2SVP EmpiricalDispersion=GD3BJ "
                     f"Opt=(MaxCycles=80) SCF=({gen.SCF_FIRST}) NoSymm", "",
                     f"{tag} Ag at {what}", "", "0 2"]
            for a, c in zip(syms, X):
                lines.append(f" {a:<2s} {-1:>2d} {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
            lines.append(f" {'Ag':<2s} {0:>2d} {pos[0]:14.8f} {pos[1]:14.8f} {pos[2]:14.8f}")
            with open(os.path.join(d, f"{name}.gjf"), "w") as f:
                f.write("\n".join(lines) + "\n\n")
            order.append(name)
            n += 1
            print(f"  {tag:<7} {label:<14} start {dmin:.2f} A from the surface"
                  f"   ({what})")
        with open(os.path.join(d, "ORDER.txt"), "w") as f:
            f.write("\n".join(order) + "\n")
    with open(os.path.join(OUT, "run_campaign.py"), "w") as f:
        f.write(gen.RUNNER)
    with open(os.path.join(OUT, "RUN_ALL.bat"), "w", newline="\r\n") as f:
        f.write(gen.LAUNCHER)
    print(f"\n{n} relaxations -> {os.path.relpath(OUT)}")
    print("Each is one Ag on a frozen cluster, 12-13 atoms: minutes, not hours.")


def reference(tag):
    """E(molecule) + E(Ag) for this cluster, from the hollow-site result.

    The substrate is frozen at one geometry across every site, so both
    reference energies are the same ones the hollow's E_b was measured with;
    their sum is E(complex) + E_b and needs no new job.
    """
    eb = json.load(open(os.path.join(RUNS, "binding_energies_pbe0.json")))[tag]
    p = os.path.join(STRUCT, f"{tag}_Ag_pbe0.xyz")
    e = float(re.search(r"E=(-?\d+\.\d+)", open(p).readlines()[1]).group(1))
    return e + eb / H2EV, e, eb


def harvest(root):
    for tag in CLUSTERS:
        ref, e_hollow, eb_hollow = reference(tag)
        syms0, X0 = read_xyz(os.path.join(STRUCT, f"{tag}.xyz"))
        print(f"\n{tag}   reference E(mol)+E(Ag) = {ref:.6f} Ha")
        print(f"{'site':<16}{'E_b (eV)':>10}{'conv':>7}{'Ag moved':>10}"
              f"{'nearest':>16}")
        print(f"{'hollow (had)':<16}{eb_hollow:>10.3f}{'yes':>7}{'-':>10}"
              f"{'-':>16}")
        rows = {}
        for p in sorted(glob.glob(os.path.join(root, tag, f"{tag}_*.out"))):
            label = os.path.basename(p)[len(tag) + 1:-4]
            txt = open(p, errors="replace").read()
            es = E_RE.findall(txt)
            g = r126.final_geometry(txt)
            if not es or g is None:
                print(f"{label:<16}  -- no energy")
                continue
            ok = "Stationary point found" in txt
            stalled = False
            if not ok and len(es) >= 8:
                tail = [float(x) for x in es[-6:]]
                steps = re.findall(r"Maximum Displacement\s+([\d.]+)", txt)
                stalled = (max(tail) - min(tail) < 2e-6 and steps
                           and all(float(x) < 1e-4 for x in steps[-5:]))
                ok = stalled
            syms, X = g
            i = syms.index("Ag")
            start = None
            for lab, pos, _ in sites(tag)[2]:
                if lab == label:
                    start = pos
            moved = float(np.linalg.norm(X[i] - start)) if start is not None else float("nan")
            others = [k for k in range(len(syms)) if k != i]
            dd = np.linalg.norm(X[others] - X[i], axis=1)
            j = others[int(dd.argmin())]
            eb = (ref - float(es[-1])) * H2EV
            rows[label] = round(eb, 4)
            conv = "yes" if ok and not stalled else ("stall" if stalled else "NO")
            print(f"{label:<16}{eb:>10.3f}{conv:>7}{moved:>9.2f} A"
                  f"{syms[j] + str(j) + ' ' + format(dd.min(), '.2f') + ' A':>16}")
        if rows:
            rows["hollow"] = round(eb_hollow, 4)
            out = os.path.join(RUNS, f"oxide_sites_{tag}.json")
            json.dump(rows, open(out, "w"), indent=1)
            print(f"wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--harvest":
        harvest(sys.argv[2])
    else:
        write()
