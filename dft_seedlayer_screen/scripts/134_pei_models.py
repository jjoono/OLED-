"""PEI as a seed layer: does an amine bind silver hard and still leave it metal?

    python scripts/134_pei_models.py             # build structures, write pei_jobs/
    python scripts/134_pei_models.py --harvest D # read it back

Solution-processed polyethyleneimine is reported again and again to wet silver
well, and it is the one candidate class this project has never touched. It also
happens to be the sharpest possible test of what package G found, because an
aliphatic amine is the opposite of an oxide in exactly the way that mattered:
MoO3 is an electron acceptor that takes 0.72 e from a silver atom and never
gives it back, while an amine nitrogen is a donor. The prediction from package G
is therefore specific and falsifiable — **PEI should bind Ag strongly through
the nitrogen lone pair while leaving q(Ag) near zero or negative**, which is the
combination package G says makes a seed layer work. If instead the amine turns
out to bind weakly, or to bind hard and oxidise the silver, the criterion is
wrong and we will know it from a candidate that was not used to build it.

Three fragments, chosen to separate the chemistry rather than to imitate the
polymer:

    ethylamine        primary amine, the PEI chain end
    trimethylamine    tertiary amine, the branch point of branched PEI
    DETA              H2N-CH2CH2-NH-CH2CH2-NH2, the linear repeat: two primary
                      and one secondary nitrogen, and a chelating pocket

Geometries are built here from standard bond lengths and angles rather than
optimised. That is not a shortcut: E_b is measured with the molecule at the
geometry it has in the complex, so it appears in both terms and any strain
cancels exactly — the same construction script 127 uses. What a hand-built
conformer does change is which sites the molecule presents, so DETA is built in
the gauche conformation that lets two nitrogens reach one silver atom, which is
the geometry the wetting literature implicitly invokes.

What this cannot say: PEI is spin-coated. Nothing here argues it could be used
under a thermally evaporated top electrode on a finished OLED stack, and the
wetting reports are bottom-contact and inverted-device geometries. The value
here is as a test of the criterion, and as a pointer to what an evaporable
molecule would need to look like.
"""
import os, re, sys, json, glob, shutil

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
gen = importlib.import_module("122_pregenerate_gjf")
r126 = importlib.import_module("126_relax_sites")
from pathgeom import STRUCT, contact, read_xyz

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pei_jobs")
RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "runs")
H2EV = 27.211386
E_RE = re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)")
S2_RE = re.compile(r"S\*\*2 before annihilation\s+(\d+\.\d+)")
OPT = os.environ.get("GAUSS_OPT", "Loose,MaxCycles=80")
BOND = {("C", "C"): 1.526, ("C", "N"): 1.469, ("N", "C"): 1.469,
        ("C", "H"): 1.093, ("N", "H"): 1.014}


def nerf(a, b, c, r, theta, phi):
    """Place an atom at distance r from c, angle theta to b, dihedral phi to a."""
    th, ph = np.radians(theta), np.radians(phi)
    bc = c - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d = np.array([-r * np.cos(th), r * np.sin(th) * np.cos(ph),
                  r * np.sin(th) * np.sin(ph)])
    return c + d[0] * bc + d[1] * m + d[2] * n


def skeleton(chain, dihedrals):
    """Heavy-atom backbone from element symbols and a list of dihedrals."""
    S = list(chain)
    X = [np.zeros(3)]
    X.append(np.array([BOND[(S[0], S[1])], 0.0, 0.0]))
    if len(S) > 2:
        r = BOND[(S[1], S[2])]
        th = np.radians(111.0)
        X.append(X[1] + r * np.array([-np.cos(th), np.sin(th), 0.0]))
    for k in range(3, len(S)):
        X.append(nerf(X[k - 3], X[k - 2], X[k - 1],
                      BOND[(S[k - 1], S[k])], 111.0, dihedrals[k - 3]))
    return S, np.array(X)


def add_h(S, X, bonds):
    """Fill each heavy atom's remaining valences with hydrogen.

    The missing directions are found by letting them repel the bonds already
    made and each other, which is VSEPR done literally and lands on the
    tetrahedron. Picking each one in turn as "furthest from what is there"
    does not: on a methyl the first hydrogen goes anti to the C-C bond instead
    of 109.5 degrees off it, and the three end up 76 degrees apart.
    """
    # Nitrogen asks for four directions, not three: the lone pair occupies one
    # of them. Without it the repulsion makes a primary amine trigonal planar
    # at 120 degrees instead of pyramidal at ~107, and the lone pair -- which is
    # the whole reason an amine binds silver -- has no direction at all. The
    # surplus direction is returned rather than given a hydrogen.
    want = {"C": 4, "N": 4}
    lone = {}
    S, X = list(S), list(np.asarray(X, float))
    n0 = len(bonds)
    nb = {i: [j for j in range(n0) if bonds[i][j]] for i in range(n0)}
    rng = np.random.default_rng(0)
    for i in range(n0):
        if S[i] not in want:
            continue
        fixed = [(X[j] - X[i]) / np.linalg.norm(X[j] - X[i]) for j in nb[i]]
        k = want[S[i]] - len(nb[i])
        if k <= 0:
            continue
        u = rng.normal(size=(k, 3))
        u /= np.linalg.norm(u, axis=1)[:, None]
        for _ in range(4000):
            f = np.zeros_like(u)
            for a in range(k):
                for v in fixed:
                    d = u[a] - v
                    f[a] += d / max(np.linalg.norm(d), 1e-3) ** 3
                for b in range(k):
                    if a == b:
                        continue
                    d = u[a] - u[b]
                    f[a] += d / max(np.linalg.norm(d), 1e-3) ** 3
            u = u + 0.01 * f
            u /= np.linalg.norm(u, axis=1)[:, None]
        r = BOND[(S[i], "H")]
        nh = k - 1 if S[i] == "N" else k
        if S[i] == "N":
            lone[i] = u[-1].copy()
        for a in range(nh):
            S.append("H")
            X.append(X[i] + r * u[a])
    return S, np.array(X), lone


def _dirs(n=4000):
    g = (1 + 5 ** 0.5) / 2
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = 2 * np.pi * i / g
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi),
                     np.cos(phi)], axis=1)


def bond_table(S, X, slack=0.35):
    n = len(S)
    B = np.zeros((n, n), bool)
    for i in range(n):
        for j in range(i + 1, n):
            key = (S[i], S[j])
            if key not in BOND:
                continue
            if np.linalg.norm(X[i] - X[j]) < BOND[key] + slack:
                B[i][j] = B[j][i] = True
    return B


def build():
    """The three fragments, heavy atoms first then hydrogens."""
    mols = {}

    # primary amine, the chain end
    S, X = skeleton("CCN", [])
    B = bond_table(S, X)
    mols["ethylamine"] = add_h(S, X, B)

    # tertiary amine, the branch point: N with three methyls
    S = ["N", "C", "C", "C"]
    X = [np.zeros(3)]
    r = BOND[("N", "C")]
    for k, ang in enumerate((0.0, 120.0, 240.0)):
        t = np.radians(70.5)          # 180 - 109.5, pyramidal about N
        a = np.radians(ang)
        X.append(r * np.array([np.sin(t) * np.cos(a), np.sin(t) * np.sin(a),
                               np.cos(t)]))
    X = np.array(X)
    mols["trimethylamine"] = add_h(S, X, bond_table(S, X))

    # linear PEI repeat, gauche so two nitrogens can reach one silver
    S, X = skeleton("NCCNCCN", [60.0, 180.0, -60.0, 180.0])
    mols["DETA"] = add_h(S, X, bond_table(S, X))
    return mols


def sites(S, X, lone):
    """Along each nitrogen's lone pair, and the N...N pocket of a chelate."""
    B = bond_table(S, X)
    ns = [i for i, s in enumerate(S) if s == "N"]
    out = []
    for i in ns:
        v = lone[i] / np.linalg.norm(lone[i])
        kind = {1: "primary", 2: "secondary", 3: "tertiary"}[
            sum(1 for j in range(len(S)) if B[i][j] and S[j] != "H")]
        out.append((f"N{i}_{kind}", _push(S, X, X[i] + 2.3 * v, v)))
    for a in range(len(ns)):
        for b in range(a + 1, len(ns)):
            i, j = ns[a], ns[b]
            if not 2.5 < np.linalg.norm(X[i] - X[j]) < 4.2:
                continue
            mid = (X[i] + X[j]) / 2
            v = mid - X.mean(axis=0)
            nv = np.linalg.norm(v)
            v = v / nv if nv > 1e-6 else np.array([0.0, 0.0, 1.0])
            out.append((f"pocket_N{i}N{j}", _push(S, X, mid + 2.2 * v, v)))
    return out


def _push(S, X, p, v):
    lim = np.array([contact(x) for x in S])
    p = np.array(p, float)
    for _ in range(60):
        if float((np.linalg.norm(X - p, axis=1) - lim).min()) >= 0.0:
            return p
        p = p + 0.1 * v
    return p


def write():
    nproc = int(os.environ.get("GAUSS_NPROC", "3"))
    mem = int(os.environ.get("GAUSS_MEM_GB", "24"))
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT, exist_ok=True)
    n_job = 1
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

    for tag, (S, X, lone) in build().items():
        p = os.path.join(STRUCT, f"{tag}.xyz")
        with open(p, "w") as f:
            f.write(f"{len(S)}\n{tag} built from standard bond lengths and "
                    f"angles (PEI fragment)\n")
            for a, c in zip(S, X):
                f.write(f"{a} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n")
        print(f"  {tag:<16} {len(S):>3} atoms  -> structures/{tag}.xyz")

        # the bare molecule at exactly this geometry: the E_b reference
        dd = os.path.join(OUT, f"{tag}_mol")
        os.makedirs(dd, exist_ok=True)
        lines = [f"%NProcShared={nproc}", f"%Mem={mem}GB",
                 f"#P {gen.FUNC}/Def2SVP EmpiricalDispersion=GD3BJ SP "
                 f"SCF=({gen.SCF_FIRST}) NoSymm", "",
                 f"{tag} bare molecule, as built", "", "0 1"]
        for a, c in zip(S, X):
            lines.append(f" {a:<2s} {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
        with open(os.path.join(dd, f"{tag}_mol.gjf"), "w") as f:
            f.write("\n".join(lines) + "\n\n")
        with open(os.path.join(dd, "ORDER.txt"), "w") as f:
            f.write(f"{tag}_mol\n")
        n_job += 1

        for label, pos in sites(S, X, lone):
            name = f"{tag}_{label}"
            dd = os.path.join(OUT, name)
            os.makedirs(dd, exist_ok=True)
            lines = [f"%Chk={name}.chk", f"%NProcShared={nproc}", f"%Mem={mem}GB",
                     f"#P {gen.FUNC}/Def2SVP EmpiricalDispersion=GD3BJ "
                     f"Opt=({OPT}) SCF=({gen.SCF_FIRST}) NoSymm", "",
                     f"{tag} Ag at {label}, molecule frozen", "", "0 2"]
            for a, c in zip(S, X):
                lines.append(f" {a:<2s} {-1:>2d} {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
            lines.append(f" {'Ag':<2s} {0:>2d} {pos[0]:14.8f} {pos[1]:14.8f} {pos[2]:14.8f}")
            with open(os.path.join(dd, f"{name}.gjf"), "w") as f:
                f.write("\n".join(lines) + "\n\n")
            with open(os.path.join(dd, "ORDER.txt"), "w") as f:
                f.write(name + "\n")
            n_job += 1
            print(f"      {label:<20} Ag starts "
                  f"{float(np.linalg.norm(X - pos, axis=1).min()):.2f} A out")

    with open(os.path.join(OUT, "run_campaign.py"), "w") as f:
        f.write(gen.RUNNER)
    with open(os.path.join(OUT, "RUN_ALL.bat"), "w", newline="\r\n") as f:
        f.write(gen.LAUNCHER)
    print(f"\n{n_job} jobs -> {os.path.relpath(OUT)}")


def energy(txt):
    hits = E_RE.findall(txt)
    if not hits:
        return None
    if "Stationary point found" in txt or "Normal termination" in txt:
        return float(hits[-1])
    return None


def harvest(root):
    e_ag = None
    for c in glob.glob(os.path.join(root, "Ag_atom", "*.out")):
        e_ag = energy(open(c, errors="replace").read())
    if e_ag is None:
        print("no Ag atom energy")
        return
    res = {}
    for tag in ("ethylamine", "trimethylamine", "DETA"):
        mp = os.path.join(root, f"{tag}_mol", f"{tag}_mol.out")
        if not os.path.exists(mp):
            continue
        e_mol = energy(open(mp, errors="replace").read())
        if e_mol is None:
            continue
        print(f"\n{tag}")
        print(f"{'site':<20}{'E_b (eV)':>10}{'q(Ag)':>8}{'S^2':>8}"
              f"{'nearest':>14}")
        rows = {}
        for p in sorted(glob.glob(os.path.join(root, f"{tag}_*", f"{tag}_*.out"))):
            label = os.path.basename(p)[len(tag) + 1:-4]
            if label == "mol":
                continue
            txt = open(p, errors="replace").read()
            e = energy(txt)
            g = r126.final_geometry(txt)
            if e is None or g is None:
                print(f"{label:<20}  -- not converged")
                continue
            S, X = g
            i = S.index("Ag")
            dd = np.linalg.norm(np.delete(X, i, axis=0) - X[i], axis=1)
            j = [k for k in range(len(S)) if k != i][int(dd.argmin())]
            s2 = float(S2_RE.findall(txt)[-1]) if S2_RE.findall(txt) else float("nan")
            b = txt.split("Mulliken charges")
            q = float("nan")
            if len(b) > 1:
                m = re.findall(r"\s+(\d+)\s+Ag\s+(-?\d+\.\d+)", b[-1])
                if m:
                    q = float(m[-1][1])
            eb = (e_mol + e_ag - e) * H2EV
            rows[label] = {"E_b_eV": round(eb, 4), "q_Ag": round(q, 4),
                           "S2": round(s2, 4), "nearest": S[j],
                           "d": round(float(dd.min()), 3)}
            print(f"{label:<20}{eb:>10.3f}{q:>8.3f}{s2:>8.3f}"
                  f"{S[j] + ' ' + format(dd.min(), '.2f') + ' A':>14}")
        if rows:
            best = max(rows.items(), key=lambda kv: kv[1]["E_b_eV"])
            print(f"{'':<20}deepest {best[0]} at {best[1]['E_b_eV']:.3f} eV, "
                  f"q(Ag) = {best[1]['q_Ag']:+.3f}")
            res[tag] = rows
    if res:
        out = os.path.join(RUNS, "pei_binding.json")
        json.dump(res, open(out, "w"), indent=1)
        print(f"\nwrote {os.path.relpath(out)}")
        print("\nfor comparison, at each substrate's deepest site:")
        print("  MoOx  2.164 eV  q=+0.722    HATCN 1.631 eV  q=+0.572")
        print("  F4TCNQ 1.208 eV q=+0.545    benzene 0.205 eV q=-0.061")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--harvest":
        harvest(sys.argv[2])
    else:
        write()
