"""HATCN derivatives: one scaffold, one electronic axis, an interior optimum.

    python scripts/135_hatcn_series.py              # build structures, write hat_series/
    python scripts/135_hatcn_series.py --harvest D  # read it back

Package F found that HATCN binds silver in the aza pocket of the
hexaazatriphenylene core -- two ring nitrogens chelating one Ag at 2.27 A,
1.631 eV -- and not on the nitriles, which are worth only 0.508. That splits the
molecule into two independent handles:

    the CORE's ring nitrogens   provide the binding pocket
    the PERIPHERAL substituents set how much charge is pulled off the silver

Those are exactly the two quantities the screening criterion needs to separate,
and they cannot be separated across unrelated materials because everything else
-- packing, surface energy, roughness, sublimation temperature -- moves with
them. Within one scaffold the film physics is held roughly fixed and only the
electronics move. That is as close to a controlled experiment as an evaporated
film gets.

Two axes:

  A (periphery, pocket kept):  HAT-H6 < HAT-F6 < HATCN in acceptor strength
  B (core, periphery kept):    triphenylene-hexacarbonitrile -- the six ring
                               nitrogens replaced by C-H, so the nitriles remain
                               and the pocket is gone

Axis B is the direct test of package F's mechanism: if the aza pocket is where
silver binds, removing it should cost most of the 1.631 eV while leaving the
nitrile tier intact. Axis A gives the criterion something to be wrong about --
E_b and q(Ag) both rise toward the nitrile, so if the good seed is the one that
binds hard and still leaves silver metallic, the best member of the series is an
INTERIOR one, not HATCN. A descriptor that is merely monotone in E_b cannot
produce an interior optimum, which is why that prediction is worth more than any
number of additional agreeing rankings.

Per molecule: the bare molecule, Ag at the pocket, and Ag2 at the pocket -- the
charge hand-back between Ag1 and Ag2 is the observable that separated HATCN from
MoOx, so it is the one worth paying for on every member.
"""
import os, re, sys, json, glob, shutil

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
gen = importlib.import_module("122_pregenerate_gjf")
r126 = importlib.import_module("126_relax_sites")
from pathgeom import STRUCT, contact, read_xyz, _bonded

OUT = os.path.join(osec := os.path.dirname(os.path.abspath(__file__)), "..", "hat_series")
RUNS = os.path.join(osec, "..", "runs")
H2EV = 27.211386
E_RE = re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)")
S2_RE = re.compile(r"S\*\*2 before annihilation\s+(\d+\.\d+)")
OPT = os.environ.get("GAUSS_OPT", "Loose,MaxCycles=80")
D_AG = 2.60
SUB_BOND = {"H": 1.086, "F": 1.340}


def parts(S, X):
    """Ring nitrogens, nitrile nitrogens, and the ring carbon each nitrile hangs
    off, read from the connectivity rather than from atom order."""
    B = _bonded(S, X)
    ring_n = [i for i, s in enumerate(S)
              if s == "N" and sum(B[i]) == 2]
    nitrile = []
    for i, s in enumerate(S):
        if s != "N" or sum(B[i]) != 1:
            continue
        c = [j for j in range(len(S)) if B[i][j]][0]
        ring_c = [j for j in range(len(S)) if B[c][j] and j != i]
        nitrile.append((i, c, ring_c[0]))
    return ring_n, nitrile


def swap_periphery(S, X, sym):
    """Replace every C#N with a single atom on the ring carbon."""
    ring_n, nitrile = parts(S, X)
    drop = set()
    keep_s, keep_x = [], []
    add = []
    for n_i, c_i, rc in nitrile:
        drop |= {n_i, c_i}
        u = X[c_i] - X[rc]
        u = u / np.linalg.norm(u)
        add.append((sym, X[rc] + SUB_BOND[sym] * u))
    for i, s in enumerate(S):
        if i in drop:
            continue
        keep_s.append(s)
        keep_x.append(X[i])
    for s, p in add:
        keep_s.append(s)
        keep_x.append(p)
    return keep_s, np.array(keep_x)


def swap_core(S, X):
    """Ring nitrogen -> C-H: the pocket removed, the nitriles left alone."""
    B = _bonded(S, X)
    ring_n, _ = parts(S, X)
    S2, X2 = list(S), [x for x in X]
    cen = X.mean(axis=0)
    add = []
    for i in ring_n:
        S2[i] = "C"
        nb = [j for j in range(len(S)) if B[i][j]]
        u = -sum((X[j] - X[i]) / np.linalg.norm(X[j] - X[i]) for j in nb)
        u = u / np.linalg.norm(u)
        add.append(("H", X[i] + 1.086 * u))
    for s, p in add:
        S2.append(s)
        X2.append(p)
    return S2, np.array(X2)


def pocket(S, X):
    """The chelating pair of ring nitrogens, and the outward normal there.

    Two ring nitrogens of the same pyrazine ring sit about 2.7 A apart and take
    one silver atom between them; that is the site package F measured at 1.631
    eV. On a molecule with no ring nitrogen left there is no pocket, and the
    fallback is the nitrile pair, which is the site worth 0.508 on HATCN.
    """
    ring_n, nitrile = parts(S, X)
    cen = X.mean(axis=0)
    nrm = np.linalg.svd(X - cen)[2][-1]
    best = None
    pool = ring_n if len(ring_n) >= 2 else [n for n, _, _ in nitrile]
    # a ring-nitrogen pair sits at 2.75 A; two nitriles on the same ring are
    # much further apart because the C#N arms point outward, and package F's
    # 0.508 eV site had Ag 2.39 A from each of a pair 4.3 A apart
    lo, hi = (2.2, 3.4) if pool is ring_n else (3.6, 4.8)
    for a in range(len(pool)):
        for b in range(a + 1, len(pool)):
            i, j = pool[a], pool[b]
            d = float(np.linalg.norm(X[i] - X[j]))
            if not lo < d < hi:
                continue
            mid = (X[i] + X[j]) / 2
            # furthest from the ring centre wins: the outer pocket, not an
            # interior pair that silver could not reach on a packed film
            r = float(np.linalg.norm(mid - cen))
            if best is None or r > best[0]:
                best = (r, mid, (i, j), d)
    if best is None:
        return None
    _, mid, pair, d = best
    lim = np.array([contact(x) for x in S])
    p = mid + 2.2 * nrm
    for _ in range(60):
        if float((np.linalg.norm(X - p, axis=1) - lim).min()) >= 0.0:
            break
        p = p + 0.1 * nrm
    return p, nrm, pair, d


def build():
    S, X = read_xyz(os.path.join(STRUCT, "HATCN.xyz"))
    mols = {"HATCN": (S, X)}
    mols["HAT-H6"] = swap_periphery(S, X, "H")
    mols["HAT-F6"] = swap_periphery(S, X, "F")
    mols["TP-CN6"] = swap_core(S, X)
    return mols


def write():
    nproc = int(os.environ.get("GAUSS_NPROC", "4"))
    mem = int(os.environ.get("GAUSS_MEM_GB", "32"))
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT, exist_ok=True)
    n = 0
    for tag, (S, X) in build().items():
        p = pocket(S, X)
        if p is None:
            print(f"  {tag:<10} no chelating pair found; skipped")
            continue
        ag, nrm, pair, dnn = p
        f = os.path.join(STRUCT, f"{tag}.xyz")
        if tag != "HATCN":
            with open(f, "w") as fh:
                fh.write(f"{len(S)}\n{tag} from HATCN by substitution, planar "
                         f"core kept\n")
                for a, c in zip(S, X):
                    fh.write(f"{a} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n")
        kind = "aza pocket" if S[pair[0]] == "N" and sum(_bonded(S, X)[pair[0]]) == 2 \
            else "nitrile pair"
        print(f"  {tag:<10} {len(S):>3} atoms  "
              + "".join(f"{e}{S.count(e)}" for e in ("C", "N", "F", "H") if S.count(e))
              + f"   site: {kind} {S[pair[0]]}{pair[0]}-{S[pair[1]]}{pair[1]} "
                f"at {dnn:.2f} A")

        # the bare molecule at this geometry: the E_b reference
        for label, extra, mult in (("mol", [], 1), ("Ag1", [ag], 2),
                                   ("Ag2", None, 1)):
            if extra is None:
                u = np.cross(nrm, [0.0, 0.0, 1.0])
                if np.linalg.norm(u) < 1e-6:
                    u = np.cross(nrm, [1.0, 0.0, 0.0])
                u /= np.linalg.norm(u)
                extra = [ag, ag + D_AG * u]
            name = f"{tag}_{label}"
            d = os.path.join(OUT, name)
            os.makedirs(d, exist_ok=True)
            opt = (f"Opt=({OPT}) " if extra else "SP ")
            lines = [f"%Chk={name}.chk", f"%NProcShared={nproc}", f"%Mem={mem}GB",
                     f"#P {gen.FUNC}/Def2SVP EmpiricalDispersion=GD3BJ {opt}"
                     f"SCF=({gen.SCF_FIRST}) NoSymm", "",
                     f"{tag} {label}, molecule frozen", "", f"0 {mult}"]
            for a, c in zip(S, X):
                code = f" {-1:>2d}" if extra else ""
                lines.append(f" {a:<2s}{code} {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
            for c in extra:
                lines.append(f" {'Ag':<2s} {0:>2d} {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
            with open(os.path.join(d, f"{name}.gjf"), "w") as fh:
                fh.write("\n".join(lines) + "\n\n")
            with open(os.path.join(d, "ORDER.txt"), "w") as fh:
                fh.write(name + "\n")
            n += 1
    d = os.path.join(OUT, "Ag_atom")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "Ag_atom.gjf"), "w") as fh:
        fh.write("\n".join([f"%NProcShared={nproc}", f"%Mem={mem}GB",
                            f"#P {gen.FUNC}/Def2SVP EmpiricalDispersion=GD3BJ SP "
                            f"SCF=({gen.SCF_FIRST}) NoSymm", "",
                            "Ag atom, doublet", "", "0 2",
                            " Ag   0.00000000   0.00000000   0.00000000"]) + "\n\n")
    with open(os.path.join(d, "ORDER.txt"), "w") as fh:
        fh.write("Ag_atom\n")
    n += 1
    with open(os.path.join(OUT, "run_campaign.py"), "w") as fh:
        fh.write(gen.RUNNER)
    with open(os.path.join(OUT, "RUN_ALL.bat"), "w", newline="\r\n") as fh:
        fh.write(gen.LAUNCHER)
    print(f"\n{n} jobs -> {os.path.relpath(OUT)}")


def energy(txt):
    hits = E_RE.findall(txt)
    if not hits:
        return None
    return float(hits[-1]) if ("Stationary point found" in txt
                               or "Normal termination" in txt) else None


def q_sum(txt, n_ag):
    b = txt.split("Mulliken charges")
    if len(b) < 2:
        return float("nan")
    q = [float(x[1]) for x in re.findall(r"\s+(\d+)\s+Ag\s+(-?\d+\.\d+)", b[-1])]
    return sum(q[-n_ag:]) if len(q) >= n_ag else float("nan")


def harvest(root):
    e_ag = None
    for c in glob.glob(os.path.join(root, "Ag_atom", "*.out")):
        e_ag = energy(open(c, errors="replace").read())
    if e_ag is None:
        print("no Ag atom energy")
        return
    print(f"{'molecule':<10}{'E_b(1)':>9}{'q(Ag1)':>9}{'E_add(2)':>10}"
          f"{'q(Ag2)':>9}{'q/atom':>8}{'hand-back':>11}")
    res = {}
    for tag in build():
        def out(lab):
            p = os.path.join(root, f"{tag}_{lab}", f"{tag}_{lab}.out")
            return open(p, errors="replace").read() if os.path.exists(p) else None
        tm, t1, t2 = out("mol"), out("Ag1"), out("Ag2")
        if tm is None or t1 is None:
            continue
        e_m, e_1 = energy(tm), energy(t1)
        if e_m is None or e_1 is None:
            continue
        eb = (e_m + e_ag - e_1) * H2EV
        q1 = q_sum(t1, 1)
        row = {"E_b1": round(eb, 4), "q1": round(float(q1), 4)}
        line = f"{tag:<10}{eb:>9.3f}{q1:>9.3f}"
        if t2 is not None and energy(t2) is not None:
            e_2 = energy(t2)
            add = (e_1 + e_ag - e_2) * H2EV
            q2 = q_sum(t2, 2)
            # the quantity that separated HATCN from MoOx: how much of the
            # adatom's charge comes back when a second atom arrives
            back = q1 - q2 / 2
            row.update({"E_add2": round(add, 4), "q2": round(float(q2), 4),
                        "q2_per_atom": round(float(q2) / 2, 4),
                        "hand_back": round(float(back), 4)})
            line += f"{add:>10.3f}{q2:>9.3f}{q2 / 2:>8.3f}{back:>+11.3f}"
        print(line)
        res[tag] = row
    if res:
        json.dump(res, open(os.path.join(RUNS, "hat_series.json"), "w"), indent=1)
        print("\nreference, measured earlier at each deepest site:")
        print("  HATCN 1.631 q=+0.572 (hand-back +0.471)   "
              "Mo3O9 2.164 q=+0.722 (hand-back +0.311)")
        print("wrote runs/hat_series.json")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--harvest":
        harvest(sys.argv[2])
    else:
        write()
