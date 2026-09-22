"""How uniform is the binding? One Ag relaxation per inequivalent site.

    python scripts/131_site_spread.py                 # write site_spread/
    python scripts/131_site_spread.py --harvest DIR   # read it back

The kMC study said the quantity that delays percolation is not the mean barrier
but its spread: at the same mean and the same flux, a patchy landscape made 4x
more islands and closed 21% later, because deep traps freeze adatoms as islands
of one that cannot coalesce. Package E then measured that spread on Mo3O9 by
accident -- four sites, 1.62 to 2.16 eV, a 0.54 eV range on one 12-atom cluster.

Nothing in the project has ever measured it on an organic. Every candidate got
exactly one Ag relaxation from one starting site, so the whole E_b column is a
single sample of a distribution whose width is the thing that matters. This
script takes it properly: one Ag-only relaxation above each symmetry-
inequivalent heavy atom, substrate frozen, at the campaign's level.

Two things come out of it. The obvious one is a better E_b -- the deepest site
found, rather than the first one tried. The one worth having is the width, per
candidate, which is what the kMC wants as input and what the screening table
has been silently assuming is zero.
"""
import os, re, sys, json, glob, shutil

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
gen = importlib.import_module("122_pregenerate_gjf")
r126 = importlib.import_module("126_relax_sites")
from pathgeom import (STRUCT, contact, read_xyz, relaxed_file, fingerprint,
                      outward, place, _bonded)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "site_spread")
RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "runs")
H2EV = 27.211386
E_RE = re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)")
S2_RE = re.compile(r"S\*\*2 before annihilation\s+(\d+\.\d+)")
S2_MAX = 0.90                 # a doublet is 0.75; 1.36 is a different state
MAX_SITES = 8
# Ag sits on a donor atom at its contact distance and on a carbon face well
# outside it: the Ag2 geometries put Ag-N at 2.22 A on HATCN but Ag-C at 2.67
# on benzene, against a 2.21 A contact distance for carbon. Starting half an
# Angstrom inside the wall on every carbon site is what made the first run
# crawl -- the optimiser walks out of a steep wall and then across a flat
# basin. Donor sites start where they belong; carbon faces start where the
# relaxed dimer says they end up.
START_SCALE = {"C": 1.20, "S": 1.15}
# The Ag coordinate on a physisorbed organic is soft, so the default thresholds
# ask for a position the energy cannot tell apart. On a flat basin Loose costs
# well under a meV and saves most of the steps.
OPT = os.environ.get("GAUSS_OPT", "Loose,MaxCycles=60")

# The organics the screening table ranks, plus the two extremes as controls.
TARGETS = [("HATCN", "HATCN_Ag_CN.xyz"), ("F4TCNQ", "F4TCNQ_Ag.xyz"),
           ("Bphen", "Bphen_Ag.xyz"), ("benzene", "benzene_Ag.xyz")]


def substrate(fn):
    syms, X = read_xyz(os.path.join(STRUCT, relaxed_file(fn)))
    i = syms.index("Ag")
    return [s for k, s in enumerate(syms) if k != i], np.delete(X, i, axis=0)


def rings(syms, X, nmax=6):
    """Simple cycles of up to nmax heavy atoms, as vertex sets."""
    B = _bonded(syms, X)
    heavy = [i for i, s in enumerate(syms) if s != "H"]
    adj = {i: [j for j in heavy if j != i and B[i][j]] for i in heavy}
    found = set()

    def walk(start, at, path):
        if len(path) > nmax:
            return
        for nxt in adj[at]:
            if nxt == start and len(path) >= 3:
                found.add(frozenset(path))
            elif nxt not in path and nxt > start:
                walk(start, nxt, path + [nxt])

    for i in heavy:
        walk(i, i, [i])
    return [sorted(r) for r in found]


def unique_sites(syms, X):
    """Atop, bridge and hollow sites, one per equivalence class, most open first.

    Atop alone is not a sample of the landscape -- benzene has exactly one
    inequivalent carbon, so an atop-only survey would report a spread of zero
    for a surface whose ring centre and C-C bridge are different sites
    entirely. Equivalence is Morgan relabelling over bonded neighbours, the
    same test the path builder uses; a bridge inherits the pair of labels and a
    hollow the multiset of its ring's.
    """
    fp = fingerprint(syms, X)
    cen = X.mean(axis=0)
    u, sv, vt = np.linalg.svd(X - cen)
    planar = sv[-1] / sv[0] < 0.10
    nrm = vt[-1]
    B = _bonded(syms, X)
    cand, seen = [], set()

    def add(key, anchor, h):
        if key in seen:
            return
        seen.add(key)
        cand.append((key, np.asarray(anchor, float), h))

    for i, s in enumerate(syms):
        if s != "H":
            add(("atop", fp[i]), X[i], contact(s))
    for i in range(len(syms)):
        for j in range(i + 1, len(syms)):
            if syms[i] == "H" or syms[j] == "H" or not B[i][j]:
                continue
            add(("bridge", tuple(sorted((fp[i], fp[j])))), (X[i] + X[j]) / 2,
                max(contact(syms[i]), contact(syms[j])))
    for r in rings(syms, X):
        add(("hollow", tuple(sorted(fp[k] for k in r))), X[r].mean(axis=0),
            max(contact(syms[k]) for k in r))

    out = []
    for key, anchor, h in cand:
        if planar:
            hint = nrm
        else:
            v = anchor - cen
            n = np.linalg.norm(v)
            hint = v / n if n > 1e-6 else nrm
        d = outward(syms, X, anchor, h, hint)
        pos = place(X, anchor, d, syms)
        # Ranking (and so the job names) is decided on the contact-distance
        # position, before the carbon back-off: the names have to stay the same
        # across reruns or finished jobs cannot be merged with new ones.
        clear = float((np.linalg.norm(X - pos, axis=1)
                       - [contact(x) for x in syms]).min())
        near = int(np.linalg.norm(X - pos, axis=1).argmin())
        near0 = int(np.linalg.norm(X - anchor, axis=1).argmin())
        pos = pos + d * (START_SCALE.get(syms[near0], 1.0) - 1.0) * h
        out.append((key[0], near, pos, clear))
    out.sort(key=lambda r: -r[3])
    return out[:MAX_SITES]


def write(pending=None, split=False):
    """pending: only these job names. split: one job per folder.

    The runner fills folders in parallel and runs one job at a time inside
    each, so four folders means four jobs at once whatever the core count. With
    a handful of long jobs left, one folder each is the difference between
    running them side by side and running them end to end.
    """
    nproc = int(os.environ.get("GAUSS_NPROC", "8"))
    mem = int(os.environ.get("GAUSS_MEM_GB", "48"))
    big = int(os.environ.get("GAUSS_BIG_ATOMS", "60"))
    os.makedirs(OUT, exist_ok=True)
    n = 0
    for tag, fn in TARGETS:
        if not os.path.exists(os.path.join(STRUCT, relaxed_file(fn))):
            print(f"  {tag:<9} no relaxed complex; skipped")
            continue
        syms, X = substrate(fn)
        d = os.path.join(OUT, tag)
        # A rerun with different sites must not leave the old ones behind: the
        # runner takes every .gjf in the folder, not only the ones in ORDER.
        shutil.rmtree(d, ignore_errors=True)
        if not split:
            os.makedirs(d, exist_ok=True)
        order = []
        np_ = nproc * 2 if len(syms) >= big else nproc
        mem_ = mem * 2 if len(syms) >= big else mem
        for k, (kind, i, pos, clear) in enumerate(unique_sites(syms, X)):
            # Two inequivalent sites can share a nearest atom, so the name is
            # numbered and the atom named in the title instead.
            name = f"{tag}_{k}{kind}"
            if pending is not None and name not in pending:
                continue
            dj = os.path.join(OUT, name) if split else d
            os.makedirs(dj, exist_ok=True)
            lines = [f"%Chk={name}.chk", f"%NProcShared={np_}", f"%Mem={mem_}GB",
                     f"#P {gen.FUNC}/Def2SVP EmpiricalDispersion=GD3BJ "
                     f"Opt=({OPT}) SCF=({gen.SCF_FIRST}) NoSymm", "",
                     f"{tag} Ag at the {kind} site nearest {syms[i]}{i}, substrate frozen", "", "0 2"]
            for a, c in zip(syms, X):
                lines.append(f" {a:<2s} {-1:>2d} {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
            lines.append(f" {'Ag':<2s} {0:>2d} {pos[0]:14.8f} {pos[1]:14.8f} {pos[2]:14.8f}")
            with open(os.path.join(dj, f"{name}.gjf"), "w") as f:
                f.write("\n".join(lines) + "\n\n")
            if split:
                with open(os.path.join(dj, "ORDER.txt"), "w") as f:
                    f.write(name + "\n")
            order.append(name)
            n += 1
        # The molecule and the Ag atom are already computed for this geometry;
        # only the complexes are new.
        if not split:
            with open(os.path.join(d, "ORDER.txt"), "w") as f:
                f.write("\n".join(order) + "\n")
        elif not order:
            shutil.rmtree(d, ignore_errors=True)
        print(f"  {tag:<9} {len(syms):>3} atoms  {len(order)} sites  {np_:>2} threads")
    with open(os.path.join(OUT, "run_campaign.py"), "w") as f:
        f.write(gen.RUNNER)
    with open(os.path.join(OUT, "RUN_ALL.bat"), "w", newline="\r\n") as f:
        f.write(gen.LAUNCHER)
    print(f"\n{n} relaxations -> {os.path.relpath(OUT)}")


def reference(tag, fn):
    """E(molecule) + E(Ag), from this candidate's own E_b and complex energy.

    The substrate is frozen at the geometry both were measured with, so the sum
    carries over unchanged and no reference job has to be repeated.
    """
    eb = json.load(open(os.path.join(RUNS, "binding_energies_pbe0.json")))[tag]
    p = os.path.join(STRUCT, relaxed_file(fn))
    e = float(re.search(r"E=(-?\d+\.\d+)", open(p).readlines()[1]).group(1))
    return e + eb / H2EV, eb


def harvest(root):
    table = {}
    for tag, fn in TARGETS:
        outs = sorted(glob.glob(os.path.join(root, tag, f"{tag}_*.out")))
        if not outs:
            continue
        ref, eb0 = reference(tag, fn)
        print(f"\n{tag}   stage-1 site {eb0:.3f} eV")
        print(f"{'site':<10}{'E_b (eV)':>10}{'S^2':>8}{'conv':>7}{'nearest':>14}")
        got = {}
        for p in outs:
            label = os.path.basename(p)[len(tag) + 1:-4]
            txt = open(p, errors="replace").read()
            es = E_RE.findall(txt)
            g = r126.final_geometry(txt)
            if not es or g is None:
                print(f"{label:<10}  -- no energy")
                continue
            s2 = float(S2_RE.findall(txt)[-1]) if S2_RE.findall(txt) else float("nan")
            ok = "Stationary point found" in txt
            syms, X = g
            i = syms.index("Ag")
            dd = np.linalg.norm(np.delete(X, i, axis=0) - X[i], axis=1)
            j = int(dd.argmin())
            js = [k for k in range(len(syms)) if k != i][j]
            eb = (ref - float(es[-1])) * H2EV
            bad = s2 == s2 and s2 > S2_MAX
            if ok and not bad:
                got[label] = round(eb, 4)
            print(f"{label:<10}{eb:>10.3f}{s2:>8.3f}"
                  f"{('yes' if ok else 'NO'):>7}"
                  f"{syms[js] + ' ' + format(dd.min(), '.2f') + ' A':>14}"
                  + ("   rejected: spin" if bad else ""))
        if got:
            v = list(got.values())
            table[tag] = {"sites": got, "deepest": max(v), "shallowest": min(v),
                          "spread": round(max(v) - min(v), 4),
                          "mean": round(float(np.mean(v)), 4)}
            print(f"{'':<10}deepest {max(v):.3f}   spread {max(v)-min(v):.3f} eV"
                  f"   over {len(v)} sites")
    if table:
        out = os.path.join(RUNS, "site_spread.json")
        json.dump(table, open(out, "w"), indent=1)
        print(f"\nwrote {os.path.relpath(out)}")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--harvest":
        harvest(sys.argv[2])
    else:
        a = sys.argv[1:]
        pend = None
        if "--pending" in a:
            pend = set(a[a.index("--pending") + 1].split(","))
        write(pending=pend, split="--split" in a)
