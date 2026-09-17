"""Stage 1 of the final geometry fix: relax Ag on every input complex.

    python scripts/126_relax_sites.py                 # write relax_sites/
    python scripts/126_relax_sites.py --harvest DIR   # read it back

The v14 probe found the cause of the last systematic error upstream of every
path fix: the input complexes put Ag where an earlier, different level of theory
left it, and on HATCN that is 1.5 A from where PBE0-D3/def2-SVP puts it. The
path starts there, ends at its symmetry image, and follows the wrong groove in
between. A physisorbed site forgives the offset -- Cu4I4's relaxed by 0.02 A --
but a chemisorbed one does not: 1.5 A cost HATCN 0.165 eV and inverted which
point on its path is the well.

So every complex gets one Ag-only optimisation, substrate frozen, at the level
the barriers are computed at. The four complexes refused so far for having Ag
inside a bond length are included; a relaxation is exactly what they needed.
Ag is first pushed out to a bonding distance so the optimiser does not start
inside a wall.

--harvest reads the finished outputs, writes structures/<name>_Ag_pbe0.xyz for
stage 2, and reports how far each site moved. That displacement is the size of
the error the campaign has been carrying for that candidate.
"""
import os, re, sys, glob

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
gen = importlib.import_module("122_pregenerate_gjf")
from pathgeom import CANDIDATES, STRUCT, contact, read_xyz

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "relax_sites")
SKIP = {"TPBi"}                       # 82 atoms; never in the campaign
E_RE = re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)")
Z2S = {1: "H", 3: "Li", 6: "C", 7: "N", 8: "O", 9: "F", 13: "Al", 15: "P",
       16: "S", 17: "Cl", 29: "Cu", 42: "Mo", 47: "Ag", 53: "I", 55: "Cs"}


def push_out(syms, xyz):
    """Move Ag out along the line from its nearest atom until it is at that
    atom's contact distance. Only acts when it starts closer than that."""
    i = syms.index("Ag")
    others = [k for k in range(len(syms)) if k != i]
    d = np.linalg.norm(xyz[others] - xyz[i], axis=1)
    j = others[int(d.argmin())]
    need = contact(syms[j])
    if d.min() >= need:
        return xyz, 0.0
    u = xyz[i] - xyz[j]
    u = u / np.linalg.norm(u)
    new = xyz.copy()
    new[i] = xyz[j] + need * u
    return new, float(need - d.min())


def final_geometry(txt):
    lines = txt.splitlines()
    last = None
    for k, l in enumerate(lines):
        if "orientation:" in l:
            j = k + 5
            z, x = [], []
            while j < len(lines) and not lines[j].startswith(" ----"):
                q = lines[j].split()
                z.append(int(q[1]))
                x.append([float(v) for v in q[3:6]])
                j += 1
            last = ([Z2S.get(a, str(a)) for a in z], np.array(x))
    return last


def write():
    nproc = int(os.environ.get("GAUSS_NPROC", "8"))
    mem = int(os.environ.get("GAUSS_MEM_GB", "48"))
    big = int(os.environ.get("GAUSS_BIG_ATOMS", "60"))
    os.makedirs(OUT, exist_ok=True)
    n = 0
    for tag, fn, rule, mult in CANDIDATES:
        if tag in SKIP:
            continue
        p = os.path.join(STRUCT, fn)
        if not os.path.exists(p):
            continue
        syms, xyz = read_xyz(p)
        xyz, pushed = push_out(syms, xyz)
        safe = tag.replace("=", "").replace("-", "")
        d = os.path.join(OUT, safe)
        os.makedirs(d, exist_ok=True)
        np_ = nproc * 2 if len(syms) >= big else nproc
        mem_ = mem * 2 if len(syms) >= big else mem
        # This is a first job with nothing to read, so it keeps the level shift
        # that finds the solution; the optimiser reuses orbitals between steps
        # on its own after that.
        lines = [f"%Chk={safe}_site.chk", f"%NProcShared={np_}", f"%Mem={mem_}GB",
                 f"#P {gen.FUNC}/Def2SVP EmpiricalDispersion=GD3BJ "
                 f"Opt=(MaxCycles=80) SCF=({gen.SCF_FIRST}) NoSymm", "",
                 f"{tag} Ag-only relaxation, substrate frozen"
                 + (f", Ag pushed out {pushed:.2f} A first" if pushed else ""),
                 "", f"0 {mult}"]
        for a, c in zip(syms, xyz):
            code = 0 if a == "Ag" else -1
            lines.append(f" {a:<2s} {code:>2d} {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
        with open(os.path.join(d, f"{safe}_site.gjf"), "w") as f:
            f.write("\n".join(lines) + "\n\n")
        with open(os.path.join(d, "ORDER.txt"), "w") as f:
            f.write(f"{safe}_site\n")
        n += 1
        print(f"  {tag:<12} {len(syms):>3} atoms  {np_:>2} threads"
              + (f"   Ag pushed out {pushed:.2f} A" if pushed else ""))
    with open(os.path.join(OUT, "run_campaign.py"), "w") as f:
        f.write(gen.RUNNER)
    with open(os.path.join(OUT, "RUN_ALL.bat"), "w", newline="\r\n") as f:
        f.write(gen.LAUNCHER)
    print(f"\n{n} relaxations -> {os.path.relpath(OUT)}")


def harvest(root):
    print(f"{'candidate':<12}{'steps':>6}{'conv':>6}{'E (Ha)':>16}{'Ag moved':>10}"
          f"{'nearest':>14}  written")
    for tag, fn, rule, mult in CANDIDATES:
        safe = tag.replace("=", "").replace("-", "")
        outs = glob.glob(os.path.join(root, safe, f"{safe}_site.out")) \
            + glob.glob(os.path.join(root, f"{safe}_site.out"))
        if not outs:
            continue
        txt = open(outs[0], errors="replace").read()
        ok = "Stationary point found" in txt
        es = E_RE.findall(txt)
        g = final_geometry(txt)
        if g is None or not es:
            print(f"{tag:<12}  -- no geometry")
            continue
        syms, xyz = g
        s0, x0 = read_xyz(os.path.join(STRUCT, fn))
        i = syms.index("Ag")
        moved = float(np.linalg.norm(xyz[i] - x0[s0.index("Ag")]))
        others = [k for k in range(len(syms)) if k != i]
        dd = np.linalg.norm(xyz[others] - xyz[i], axis=1)
        j = others[int(dd.argmin())]
        dest = os.path.join(STRUCT, fn.replace(".xyz", "_pbe0.xyz"))
        if ok:
            with open(dest, "w") as f:
                f.write(f"{len(syms)}\n{tag} Ag relaxed PBE0-D3/def2-SVP, "
                        f"substrate frozen, E={float(es[-1]):.8f}\n")
                for a, c in zip(syms, xyz):
                    f.write(f"{a} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n")
        print(f"{tag:<12}{len(es):>6}{('yes' if ok else 'NO'):>6}{float(es[-1]):>16.6f}"
              f"{moved:>9.2f} A   {syms[j]}{j} {dd.min():.2f} A"
              f"  {os.path.basename(dest) if ok else '-'}")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--harvest":
        harvest(sys.argv[2])
    else:
        write()
