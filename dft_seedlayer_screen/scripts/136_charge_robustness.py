"""Is the charge hand-back a property of the silver, or of the Mulliken scheme?

    python scripts/136_charge_robustness.py                 # write charge_jobs/
    python scripts/136_charge_robustness.py --harvest DIR   # read it back

Package G left exactly one observable standing: the charge on a growing silver
cluster. HATCN's adatom carries +0.572 and the pair hands most of it back
(+0.10 to +0.16 per atom from n = 2 on), while Mo3O9 keeps every cluster it was
measured at ionised (+0.25 to +0.41 per atom, +1.22 in total at n = 4). The
incremental energy E_add was dropped because singlet/doublet alternation is
larger than any difference between substrates; the Ag-Ag distance separates the
oxide only at n = 2. So the whole argument now rests on charges -- and Mulliken
charges in def2-SVP are the least trustworthy number in the project as absolute
values, because they divide overlap density by a rule that depends on the basis.

This package asks the same question three other ways, on the geometries already
converged, so nothing is re-optimised:

    Hirshfeld and CM5   partition the density in real space, not the basis
    NBO natural charge  plus the Ag 5s natural occupancy, which is the direct
                        statement of "metallic": Ag(0) keeps ~1 electron in 5s,
                        Ag(I) has given it away
    def2-TZVP           the same single points in a larger basis, for HATCN and
                        Mo3O9 at n = 1, 2 -- if the hand-back is a basis-set
                        artefact it changes here

The claim survives if every scheme orders the substrates the same way at every
n and the hand-back on HATCN appears in all of them. It does not need the
absolute numbers to agree; no two charge schemes ever do.
"""
import os, re, sys, json, glob, shutil

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
gen = importlib.import_module("122_pregenerate_gjf")
r126 = importlib.import_module("126_relax_sites")
from pathgeom import STRUCT, read_xyz

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "charge_jobs")
RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "runs")
CLU = os.path.join(STRUCT, "clusters")
TAGS = ["HATCN", "Mo3O9", "F4TCNQ", "benzene"]
TZVP = {("HATCN", 1), ("HATCN", 2), ("Mo3O9", 1), ("Mo3O9", 2)}


def save_geometries(root):
    """Copy the converged Ag_n geometries into the repository, once.

    The package G outputs live only on the workstation and in a scratch copy;
    the charge jobs, and anyone checking them, need the geometries durable.
    """
    os.makedirs(CLU, exist_ok=True)
    for tag in TAGS:
        S, X = read_xyz(os.path.join(STRUCT, f"{tag}_Ag1_deep.xyz"))
        _write(os.path.join(CLU, f"{tag}_Ag1.xyz"), S, X, f"{tag} Ag1, package F/E")
        for n in (2, 3, 4):
            p = os.path.join(root, f"{tag}_Ag{n}", f"{tag}_Ag{n}.out")
            if not os.path.exists(p):
                continue
            g = r126.final_geometry(open(p, errors="replace").read())
            if g is not None:
                _write(os.path.join(CLU, f"{tag}_Ag{n}.xyz"), g[0], g[1],
                       f"{tag} Ag{n}, package G, PBE0-D3/def2-SVP")


def _write(p, S, X, title):
    with open(p, "w") as f:
        f.write(f"{len(S)}\n{title}\n")
        for a, c in zip(S, X):
            f.write(f"{a} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n")


def write():
    nproc = int(os.environ.get("GAUSS_NPROC", "4"))
    mem = int(os.environ.get("GAUSS_MEM_GB", "32"))
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT, exist_ok=True)
    jobs = 0
    for tag in TAGS:
        for n in (1, 2, 3, 4):
            p = os.path.join(CLU, f"{tag}_Ag{n}.xyz")
            if not os.path.exists(p):
                continue
            S, X = read_xyz(p)
            mult = 2 if n % 2 else 1
            for basis in ("Def2SVP", "Def2TZVP"):
                if basis == "Def2TZVP" and (tag, n) not in TZVP:
                    continue
                name = f"{tag}_Ag{n}_{'svp' if basis == 'Def2SVP' else 'tzvp'}"
                d = os.path.join(OUT, name)
                os.makedirs(d, exist_ok=True)
                lines = [f"%NProcShared={nproc}", f"%Mem={mem}GB",
                         f"#P {gen.FUNC}/{basis} EmpiricalDispersion=GD3BJ SP "
                         f"SCF=({gen.SCF_FIRST}) Pop=(NBO,Hirshfeld) NoSymm", "",
                         f"{tag} Ag{n} charges, {basis}, package G geometry", "",
                         f"0 {mult}"]
                for a, c in zip(S, X):
                    lines.append(f" {a:<2s} {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
                with open(os.path.join(d, f"{name}.gjf"), "w") as f:
                    f.write("\n".join(lines) + "\n\n")
                with open(os.path.join(d, "ORDER.txt"), "w") as f:
                    f.write(name + "\n")
                jobs += 1
    with open(os.path.join(OUT, "run_campaign.py"), "w") as f:
        f.write(gen.RUNNER)
    with open(os.path.join(OUT, "RUN_ALL.bat"), "w", newline="\r\n") as f:
        f.write(gen.LAUNCHER)
    print(f"{jobs} single points -> {os.path.relpath(OUT)}")


def _block(txt, head, stop):
    i = txt.rfind(head)
    if i < 0:
        return ""
    j = txt.find(stop, i + len(head))
    return txt[i:j if j > 0 else None]


def charges(txt, n_ag):
    """Summed Ag charge under each scheme, and the mean Ag 5s occupation."""
    out = {}
    m = _block(txt, "Mulliken charges", "Sum of Mulliken")
    q = [float(x) for x in re.findall(r"\d+\s+Ag\s+(-?\d+\.\d+)", m)]
    out["mulliken"] = sum(q[-n_ag:]) if len(q) >= n_ag else None
    h = _block(txt, "Hirshfeld charges, spin densities", "Tot ")
    rows = re.findall(r"\d+\s+Ag\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
                      r"\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", h)
    if len(rows) >= n_ag:
        out["hirshfeld"] = sum(float(r[0]) for r in rows[-n_ag:])
        out["cm5"] = sum(float(r[5]) for r in rows[-n_ag:])
    nb = _block(txt, "Summary of Natural Population Analysis", "* Total *")
    q = [float(x) for x in re.findall(r"Ag\s+\d+\s+(-?\d+\.\d+)", nb)]
    out["nbo"] = sum(q[-n_ag:]) if len(q) >= n_ag else None
    conf = re.findall(r"Ag\s+\d+\s+\[core\]5S\(\s*(\d+\.\d+)\)", txt)
    out["ag5s"] = float(np.mean([float(c) for c in conf[-n_ag:]])) \
        if len(conf) >= n_ag else None
    return out


def harvest(root):
    res = {}
    schemes = ("mulliken", "hirshfeld", "cm5", "nbo")
    print(f"{'job':<20}" + "".join(f"{s:>11}" for s in schemes) + f"{'Ag 5s':>8}")
    for tag in TAGS:
        for n in (1, 2, 3, 4):
            for b in ("svp", "tzvp"):
                p = os.path.join(root, f"{tag}_Ag{n}_{b}", f"{tag}_Ag{n}_{b}.out")
                if not os.path.exists(p):
                    continue
                txt = open(p, errors="replace").read()
                if "Normal termination" not in txt:
                    print(f"{tag}_Ag{n}_{b:<10} not finished")
                    continue
                c = charges(txt, n)
                res[f"{tag}_Ag{n}_{b}"] = c
                print(f"{tag + '_Ag' + str(n) + '_' + b:<20}"
                      + "".join(f"{c[s] / n:>+11.3f}" if c.get(s) is not None
                                else f"{'--':>11}" for s in schemes)
                      + (f"{c['ag5s']:>8.3f}" if c.get("ag5s") is not None else f"{'--':>8}"))
    print("(charges are per silver atom; Ag 5s is the mean natural occupation)")
    if res:
        json.dump(res, open(os.path.join(RUNS, "charge_robustness.json"), "w"), indent=1)
        # the test the claim has to pass: same substrate order in every scheme
        for s in schemes:
            for n in (1, 2, 3, 4):
                v = {t: res.get(f"{t}_Ag{n}_svp", {}).get(s) for t in TAGS}
                v = {t: x / n for t, x in v.items() if x is not None}
                if len(v) == len(TAGS):
                    print(f"  {s:<10} n={n}: " + " < ".join(
                        t for t, _ in sorted(v.items(), key=lambda kv: kv[1])))
        print("wrote runs/charge_robustness.json")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--harvest":
        harvest(sys.argv[2])
    else:
        save_geometries(os.environ.get(
            "PKG_G_DIR", "/tmp/claude-0/-home-user-OLED-/"
            "965acd1b-0fc7-5d29-a851-865963011861/scratchpad/cg"))
        write()
