"""Gaussian package for Ag on n-doped Bphen (Bphen(-)Cs(+)), and its harvester.

    python scripts/154_bphencs_gaussian.py            # write bphencs_jobs/
    python scripts/154_bphencs_gaussian.py --harvest  # read the returned .out/.log

Same level and runner as the earlier packages: PBE1PBE/Def2SVP GD3BJ
(Gaussian folds D3 into "SCF Done"), one folder per chain, run_campaign.py.
The whole molecule is used; the PySCF run of the same system had not
finished one SCF in 2 h on the 4-core container.

Folders
  ref      BphenCs (doublet) and the Ag atom, single points at the frozen geometry
  site_Cs  Ag next to Cs        -- singlet (UKS, Guess=Mix: broken symmetry allowed)
  site_N   Ag at the N side        and triplet; Ag position optimised, molecule
  site_C   Ag on the pi face       frozen (-1 flags), as for every earlier E_b
  ag2      Ag2 at the Cs site (doublet), both Ag optimised: charge with two atoms

E_b = E(BphenCs) + E(Ag) - min(E_singlet, E_triplet); the D3 part is read
from the "Dispersion energy" lines. Sites are the deepest of each kind from
the GFN2 survey (runs/bphen_cs_seed.json, script 153). The undoped baseline is
the existing Gaussian value for Bphen (chelate, 0.621 eV, same level).
"""
import importlib
import json
import os
import re
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pathgeom import read_xyz  # noqa: E402

ROOT = os.path.join(HERE, "..")
PKG = os.path.join(ROOT, "bphencs_jobs")
SURVEY = os.path.join(ROOT, "runs", "bphen_cs_seed.json")
XYZ = os.path.join(ROOT, "structures", "common", "BphenCs_opt.xyz")
OUT = os.path.join(ROOT, "runs", "bphencs_gaussian.json")
EV = 27.211386

HEAD = "%Chk={name}.chk\n%NProcShared={nproc}\n%Mem=24GB\n"
NPROC = 3
SCF = "SCF=(XQC,MaxCycle=128,VShift=300) NoSymm Pop=Hirshfeld"
R_SP = "#P PBE1PBE/Def2SVP EmpiricalDispersion=GD3BJ SP " + SCF
R_OPT = "#P UPBE1PBE/Def2SVP EmpiricalDispersion=GD3BJ Opt=(Loose,MaxCycles=80) " + SCF
R_OPT_BS = R_OPT + " Guess=Mix"

# v2 (2026-10-29): the v1 site jobs optimised Ag at PBE0 with 3 cores each and
# were still running hours after the 52-min molecule single point. As in
# script 151: locate Ag at PBE-D3 with density fitting (no exact exchange),
# then score at PBE0 in a --Link1-- step from the same checkpoint; 8 cores
# per job (the runner fills 32 physical cores with 4 such jobs).
V2 = "--v2" in sys.argv
SCF2 = "SCF=(XQC,MaxCycle=128,VShift=300) NoSymm"
R_LOC = "#P UPBEPBE/Def2SVP/W06 EmpiricalDispersion=GD3BJ Opt=(Loose,MaxCycles=60) " + SCF2
R_LOC_BS = R_LOC + " Guess=Mix"
R_SCORE = ("#P UPBE1PBE/Def2SVP EmpiricalDispersion=GD3BJ SP Geom=Check Guess=Read "
           + SCF2 + " Pop=Hirshfeld")


def gjf(folder, name, route, title, charge, mult, S, X, free=()):
    os.makedirs(folder, exist_ok=True)
    lines = [HEAD.format(name=name, nproc=NPROC) + route, "", title, "", f"{charge} {mult}"]
    for i, (s, x) in enumerate(zip(S, X)):
        flag = "" if free is None else (" 0" if i in free else " -1")
        lines.append(f" {s:<2}{flag}  {x[0]:13.8f} {x[1]:13.8f} {x[2]:13.8f}")
    lines += [""]
    if V2 and route.startswith("#P UPBEPBE"):
        # --Link1--: PBE0 single point at the PBE-D3 geometry, from its orbitals
        lines += ["--Link1--", HEAD.format(name=name, nproc=NPROC) + R_SCORE, "",
                  title + " -- PBE0 at the PBE-D3 geometry", "", f"{charge} {mult}", ""]
    lines += [""]
    with open(os.path.join(folder, name + ".gjf"), "w", newline="\r\n") as f:
        f.write("\n".join(lines))
    with open(os.path.join(folder, "ORDER.txt"), "a", newline="\r\n") as f:
        f.write(name + "\n")


def site_picks(d):
    by_el = {}
    for lab, v in sorted(d["sites"].items(), key=lambda kv: -kv[1]["E_b"]):
        el = "".join(c for c in v["near"] if c.isalpha())
        by_el.setdefault(el, lab)
    return by_el


def ag2_start(S, X, d, deep):
    """GFN2 Ag2 at the deepest site (as in 153 stage 3), positions returned."""
    s153 = importlib.import_module("153_bphen_cs_seed")
    n = len(S)
    ag1 = np.array(d["sites"][deep]["pos"])
    dd = np.linalg.norm(X - ag1, axis=1)
    nv = ag1 - X[int(dd.argmin())]
    nv /= np.linalg.norm(nv)
    t1 = np.cross(nv, [0, 0, 1.0])
    t1 /= np.linalg.norm(t1)
    t2 = np.cross(nv, t1)
    best = None
    for th in (0, 2 * np.pi / 3, 4 * np.pi / 3):
        p2 = ag1 + 2.6 * (np.cos(th) * t1 + np.sin(th) * t2)
        r = s153.relax(S + ["Ag", "Ag"], np.vstack([X, ag1, p2]), n, 2)
        if r and (best is None or r[0] < best[0]):
            best = r
    return best[2][n:]


def write():
    d = json.load(open(SURVEY))
    S, X = read_xyz(XYZ)
    X = np.asarray(X, float)
    n = len(S)
    global NPROC
    pkg = PKG + "_v2" if V2 else PKG
    if os.path.exists(pkg):
        shutil.rmtree(pkg)
    os.makedirs(pkg)
    if V2:
        NPROC = 8
        r_s, r_t = R_LOC_BS, R_LOC
    else:
        r_s, r_t = R_OPT_BS, R_OPT
        for f in ("RUN_ALL.bat", "run_campaign.py"):
            shutil.copy(os.path.join(ROOT, "cluster_growth", f), pkg)
        ref = os.path.join(pkg, "ref")
        gjf(ref, "BphenCs_mol", R_SP, "Bphen(-)Cs(+), doublet, GFN2 geometry", 0, 2, S, X, None)
        gjf(ref, "Ag_atom", R_SP, "Ag atom, doublet", 0, 2, ["Ag"], [[0, 0, 0]], None)
    picks = site_picks(d)
    for el, lab in picks.items():
        folder = os.path.join(pkg, f"site_{el}")
        ag = np.array(d["sites"][lab]["pos"])
        S2, X2 = S + ["Ag"], np.vstack([X, ag])
        gjf(folder, f"{el}_singlet", r_s, f"Ag at {lab}, singlet (broken symmetry allowed),"
            " molecule frozen", 0, 1, S2, X2, {n})
        gjf(folder, f"{el}_triplet", r_t, f"Ag at {lab}, triplet, molecule frozen",
            0, 3, S2, X2, {n})
    deep = sorted(d["sites"], key=lambda k: -d["sites"][k]["E_b"])[0]
    P = ag2_start(S, X, d, deep)
    gjf(os.path.join(pkg, "ag2"), "Ag2_Cs", r_t, f"Ag2 at {deep}, doublet, molecule frozen",
        0, 2, S + ["Ag", "Ag"], np.vstack([X, P]), {n, n + 1})
    if V2:
        with open(os.path.join(pkg, "README_v2.txt"), "w", newline="\r\n") as f:
            f.write(
                "v2: replaces site_C, site_Cs, site_N and ag2 of bphencs_jobs (ref is kept).\n"
                "Ag is located at PBE-D3 (density fitting, fast) and scored at PBE0 in a\n"
                "--Link1-- step; 8 cores per job.\n\n"
                "1. Close the running window and end any g16 process (Task Manager).\n"
                "2. In bphencs_jobs, delete site_C, site_Cs, site_N and ag2.\n"
                "3. Copy the four folders of this zip into bphencs_jobs.\n"
                "4. Run RUN_ALL.bat again: ref is already done and is skipped.\n")
        print("wrote", pkg, "sites", picks)
        return
    with open(os.path.join(pkg, "README.txt"), "w", newline="\r\n") as f:
        f.write(
            "Ag on n-doped Bphen (Bphen(-)Cs(+)) -- PBE1PBE/Def2SVP GD3BJ\n\n"
            "Run: double-click RUN_ALL.bat (or: python run_campaign.py --workers 3).\n"
            "Finished jobs are skipped on a restart. Return the whole folder.\n\n"
            "ref/      molecule and Ag atom single points (fast)\n"
            "site_*/   Ag position optimised, molecule frozen; singlet and triplet\n"
            "ag2/      two Ag atoms at the deepest (Cs) site\n\n"
            f"Sites: {json.dumps(picks)}  (deepest of each kind, GFN2 survey)\n")
    print("wrote", pkg, "sites", picks)


E_RE = re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)")
D_RE = re.compile(r"Dispersion energy=\s*(-?\d+\.\d+)\s*Hartrees")
S2_RE = re.compile(r"S\*\*2 before annihilation\s+(\d+\.\d+)")


def read_out(folder, stem):
    for ext in (".out", ".log"):
        p = os.path.join(folder, stem + ext)
        if os.path.exists(p):
            t = open(p, errors="ignore").read()
            break
    else:
        return None
    r = {"normal": "Normal termination" in t[-4000:]}
    e = E_RE.findall(t)
    if e:
        r["E"] = float(e[-1])
    dsp = D_RE.findall(t)
    if dsp:
        r["D3"] = float(dsp[-1])
    s2 = S2_RE.findall(t)
    if s2:
        r["S2"] = float(s2[-1])
    blocks = t.split("Mulliken charges")
    if len(blocks) > 1:
        q = []
        for line in blocks[-1].splitlines()[2:]:
            p = line.split()
            if len(p) < 3 or not p[0].isdigit():
                break
            q.append((p[1], float(p[2])))
        r["q_Ag"] = [v for s, v in q if s == "Ag"]
        r["q_Cs"] = [v for s, v in q if s == "Cs"]
    return r


def harvest():
    res = {}
    ref = {k: read_out(os.path.join(PKG, "ref"), k) for k in ("BphenCs_mol", "Ag_atom")}
    res["ref"] = ref
    if not all(ref.values()) or "E" not in ref["BphenCs_mol"] or "E" not in ref["Ag_atom"]:
        print("ref not finished")
        json.dump(res, open(OUT, "w"), indent=1)
        return
    for folder in sorted(os.listdir(PKG)):
        if not folder.startswith("site_"):
            continue
        el = folder[5:]
        s = read_out(os.path.join(PKG, folder), f"{el}_singlet")
        t = read_out(os.path.join(PKG, folder), f"{el}_triplet")
        row = {"singlet": s, "triplet": t}
        done = [x for x in (s, t) if x and x.get("normal") and "E" in x]
        if done:
            best = min(done, key=lambda x: x["E"])
            eb = (ref["BphenCs_mol"]["E"] + ref["Ag_atom"]["E"] - best["E"]) * EV
            disp = (best.get("D3", 0) - ref["BphenCs_mol"].get("D3", 0)) * EV
            row.update(E_b=round(eb, 4), E_disp=round(-disp, 4), E_elec=round(eb + disp, 4),
                       state="singlet" if best is s else "triplet",
                       q_Ag=best.get("q_Ag"), q_Cs=best.get("q_Cs"),
                       complete=len(done) == 2)
        res[folder] = row
        print(f"{folder:<8} E_b {row.get('E_b')}  state {row.get('state')}  q_Ag {row.get('q_Ag')}")
    a2 = read_out(os.path.join(PKG, "ag2"), "Ag2_Cs")
    res["ag2"] = a2
    if a2 and a2.get("q_Ag"):
        print("Ag2 total charge", round(sum(a2["q_Ag"]), 3), "q_Cs", a2.get("q_Cs"))
    json.dump(res, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    harvest() if "--harvest" in sys.argv else write()
