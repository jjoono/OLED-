"""Beyond-HATCN screen: charge-transport-compatible materials with potentially
stronger Ag anchoring than HATCN (E_b = 1.03 eV, PBE-D3BJ/def2-SVP CP).

Candidates (cluster/site models, same protocol as 03/05/06/10):
  - CuSCN  (p-type transparent HTL; (CuSCN)3 ring; Ag@S and Ag@N sites)
  - CuI    (p-type HTL; Cu4I4 cubane; Ag@I site)
  - Ph3PO  (model for phosphine-oxide ETL hosts DPEPO/TSPO1; Ag@O=P)
  - BTD    (2,1,3-benzothiadiazole, n-type acceptor unit; Ag@N)
  - thiophene (PEDOT/thiophene HTL model; Ag@S)
  - s-triazine (triazine ETL core; Ag@N)

Pipeline: build (RDKit or hand-coded) -> GFN2-xTB opt (geometry only)
-> DFT rigid scan of Ag along Ag-anchor axis (xTB overbinds Ag)
-> PBE-D3BJ/def2-SVP UKS CP binding energy at scan minimum.
"""
import numpy as np, os, subprocess, shutil, re, json, sys
from rdkit import Chem
from rdkit.Chem import AllChem

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STR, RUNS = os.path.join(BASE, "structures"), os.path.join(BASE, "runs")
os.makedirs(STR, exist_ok=True); os.makedirs(RUNS, exist_ok=True)

def write_xyz(path, syms, xyz, cm=""):
    with open(path, "w") as f:
        f.write(f"{len(syms)}\n{cm}\n")
        for s, c in zip(syms, xyz):
            f.write(f"{s} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n")

def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0]); syms, xyz = [], []
    for l in lines[2:2+n]:
        p = l.split(); syms.append(p[0]); xyz.append(list(map(float, p[1:4])))
    return syms, np.array(xyz)

def run_xtb(xyz_path, tag, uhf=0, chrg=0):
    wd = os.path.join(RUNS, tag); os.makedirs(wd, exist_ok=True)
    if os.path.exists(os.path.join(wd, "xtbopt.xyz")):
        print(f"[xtb] {tag}: cached", flush=True); return True
    shutil.copy(xyz_path, os.path.join(wd, "in.xyz"))
    with open(os.path.join(wd, "xtb.log"), "w") as log:
        subprocess.run(["xtb", "in.xyz", "--gfn", "2", "--chrg", str(chrg),
                        "--uhf", str(uhf), "--opt", "tight"],
                       cwd=wd, stdout=log, stderr=subprocess.STDOUT)
    ok = os.path.exists(os.path.join(wd, "xtbopt.xyz"))
    print(f"[xtb] {tag}: {'ok' if ok else 'FAILED'}", flush=True)
    return ok

def rdkit_mol(smi, seed=7):
    mol = Chem.MolFromSmiles(smi); mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    try: AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    except Exception: pass
    conf = mol.GetConformer()
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                     conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
    return syms, xyz

# ---------------- 1. Build substrates ----------------
def build_all():
    # organics via RDKit
    for name, smi in {
        "thiophene": "c1ccsc1",
        "BTD": "c1ccc2nsnc2c1",              # 2,1,3-benzothiadiazole
        "triazine": "c1ncncn1",
        "Ph3PO": "O=P(c1ccccc1)(c1ccccc1)c1ccccc1",
        "Me3PO": "O=P(C)(C)C",
    }.items():
        p = os.path.join(STR, f"{name}.xyz")
        if not os.path.exists(p.replace(".xyz", "_DONExtb")):
            syms, xyz = rdkit_mol(smi)
            write_xyz(p, syms, xyz, name)
            run_xtb(p, name)

    # (CuSCN)3 ring, hand-built planar guess -> xtb relax
    if True:
        units, R = 3, 3.4
        syms, xyz = [], []
        for k in range(units):
            th = 2*np.pi*k/units
            c, s = np.cos(th), np.sin(th)
            rot = np.array([[c, -s], [s, c]])
            # local fragment (x,y): Cu(0,0) S(2.25,0.6) C(3.4,1.5) N(4.2,2.3)->next Cu
            frag = {"Cu": (0.0, 0.0), "S": (2.25, 0.8), "C": (3.3, 1.9), "N": (3.9, 3.0)}
            for sym, (fx, fy) in frag.items():
                v = rot @ np.array([fx + R, fy])
                syms.append(sym); xyz.append([v[0], v[1], 0.0])
        write_xyz(os.path.join(STR, "CuSCN3.xyz"), syms, np.array(xyz), "CuSCN trimer ring guess")
        run_xtb(os.path.join(STR, "CuSCN3.xyz"), "CuSCN3")

    # Cu4I4 cubane
    a = 2.62
    cu = [(0,0,0),(a,a,0),(a,0,a),(0,a,a)]
    ii = [(a,0,0),(0,a,0),(0,0,a),(a,a,a)]
    syms = ["Cu"]*4 + ["I"]*4
    xyz = np.array(cu + ii, float)
    write_xyz(os.path.join(STR, "Cu4I4.xyz"), syms, xyz, "Cu4I4 cubane guess")
    run_xtb(os.path.join(STR, "Cu4I4.xyz"), "Cu4I4")

# ---------------- 2. Place Ag & xtb opt complex ----------------
def unit(v): return v / np.linalg.norm(v)

def place_ag(tag):
    """Return list of (complex_tag, syms, xyz_with_Ag, anchor_index)."""
    syms, x = read_xyz(os.path.join(RUNS, tag, "xtbopt.xyz"))
    cen = x.mean(axis=0)
    out = []
    def plane_normal(pts):
        u, s, vt = np.linalg.svd(pts - pts.mean(axis=0)); return vt[2]
    if tag == "thiophene":
        si = syms.index("S")
        n = plane_normal(x)
        d = unit(unit(x[si]-cen) + 0.9*n)
        out.append((f"{tag}_Ag", syms+["Ag"], np.vstack([x, x[si]+2.55*d]), si))
    elif tag == "BTD":
        ni = [i for i,s in enumerate(syms) if s == "N"][0]
        d = unit(x[ni]-cen)
        out.append((f"{tag}_Ag", syms+["Ag"], np.vstack([x, x[ni]+2.30*d]), ni))
    elif tag == "triazine":
        ni = [i for i,s in enumerate(syms) if s == "N"][0]
        d = unit(x[ni]-cen)
        out.append((f"{tag}_Ag", syms+["Ag"], np.vstack([x, x[ni]+2.30*d]), ni))
    elif tag in ("Ph3PO", "Me3PO"):
        oi = syms.index("O"); pi = syms.index("P")
        d = unit(x[oi]-x[pi])
        out.append((f"{tag}_Ag", syms+["Ag"], np.vstack([x, x[oi]+2.20*d]), oi))
    elif tag == "CuSCN3":
        n = plane_normal(x)
        si = [i for i,s in enumerate(syms) if s == "S"][0]
        ni = [i for i,s in enumerate(syms) if s == "N"][0]
        out.append((f"{tag}_AgS", syms+["Ag"], np.vstack([x, x[si]+2.55*unit(0.6*unit(x[si]-cen)+n)]), si))
        out.append((f"{tag}_AgN", syms+["Ag"], np.vstack([x, x[ni]+2.30*unit(0.6*unit(x[ni]-cen)+n)]), ni))
    elif tag == "Cu4I4":
        idx = [i for i,s in enumerate(syms) if s == "I"]
        # pick iodine farthest from centroid, go outward
        di = max(idx, key=lambda i: np.linalg.norm(x[i]-cen))
        d = unit(x[di]-cen)
        out.append((f"{tag}_Ag", syms+["Ag"], np.vstack([x, x[di]+2.85*d]), di))
    return out

# ---------------- 3. psi4 CP with DFT distance refinement ----------------
def psi4_stage(jobs):
    import psi4
    psi4.set_memory("6 GB")
    psi4.set_num_threads(os.cpu_count() or 4)
    psi4.core.set_output_file(os.path.join(RUNS, "psi4_beyond_hatcn.out"), False)
    base_opts = {"basis": "def2-svp", "scf_type": "df", "reference": "uks",
                 "maxiter": 150, "guess": "sad"}
    psi4.set_options(base_opts)
    H2EV = 27.211386
    M = "pbe-d3bj"

    def gstr(syms, xyz, ghost=None, mult=1):
        s = f"0 {mult}\n"
        for i, (sym, c) in enumerate(zip(syms, xyz)):
            t = f"Gh({sym})" if ghost and i in ghost else sym
            s += f"{t} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n"
        return s + "symmetry c1\nno_reorient\nno_com\n"

    def robust_energy(geom_str):
        for attempt, extra in enumerate([{}, {"level_shift": 1.0, "level_shift_cutoff": 0.01,
                                              "damping_percentage": 15},
                                         {"soscf": True, "soscf_max_iter": 40}]):
            try:
                psi4.set_options({**base_opts, **extra})
                e = psi4.energy(M, molecule=psi4.geometry(geom_str))
                psi4.set_options(base_opts)
                return e
            except Exception as ex:
                print(f"    scf attempt {attempt} failed: {type(ex).__name__}", flush=True)
                psi4.core.clean()
        return None

    res_path = os.path.join(RUNS, "beyond_hatcn_binding_eV.json")
    res = json.load(open(res_path)) if os.path.exists(res_path) else {}

    for ctag, anchor in jobs:
        if ctag in res and isinstance(res[ctag], dict) and "Eb_eV" in res[ctag]:
            print(f"[skip] {ctag} done: {res[ctag]['Eb_eV']:.3f} eV", flush=True); continue
        syms, x = read_xyz(os.path.join(RUNS, ctag, "xtbopt.xyz"))
        agi = len(syms) - 1
        # anchor = nearest heavy atom to Ag after xtb opt (Ag may have migrated)
        heavy = [i for i in range(agi) if syms[i] != "H"]
        anchor = min(heavy, key=lambda i: np.linalg.norm(x[agi] - x[i]))
        print(f"  [{ctag}] anchor after xtb: {syms[anchor]}{anchor}", flush=True)
        ax = unit(x[agi] - x[anchor])
        r0 = np.linalg.norm(x[agi] - x[anchor])
        # rigid scan: complex energy vs Ag displacement along anchor axis
        best = (None, None)
        scan = {}
        dts = [-0.2, -0.1, 0.0, 0.1, 0.25, 0.45, 0.7]
        if len(syms) > 25: dts = [-0.1, 0.0, 0.15, 0.35]   # big systems: fewer points
        for dt in dts:
            xs = x.copy(); xs[agi] = x[anchor] + (r0 + dt) * ax
            e = robust_energy(gstr(syms, xs, None, 2))
            if e is None: continue
            scan[round(r0+dt, 3)] = e
            if best[1] is None or e < best[1]: best = (dt, e)
            print(f"  [{ctag}] r={r0+dt:.2f} E={e:.6f}", flush=True)
        if best[1] is None:
            res[ctag] = {"error": "scan_failed"}; json.dump(res, open(res_path, "w"), indent=2); continue
        # parabolic refine around best
        dts = sorted(scan.keys()); rbest = r0 + best[0]
        xs = x.copy(); xs[agi] = x[anchor] + rbest * ax
        write_xyz(os.path.join(STR, f"{ctag}_refined.xyz"), syms, xs, f"{ctag} DFT-scan min")
        e_cx = scan[round(rbest, 3)]
        e_sub = robust_energy(gstr(syms, xs, {agi}, 1))
        e_ag = robust_energy(gstr(syms, xs, set(range(agi)), 2))
        if None in (e_sub, e_ag):
            res[ctag] = {"error": "cp_failed"}; json.dump(res, open(res_path, "w"), indent=2); continue
        eb = (e_sub + e_ag - e_cx) * H2EV
        res[ctag] = {"Eb_eV": eb, "r_A": rbest,
                     "scan": {str(k): v for k, v in scan.items()}}
        print(f"[RESULT] {ctag}: E_b(CP) = {eb:.3f} eV at r = {rbest:.2f} A", flush=True)
        json.dump(res, open(res_path, "w"), indent=2)
        psi4.core.clean()

SKIP_XTB = {"CuSCN3_AgS", "CuSCN3_AgN", "Cu4I4_Ag"}   # GFN2 diverges on Cu+Ag open shell;
# use frozen substrate + hand-placed Ag, refined by the DFT distance scan instead.

if __name__ == "__main__":
    build_all()
    jobs = []
    for tag in ["CuSCN3", "Cu4I4", "Me3PO", "thiophene", "BTD", "triazine", "Ph3PO"]:
        if not os.path.exists(os.path.join(RUNS, tag, "xtbopt.xyz")):
            print(f"[warn] substrate {tag} xtb failed, skipping", flush=True); continue
        for ctag, syms, xyz, anchor in place_ag(tag):
            p = os.path.join(STR, f"{ctag}.xyz")
            write_xyz(p, syms, xyz, ctag)
            if ctag in SKIP_XTB:
                wd = os.path.join(RUNS, ctag); os.makedirs(wd, exist_ok=True)
                shutil.copy(p, os.path.join(wd, "xtbopt.xyz"))
                print(f"[xtb] {ctag}: skipped (frozen substrate + DFT scan)", flush=True)
                jobs.append((ctag, anchor))
            elif run_xtb(p, ctag, uhf=1):
                jobs.append((ctag, anchor))
    print("jobs:", jobs, flush=True)
    psi4_stage(jobs)
    print("ALL DONE", flush=True)
