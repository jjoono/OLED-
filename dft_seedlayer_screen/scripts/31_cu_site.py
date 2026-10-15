"""Addendum to 30: Ag approaching the METAL (Cu) site of Cu4I4 and (CuSCN)3.
Rationale: JMCA 2023 reports Ag-induced Cu+ disproportionation / Cu-Ag alloying at
CuSCN|Ag interfaces; Ag-Cu dimer BE is 1.93 eV -> the strong anchoring channel on
Cu(I) compounds may be metallophilic Ag-Cu, not Ag-anion. Frozen substrate +
DFT distance scan + CP, same level (PBE-D3BJ/def2-SVP, UKS).
"""
import numpy as np, os, json
BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STR, RUNS = os.path.join(BASE, "structures"), os.path.join(BASE, "runs")

def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0]); syms, xyz = [], []
    for l in lines[2:2+n]:
        p = l.split(); syms.append(p[0]); xyz.append(list(map(float, p[1:4])))
    return syms, np.array(xyz)

def unit(v): return v / np.linalg.norm(v)

import psi4
psi4.set_memory("6 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_cu_site.out"), False)
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

def robust_energy(geom_str, mult):
    for attempt, extra in enumerate([{}, {"level_shift": 1.0, "level_shift_cutoff": 0.01,
                                          "damping_percentage": 15}]):
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

for sub, ctag in [("Cu4I4", "Cu4I4_AgCu"), ("CuSCN3", "CuSCN3_AgCu")]:
    if ctag in res and "Eb_eV" in res.get(ctag, {}):
        print(f"[skip] {ctag}", flush=True); continue
    syms, x = read_xyz(os.path.join(RUNS, sub, "xtbopt.xyz"))
    cen = x.mean(axis=0)
    cus = [i for i, s in enumerate(syms) if s == "Cu"]
    # least-coordinated / most-exposed Cu: max distance from centroid
    ci = max(cus, key=lambda i: np.linalg.norm(x[i] - cen))
    d = unit(x[ci] - cen)
    scan = {}
    best = (None, None)
    for r in [2.35, 2.45, 2.55, 2.70, 2.90, 3.15, 3.5]:
        xs = np.vstack([x, x[ci] + r * d])
        e = robust_energy(gstr(list(syms) + ["Ag"], xs, None, 2), 2)
        if e is None: continue
        scan[r] = e
        if best[1] is None or e < best[1]: best = (r, e)
        print(f"  [{ctag}] r={r:.2f} E={e:.6f}", flush=True)
    if best[1] is None:
        res[ctag] = {"error": "scan_failed"}; json.dump(res, open(res_path, "w"), indent=2); continue
    rb = best[0]
    xs = np.vstack([x, x[ci] + rb * d])
    allsyms = list(syms) + ["Ag"]; agi = len(allsyms) - 1
    e_cx = scan[rb]
    e_sub = robust_energy(gstr(allsyms, xs, {agi}, 1), 1)
    e_ag = robust_energy(gstr(allsyms, xs, set(range(agi)), 2), 2)
    if None in (e_sub, e_ag):
        res[ctag] = {"error": "cp_failed"}
    else:
        eb = (e_sub + e_ag - e_cx) * H2EV
        res[ctag] = {"Eb_eV": eb, "r_A": rb, "scan": {str(k): v for k, v in scan.items()}}
        print(f"[RESULT] {ctag}: E_b(CP) = {eb:.3f} eV at r = {rb:.2f} A", flush=True)
    json.dump(res, open(res_path, "w"), indent=2)
    psi4.core.clean()
print("ALL DONE", flush=True)
