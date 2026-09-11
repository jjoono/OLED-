"""Finish remaining pieces:
1. HATCN_Ag_face: CP binding at scan point +0.15 (complex E known: -1487.27468760)
2. LiF32_Ag: short outward scan + CP at min
3. Mo3O9_Ag: extended outward scan (+0.8..+1.8) + CP at min
maxiter capped to avoid pathological grinds.
"""
import psi4, os, json, sys
import numpy as np

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")
psi4.set_memory("6 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_finish.out"), False)
METHOD = "pbe-d3bj"
H2EV = 27.211386
BASEOPTS = {"basis": "def2-svp", "scf_type": "df", "reference": "uks",
            "maxiter": 150, "guess": "sad"}

def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0]); syms, xyz = [], []
    for l in lines[2:2+n]:
        p = l.split(); syms.append(p[0]); xyz.append(list(map(float, p[1:4])))
    return syms, np.array(xyz)

def gstr(syms, xyz, ghost=None, charge=0, mult=1):
    s = f"{charge} {mult}\n"
    for i, (sym, c) in enumerate(zip(syms, xyz)):
        t = f"Gh({sym})" if ghost and i in ghost else sym
        s += f"{t} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}\n"
    return s + "symmetry c1\nno_reorient\nno_com\n"

def energy(syms, xyz, ghost=None, charge=0, mult=1, extra=None):
    o = dict(BASEOPTS)
    if extra: o.update(extra)
    psi4.set_options(o)
    mol = psi4.geometry(gstr(syms, xyz, ghost, charge, mult))
    e = psi4.energy(METHOD, molecule=mol)
    psi4.core.clean()
    return e

out_json = os.path.join(RUNS, "psi4_binding_refined_eV.json")
results = json.load(open(out_json)) if os.path.exists(out_json) else {}
def save(): json.dump(results, open(out_json, "w"), indent=2)

MOOPTS = {"level_shift": 1.0, "level_shift_cutoff": 1e-3,
          "damping_percentage": 20.0, "maxiter": 400}

def face_axis(syms, xyz, agi):
    heavy = [i for i in range(agi) if syms[i] != "H"]
    c = xyz[heavy].mean(axis=0)
    u, s, vt = np.linalg.svd(xyz[heavy] - c)
    ax = vt[2]
    if np.dot(xyz[agi] - c, ax) < 0: ax = -ax
    return ax

# ---- 1. HATCN face CP at +0.15 ----
if "HATCN_Ag_face" not in results:
    syms, xyz = read_xyz(os.path.join(RUNS, "HATCN_Ag_face", "xtbopt.xyz"))
    agi = len(syms) - 1
    ax = face_axis(syms, xyz[:agi], agi) if False else None
    # recompute axis over substrate only
    heavy = [i for i in range(agi) if syms[i] != "H"]
    c = xyz[heavy].mean(axis=0)
    u, s, vt = np.linalg.svd(xyz[heavy] - c)
    ax = vt[2]
    if np.dot(xyz[agi] - c, ax) < 0: ax = -ax
    xb = xyz.copy(); xb[agi] = xyz[agi] + 0.15 * ax
    e_cx = -1487.27468760
    e_sub = energy(syms, xb, {agi}, 0, 1)
    print("HATCN face e_sub done", flush=True)
    e_ag = energy(syms, xb, set(range(agi)), 0, 2)
    eb = (e_sub + e_ag - e_cx) * H2EV
    results["HATCN_Ag_face"] = eb; results["HATCN_Ag_face_offset"] = 0.15
    print(f"HATCN_Ag_face: E_b(CP) = {eb:.3f} eV", flush=True); save()

# ---- 2. LiF scan ----
if "LiF32_Ag" not in results:
    syms, xyz = read_xyz(os.path.join(RUNS, "LiF32_Ag", "xtbopt.xyz"))
    agi = len(syms) - 1
    sub = xyz[:agi]
    d = np.linalg.norm(sub - xyz[agi], axis=1); j = int(np.argmin(d))
    ax = xyz[agi] - sub[j]; ax /= np.linalg.norm(ax)
    es = {}
    for off in [0.0, 0.15, 0.3, 0.5]:
        x2 = xyz.copy(); x2[agi] = xyz[agi] + off * ax
        try:
            es[off] = energy(syms, x2, None, 0, 2)
            print(f"LiF off {off:+.2f}: {es[off]:.8f}", flush=True)
        except Exception as ex:
            print(f"LiF off {off:+.2f} FAILED {type(ex).__name__}", flush=True)
    best = min(es, key=es.get)
    xb = xyz.copy(); xb[agi] = xyz[agi] + best * ax
    e_sub = energy(syms, xb, {agi}, 0, 1)
    e_ag = energy(syms, xb, set(range(agi)), 0, 2)
    eb = (e_sub + e_ag - es[best]) * H2EV
    results["LiF32_Ag"] = eb; results["LiF32_Ag_offset"] = best
    print(f"LiF32_Ag: E_b(CP) at {best:+.2f} = {eb:.3f} eV", flush=True); save()

# ---- 3. Mo3O9 extended ----
if results.get("Mo3O9_Ag_offset") == 0.8:
    syms, xyz = read_xyz(os.path.join(RUNS, "Mo3O9_Ag", "xtbopt.xyz"))
    agi = len(syms) - 1
    heavy = list(range(agi))
    c = xyz[heavy].mean(axis=0)
    u, s, vt = np.linalg.svd(xyz[heavy] - c)
    ax = vt[2]
    if np.dot(xyz[agi] - c, ax) < 0: ax = -ax
    es = {0.8: None}
    for off in [0.8, 1.1, 1.4, 1.8]:
        x2 = xyz.copy(); x2[agi] = xyz[agi] + off * ax
        try:
            es[off] = energy(syms, x2, None, 0, 2, MOOPTS)
            print(f"Mo3O9 off {off:+.2f}: {es[off]:.8f}", flush=True)
        except Exception as ex:
            print(f"Mo3O9 off {off:+.2f} FAILED {type(ex).__name__}", flush=True)
    es = {k: v for k, v in es.items() if v is not None}
    if es:
        best = min(es, key=es.get)
        xb = xyz.copy(); xb[agi] = xyz[agi] + best * ax
        e_sub = energy(syms, xb, {agi}, 0, 1, MOOPTS)
        e_ag = energy(syms, xb, set(range(agi)), 0, 2, MOOPTS)
        eb = (e_sub + e_ag - es[best]) * H2EV
        results["Mo3O9_Ag"] = eb; results["Mo3O9_Ag_offset"] = best
        print(f"Mo3O9_Ag: E_b(CP) at {best:+.2f} = {eb:.3f} eV", flush=True); save()

print(json.dumps(results, indent=2))
