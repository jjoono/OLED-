"""Robust final pass: each system in try/except; Bphen complex with SOSCF."""
import os, json, traceback, psi4

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")
psi4.set_memory("5 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_bphen5.out"), False)
H2EV = 27.211386

def read_xyz(p):
    L = open(p).read().strip().splitlines(); n = int(L[0])
    return [l.split()[:4] for l in L[2:2+n]]

def gstr(at, ghost=None, charge=0, mult=1):
    s = f"{charge} {mult}\n"
    for i, (sym, x, y, z) in enumerate(at):
        t = f"Gh({sym})" if ghost and i in ghost else sym
        s += f"{t} {x} {y} {z}\n"
    return s + "symmetry c1\nno_reorient\nno_com\n"

PLAIN = {"basis": "def2-svp", "scf_type": "df", "reference": "uks",
         "maxiter": 300, "guess": "sad", "guess_mix": False, "soscf": False}
SOSCF = {"basis": "def2-svp", "scf_type": "df", "reference": "uks",
         "maxiter": 400, "guess": "gwh", "guess_mix": True,
         "soscf": True, "soscf_start_convergence": 1e-2}

def en(at, ghost, charge, mult, opts):
    psi4.set_options(opts)
    e = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(at, ghost, charge, mult)))
    psi4.core.clean()
    return e

jp = os.path.join(RUNS, "bphen_cs_binding_eV.json")
res = json.load(open(jp)) if os.path.exists(jp) else {}
def save(): json.dump(res, open(jp, "w"), indent=2)

# 1. Bphen_Ag (complex with SOSCF)
if "Bphen_Ag" not in res:
    try:
        at = read_xyz(os.path.join(RUNS, "Bphen_Ag", "xtbopt.xyz")); agi = len(at)-1
        e_cx = en(at, None, 0, 2, SOSCF); print("Bphen complex ok", flush=True)
        e_sub = en(at, {agi}, 0, 1, PLAIN)
        e_ag = en(at, set(range(agi)), 0, 2, PLAIN)
        res["Bphen_Ag"] = (e_sub + e_ag - e_cx) * H2EV
        print(f"Bphen_Ag: E_b(CP) = {res['Bphen_Ag']:.3f} eV", flush=True); save()
    except Exception:
        print("Bphen_Ag FAILED:", traceback.format_exc().splitlines()[-1], flush=True)

# 2. Cs2CO3_Ag
if "Cs2CO3_Ag" not in res:
    try:
        at = read_xyz(os.path.join(RUNS, "Cs2CO3_Ag", "xtbopt.xyz")); agi = len(at)-1
        e_cx = en(at, None, 0, 2, PLAIN)
        e_sub = en(at, {agi}, 0, 1, PLAIN)
        e_ag = en(at, set(range(agi)), 0, 2, PLAIN)
        res["Cs2CO3_Ag"] = (e_sub + e_ag - e_cx) * H2EV
        print(f"Cs2CO3_Ag: E_b(CP) = {res['Cs2CO3_Ag']:.3f} eV", flush=True); save()
    except Exception:
        print("Cs2CO3_Ag FAILED:", traceback.format_exc().splitlines()[-1], flush=True)

# 3. Bphen anion vertical
if "Bphen_anion_Ag_vertical" not in res:
    try:
        at = read_xyz(os.path.join(RUNS, "Bphen_Ag", "xtbopt.xyz")); agi = len(at)-1
        e_cx = en(at, None, -1, 1, SOSCF)
        e_sub = en(at, {agi}, -1, 2, PLAIN)
        e_ag = en(at, set(range(agi)), 0, 2, PLAIN)
        res["Bphen_anion_Ag_vertical"] = (e_sub + e_ag - e_cx) * H2EV
        print(f"Bphen_anion_Ag_vertical: E_b(CP) = {res['Bphen_anion_Ag_vertical']:.3f} eV", flush=True); save()
    except Exception:
        print("Bphen_anion FAILED:", traceback.format_exc().splitlines()[-1], flush=True)

print(json.dumps(res, indent=2))
