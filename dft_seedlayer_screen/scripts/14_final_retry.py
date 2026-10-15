"""Final retry: Bphen_Ag, Cs2CO3_Ag (neutral) + Bphen anion effect (vertical, at
neutral Bphen_Ag geometry). ADIIS disabled (conflicts with damping), level shift only.
"""
import os, json, psi4

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")
psi4.set_memory("5 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_bphen3.out"), False)
H2EV = 27.211386

OPTS = {"basis": "def2-svp", "scf_type": "df", "reference": "uks",
        "maxiter": 500, "guess": "sad",
        "scf_initial_accelerator": "none",
        "level_shift": 2.0, "level_shift_cutoff": 1e-4}

def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0]); syms, xyz = [], []
    for l in lines[2:2+n]:
        p = l.split(); syms.append(p[0]); xyz.append(p[1:4])
    return syms, xyz

def gstr(syms, xyz, ghost=None, charge=0, mult=1):
    st = f"{charge} {mult}\n"
    for i, (sym, c) in enumerate(zip(syms, xyz)):
        t = f"Gh({sym})" if ghost and i in ghost else sym
        st += f"{t} {c[0]} {c[1]} {c[2]}\n"
    return st + "symmetry c1\nno_reorient\nno_com\n"

def en(syms, xyz, ghost, charge, mult):
    psi4.set_options(OPTS)
    e = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, ghost, charge, mult)))
    psi4.core.clean()
    return e

res_path = os.path.join(RUNS, "bphen_cs_binding_eV.json")
res = json.load(open(res_path)) if os.path.exists(res_path) else {}

# neutral Bphen_Ag and Cs2CO3_Ag
for tag in ["Bphen_Ag", "Cs2CO3_Ag"]:
    if tag in res:
        print(f"{tag}: cached {res[tag]:.3f}", flush=True); continue
    syms, xyz = read_xyz(os.path.join(RUNS, tag, "xtbopt.xyz"))
    agi = len(syms) - 1
    e_cx = en(syms, xyz, None, 0, 2)
    e_sub = en(syms, xyz, {agi}, 0, 1)
    e_ag = en(syms, xyz, set(range(agi)), 0, 2)
    res[tag] = (e_sub + e_ag - e_cx) * H2EV
    print(f"{tag}: E_b(CP) = {res[tag]:.3f} eV", flush=True)
    json.dump(res, open(res_path, "w"), indent=2)

# anion (vertical, at neutral geometry): [Bphen-Ag]^- singlet vs Bphen^-(doublet)+Ag(doublet)
if "Bphen_anion_Ag_vertical" not in res:
    syms, xyz = read_xyz(os.path.join(RUNS, "Bphen_Ag", "xtbopt.xyz"))
    agi = len(syms) - 1
    e_cx = en(syms, xyz, None, -1, 1)
    e_sub = en(syms, xyz, {agi}, -1, 2)
    e_ag = en(syms, xyz, set(range(agi)), 0, 2)
    res["Bphen_anion_Ag_vertical"] = (e_sub + e_ag - e_cx) * H2EV
    print(f"Bphen_anion_Ag_vertical: E_b(CP) = {res['Bphen_anion_Ag_vertical']:.3f} eV", flush=True)
    json.dump(res, open(res_path, "w"), indent=2)

print(json.dumps(res, indent=2))
