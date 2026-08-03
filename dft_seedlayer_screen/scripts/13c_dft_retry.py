"""DFT-only retry for Bphen_Ag and Cs2CO3_Ag with robust SCF options."""
import os, json, psi4

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")
psi4.set_memory("5 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_bphen2.out"), False)
H2EV = 27.211386

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

ROBUST = {"basis": "def2-svp", "scf_type": "df", "reference": "uks",
          "maxiter": 500, "guess": "sad",
          "level_shift": 1.0, "level_shift_cutoff": 1e-3,
          "damping_percentage": 20.0}

res_path = os.path.join(RUNS, "bphen_cs_binding_eV.json")
res = json.load(open(res_path)) if os.path.exists(res_path) else {}

for tag, mults in [("Bphen_Ag", (2, 1, 2)), ("Cs2CO3_Ag", (2, 1, 2))]:
    if tag in res: continue
    syms, xyz = read_xyz(os.path.join(RUNS, tag, "xtbopt.xyz"))
    n = len(syms); agi = n - 1
    psi4.set_options(ROBUST)
    e_cx = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, None, 0, mults[0]))); psi4.core.clean()
    e_sub = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, {agi}, 0, mults[1]))); psi4.core.clean()
    e_ag = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, set(range(agi)), 0, mults[2]))); psi4.core.clean()
    eb = (e_sub + e_ag - e_cx) * H2EV
    res[tag] = eb
    print(f"{tag}: E_b(CP) = {eb:.3f} eV", flush=True)
    json.dump(res, open(res_path, "w"), indent=2)
