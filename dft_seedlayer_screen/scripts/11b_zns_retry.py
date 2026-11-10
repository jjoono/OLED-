import numpy as np, os, json, psi4

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")
psi4.set_memory("5 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_zns2.out"), False)
psi4.set_options({"basis": "def2-svp", "scf_type": "df", "reference": "rks",
                  "maxiter": 400, "guess": "sad",
                  "level_shift": 1.0, "level_shift_cutoff": 1e-3,
                  "damping_percentage": 10.0, "soscf": True})
H2EV = 27.211386

def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0]); syms, xyz = [], []
    for l in lines[2:2+n]:
        p = l.split(); syms.append(p[0]); xyz.append(list(map(float, p[1:4])))
    return syms, np.array(xyz)

def gstr(syms, xyz, ghost=None, mult=1):
    s = f"0 {mult}\n"
    for i, (sym, c) in enumerate(zip(syms, xyz)):
        t = f"Gh({sym})" if ghost and i in ghost else sym
        s += f"{t} {c[0]} {c[1]} {c[2]}\n"
    return s + "symmetry c1\nno_reorient\nno_com\n"

syms, xyz = read_xyz(os.path.join(RUNS, "HATCN_ZnS", "xtbopt.xyz"))
n = len(syms)
e_cx = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, None, 1)))
print("complex ok", flush=True); psi4.core.clean()
e_sub = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, {n-2, n-1}, 1)))
print("sub ok", flush=True); psi4.core.clean()
e_ad = psi4.energy("pbe-d3bj", molecule=psi4.geometry(gstr(syms, xyz, set(range(n-2)), 1)))
eb = (e_sub + e_ad - e_cx) * H2EV
print(f"ZnS(molecule) on HATCN: {eb:.3f} eV", flush=True)

j = os.path.join(RUNS, "zns_binding_eV.json")
d = json.load(open(j)) if os.path.exists(j) else {}
d["ZnS_on_HATCN"] = eb
json.dump(d, open(j, "w"), indent=2)
