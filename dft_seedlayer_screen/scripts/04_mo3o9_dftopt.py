"""DFT re-optimization of Mo3O9 + Ag (and bare Mo3O9) — GFN2 geometry was suspect.
UKS PBE-D3BJ/def2-SVP, then CP-corrected binding at the DFT minimum.
"""
import psi4, os, json

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")

psi4.set_memory("6 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_mo3o9_opt.out"), False)
psi4.set_options({
    "basis": "def2-svp", "scf_type": "df", "reference": "uks",
    "maxiter": 300, "guess": "sad",
    "geom_maxiter": 100,
})
METHOD = "pbe-d3bj"
H2EV = 27.211386

def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    return [l.split()[:4] for l in lines[2:2 + int(lines[0])]]

def gstr(atoms, ghost=None, charge=0, mult=1):
    s = f"{charge} {mult}\n"
    for i, (sym, x, y, z) in enumerate(atoms):
        t = f"Gh({sym})" if ghost and i in ghost else sym
        s += f"{t} {x} {y} {z}\n"
    return s + "symmetry c1\nno_reorient\nno_com\n"

# optimize complex (doublet)
atoms = read_xyz(os.path.join(RUNS, "Mo3O9_Ag", "xtbopt.xyz"))
mol = psi4.geometry(gstr(atoms, None, 0, 2))
e_cx = psi4.optimize(METHOD, molecule=mol)
print(f"complex opt E = {e_cx:.8f}", flush=True)
geom = mol.geometry().np * 0.52917721067
syms = [mol.symbol(i) for i in range(mol.natom())]
opt_atoms = [[syms[i], *map(str, geom[i])] for i in range(len(syms))]
with open(os.path.join(RUNS, "Mo3O9_Ag", "dftopt.xyz"), "w") as f:
    f.write(f"{len(opt_atoms)}\nPBE-D3BJ/def2-SVP optimized\n")
    for a in opt_atoms:
        f.write(" ".join(a) + "\n")

# optimize bare Mo3O9 (singlet)
atoms0 = read_xyz(os.path.join(RUNS, "Mo3O9", "xtbopt.xyz"))
mol0 = psi4.geometry(gstr(atoms0, None, 0, 1))
e_sub_rel = psi4.optimize(METHOD, molecule=mol0)
print(f"bare Mo3O9 opt E = {e_sub_rel:.8f}", flush=True)

# CP pieces at complex geometry
agi = len(opt_atoms) - 1
m = psi4.geometry(gstr(opt_atoms, {agi}, 0, 1)); e_sub_cp = psi4.energy(METHOD, molecule=m)
m = psi4.geometry(gstr(opt_atoms, set(range(agi)), 0, 2)); e_ag_cp = psi4.energy(METHOD, molecule=m)
eb_cp = (e_sub_cp + e_ag_cp - e_cx) * H2EV          # CP, no deformation
# also fully relaxed reference (with free Ag atom energy from CP piece basis ~ ok)
eb_rel = (e_sub_rel + (e_ag_cp) - e_cx) * H2EV       # approx: relaxed sub + Ag(ghost basis)

out = {"Mo3O9_Ag_dftopt_CP": eb_cp, "Mo3O9_Ag_dftopt_relaxedSub": eb_rel}
print(json.dumps(out, indent=2), flush=True)
j = os.path.join(RUNS, "psi4_binding_eV.json")
d = json.load(open(j)) if os.path.exists(j) else {}
d.update(out)
json.dump(d, open(j, "w"), indent=2)
