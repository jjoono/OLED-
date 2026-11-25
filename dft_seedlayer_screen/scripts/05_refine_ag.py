"""Refine Ag position at DFT level (substrate frozen), then CP binding energy.
GFN2 geometries are systematically too tight for PBE-D3 -> single points underbind.
Optimizes only the Ag atom (frozen_cartesian on all substrate atoms).
"""
import psi4, os, json, sys

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")

psi4.set_memory("6 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_refine.out"), False)
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

systems = [  # (tag, use dftopt if exists)
    "HATCN_Ag_face", "HATCN_Ag_CN", "pbPPhenB_Ag", "TPBi_Ag", "LiF32_Ag",
]
only = sys.argv[1:] if len(sys.argv) > 1 else None
out_json = os.path.join(RUNS, "psi4_binding_refined_eV.json")
results = json.load(open(out_json)) if os.path.exists(out_json) else {}

for tag in systems:
    if only and tag not in only:
        continue
    if tag in results:
        print(f"{tag}: cached {results[tag]:.3f} eV", flush=True); continue
    p = os.path.join(RUNS, tag, "xtbopt.xyz")
    if not os.path.exists(p):
        p = os.path.join(RUNS, tag, "in.xyz")
    atoms = read_xyz(p)
    agi = len(atoms) - 1
    assert atoms[agi][0] == "Ag"
    # frozen cartesian list (1-indexed, "atom xyz" strings) for all substrate atoms
    frozen = " ".join(f"{i+1} xyz" for i in range(agi))
    psi4.set_options({
        "basis": "def2-svp", "scf_type": "df", "reference": "uks",
        "maxiter": 300, "guess": "sad", "geom_maxiter": 60,
        "optking__frozen_cartesian": frozen,
    })
    mol = psi4.geometry(gstr(atoms, None, 0, 2))
    print(f"== {tag}: optimizing Ag position ==", flush=True)
    e_cx = psi4.optimize(METHOD, molecule=mol)
    geom = mol.geometry().np * 0.52917721067
    syms = [mol.symbol(i) for i in range(mol.natom())]
    opt_atoms = [[syms[i], *map(str, geom[i])] for i in range(len(syms))]
    with open(os.path.join(RUNS, tag, "dft_ag_refined.xyz"), "w") as f:
        f.write(f"{len(opt_atoms)}\nAg-refined PBE-D3BJ/def2-SVP\n")
        for a in opt_atoms:
            f.write(" ".join(a) + "\n")
    psi4.set_options({"optking__frozen_cartesian": ""})
    m = psi4.geometry(gstr(opt_atoms, {agi}, 0, 1)); e_sub = psi4.energy(METHOD, molecule=m)
    m = psi4.geometry(gstr(opt_atoms, set(range(agi)), 0, 2)); e_ag = psi4.energy(METHOD, molecule=m)
    eb = (e_sub + e_ag - e_cx) * H2EV
    results[tag] = eb
    print(f"{tag}: refined E_b(CP) = {eb:.3f} eV", flush=True)
    json.dump(results, open(out_json, "w"), indent=2)

print(json.dumps(results, indent=2))
