"""DFT (Psi4) counterpoise-corrected Ag binding energies at GFN2-xTB geometries.
Level: UKS PBE-D3BJ / def2-SVP (def2-ECP on Ag, Mo via basis set defaults).
"""
import psi4, os, json, sys

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")

psi4.set_memory("6 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_sp.out"), False)
psi4.set_options({
    "basis": "def2-svp",
    "scf_type": "df",
    "reference": "uks",
    "dft_spherical_points": 434,
    "dft_radial_points": 75,
    "maxiter": 300,
    "fail_on_maxiter": True,
    "guess": "sad",
})

METHOD = "pbe-d3bj"

def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0])
    atoms = [l.split()[:4] for l in lines[2:2 + n]]
    return atoms  # [sym, x, y, z]

def geom_string(atoms, ghost_idx=None, charge=0, mult=1):
    s = f"{charge} {mult}\n"
    for i, (sym, x, y, z) in enumerate(atoms):
        tag = f"Gh({sym})" if ghost_idx and i in ghost_idx else sym
        s += f"{tag} {x} {y} {z}\n"
    s += "symmetry c1\nno_reorient\nno_com\n"
    return s

def energy(atoms, ghost_idx=None, charge=0, mult=1, label=""):
    mol = psi4.geometry(geom_string(atoms, ghost_idx, charge, mult))
    e = psi4.energy(METHOD, molecule=mol)
    print(f"  {label}: {e:.8f} Eh", flush=True)
    psi4.core.clean()
    return e

H2EV = 27.211386
results = {}
out_json = os.path.join(RUNS, "psi4_binding_eV.json")
if os.path.exists(out_json):
    results = json.load(open(out_json))

# systems: (tag of complex run dir, n_substrate_atoms) — Ag is last atom
systems = [
    ("Mo3O9_Ag", 12),
    ("LiF32_Ag", 32),
    ("HATCN_Ag_face", 30),
    ("HATCN_Ag_CN", 30),
    ("pbPPhenB_Ag", 72),
    ("TPBi_Ag", 81),
]
only = sys.argv[1:] if len(sys.argv) > 1 else None

for tag, nsub in systems:
    if only and tag not in only:
        continue
    if tag in results:
        print(f"{tag}: cached {results[tag]:.3f} eV"); continue
    xyz = os.path.join(RUNS, tag, "xtbopt.xyz")
    if not os.path.exists(xyz):
        xyz = os.path.join(RUNS, tag, "in.xyz")  # LiF frozen case fallback
    atoms = read_xyz(xyz)
    assert atoms[-1][0] == "Ag", f"last atom not Ag in {tag}"
    agi = len(atoms) - 1
    print(f"== {tag} ({len(atoms)} atoms) ==", flush=True)
    e_cx  = energy(atoms, None, 0, 2, "complex(doublet)")
    e_sub = energy(atoms, {agi}, 0, 1, "substrate+ghostAg")
    e_ag  = energy(atoms, set(range(agi)), 0, 2, "Ag+ghostSub")
    eb = (e_sub + e_ag - e_cx) * H2EV
    results[tag] = eb
    print(f"{tag}: E_b(CP) = {eb:.3f} eV", flush=True)
    json.dump(results, open(out_json, "w"), indent=2)

# Ag2 reference at same level
if "Ag2" not in results and (not only or "Ag2" in only):
    xyz = os.path.join(RUNS, "Ag2", "xtbopt.xyz")
    atoms = read_xyz(xyz)
    e2 = energy(atoms, None, 0, 1, "Ag2")
    ea = energy(atoms, {1}, 0, 2, "Ag+ghost")
    eb2 = energy(atoms, {0}, 0, 2, "ghost+Ag")
    results["Ag2"] = (ea + eb2 - e2) * H2EV
    json.dump(results, open(out_json, "w"), indent=2)

print(json.dumps(results, indent=2))
