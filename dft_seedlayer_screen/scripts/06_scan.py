"""Rigid 1-D scan of Ag along its approach axis at DFT (PBE-D3BJ/def2-SVP),
then CP-corrected binding at the scan minimum. Robust replacement for optking
freezing. Axis: from nearest substrate heavy atom (or face normal) to Ag.
"""
import psi4, os, json, sys
import numpy as np

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")

psi4.set_memory("6 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_scan.out"), False)
METHOD = "pbe-d3bj"
H2EV = 27.211386

BASEOPTS = {"basis": "def2-svp", "scf_type": "df", "reference": "uks",
            "maxiter": 300, "guess": "sad"}

def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0])
    syms, xyz = [], []
    for l in lines[2:2 + n]:
        p = l.split(); syms.append(p[0]); xyz.append(list(map(float, p[1:4])))
    return syms, np.array(xyz)

def gstr(syms, xyz, ghost=None, charge=0, mult=1):
    s = f"{charge} {mult}\n"
    for i, (sym, c) in enumerate(zip(syms, xyz)):
        t = f"Gh({sym})" if ghost and i in ghost else sym
        s += f"{t} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}\n"
    return s + "symmetry c1\nno_reorient\nno_com\n"

def energy(syms, xyz, ghost=None, charge=0, mult=1, extra=None):
    opts = dict(BASEOPTS)
    if extra: opts.update(extra)
    psi4.set_options(opts)
    mol = psi4.geometry(gstr(syms, xyz, ghost, charge, mult))
    e = psi4.energy(METHOD, molecule=mol)
    psi4.core.clean()
    return e

SYSTEMS = {
    # tag: (geometry file rel to runs, special scf opts)
    "Mo3O9_Ag":      ("Mo3O9_Ag/xtbopt.xyz", {"level_shift": 1.0, "level_shift_cutoff": 1e-3,
                                              "damping_percentage": 20.0, "maxiter": 500}),
    "Mo3O8_Ag":      ("Mo3O8_Ag/xtbopt.xyz", {"level_shift": 1.0, "level_shift_cutoff": 1e-3,
                                              "damping_percentage": 20.0, "maxiter": 500}),
    "HATCN_Ag_face": ("HATCN_Ag_face/xtbopt.xyz", None),
    "HATCN_Ag_CN":   ("HATCN_Ag_CN/xtbopt.xyz", None),
    "pbPPhenB_Ag":   ("pbPPhenB_Ag/xtbopt.xyz", None),
    "TPBi_Ag":       ("TPBi_Ag/xtbopt.xyz", None),
    "LiF32_Ag":      ("LiF32_Ag/xtbopt.xyz", None),
}
OFFSETS = json.loads(os.environ.get("SCAN_OFFSETS", "[-0.15, 0.0, 0.15, 0.3, 0.5, 0.8]"))

only = sys.argv[1:] if len(sys.argv) > 1 else None
out_json = os.path.join(RUNS, "psi4_binding_refined_eV.json")
results = json.load(open(out_json)) if os.path.exists(out_json) else {}

for tag, (rel, extra) in SYSTEMS.items():
    if only and tag not in only:
        continue
    if tag in results:
        print(f"{tag}: cached {results[tag]:.3f} eV", flush=True); continue
    p = os.path.join(RUNS, rel)
    if not os.path.exists(p):
        p = os.path.join(RUNS, tag, "in.xyz")
    syms, xyz = read_xyz(p)
    agi = len(syms) - 1
    assert syms[agi] == "Ag"
    # approach axis
    sub = xyz[:agi]
    if tag in ("HATCN_Ag_face", "Mo3O9_Ag"):
        heavy = [i for i in range(agi) if syms[i] != "H"]
        c = sub[heavy].mean(axis=0)
        u, s, vt = np.linalg.svd(sub[heavy] - c)
        axis = vt[2]
        if np.dot(xyz[agi] - c, axis) < 0: axis = -axis
    else:
        d = np.linalg.norm(sub - xyz[agi], axis=1)
        j = int(np.argmin(d))
        axis = xyz[agi] - sub[j]; axis /= np.linalg.norm(axis)
    print(f"== {tag}: scanning ==", flush=True)
    energies = {}
    for off in OFFSETS:
        x2 = xyz.copy(); x2[agi] = xyz[agi] + off * axis
        try:
            e = energy(syms, x2, None, 0, 2, extra)
        except Exception as ex:
            print(f"  off {off:+.2f}: FAILED {type(ex).__name__}", flush=True)
            continue
        energies[off] = e
        print(f"  off {off:+.2f}: {e:.8f}", flush=True)
    if not energies:
        print(f"{tag}: all scan points failed", flush=True); continue
    best = min(energies, key=energies.get)
    xb = xyz.copy(); xb[agi] = xyz[agi] + best * axis
    e_cx = energies[best]
    e_sub = energy(syms, xb, {agi}, 0, 1, extra)
    e_ag  = energy(syms, xb, set(range(agi)), 0, 2, extra)
    eb = (e_sub + e_ag - e_cx) * H2EV
    results[tag] = eb
    results[tag + "_offset"] = best
    print(f"{tag}: E_b(CP) at off {best:+.2f} = {eb:.3f} eV", flush=True)
    json.dump(results, open(out_json, "w"), indent=2)

print(json.dumps(results, indent=2))
