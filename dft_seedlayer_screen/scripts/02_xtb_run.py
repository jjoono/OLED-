"""Run GFN2-xTB optimizations and compute Ag binding energies."""
import subprocess, os, re, json, shutil, sys

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STR = os.path.join(BASE, "structures")
RUNS = os.path.join(BASE, "runs")
os.makedirs(RUNS, exist_ok=True)

def run_xtb(xyz, tag, uhf=0, opt=True, fix_atoms=None, gfn=2):
    wd = os.path.join(RUNS, tag)
    os.makedirs(wd, exist_ok=True)
    shutil.copy(xyz, os.path.join(wd, "in.xyz"))
    args = ["xtb", "in.xyz", "--gfn", str(gfn), "--chrg", "0", "--uhf", str(uhf)]
    if opt:
        args += ["--opt", "tight"]
    if fix_atoms:
        with open(os.path.join(wd, "xtb.inp"), "w") as f:
            f.write("$fix\n  atoms: " + ",".join(map(str, fix_atoms)) + "\n$end\n")
        args += ["--input", "xtb.inp"]
    with open(os.path.join(wd, "xtb.log"), "w") as log:
        r = subprocess.run(args, cwd=wd, stdout=log, stderr=subprocess.STDOUT)
    txt = open(os.path.join(wd, "xtb.log"), encoding="utf-8", errors="ignore").read()
    m = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", txt)
    if not m:
        print(f"FAILED: {tag} (exit {r.returncode})"); sys.exit(1)
    e = float(m[-1])
    print(f"{tag}: E = {e:.6f} Eh")
    return e

H2EV = 27.211386
E = {}

# references
E["Ag1"] = run_xtb(os.path.join(STR, "Ag1.xyz"), "Ag1", uhf=1, opt=False)
E["Ag2"] = run_xtb(os.path.join(STR, "Ag2.xyz"), "Ag2", uhf=0, opt=True)

# bare substrates
E["HATCN"]    = run_xtb(os.path.join(STR, "HATCN.xyz"), "HATCN")
E["TPBi"]     = run_xtb(os.path.join(STR, "TPBi.xyz"), "TPBi")
E["pbPPhenB"] = run_xtb(os.path.join(STR, "pbPPhenB.xyz"), "pbPPhenB")
E["Mo3O9"]    = run_xtb(os.path.join(STR, "Mo3O9.xyz"), "Mo3O9")
E["LiF32"]    = run_xtb(os.path.join(STR, "LiF32.xyz"), "LiF32", opt=False)  # frozen bulk

# complexes (doublets)
E["HATCN_Ag_face"] = run_xtb(os.path.join(STR, "HATCN_Ag_face.xyz"), "HATCN_Ag_face", uhf=1)
E["HATCN_Ag_CN"]   = run_xtb(os.path.join(STR, "HATCN_Ag_CN.xyz"), "HATCN_Ag_CN", uhf=1)
E["TPBi_Ag"]       = run_xtb(os.path.join(STR, "TPBi_Ag_N.xyz"), "TPBi_Ag", uhf=1)
E["pbPPhenB_Ag"]   = run_xtb(os.path.join(STR, "pbPPhenB_Ag_chelate.xyz"), "pbPPhenB_Ag", uhf=1)
E["Mo3O9_Ag"]      = run_xtb(os.path.join(STR, "Mo3O9_Ag.xyz"), "Mo3O9_Ag", uhf=1)
E["LiF32_Ag"]      = run_xtb(os.path.join(STR, "LiF32_Ag.xyz"), "LiF32_Ag", uhf=1,
                             fix_atoms=list(range(1, 33)))  # freeze LiF, relax Ag

res = {}
def eb(cx, sub):
    return (E[sub] + E["Ag1"] - E[cx]) * H2EV

res["Ag2_BE_eV"] = (2 * E["Ag1"] - E["Ag2"]) * H2EV
res["HATCN_face"] = eb("HATCN_Ag_face", "HATCN")
res["HATCN_CN"]   = eb("HATCN_Ag_CN", "HATCN")
res["TPBi"]       = eb("TPBi_Ag", "TPBi")
res["pbPPhenB"]   = eb("pbPPhenB_Ag", "pbPPhenB")
res["Mo3O9"]      = eb("Mo3O9_Ag", "Mo3O9")
res["LiF"]        = eb("LiF32_Ag", "LiF32")

print(json.dumps(res, indent=2))
with open(os.path.join(RUNS, "xtb_binding_eV.json"), "w") as f:
    json.dump(res, f, indent=2)
