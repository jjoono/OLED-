"""F6TCNNQ: the candidate that was dismissed on reasoning that has since collapsed.

WHY IT IS BACK. F6TCNNQ appears twice in the project notes, both times as a
rejection:

  LIT_NOTES: "nitrile chemistry = same family as HATCN (F4TCNQ 0.97 < 1.03),
              unlikely to be better"
  REPORT.md: "F4TCNQ/F6-TCNNQ (CT band; ... optically excluded)"

Both arguments are now known to be wrong. The first rests on the cluster ordering
F4TCNQ 0.97 < HATCN 1.03, which the periodic slab REVERSED (F4TCNQ 1.556 vs HATCN
1.346, scripts/43) -- so the premise the extrapolation was built on is false. The
second rests on the anion's CT band disqualifying this family optically, which
TD-DFT contradicted: HATCN's anion absorbs MORE in the visible than F4TCNQ's
(0.299 vs 0.197, scripts/44). Neither rejection survives, and F6TCNNQ has never
actually been calculated.

WHY IT MATTERS. F6TCNNQ is the most widely used evaporable p-dopant in OLED
practice, preferred over F4TCNQ precisely because its lower vapour pressure makes
thermal evaporation controllable -- the same process constraint that applies here.
A reviewer of a paper about nitrile acceptors as Ag seeds will ask about it first.

WHAT THE DESIGN RULE PREDICTS. Extending the naphthoquinodimethane core and adding
two more fluorines should raise the electron affinity above F4TCNQ's, and if E_b
and EA are as coupled as Fig. 1(b) claims, E_b should rise too. That puts F6TCNNQ
further along the same trade-off line: better anchoring, worse oxidation risk. The
prediction is falsifiable -- if EA rises and E_b does not, the correlation the
figure is built on is not there, and the design rule needs rethinking. Either
outcome is informative, which is the point of running it.

PROTOCOL. Identical to scripts/10 and 30 so the number is comparable:
RDKit -> GFN2-xTB opt -> PBE-D3BJ/def2-SVP UKS rigid scan of Ag along the C-N axis
-> counterpoise-corrected E_b at the scan minimum. EA then goes through scripts/45
unchanged, so it lands at the same level as HATCN/F4TCNQ/TCNQ and in the same cache.
"""
import os, sys, json, shutil, subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STR, RUNS = os.path.join(BASE, "structures"), os.path.join(BASE, "runs")
H2EV = 27.211386

# 2,2'-(perfluoronaphthalene-2,6-diylidene)dimalononitrile.
# Built as F4TCNQ's SMILES (scripts/10) extended from a benzene to a naphthalene
# core: the 2,6 positions carry the dicyanomethylene groups and all six remaining
# ring positions are fluorinated. Expected formula C16F6N4 (MW 362.2) -- checked
# at runtime rather than trusted, because a mis-drawn SMILES here would silently
# produce a different molecule.
SMILES = "N#CC(C#N)=C1C(F)=C(F)C2=C(F)C(=C(C#N)C#N)C(F)=C(F)C2=C1F"
EXPECT_FORMULA = "C16F6N4"
TAG = "F6TCNNQ"

REF = {"HATCN": 1.029, "F4TCNQ": 0.966}          # cluster E_b, same protocol
REF_SLAB = {"HATCN": 1.346, "F4TCNQ": 1.556}     # monolayer E_b
REF_EA = {"TCNQ": 3.536, "F4TCNQ": 4.020, "HATCN": 3.378}


def write_xyz(path, syms, xyz, cm=""):
    with open(path, "w") as f:
        f.write(f"{len(syms)}\n{cm}\n")
        for s, c in zip(syms, xyz):
            f.write(f"{s} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n")


def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0]); syms, xyz = [], []
    for l in lines[2:2 + n]:
        p = l.split(); syms.append(p[0]); xyz.append(list(map(float, p[1:4])))
    return syms, np.array(xyz)


def run_xtb(xyz_path, tag, uhf=0):
    wd = os.path.join(RUNS, tag)
    os.makedirs(wd, exist_ok=True)
    if os.path.exists(os.path.join(wd, "xtbopt.xyz")):
        print(f"[xtb] {tag}: cached", flush=True)
        return True
    shutil.copy(xyz_path, os.path.join(wd, "in.xyz"))
    with open(os.path.join(wd, "xtb.log"), "w") as log:
        subprocess.run(["xtb", "in.xyz", "--gfn", "2", "--uhf", str(uhf), "--opt", "tight"],
                       cwd=wd, stdout=log, stderr=subprocess.STDOUT)
    ok = os.path.exists(os.path.join(wd, "xtbopt.xyz"))
    print(f"[xtb] {tag}: {'ok' if ok else 'FAILED'}", flush=True)
    return ok


def build():
    mol = Chem.MolFromSmiles(SMILES)
    if mol is None:
        raise SystemExit("SMILES did not parse")
    formula = CalcMolFormula(mol)
    print(f"{TAG}: {formula}", flush=True)
    if formula.replace("+", "").replace("-", "") != EXPECT_FORMULA:
        raise SystemExit(f"formula {formula} != expected {EXPECT_FORMULA}; "
                         "the SMILES does not describe F6TCNNQ -- fix before running")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    conf = mol.GetConformer()
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                     conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
    p = os.path.join(STR, f"{TAG}.xyz")
    write_xyz(p, syms, xyz, TAG)
    run_xtb(p, TAG)

    # copy the relaxed geometry to structures/ so scripts/45 can find it by name
    s, x = read_xyz(os.path.join(RUNS, TAG, "xtbopt.xyz"))
    write_xyz(p, s, x, f"{TAG} GFN2-xTB optimised")
    return s, x


def place_ag(syms, x):
    """Ag on the outermost nitrile N, along the C-N axis -- same site as F4TCNQ."""
    cen = x.mean(axis=0)
    ns = [i for i, s in enumerate(syms) if s == "N"]
    i_n = max(ns, key=lambda k: np.linalg.norm(x[k] - cen))
    d = np.linalg.norm(x - x[i_n], axis=1); d[i_n] = 9e9
    i_c = int(np.argmin(d))
    ax = x[i_n] - x[i_c]; ax /= np.linalg.norm(ax)
    return np.vstack([x, x[i_n] + 2.30 * ax]), i_n


def cluster_eb(syms, x):
    import psi4
    psi4.set_memory(os.environ.get("PSI4_MEM", "6 GB"))
    psi4.set_num_threads(int(os.environ.get("PSI4_THREADS", "4")))
    psi4.core.set_output_file(os.path.join(RUNS, "psi4_f6tcnnq.out"), False)
    # WHY THIS DIFFERS FROM scripts/30. The first attempt used that script's
    # settings verbatim (maxiter 150, SAD guess, near-to-far scan) and every single
    # scan point failed to converge -- 70 minutes, zero usable energies. That is
    # not a numerical accident, and it is worth recording because it is the same
    # chemistry the rest of the project is about: F6TCNNQ's electron affinity is
    # high enough that Ag(0) + neutral molecule and Ag(+) + radical anion are
    # nearly degenerate, so the SCF oscillates between two charge states. HATCN and
    # F4TCNQ converged without complaint at this level; F6TCNNQ does not.
    #
    # The fix is the one used for the slabs in scripts/39: start FAR, where the two
    # configurations are cleanly separated and the neutral solution is unambiguous,
    # then walk inward carrying the converged orbitals forward as the guess. Also
    # raises maxiter, since a hard case needs iterations more than it needs tricks.
    base = {"basis": "def2-svp", "scf_type": "df", "reference": "uks",
            "maxiter": 500, "guess": "sad"}
    psi4.set_options(base)
    M = "pbe-d3bj"

    xs_all, i_n = place_ag(syms, x)
    s_all = syms + ["Ag"]
    agi = len(s_all) - 1

    def gstr(sy, xyz, ghost=None, mult=1):
        s = f"0 {mult}\n"
        for i, (sym, c) in enumerate(zip(sy, xyz)):
            t = f"Gh({sym})" if ghost and i in ghost else sym
            s += f"{t} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n"
        return s + "symmetry c1\nno_reorient\nno_com\n"

    def energy(g, carry=False):
        """SCF with orbital carry-over. `carry` reads the previous point's orbitals.

        Returns (energy, converged_ok). The guess is only reset to SAD if reading
        fails, because a converged neighbour is a far better starting point than
        any atomic guess for a system this close to a charge-transfer crossing.
        """
        # SOSCF FIRST, and this is the whole fix. Diagnosis of the failed runs:
        # near convergence the SCF settles at -1571.9958, which IS the answer (the
        # fragment sum at 20 A is -1571.9962), and then ADIIS throws it to
        # -1571.8486 and back, forever. The two states are 0.147 Ha = 4.0 eV apart,
        # the scale of Ag(0)+neutral versus Ag(+)+anion, so DIIS extrapolation keeps
        # mixing two charge states instead of picking one.
        #
        # More iterations do not help -- the oscillation is stable, and 500 of them
        # bought nothing. Damping and level shifts do not help either; both failed.
        # Second-order SCF does, because it stops extrapolating and follows the
        # gradient once DIIS has got close. soscf_start_convergence hands over at
        # 1e-2, well before the oscillation sets in.
        #
        # (The hundreds-of-Hartree swings visible early in the log are ordinary
        # start-up transients from the atomic guess, not the failure mode. Reading
        # them as the failure sent this diagnosis down the wrong path once already.)
        soscf = {"soscf": True, "soscf_start_convergence": 1.0e-2,
                 "soscf_max_iter": 60, "soscf_conv": 1.0e-6}
        ladder = [
            {**soscf, **({"guess": "read"} if carry else {})},
            {**soscf, "guess": "sad", "damping_percentage": 30},
            {**soscf, "guess": "gwh", "damping_percentage": 50,
             "level_shift": 0.5, "level_shift_cutoff": 1e-3},
            {"guess": "sad", "damping_percentage": 60},   # last resort, plain DIIS
        ]
        for i, extra in enumerate(ladder):
            try:
                psi4.set_options({**base, **extra})
                e, wfn = psi4.energy(M, molecule=psi4.geometry(g), return_wfn=True)
                psi4.set_options(base)
                return e, wfn
            except Exception as ex:
                print(f"    scf attempt {i} failed ({type(ex).__name__})", flush=True)
                if i < len(ladder) - 1:
                    psi4.core.clean()      # keep scratch on the last try for diagnosis
        return None, None

    ax = xs_all[agi] - xs_all[i_n]
    r0 = np.linalg.norm(ax); ax = ax / r0
    scan = {}
    # FAR -> NEAR. At 6 A the Ag atom is a spectator and the neutral solution is
    # unambiguous; each converged point then seeds the next. Scanning the other way
    # (as the first attempt did) starts in the region where the two charge states
    # compete, with only an atomic guess to pick between them.
    radii = [6.0, 5.0, 4.0, 3.4, 3.0, 2.7, 2.45, 2.3, 2.2]
    first = True
    for r in radii:
        xs = xs_all.copy()
        xs[agi] = xs_all[i_n] + r * ax
        e, _ = energy(gstr(s_all, xs, None, 2), carry=not first)
        if e is None:
            print(f"  r = {r:.2f} A   FAILED", flush=True)
            continue
        first = False
        scan[round(r, 3)] = e
        print(f"  r = {r:.2f} A   E = {e:.6f}", flush=True)
    if not scan:
        return {"error": "scan_failed"}

    rb = min(scan, key=scan.get)
    if rb == max(scan):
        # Minimum at the outermost point means the curve is still descending
        # outward, i.e. there is no bound minimum on this grid -- not a binding
        # energy. Report it rather than quoting the endpoint as if it were one.
        print("  WARNING: minimum at the largest radius scanned", flush=True)
    xs = xs_all.copy(); xs[agi] = xs_all[i_n] + rb * ax
    write_xyz(os.path.join(STR, f"{TAG}_Ag.xyz"), s_all, xs, f"{TAG}+Ag, DFT scan min")
    e_cx = scan[rb]
    e_sub, _ = energy(gstr(s_all, xs, {agi}, 1))          # molecule + Ag ghost
    e_ag, _ = energy(gstr(s_all, xs, set(range(agi)), 2))  # Ag + molecule ghost
    if None in (e_sub, e_ag):
        return {"error": "cp_failed", "scan": {str(k): v for k, v in scan.items()}}
    return {"Eb_eV": (e_sub + e_ag - e_cx) * H2EV, "r_A": rb,
            "n_scan_ok": len(scan), "n_scan_tried": len(radii),
            "scan": {str(k): v for k, v in scan.items()}}


def main():
    syms, x = build()

    res_path = os.path.join(RUNS, "f6tcnnq.json")
    res = json.load(open(res_path)) if os.path.exists(res_path) else {}

    if "Eb_eV" not in res:
        print("\n=== cluster E_b (PBE-D3BJ/def2-SVP, counterpoise) ===", flush=True)
        res.update(cluster_eb(syms, x))
        json.dump(res, open(res_path, "w"), indent=2)

    # EA through scripts/45, unchanged, so the level matches the other three
    print("\n=== electron affinity (wB97X/def2-TZVP, vertical dSCF) ===", flush=True)
    import importlib.util
    _s = importlib.util.spec_from_file_location(
        "s45", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "45_ag_extraction.py"))
    S45 = importlib.util.module_from_spec(_s)
    _s.loader.exec_module(S45)
    try:
        res["EA_eV"] = S45.vertical_ea(TAG, f"{TAG}.xyz", S45.BASIS)
    except Exception as ex:
        print(f"  EA FAILED: {type(ex).__name__}: {ex}", flush=True)
        res["EA_eV"] = None
    json.dump(res, open(res_path, "w"), indent=2)

    # ------------------------------------------------------------------ report
    print("\n" + "=" * 66)
    print("F6TCNNQ vs the family")
    print("=" * 66)
    eb, ea = res.get("Eb_eV"), res.get("EA_eV")
    print(f"{'molecule':<10}{'cluster E_b':>13}{'slab E_b':>11}{'EA':>9}")
    for m in ("HATCN", "F4TCNQ"):
        print(f"{m:<10}{REF[m]:>13.3f}{REF_SLAB[m]:>11.3f}{REF_EA[m]:>9.3f}")
    print(f"{TAG:<10}{eb if eb is None else f'{eb:>13.3f}'}"
          f"{'not run':>11}{ea if ea is None else f'{ea:>9.3f}'}")

    if eb is None or ea is None:
        print("\nincomplete -- see errors above")
        return

    print("\nTEST OF THE DESIGN RULE")
    d_ea = ea - REF_EA["F4TCNQ"]
    d_eb = eb - REF["F4TCNQ"]
    print(f"  vs F4TCNQ:  dEA = {d_ea:+.3f} eV,  dE_b = {d_eb:+.3f} eV")
    if d_ea > 0.05 and d_eb > 0.02:
        print("  -> Both rise together, as Fig. 1(b) assumes. F6TCNNQ sits further")
        print("     along the trade-off line: stronger anchoring, higher oxidation")
        print("     driving force. The correlation gains a third point.")
    elif d_ea > 0.05 and d_eb <= 0.02:
        print("  -> EA rises but E_b does NOT. The two axes are then not coupled the")
        print("     way Fig. 1(b) claims, and the trade-off framing needs rethinking")
        print("     before it goes in a manuscript. This is the outcome that would")
        print("     falsify the design rule.")
    else:
        print("  -> EA did not rise above F4TCNQ. Check the geometry and the")
        print("     assumption that a larger fluorinated core is a stronger acceptor.")

    print("\n  NOTE: cluster E_b is the comparable number here. The slab value is")
    print("  what actually ranks the family (it reversed HATCN vs F4TCNQ), and it")
    print("  has NOT been computed for F6TCNNQ. Do not rank on the cluster number")
    print("  alone -- that is the mistake that got F6TCNNQ dismissed in the first")
    print("  place. scripts/43 is the template if this looks worth pursuing.")


if __name__ == "__main__":
    main()
