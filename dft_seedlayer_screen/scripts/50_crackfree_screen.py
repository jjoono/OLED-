"""Screen for organics that anchor Ag like HATCN but resist the cracking.

WHY. The HATCN 30 nm films crack, and the working diagnosis (see the crack-optics
discussion and REPORT threads) is NOT classical stress fracture -- channel-cracking
h_c comes out in micrometres for any sane parameters -- but crystallisation:
HATCN is small, rigid, perfectly planar and highly symmetric, i.e. everything a
molecule needs to be a BAD glass former. Grain growth through the film thickness
then reads as cracks/voids in SEM. So "seed like HATCN, without the cracks" means:

    keep the anchor chemistry (nitriles / chelating N),
    change the molecular SHAPE (3D, twisted, asymmetric, bigger).

TWO AXES, TWO METHODS.

  (1) Crystallisation risk -- cheap shape descriptors on the full molecule.
      * 3D-ness: NPR1+NPR2 (normalised principal moments of inertia). A flat disc
        or a rod gives ~1.0 (both pack and crystallise well); a globular propeller
        tends to 2.0. HATCN is the disc limit.
      * planarity: RMS out-of-plane extent (3rd singular value of the centred
        coordinates), in A.
      * rigidity: rotatable-bond count (rigid -> crystallises; some floppiness
        frustrates packing).
      * symmetry: fraction of symmetry-distinct heavy atoms (RDKit canonical
        ranks). Fewer distinct classes = more symmetric = better packing.
      * size: MW. Bigger molecules nucleate slower.
      The combined score is a HEURISTIC with transparent weights, and it is
      validated on knowns before being trusted: HATCN must come out worst,
      TPBi (Tg ~122 C, famously amorphous) among the best. If that fails, the
      score prints a warning telling you not to use it.

  (2) Anchoring -- fragment DFT at the established level (GFN2-xTB geometry,
      PBE-D3BJ/def2-SVP counterpoise, the scripts/30 pipeline), on the anchor
      moiety rather than the full molecule. Full 4CzIPN (88 atoms) is beyond this
      container, but E_b is local to the Ag-N contact; what the big substituents
      change is (a) the acceptor strength of the ring the nitrile hangs on and
      (b) the film-level pocket geometry. (a) is probed directly and cheaply:

        benzonitrile          PhCN      the neutral reference
        terephthalonitrile    p-(CN)2   acceptor-substituted ring  -> E_b up?
        4-(NMe2)benzonitrile  DMABN     donor-substituted ring     -> E_b down?
        phthalonitrile        o-(CN)2   can two adjacent nitriles chelate one Ag?

      The donor/acceptor pair brackets what carbazole substitution (4CzIPN-type
      donors) does to nitrile binding. The phthalonitrile case tests whether an
      o-dinitrile chelate can substitute for HATCN's inter-molecular pocket
      (+0.53 eV, scripts/41) WITHIN one molecule -- if yes, that is the design
      handle: a twisted glass former carrying an o-dicyanoarene anchor.

WHAT THIS DOES NOT DO. No full-molecule E_b for the big candidates, no monolayer
pocket correction, no Tg prediction (the quoted Tg values are from memory and the
session cannot reach literature -- they are labelled UNVERIFIED and used only as a
sanity check on the descriptor ranking, never as data). Transparency flags for the
CT emitters (4CzIPN, 2CzPN absorb in the visible -- that is what makes them
emitters) are chemistry knowledge, flagged for verification the same way.
"""
import os, sys, json, shutil, subprocess
import numpy as np

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STR, RUNS = os.path.join(BASE, "structures"), os.path.join(BASE, "runs")
H2EV = 27.211386

MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
        "P": 30.974, "S": 32.06}

# ---------------------------------------------------------------- candidates
# name: (SMILES, expected formula, anchor description, transparency note)
FULL = {
    "HATCN":    ("N#Cc1nc2c3nc(C#N)c(C#N)nc3c3nc(C#N)c(C#N)nc3c2nc1C#N",
                 "C18N12", "6x CN, pocket", "OK (<450 nm)"),
    "F4TCNQ":   ("N#CC(C#N)=C1C(F)=C(F)C(=C(C#N)C#N)C(F)=C1F",
                 "C12F4N4", "4x CN", "OK neutral; anion issue"),
    "F6TCNNQ":  ("N#CC(C#N)=C1C(F)=C(F)C2=C(F)C(=C(C#N)C#N)C(F)=C(F)C2=C1F",
                 "C16F6N4", "4x CN", "check"),
    "4CzIPN":   ("N#Cc1c(-n2c3ccccc3c3ccccc32)c(-n2c3ccccc3c3ccccc32)c(C#N)"
                 "c(-n2c3ccccc3c3ccccc32)c1-n1c2ccccc2c2ccccc21",
                 "C56H32N6", "2x CN (donor-flanked)",
                 "ABSORBS VISIBLE (green TADF emitter) -- likely disqualifying"),
    "2CzPN":    ("N#Cc1cc(-n2c3ccccc3c3ccccc32)c(-n2c3ccccc3c3ccccc32)cc1C#N",
                 "C32H18N4", "o-di-CN (donor-flanked)",
                 "ABSORBS VISIBLE (sky-blue emitter) -- likely disqualifying"),
    "mCPCN":    ("N#Cc1ccc2c(c1)c1ccccc1n2-c1cccc(-n2c3ccccc3c3ccccc32)c1",
                 "C31H19N3", "1x CN on carbazole", "OK (UV host)"),
    "TPBi":     ("c1ccc(cc1)n1c(nc2ccccc21)-c1cc(-c2nc3ccccc3n2-c2ccccc2)cc"
                 "(-c2nc3ccccc3n2-c2ccccc2)c1",
                 "C45H30N6", "3x benzimidazole N", "OK"),
    "B3PyMPM":  ("Cc1nc(-c2cc(-c3cccnc3)cc(-c3cccnc3)c2)cc(-c2cc(-c3cccnc3)cc"
                 "(-c3cccnc3)c2)n1",
                 "C37H26N6", "4x pyridyl N", "OK"),
    "Bphen":    ("c1ccc(cc1)-c1ccnc2c1ccc1c(-c3ccccc3)ccnc21",
                 "C24H16N2", "phen N,N chelate", "OK"),
    "TmPyPB":   ("c1cc(cnc1)-c1cccc(c1)-c1cc(cc(c1)-c1cccc(c1)-c1cccnc1)"
                 "-c1cccc(c1)-c1cccnc1",
                 "C39H27N3", "3x pyridyl N", "OK"),
    "PO-T2T":   ("O=P(c1ccccc1)(c1ccccc1)c1cccc(-c2nc(-c3cccc(P(=O)(c4ccccc4)"
                 "c4ccccc4)c3)nc(-c3cccc(P(=O)(c4ccccc4)c4ccccc4)c3)n2)c1",
                 "C57H42N3O3P3", "triazine N + 3x P=O", "OK"),
}

# Tg from memory. UNVERIFIED -- the session cannot reach the literature. Used only
# to sanity-check the descriptor ranking, never reported as data.
TG_MEMORY_C = {"TPBi": 122, "Bphen": 62, "TmPyPB": 79, "B3PyMPM": 108}

FRAGMENTS = {
    "PhCN":   ("N#Cc1ccccc1", "C7H5N"),
    "pDCNB":  ("N#Cc1ccc(C#N)cc1", "C8H4N2"),      # terephthalonitrile
    "DMABN":  ("N#Cc1ccc(N(C)C)cc1", "C9H10N2"),
    "oDCNB":  ("N#Cc1ccccc1C#N", "C8H4N2"),        # phthalonitrile, chelate test
}

EB_KNOWN = {  # cluster values already computed at this level elsewhere in the repo
    "HATCN": 1.029, "F4TCNQ": 0.966,
}


# ------------------------------------------------------------------ helpers
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


def embed(smiles, expect, name):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcNumRotatableBonds
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, f"SMILES parse failed"
    formula = CalcMolFormula(mol).replace("+", "").replace("-", "")
    if formula != expect:
        return None, f"formula {formula} != {expect} (mis-drawn SMILES, fix it)"
    molH = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(molH, randomSeed=7) != 0:
        AllChem.EmbedMolecule(molH, randomSeed=7, useRandomCoords=True)
    try:
        AllChem.MMFFOptimizeMolecule(molH, maxIters=4000)
    except Exception:
        pass
    conf = molH.GetConformer()
    syms = [a.GetSymbol() for a in molH.GetAtoms()]
    xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                     conf.GetAtomPosition(i).z] for i in range(molH.GetNumAtoms())])
    # symmetry: fraction of symmetry-distinct HEAVY atoms
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    heavy = mol.GetNumAtoms()
    classes = len(set(ranks))
    nrot = CalcNumRotatableBonds(mol)
    from rdkit.Chem import Descriptors
    return {"syms": syms, "xyz": xyz, "mw": Descriptors.MolWt(mol),
            "nrot": nrot, "sym_frac": classes / heavy}, None


def shape(syms, xyz):
    m = np.array([MASS.get(s, 12.0) for s in syms])
    x = xyz - (xyz * m[:, None]).sum(0) / m.sum()
    # inertia tensor
    I = np.zeros((3, 3))
    for mi, xi in zip(m, x):
        I += mi * (np.dot(xi, xi) * np.eye(3) - np.outer(xi, xi))
    ev = np.sort(np.linalg.eigvalsh(I))
    npr1, npr2 = ev[0] / ev[2], ev[1] / ev[2]
    # planarity: RMS out-of-plane thickness
    _, sv, _ = np.linalg.svd(x - x.mean(0))
    thick = sv[2] / np.sqrt(len(syms))
    return npr1 + npr2, thick


def risk_score(threeD, nrot, sym_frac, mw):
    """0 (glassy, safe) .. 100 (HATCN-like crystalliser). Heuristic; weights are
    transparent and the ranking is validated on knowns before use."""
    flat = np.clip((1.20 - threeD) / 0.20, 0, 1)         # 1.0 = disc/rod limit
    rigid = np.clip((6 - nrot) / 6, 0, 1)                # 0 rot bonds = rigid
    sym = np.clip((0.55 - sym_frac) / 0.45, 0, 1)        # few classes = symmetric
    small = np.clip((750 - mw) / 500, 0, 1)              # small = fast nucleation
    return 100 * (0.40 * flat + 0.20 * rigid + 0.25 * sym + 0.15 * small)


# ------------------------------------------------------- descriptor screen
def descriptor_screen():
    rows = []
    for name, (smi, expect, anchor, transp) in FULL.items():
        d, err = embed(smi, expect, name)
        if err:
            print(f"  [!] {name}: {err}")
            continue
        threeD, thick = shape(d["syms"], d["xyz"])
        r = risk_score(threeD, d["nrot"], d["sym_frac"], d["mw"])
        rows.append({"name": name, "mw": d["mw"], "threeD": threeD,
                     "thick": thick, "nrot": d["nrot"],
                     "sym_frac": d["sym_frac"], "risk": r,
                     "anchor": anchor, "transp": transp})
    rows.sort(key=lambda r: r["risk"])

    print(f"\n{'':24}{'MW':>6}{'3D-ness':>9}{'thick':>7}{'Nrot':>6}"
          f"{'symcls':>8}{'RISK':>6}   anchor / transparency")
    print("-" * 110)
    for r in rows:
        tg = f"  [Tg~{TG_MEMORY_C[r['name']]}C mem]" if r["name"] in TG_MEMORY_C else ""
        print(f"{r['name']:<24}{r['mw']:>6.0f}{r['threeD']:>9.2f}{r['thick']:>7.2f}"
              f"{r['nrot']:>6}{r['sym_frac']:>8.2f}{r['risk']:>6.0f}   "
              f"{r['anchor']} | {r['transp']}{tg}")

    # validation gate on the heuristic itself
    byname = {r["name"]: r for r in rows}
    ok = True
    if "HATCN" in byname and byname["HATCN"]["risk"] < max(r["risk"] for r in rows) - 1:
        print("\n  WARNING: HATCN is not the top-risk molecule; heuristic FAILED")
        print("  its known-bad anchor point. Do not use the ranking.")
        ok = False
    if "TPBi" in byname and byname["TPBi"]["risk"] > np.median([r["risk"] for r in rows]):
        print("\n  WARNING: TPBi (famously amorphous) scores above median risk;")
        print("  heuristic FAILED its known-good anchor point. Do not use it.")
        ok = False
    if ok:
        print("\n  validation: HATCN worst, TPBi in the glassy half -- heuristic holds.")
    json.dump({"rows": rows, "validated": ok},
              open(os.path.join(RUNS, "crackfree_descriptors.json"), "w"), indent=2)
    return rows, ok


# ------------------------------------------------------- fragment DFT stage
def run_xtb(tag, uhf=0):
    wd = os.path.join(RUNS, tag)
    os.makedirs(wd, exist_ok=True)
    if os.path.exists(os.path.join(wd, "xtbopt.xyz")):
        return True
    shutil.copy(os.path.join(STR, f"{tag}.xyz"), os.path.join(wd, "in.xyz"))
    with open(os.path.join(wd, "xtb.log"), "w") as log:
        subprocess.run(["xtb", "in.xyz", "--gfn", "2", "--uhf", str(uhf),
                        "--opt", "tight"],
                       cwd=wd, stdout=log, stderr=subprocess.STDOUT)
    return os.path.exists(os.path.join(wd, "xtbopt.xyz"))


def build_fragment_complexes():
    """Fragment + Ag, xtb-relaxed. Returns list of complex tags."""
    jobs = []
    for tag, (smi, expect) in FRAGMENTS.items():
        d, err = embed(smi, expect, tag)
        if err:
            print(f"  [!] {tag}: {err}")
            continue
        # Never clobber an existing structure file: PhCN.xyz predates this script
        # and is the committed input of the PhCN_Zn/PhCN_ZnS results, so silently
        # replacing it would detach those results from their input. The first run
        # of this script did exactly that; the file was restored from git.
        frag_path = os.path.join(STR, f"{tag}.xyz")
        if not os.path.exists(frag_path):
            write_xyz(frag_path, d["syms"], d["xyz"], tag)
        if not run_xtb(tag):
            print(f"  [!] {tag}: xtb failed")
            continue
        syms, x = read_xyz(os.path.join(RUNS, tag, "xtbopt.xyz"))
        cen = x.mean(0)
        ns = [i for i, s in enumerate(syms) if s == "N"
              and min(np.linalg.norm(x - x[i], axis=1)[np.arange(len(syms)) != i]) < 1.25]
        # ^ nitrile N only: its nearest neighbour is the C#N carbon at ~1.16 A.
        #   DMABN's amine N has neighbours at ~1.38 A and must not be the anchor.
        if tag == "oDCNB" and len(ns) == 2:
            # chelate guess: Ag on the bisector of the two nitrile N, 2.4 A out
            mid = 0.5 * (x[ns[0]] + x[ns[1]])
            d_out = mid - cen; d_out /= np.linalg.norm(d_out)
            ag = mid + 2.4 * d_out
            ctag = f"{tag}_Ag_chel"
            write_xyz(os.path.join(STR, f"{ctag}.xyz"), syms + ["Ag"],
                      np.vstack([x, ag]), "chelate guess")
            if run_xtb(ctag, uhf=1):
                jobs.append(ctag)
        i_n = max(ns, key=lambda k: np.linalg.norm(x[k] - cen))
        dd = np.linalg.norm(x - x[i_n], axis=1); dd[i_n] = 9e9
        i_c = int(np.argmin(dd))
        ax = x[i_n] - x[i_c]; ax /= np.linalg.norm(ax)
        ctag = f"{tag}_Ag"
        write_xyz(os.path.join(STR, f"{ctag}.xyz"), syms + ["Ag"],
                  np.vstack([x, x[i_n] + 2.30 * ax]), "end-on guess")
        if run_xtb(ctag, uhf=1):
            jobs.append(ctag)
    return jobs


def psi4_stage(jobs):
    """scripts/30's scan+CP, unchanged, on the xtb-relaxed complexes."""
    import psi4
    psi4.set_memory(os.environ.get("PSI4_MEM", "6 GB"))
    psi4.set_num_threads(int(os.environ.get("PSI4_THREADS", "4")))
    psi4.core.set_output_file(os.path.join(RUNS, "psi4_crackfree.out"), False)
    base = {"basis": "def2-svp", "scf_type": "df", "reference": "uks",
            "maxiter": 200, "guess": "sad"}
    psi4.set_options(base)
    M = "pbe-d3bj"

    def gstr(syms, xyz, ghost=None, mult=1):
        s = f"0 {mult}\n"
        for i, (sym, c) in enumerate(zip(syms, xyz)):
            t = f"Gh({sym})" if ghost and i in ghost else sym
            s += f"{t} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n"
        return s + "symmetry c1\nno_reorient\nno_com\n"

    def energy(g):
        for extra in ({}, {"damping_percentage": 20, "level_shift": 0.5,
                           "level_shift_cutoff": 1e-4}, {"soscf": True}):
            try:
                psi4.set_options({**base, **extra})
                e = psi4.energy(M, molecule=psi4.geometry(g))
                psi4.set_options(base)
                return e
            except Exception as ex:
                print(f"    scf retry ({type(ex).__name__})", flush=True)
                psi4.core.clean()
        return None

    res_path = os.path.join(RUNS, "crackfree_fragments_eV.json")
    res = json.load(open(res_path)) if os.path.exists(res_path) else {}
    for ctag in jobs:
        if ctag in res and "Eb_eV" in res.get(ctag, {}):
            print(f"[skip] {ctag}: {res[ctag]['Eb_eV']:.3f} eV", flush=True)
            continue
        syms, x = read_xyz(os.path.join(RUNS, ctag, "xtbopt.xyz"))
        agi = len(syms) - 1
        heavy = [i for i in range(agi) if syms[i] != "H"]
        anchor = min(heavy, key=lambda i: np.linalg.norm(x[agi] - x[i]))
        ax = x[agi] - x[anchor]
        r0 = np.linalg.norm(ax); ax /= r0
        print(f"[{ctag}] anchor {syms[anchor]}{anchor}, r0 = {r0:.2f}", flush=True)
        scan, best = {}, (None, None)
        for dt in (-0.2, -0.1, 0.0, 0.1, 0.25, 0.45):
            xs = x.copy(); xs[agi] = x[anchor] + (r0 + dt) * ax
            e = energy(gstr(syms, xs, None, 2))
            if e is None:
                continue
            scan[round(r0 + dt, 3)] = e
            if best[1] is None or e < best[1]:
                best = (r0 + dt, e)
            print(f"    r = {r0+dt:.2f}  E = {e:.6f}", flush=True)
        if best[1] is None:
            res[ctag] = {"error": "scan_failed"}
            json.dump(res, open(res_path, "w"), indent=2)
            continue
        xs = x.copy(); xs[agi] = x[anchor] + best[0] * ax
        e_sub = energy(gstr(syms, xs, {agi}, 1))
        e_ag = energy(gstr(syms, xs, set(range(agi)), 2))
        if None in (e_sub, e_ag):
            res[ctag] = {"error": "cp_failed"}
        else:
            res[ctag] = {"Eb_eV": (e_sub + e_ag - best[1]) * H2EV, "r_A": best[0]}
            print(f"    E_b = {res[ctag]['Eb_eV']:.3f} eV", flush=True)
        json.dump(res, open(res_path, "w"), indent=2)

    print("\n=== fragment summary ===")
    print(f"{'fragment':<16}{'E_b (eV)':>9}   meaning")
    notes = {
        "PhCN_Ag": "single CN, neutral ring -- the baseline",
        "pDCNB_Ag": "acceptor-substituted ring: does EA pull raise E_b?",
        "DMABN_Ag": "donor-substituted ring: the 4CzIPN direction",
        "oDCNB_Ag": "o-dinitrile, end-on",
        "oDCNB_Ag_chel": "o-dinitrile, chelate -- in-molecule pocket?",
    }
    for k in ("PhCN_Ag", "pDCNB_Ag", "DMABN_Ag", "oDCNB_Ag", "oDCNB_Ag_chel"):
        v = res.get(k, {})
        s = f"{v['Eb_eV']:>9.3f}" if "Eb_eV" in v else f"{v.get('error','--'):>9}"
        print(f"{k:<16}{s}   {notes.get(k,'')}")
    print(f"\n  reference: HATCN cluster 1.029, F4TCNQ 0.966 (same protocol)")


if __name__ == "__main__":
    print("=" * 110)
    print("PART 1 -- crystallisation-risk descriptors (full molecules)")
    print("=" * 110)
    rows, ok = descriptor_screen()
    if "--descriptors-only" in sys.argv:
        sys.exit(0)
    print("\n" + "=" * 110)
    print("PART 2 -- anchor-fragment DFT (substituent + chelate effects)")
    print("=" * 110)
    jobs = build_fragment_complexes()
    print(f"fragment complexes ready: {jobs}", flush=True)
    psi4_stage(jobs)
