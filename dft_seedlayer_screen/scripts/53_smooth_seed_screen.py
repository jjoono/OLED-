"""Smooth, amorphous, and still binds Ag: screen 3D acceptor architectures.

WHERE THIS COMES FROM. Four results constrain the design:

  1. The cracking is crystallisation, not fracture (channel-crack h_c is in
     micrometres for any sane parameter). HATCN is small, rigid, perfectly planar
     and highly symmetric -- the worst case for glass formation, risk score 92.
  2. An intramolecular chelate does NOT work: the o-dinitrile bridge binds Ag at
     0.130 eV against 0.292 eV for a single end-on nitrile (scripts/51). Binding
     cannot be bought by putting two nitriles next to each other.
  3. Donor substitution is FREE on the anchor axis: DMABN 0.287 vs benzonitrile
     0.292 eV. Twisting a molecule with carbazole-type donors costs no binding.
  4. Acceptor substitution HELPS: terephthalonitrile 0.379 eV, +0.087 over
     benzonitrile. Binding tracks how strongly the ring pulls charge off Ag.

So binding has to come from nitrile COUNT and from an electron-poor core, while
the SHAPE has to stop the molecule packing. Those pull against each other, because
what makes a good acceptor is an extended flat pi system, and flat pi systems
crystallise. The way out is an sp3 centre holding intact pi units in
non-coplanar orientations -- spiro, tetrahedral, propeller -- which keeps each
acceptor unit's electronics while destroying the packing.

TWO PARTS.

  PART A validates a cheap proxy. Running E_b with an Ag atom is expensive and,
  for the strongest acceptors, does not converge at all. If E_b tracks the
  molecule's electron affinity across the fragment series already computed, EA
  becomes a screening stand-in and the whole candidate list can be ranked without
  ever placing an Ag atom. That is a REAL test on four points spanning donor to
  acceptor substitution, and it is allowed to fail -- Fig. 1(b) asserts this
  correlation on two points, which is not enough to assert anything.

  PART B scores an expanded candidate list on crystallisation risk, adding the
  3D-acceptor architectures that parts 1-4 point to and that the earlier list
  lacked: a tetrahedral sp3 core carrying four benzonitrile arms, its silicon
  analogue, a spiro-linked dinitrile, a propeller borane, and an sp3
  phosphine oxide.

WHAT IS NOT DONE HERE. No full-molecule E_b (the container cannot), no monolayer
pocket correction, no Tg. Sublimability is not assessed and matters: a tetrahedral
tetranitrile at MW ~460 should evaporate, but that is chemistry judgement, not a
computed result, and it is flagged rather than assumed.
"""
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _ckpt import Checkpoint

BASE = os.path.abspath(os.path.join(HERE, ".."))
STR, RUNS = os.path.join(BASE, "structures"), os.path.join(BASE, "runs")
HA = 27.211386

MASS = {"H": 1.008, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999,
        "F": 18.998, "Si": 28.086, "P": 30.974, "S": 32.06}

# ---------------------------------------------------------------- PART A
# Fragments whose E_b is already known at PBE-D3BJ/def2-SVP + counterpoise.
# If EA correlates with these, EA can screen the big molecules.
FRAG_EB = {
    "DMABN": 0.287,     # donor-substituted ring
    "PhCN": 0.292,      # neutral reference
    "pDCNB": 0.379,     # acceptor-substituted ring
    "oDCNB": 0.130,     # o-dinitrile, but measured as a CHELATE -- different
                        # binding mode, so it is plotted and excluded from the fit
}
FIT_EXCLUDE = {"oDCNB"}

# ---------------------------------------------------------------- PART B
# name: (SMILES, expected formula, anchors, note)
CANDIDATES = {
    # --- the planar nitrile incumbents, as the known-bad anchor points ---
    "HATCN": ("N#Cc1nc2c3nc(C#N)c(C#N)nc3c3nc(C#N)c(C#N)nc3c2nc1C#N",
              "C18N12", "6x CN", "incumbent; risk anchor point"),
    "F4TCNQ": ("N#CC(C#N)=C1C(F)=C(F)C(=C(C#N)C#N)C(F)=C1F",
               "C12F4N4", "4x CN", "planar"),

    # --- sp3 centre, pi units forced non-coplanar: the design target ---
    "TCPM": ("N#Cc1ccc(cc1)C(c1ccc(C#N)cc1)(c1ccc(C#N)cc1)c1ccc(C#N)cc1",
             "C29H16N4", "4x CN, tetrahedral",
             "tetraphenylmethane core; 4 nitriles, maximally non-planar"),
    "TCPSi": ("N#Cc1ccc(cc1)[Si](c1ccc(C#N)cc1)(c1ccc(C#N)cc1)c1ccc(C#N)cc1",
              "C28H16N4Si", "4x CN, tetrahedral",
              "Si analogue; longer bonds, even less able to stack"),
    "SBF2CN": ("N#Cc1ccc2c(c1)-c1ccccc1C21c2ccccc2-c2ccc(C#N)cc21",
               "C27H14N2", "2x CN, spiro",
               "spirobifluorene: two orthogonal pi units, classic glass former"),
    "TCPB": ("N#Cc1ccc(cc1)-c1cc(-c2ccc(C#N)cc2)cc(-c2ccc(C#N)cc2)c1",
             "C27H15N3", "3x CN, star",
             "planar-ish star -- control for the sp3 cores"),

    # --- 3D non-nitrile acceptors already in the project's vocabulary ---
    "TPBi": ("c1ccc(cc1)n1c(nc2ccccc21)-c1cc(-c2nc3ccccc3n2-c2ccccc2)cc"
             "(-c2nc3ccccc3n2-c2ccccc2)c1",
             "C45H30N6", "3x benzimidazole N", "glassy anchor point"),
    "DPEPO": ("O=P(c1ccccc1)(c1ccccc1)c1ccccc1Oc1ccccc1P(=O)(c1ccccc1)c1ccccc1",
              "C36H28O3P2", "2x P=O", "sp3 P centres, very 3D"),
    "mCPCN": ("N#Cc1ccc2c(c1)c1ccccc1n2-c1cccc(-n2c3ccccc3c3ccccc32)c1",
              "C31H19N3", "1x CN", "twisted carbazole-nitrile"),
    "TCTA": ("c1ccc2c(c1)c1ccccc1n2-c1ccc(cc1)N(c1ccc(cc1)-n1c2ccccc2c2ccccc21)"
             "c1ccc(cc1)-n1c2ccccc2c2ccccc21",
             "C54H36N4", "amine N only", "3D propeller, WEAK acceptor -- control"),
}


def write_xyz(path, syms, xyz, cm=""):
    with open(path, "w") as f:
        f.write(f"{len(syms)}\n{cm}\n")
        for s, c in zip(syms, xyz):
            f.write(f"{s} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n")


def embed(smiles, expect):
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcNumRotatableBonds
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "SMILES parse failed"
    formula = CalcMolFormula(mol).replace("+", "").replace("-", "")
    if formula != expect:
        return None, f"formula {formula} != {expect}"
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
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    return {"syms": syms, "xyz": xyz, "mw": Descriptors.MolWt(mol),
            "nrot": CalcNumRotatableBonds(mol),
            "sym_frac": len(set(ranks)) / mol.GetNumAtoms()}, None


def shape(syms, xyz):
    m = np.array([MASS.get(s, 12.0) for s in syms])
    x = xyz - (xyz * m[:, None]).sum(0) / m.sum()
    I = np.zeros((3, 3))
    for mi, xi in zip(m, x):
        I += mi * (np.dot(xi, xi) * np.eye(3) - np.outer(xi, xi))
    ev = np.sort(np.linalg.eigvalsh(I))
    thick = np.linalg.svd(x - x.mean(0))[1][2] / np.sqrt(len(syms))
    return ev[0] / ev[2] + ev[1] / ev[2], thick


def risk(threeD, nrot, sym_frac, mw):
    flat = np.clip((1.20 - threeD) / 0.20, 0, 1)
    rigid = np.clip((6 - nrot) / 6, 0, 1)
    sym = np.clip((0.55 - sym_frac) / 0.45, 0, 1)
    small = np.clip((750 - mw) / 500, 0, 1)
    return 100 * (0.40 * flat + 0.20 * rigid + 0.25 * sym + 0.15 * small)


# ============================================================= PART A
def part_a():
    import psi4
    psi4.set_memory(os.environ.get("PSI4_MEM", "5 GB"))
    psi4.set_num_threads(int(os.environ.get("PSI4_THREADS", "4")))
    psi4.core.set_output_file(os.path.join(RUNS, "psi4_proxy.out"), False)
    ck = Checkpoint("dft_seedlayer_screen/runs/eb_ea_proxy.json", label="proxy EA")

    def sp(fname, charge, mult):
        psi4.core.clean()
        lines = open(os.path.join(RUNS, fname, "xtbopt.xyz")).read().splitlines()
        n = int(lines[0])
        body = "\n".join(lines[2:2 + n])
        mol = psi4.geometry(f"{charge} {mult}\n{body}\n"
                            "symmetry c1\nno_reorient\nno_com\n")
        psi4.set_options({"basis": "def2-svp", "scf_type": "df",
                          "reference": "uks" if mult > 1 else "rks",
                          "maxiter": 300, "guess": "sad",
                          # threshold set from the measured DF noise floor, as in
                          # scripts/52 -- 1e-8 cannot be reached and hangs the run
                          "d_convergence": 1e-7, "e_convergence": 1e-5,
                          "damping_percentage": 20.0})
        return psi4.energy("wb97x", molecule=mol)

    ea = {}
    for tag in FRAG_EB:
        for suffix, chg, mult in (("neutral", 0, 1), ("anion", -1, 2)):
            key = f"{tag}_{suffix}"
            if ck.has(key):
                continue
            print(f"  {tag:<7} {suffix:<8} running...", flush=True)
            ck.put(key, sp(tag, chg, mult))
        e0, em = ck.get(f"{tag}_neutral"), ck.get(f"{tag}_anion")
        if e0 is None or em is None:
            print(f"  {tag}: incomplete"); return None
        ea[tag] = (e0 - em) * HA
        print(f"  {tag:<7} EA = {ea[tag]:>6.3f} eV,  E_b = {FRAG_EB[tag]:.3f} eV",
              flush=True)

    fit = [(ea[t], FRAG_EB[t]) for t in ea if t not in FIT_EXCLUDE]
    print("\n  " + "-" * 56)
    print(f"  {'fragment':<9}{'EA (eV)':>9}{'E_b (eV)':>10}   note")
    for t in sorted(ea, key=lambda k: ea[k]):
        note = "CHELATE -- different binding mode, excluded from fit" \
            if t in FIT_EXCLUDE else ""
        print(f"  {t:<9}{ea[t]:>9.3f}{FRAG_EB[t]:>10.3f}   {note}")

    ok = False
    if len(fit) >= 3:
        xs = np.array([f[0] for f in fit]); ys = np.array([f[1] for f in fit])
        r = float(np.corrcoef(xs, ys)[0, 1])
        slope = float(np.polyfit(xs, ys, 1)[0])
        print(f"\n  correlation over {len(fit)} end-on fragments: r = {r:.3f}, "
              f"slope = {slope:.3f} eV(E_b)/eV(EA)")
        ok = r > 0.9
        if ok:
            print("  -> EA tracks E_b across donor-to-acceptor substitution.")
            print("     EA is usable as a screening proxy for anchor strength,")
            print("     which is what makes the candidate list rankable at all:")
            print("     EA needs no Ag atom, and it converges where E_b does not.")
        else:
            print("  -> EA does NOT track E_b on this series. The proxy fails, and")
            print("     with it the Fig. 1(b) claim that the two axes are coupled.")
            print("     Rank the candidates on shape only, and say so.")
    ck.put("EA_eV", ea); ck.put("proxy_valid", ok)
    return ea, ok


# ============================================================= PART B
def part_b():
    rows = []
    for name, (smi, expect, anchor, note) in CANDIDATES.items():
        d, err = embed(smi, expect)
        if err:
            print(f"  [!] {name}: {err}")
            continue
        t3, thick = shape(d["syms"], d["xyz"])
        write_xyz(os.path.join(STR, f"{name}.xyz"), d["syms"], d["xyz"], name) \
            if not os.path.exists(os.path.join(STR, f"{name}.xyz")) else None
        rows.append({"name": name, "mw": d["mw"], "threeD": t3, "thick": thick,
                     "nrot": d["nrot"], "sym_frac": d["sym_frac"],
                     "risk": risk(t3, d["nrot"], d["sym_frac"], d["mw"]),
                     "anchor": anchor, "note": note})
    rows.sort(key=lambda r: r["risk"])
    print(f"\n{'':10}{'MW':>6}{'3D':>7}{'thick':>7}{'Nrot':>6}{'RISK':>6}"
          f"   anchors | note")
    print("-" * 104)
    for r in rows:
        print(f"{r['name']:<10}{r['mw']:>6.0f}{r['threeD']:>7.2f}{r['thick']:>7.2f}"
              f"{r['nrot']:>6}{r['risk']:>6.0f}   {r['anchor']} | {r['note']}")
    hat = next((r for r in rows if r["name"] == "HATCN"), None)
    if hat and hat["risk"] < max(r["risk"] for r in rows) - 1:
        print("\n  WARNING: HATCN is not top-risk; heuristic failed its anchor point")
    json.dump(rows, open(os.path.join(RUNS, "smooth_seed_descriptors.json"), "w"),
              indent=2)
    return rows


if __name__ == "__main__":
    print("=" * 104)
    print("PART B -- crystallisation risk, expanded to 3D acceptor architectures")
    print("=" * 104)
    rows = part_b()
    if "--descriptors-only" in sys.argv:
        sys.exit(0)
    print("\n" + "=" * 104)
    print("PART A -- does EA predict E_b? (validating a cheap screening proxy)")
    print("=" * 104)
    out = part_a()
    if out and rows:
        ea, ok = out
        if ok:
            print("\n" + "=" * 104)
            print("COMBINED: low crystallisation risk AND many strong-acceptor anchors")
            print("=" * 104)
            for r in rows[:5]:
                print(f"  {r['name']:<10} risk {r['risk']:>3.0f}   {r['anchor']}")
            print("\n  Next step is EA on these full molecules -- no Ag atom needed,")
            print("  so it converges and is affordable, unlike E_b.")
