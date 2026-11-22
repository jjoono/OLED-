"""Two-tier seed descriptor: D = (anchor-site areal density) x P_trap(E_b).

    python scripts/174_two_tier_descriptor.py   ->  runs/two_tier.json

Why (notes/KONICA_TABLE2_RESCORE.md): on the ten N-heterocycle seeds of
Konica US 9,923,166 (Example 2, Ag 5 nm), our PBE0 E_b ranks the sheet
resistance at rho = -0.46 and the patent's own Ag-N energy dE at rho = 0.07,
while the areal density of lone-pair N ranks it at rho = 0.89. Within one
binding class, a site above the trap threshold traps; binding harder adds
nothing, and how many such sites face the vapour is what sets the nucleus
density. So:

    D = sum over anchor sites of P_trap(E_b) / SASA
    P_trap(E_b) = 1 - Phi((E* - E_b) / sigma),   sigma = 0.10 eV

Anchor sites: N with a lone pair (2-coordinate ring N, nitrile N) and the O
of P=O. pi sites are not anchors: amine HTLs and anthracene do not seed Ag
(Konica No01 unmeasurable; lit. Ag on NPB/amine HTL closes at > 7 nm) while
their pi binding is 0.13-0.28 eV. E_b per site is the molecule's deepest
PBE0 value (screen level, PBE0-D3/def2-SVP); hindered sites are already low
in it (BCP, NBphen, T2T).

E* (screen level) is fixed from external data before any of our own
experiments:
  upper bound: Konica pyridyl sites (0.35-0.50 eV) must all trap: D ranks
               the patent at rho >= 0.86 for E* <= 0.35 and drops to 0.47-0.75
               at the kinetic estimate 0.465 (script 165);
  lower bound: pi sites (<= 0.28 eV) must not trap.
  -> E* = 0.30 eV, sensitivity 0.25 / 0.35.

Applied to: (1) the Konica set (validation), (2) the commercial set (the
pre-registered prediction), (3) the HAT design library (re-rank with D in
place of E_b as the nucleation objective; E_b becomes a gate).
"""
import json
import os
import sys

import numpy as np
from scipy.stats import norm, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "runs", "two_tier.json")
E_STAR, SIG = 0.30, 0.10
E_STAR_SENS = (0.25, 0.35)


def p_trap(eb, e_star=E_STAR, sig=SIG):
    return float(1 - norm.cdf(e_star, loc=eb, scale=sig))


def anchors_sasa(mol):
    """(anchor count, SASA A^2) of an RDKit molecule with H and 3D coordinates."""
    from rdkit.Chem import rdFreeSASA
    n = 0
    X = mol.GetConformer().GetPositions()
    metals = [i for i, a in enumerate(mol.GetAtoms()) if a.GetSymbol() in ("Al", "Li", "Cs", "Ag", "B")]
    for a in mol.GetAtoms():
        s, deg = a.GetSymbol(), a.GetDegree()
        # N already coordinated to a metal (Alq3's Al-N, not in the bond table) is not free
        if any(np.linalg.norm(X[a.GetIdx()] - X[m]) < 2.4 for m in metals):
            continue
        if s == "N" and deg <= 2:
            n += 1
        elif s == "O" and deg == 1 and a.GetNeighbors()[0].GetSymbol() == "P":
            n += 1
    sasa = rdFreeSASA.CalcSASA(mol, rdFreeSASA.classifyAtoms(mol))
    return n, float(sasa)


def mol_from_file(path):
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    if path.endswith(".mol"):
        return Chem.MolFromMolFile(path, removeHs=False)
    m = Chem.MolFromXYZFile(path)
    rdDetermineBonds.DetermineConnectivity(m)
    return m


def mol_from_smiles(smi, seed=7):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    if AllChem.EmbedMolecule(m, randomSeed=seed) != 0:
        p = AllChem.ETKDGv3()
        p.useRandomCoords = True
        AllChem.EmbedMolecule(m, p)
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except Exception:
        pass
    return m


def D(n, sasa, eb, e_star=E_STAR):
    return 1e3 * n * p_trap(eb, e_star) / sasa       # per 1000 A^2


def konica():
    s = json.load(open(os.path.join(ROOT, "runs", "konica", "summary.json")))["rows"]
    rows, y = {}, []
    ks = [k for k in s if k != "No01"]
    for k in ks:
        m = mol_from_file(os.path.join(ROOT, "structures", "konica", f"{k}.mol"))
        n, sa = anchors_sasa(m)
        rs = s[k]["patent"]["Rs"]
        rows[k] = {"n_anchor": n, "SASA": round(sa, 1), "E_b": s[k]["E_b"], "Rs": rs,
                   "D": round(D(n, sa, s[k]["E_b"]), 3), "patent_dEef": s[k]["patent"]["dEef"]}
        y.append(-(np.log10(rs) if rs else 9.0))
    rho = lambda xs: round(float(spearmanr(xs, y).correlation), 3)  # noqa: E731
    scan = {f"{e:.2f}": rho([D(rows[k]["n_anchor"], rows[k]["SASA"], rows[k]["E_b"], e) for k in ks])
            for e in np.arange(0.15, 0.56, 0.05)}
    return {"rows": rows, "rho_D": rho([rows[k]["D"] for k in ks]),
            "rho_E_b": rho([rows[k]["E_b"] for k in ks]),
            "rho_patent_dEef": rho([-rows[k]["patent_dEef"] for k in ks]),
            "rho_D_vs_E_star": scan,
            "note": "Rs 'unmeasurable' (No02) ranked worst; Spearman vs -log10 Rs"}


COMMERCIAL_FILES = {
    "HATCN": "HATCN.xyz", "Bphen": "Bphen.xyz", "B3PyMPM": "B3PyMPM.xyz", "p-bPPhenB": "pbPPhenB.xyz",
}


def commercial():
    rank = {r["name"]: r for r in json.load(open(os.path.join(ROOT, "runs", "common_ranking.json")))}
    out = {}
    for name, r in rank.items():
        f = COMMERCIAL_FILES.get(name)
        cands = ([os.path.join(ROOT, "structures", f)] if f else []) + [
            os.path.join(ROOT, "structures", "common", f"{name}.mol")]
        path = next((p for p in cands if os.path.exists(p)), None)
        if path is None or r.get("E_b") is None:
            continue
        n, sa = anchors_sasa(mol_from_file(path))
        out[name] = {"role": r.get("role"), "E_b": r["E_b"], "site": r.get("site"),
                     "n_anchor": n, "SASA": round(sa, 1),
                     "P_trap": round(p_trap(r["E_b"]), 3),
                     "D": round(D(n, sa, r["E_b"]), 3),
                     "D_E25": round(D(n, sa, r["E_b"], E_STAR_SENS[0]), 3),
                     "D_E35": round(D(n, sa, r["E_b"], E_STAR_SENS[1]), 3)}
    h = out["HATCN"]["D"]
    for v in out.values():
        v["D_rel_HATCN"] = round(v["D"] / h, 3)
    return out


def hat_library():
    """Re-rank the design library: D replaces E_b as the nucleation objective."""
    sys.path.insert(0, HERE)
    s167 = __import__("importlib").import_module("167_design_pareto")
    P = json.load(open(os.path.join(ROOT, "runs", "design", "pareto.json")))
    lib = json.load(open(os.path.join(ROOT, "runs", "design", "library.json")))
    cache_p = os.path.join(ROOT, "runs", "design", "anchor_sasa.json")
    cache = json.load(open(cache_p)) if os.path.exists(cache_p) else {}
    rows = list(P["pass"]) + list(P["HATCN"])
    for r in rows:
        n = r["name"]
        if n not in cache:
            m = mol_from_smiles(lib[n]["smiles"])
            cache[n] = anchors_sasa(m)
    json.dump(cache, open(cache_p, "w"))
    for r in rows:
        na, sa = cache[r["name"]]
        r["n_anchor"], r["SASA"] = na, round(sa, 1)
        r["D"] = round(D(na, sa, r["E_b_safe"]), 3)
    hat = [r for r in rows if r["name"] == "HATCN"][0]
    ok = [r for r in rows if r["name"] != "HATCN"]
    F = np.array([[r["D"], -r["risk"], r["gap"], -(r["synth_class"] - 1)] for r in ok])
    fr = s167.fronts(F)
    for r, k in zip(ok, fr):
        r["front_D"] = int(k)
    ok.sort(key=lambda r: (r["front_D"], -r["D"]))
    beats = [r["name"] for r in ok if r["D"] >= hat["D"] and r["gap"] >= hat["gap"]]
    keep = ("name", "E_b", "E_b_safe", "n_anchor", "SASA", "D", "gap", "risk", "MW",
            "synth_class", "pareto_front", "front_D")
    return {"HATCN": {k: hat.get(k) for k in keep},
            "n_ranked": len(ok),
            "D_max_candidate": max(r["D"] for r in ok),
            "n_with_D_ge_HATCN": sum(r["D"] >= hat["D"] for r in ok),
            "beat_HATCN_on_D_and_gap": beats,
            "front0": [{k: r.get(k) for k in keep} for r in ok if r["front_D"] == 0],
            "top20_by_D": [{k: r.get(k) for k in keep} for r in sorted(ok, key=lambda r: -r["D"])[:20]]}


def main():
    out = {"params": {"E_star_eV": E_STAR, "sigma_eV": SIG, "E_star_sensitivity": E_STAR_SENS,
                      "anchors": "N with <= 2 neighbours, O of P=O; pi sites excluded",
                      "D_units": "trapping anchors per 1000 A^2 SASA"}}
    out["konica"] = konica()
    k = out["konica"]
    print(f"Konica: rho(D) {k['rho_D']}  rho(E_b) {k['rho_E_b']}  rho(dEef) {k['rho_patent_dEef']}")
    print("  rho(D) vs E*:", k["rho_D_vs_E_star"])
    out["commercial"] = commercial()
    print(f"\n{'material':<10} {'E_b':>5} {'Ptrap':>5} {'n':>3} {'SASA':>6} {'D':>6} {'D/HATCN':>7} "
          f"{'D(E*.25)':>8} {'D(E*.35)':>8}")
    for n, v in sorted(out["commercial"].items(), key=lambda kv: -kv[1]["D"]):
        print(f"{n:<10} {v['E_b']:5.2f} {v['P_trap']:5.2f} {v['n_anchor']:3d} {v['SASA']:6.0f} "
              f"{v['D']:6.2f} {v['D_rel_HATCN']:7.3f} {v['D_E25']:8.2f} {v['D_E35']:8.2f}")
    out["hat_library"] = hat_library()
    h = out["hat_library"]
    print(f"\nHAT library: {h['n_ranked']} ranked; HATCN D {h['HATCN']['D']}; best candidate D "
          f"{h['D_max_candidate']}; {h['n_with_D_ge_HATCN']} with D >= HATCN; "
          f"beat HATCN on D and gap: {h['beat_HATCN_on_D_and_gap']}")
    for r in h["top20_by_D"][:10]:
        print(f"  {r['name']:<34} D {r['D']:.2f} n {r['n_anchor']} SASA {r['SASA']:.0f} "
              f"Eb_safe {r['E_b_safe']:.2f} gap {r['gap']:.2f} MW {r['MW']}")
    json.dump(out, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
