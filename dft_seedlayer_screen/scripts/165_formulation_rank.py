"""Step 3: rank seed FORMULATIONS (blends, doping) for the CGL p-side seed.

    python scripts/165_formulation_rank.py [--kmc]  ->  runs/design/formulations.json

Same role and thresholds as notes/DESIGN_SPEC_CGL_SEED.md. A formulation is
judged on five things, each computed, none assumed:

1. Nucleation (kMC). Script 68's trap-map kMC (scripts/36 SOS model,
   Venables-validated) gives island density N_sat vs partner fraction x and
   the partner's effective trap fraction theta_eff. What changed since script
   68: its partner E_b were single-site numbers (TPBi 0.89, NBphen 0.87). The
   multi-site PBE0 ranking (script 150) gives 0.46 for both, which puts them
   right on the trap/diluent boundary (tau ~ t_capture). So a partner is now a
   trap with a probability:
        p_trap = P[ tau(E_b) > t_capture ],  E_b ~ N(E_b,PBE0, 0.10 eV)
   and theta_eff = theta_site * p_trap. N_sat is interpolated in theta_eff
   from the kMC grid (x = 0.5 from script 68; x = 1, neat partner, run here
   with --kmc).
2. Crystallisation kill: HATCN-domain spanning dies at x >= 0.45 (script 68).
   The partner must not crystallise itself (its own risk score, script 53).
3. CT absorption with HATCN: a donor partner with a low IP forms a
   charge-transfer complex with HATCN that absorbs in the visible/NIR
   (NPB:HATCN, TAPC:HATCN are documented CT blends; TPBi/Bphen-class ETLs
   with HATCN are not). An absolute CT energy from gas-phase IP/EA needs a
   solid-state correction we cannot compute reliably, so the rule is
   calibrated on those known cases instead. IPEA1-xTB vertical IPs split
   cleanly: the amine donors (TAPC, NPB, TPD, TCTA) at 6.29-6.31 eV, Alq3/
   CBP/mCP at 6.7-7.0, every ETL at >= 7.08. So: IP < 6.5 -> CT (as for the
   documented blends); 6.5-7.0 -> weak CT possible; >= 7.0 -> none.
4. Ag charge veto: unchanged by blending (every partner passes it).
5. Electrical (CGL p-side): HATCN must stay connected through the film
   (3D site percolation, p_c ~ 0.31 for random mixing) -> x <= 0.6. So the
   window is 0.45 <= x <= 0.6: enough partner to kill 2D domains, little
   enough to keep HATCN percolating.

Bphen:Cs2CO3 as a stand-alone seed (no HATCN) is ranked in the same frame:
the chelate (0.62 eV) and the Cs sites (0.89 eV, script 155) are its traps.
It fails the electrical role (it is the n-side) but is the reference for
"no p-seed" stacks.
"""
import importlib
import json
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402
from scipy.stats import norm  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
D = os.path.join(ROOT, "runs", "design")
OUT = os.path.join(D, "formulations.json")
PCACHE = os.path.join(D, "partner_props.json")
KMC1 = os.path.join(D, "kmc_neat_partner.json")

KB, NU0, T = 8.617333262e-5, 1.0e13, 300.0
E_D, E_FACE, T_CAP = 0.29, 0.10, 0.01
SIG_EB = 0.10
X_KILL, X_PERC = 0.45, 0.60
EA_HATCN = 3.693

# fraction of a partner's surface sites that are anchors (script 68's values;
# chelate/pyridyl/P=O exposure estimates)
THETA_SITE = {"TPBi": 0.50, "NBphen": 0.35, "Bphen": 0.35, "BCP": 0.30, "p-bPPhenB": 0.35,
              "B3PyMPM": 0.40, "TmPyPB": 0.40, "3TPYMB": 0.40, "T2T": 0.30, "PO-T2T": 0.40,
              "DPEPO": 0.30, "TSPO1": 0.30, "Alq3": 0.30, "Liq": 0.40,
              "CBP": 0.15, "TCTA": 0.15, "mCP": 0.15, "TAPC": 0.15, "NPB": 0.15, "TPD": 0.15}
SMILES_EXTRA = {
    "p-bPPhenB": "c1ccc(cc1)-c1cc(-c2ccc(cc2)-c2cc(-c3ccccc3)nc3c2ccc2cccnc23)c2ccc3cccnc3c2n1",
    "B3PyMPM": "Cc1nc(-c2cc(-c3cccnc3)cc(-c3cccnc3)c2)cc(-c2cc(-c3cccnc3)cc(-c3cccnc3)c2)n1",
    "Liq": "[Li]Oc1cccc2cccnc12",
}


def tau(eb):
    return np.exp((E_D + eb - E_FACE) / (KB * T)) / NU0


def p_trap(eb):
    # tau > t_cap  <=>  E_b > E* ;  E* = kT ln(t_cap nu0) - E_D + E_FACE
    e_star = KB * T * np.log(T_CAP * NU0) - E_D + E_FACE
    return float(1 - norm.cdf(e_star, loc=eb, scale=SIG_EB)), float(e_star)


def partner_props(names):
    cache = json.load(open(PCACHE)) if os.path.exists(PCACHE) else {}
    s161 = importlib.import_module("161_design_props")
    from pathgeom import read_xyz
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    props = {}
    pj = os.path.join(D, "props.jsonl")
    if os.path.exists(pj):
        props = {json.loads(l)["name"]: json.loads(l) for l in open(pj) if l.strip()}
    for n in names:
        if n in cache:
            continue
        if n in props:
            p = props[n]
            cache[n] = {"IP": p["IP"], "EA": p["EA"], "risk": p["cryst_risk"]}
        elif n in SMILES_EXTRA:
            m = Chem.MolFromSmiles(SMILES_EXTRA[n])
            rk = list(Chem.CanonicalRankAtoms(m, breakTies=False))
            d = {"smiles": SMILES_EXTRA[n], "nrot": rdMolDescriptors.CalcNumRotatableBonds(m),
                 "sym_frac": len(set(rk)) / m.GetNumHeavyAtoms(), "MW": Descriptors.MolWt(m)}
            r = s161.work((n, d))
            cache[n] = {"IP": r["IP"], "EA": r["EA"], "risk": r["cryst_risk"]}
        else:
            base = os.path.join(ROOT, "structures", "common", n)
            m = Chem.MolFromMolFile(base + ".mol", removeHs=False)
            S, X = read_xyz(base + ".xyz")
            X = np.asarray(X, float)
            from ase import Atoms
            Z = np.array(Atoms(S).numbers)
            e0 = s161.sp("IPEA1-xTB", Z, X, 0, 0).get("energy")
            ep = s161.sp("IPEA1-xTB", Z, X, 1, 1).get("energy")
            em = s161.sp("IPEA1-xTB", Z, X, -1, 1).get("energy")
            mh = Chem.RemoveHs(m)
            rk = list(Chem.CanonicalRankAtoms(mh, breakTies=False))
            td, _ = s161.shape(S, X)
            risk = s161.risk(td, rdMolDescriptors.CalcNumRotatableBonds(mh),
                             len(set(rk)) / mh.GetNumHeavyAtoms(), Descriptors.MolWt(mh))
            cache[n] = {"IP": round(float((ep - e0) * s161.HA - s161.SHIFT), 3),
                        "EA": round(float((e0 - em) * s161.HA - s161.SHIFT), 3),
                        "risk": round(risk, 1)}
        print(n, cache[n], flush=True)
        json.dump(cache, open(PCACHE, "w"), indent=1)
    return cache


def kmc_neat(thetas=(0.0, 0.1, 0.2, 0.35, 0.5, 1.0), nseed=3):
    """N_sat on a NEAT partner film (x = 1) vs its trap fraction."""
    spec = importlib.util.spec_from_file_location("k68", os.path.join(HERE, "68_mixed_seed_kmc.py"))
    k68 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(k68)
    out = {}
    for th in thetas:
        ns, tp = [], []
        for s in range(nseed):
            if th >= 1.0:   # every site a trap = neat HATCN
                n, t = k68.run_case(128, 1.0e4, 0.0, 0.0, False, seed=s)
            else:
                n, t = k68.run_case(128, 1.0e4, 1.0, th, th > 0, seed=s)
            if n:
                ns.append(n)
            if t:
                tp.append(t)
        out[str(th)] = {"N_sat": float(np.mean(ns)) if ns else None,
                        "theta_perc": float(np.mean(tp)) if tp else None}
        print("neat partner theta", th, out[str(th)], flush=True)
    json.dump(out, open(KMC1, "w"), indent=1)
    return out


def main():
    rank = {r["name"]: r for r in json.load(open(os.path.join(ROOT, "runs", "common_ranking.json")))}
    k68 = json.load(open(os.path.join(ROOT, "runs", "mixed_seed_kmc.json")))
    names = [n for n in THETA_SITE if n in rank]
    pp = partner_props(names)
    neat = json.load(open(KMC1)) if os.path.exists(KMC1) and "--kmc" not in sys.argv else kmc_neat()
    # N_sat at x = 0.5 vs theta_eff (script 68 grid), relative to neat HATCN
    grid = sorted((float(k), [r for r in v if r["x"] == 0.5][0]["N_sat"])
                  for k, v in k68["curves_by_theta_p"].items())
    n_hat = k68["curves_by_theta_p"]["0.0"][0]["N_sat"]
    gx, gy = np.array([g[0] for g in grid]), np.array([g[1] for g in grid])
    nx = np.array(sorted(float(k) for k in neat))
    ny = np.array([neat[str(k) if str(k) in neat else repr(k)]["N_sat"] for k in nx])
    rows = []
    for n in names:
        eb = rank[n]["E_b"]
        pt, e_star = p_trap(eb)
        th = THETA_SITE[n] * pt
        ip = pp[n]["IP"]
        ct = "strong" if ip < 6.5 else ("weak possible" if ip < 7.0 else "none")
        rows.append({
            "formulation": f"HATCN:{n} 1:1", "partner": n, "E_b_partner": eb,
            "tau_s": float(f"{tau(eb):.3g}"), "p_trap": round(pt, 2), "theta_eff": round(th, 3),
            "N_sat_rel": round(float(np.interp(th, gx, gy)) / n_hat, 3),
            "hatcn_domains_killed": True, "partner_risk": pp[n]["risk"],
            "IP_partner": ip, "CT_class": ct, "CT_flag": ct != "none",
            "dq_partner": rank[n]["dq"]})
    # neat references
    rows.append({"formulation": "HATCN neat", "partner": None, "E_b_partner": 1.631, "p_trap": 1.0,
                 "theta_eff": 1.0, "N_sat_rel": 1.0, "hatcn_domains_killed": False,
                 "partner_risk": None, "CT_flag": False})
    # Bphen:Cs2CO3 alone: chelate theta 0.35 at 0.62 eV + Cs sites at 0.89 eV
    pt_b, _ = p_trap(0.621)
    pt_cs, _ = p_trap(0.886)
    for y in (0.0, 0.05, 0.10, 0.20):
        # one Cs+ per ~(1/y) Bphen; each Cs exposes ~4 surface sites of a 4-site block
        th = min(1.0, 0.35 * pt_b + y * 4 / 4 * pt_cs)
        rows.append({"formulation": f"Bphen:Cs {int(y*100)}% (no HATCN)", "partner": "Bphen",
                     "E_b_partner": 0.886 if y else 0.621, "p_trap": round(pt_b, 2),
                     "theta_eff": round(th, 3),
                     "N_sat_rel": round(float(np.interp(th, nx, ny)) / n_hat, 3),
                     "hatcn_domains_killed": None, "partner_risk": pp["Bphen"]["risk"],
                     "CT_flag": False, "electrical_role": "n-side: fails CGL p-side spec"})
    # verdict
    for r in rows:
        bad = []
        if r.get("CT_flag"):
            bad.append(f"CT with HATCN ({r.get('CT_class')})")
        if r["partner_risk"] is not None and r["partner"] and r["partner_risk"] > 58 and "Cs" not in r["formulation"]:
            bad.append("partner crystallises itself")
        if r["hatcn_domains_killed"] is False:
            bad.append("HATCN domains span the film (crystallisation, Park & Suh)")
        if r.get("electrical_role"):
            bad.append("not a p-side seed")
        r["issues"] = bad
        r["score"] = round(r["N_sat_rel"] - 0.15 * len(bad)
                           - 0.002 * max(0, (r["partner_risk"] or 0) - 30), 3)
    rows.sort(key=lambda r: -r["score"])
    out = {"params": {"sigma_Eb": SIG_EB, "E_star_trap_eV": round(e_star, 3), "t_capture_s": T_CAP,
                      "x_window": [X_KILL, X_PERC]}, "rows": rows,
           "kmc_neat_partner": neat}
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"trap threshold E* = {e_star:.3f} eV (tau = t_capture)")
    for r in rows:
        print(f"{r['formulation']:<28} Eb {r['E_b_partner']:.3f} ptrap {r['p_trap']:.2f} "
              f"Nsat {r['N_sat_rel']:.3f} risk {r['partner_risk']} CT {r.get('CT_class','-')} "
              f"score {r['score']:+.3f} {'; '.join(r['issues'])}")


if __name__ == "__main__":
    main()
