"""Rank common evaporable OLED materials as Ag seed layers.

Ranking variable: PBE0 E_b of a lone Ag atom at the deepest site (multi-site).
Gates:  condensation E_b >= 0.5 eV (below it the first atoms re-evaporate or
        diffuse to aggregate, see notes/LIT_SEED_SYNTHESIS_20261011.md section 3);
        charge-capacity veto dq(Ag1->Ag2) >= +0.3 e.

Sources, in order of trust:
  pbe0       full-molecule PBE0/def2-SVP + D3BJ, multi-site (binding_energies_pbe0 etc.)
  pbe0_frag  fragment PBE0 at the GFN2 site (script 149, validated to 0.013 eV)
  analog     not computed at PBE0; the value of the closest computed motif,
             because GFN2 E_b is not usable across elements (PBE0/GFN2 ratio
             0.19-0.68 on N sites, 0.7-0.85 on all-carbon sites)
The veto always comes from GFN2 (charge ordering validated against PBE0).
"""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# earlier candidates, full-molecule PBE0 (notes/EB_SCREEN_TABLE.md, hidden_sites_pbe0.json)
# role, E_b, site; dq(Ag1->Ag2) = 2 q_Ag2_per_atom - q_Ag1 from runs/xtb_survey.json
EARLIER = {
    "HATCN":     ("HIL (CGL)",          1.631, "aza N"),
    "F4TCNQ":    ("p-dopant",           1.208, "nitrile N"),
    "p-bPPhenB": ("ETL (n-host)",       0.785, "phen chelate"),
    "Bphen":     ("ETL",                0.621, "phen chelate"),
    "B3PyMPM":   ("ETL",                0.416, "pyrimidine/pyridyl N"),
    "Liq":       ("EIL / ETL co-dep",   0.398, "Li-O"),
}
ROLE = {
    "TPBi": "ETL / HBL", "BCP": "HBL / ETL", "NBphen": "ETL", "TmPyPB": "ETL / HBL",
    "3TPYMB": "ETL / HBL", "T2T": "ETL / host", "PO-T2T": "ETL / exciplex acceptor",
    "DPEPO": "blue host / HBL", "TSPO1": "blue host / HBL", "Alq3": "ETL / green host",
    "F6-TCNNQ": "p-dopant", "NPB": "HTL", "TPD": "HTL", "TAPC": "HTL / EBL",
    "TCTA": "HTL / host", "CBP": "host", "mCP": "host / EBL",
}
# closest computed motif for anything without a PBE0 value yet
ANALOG = {
    "TmPyPB": ("B3PyMPM 3-pyridyl N", 0.416),
    "3TPYMB": ("B3PyMPM 3-pyridyl N", 0.416),
    "TPBi":   ("pyridine-type N (benzimidazole N3)", 0.373),
    "NBphen": ("Bphen chelate, minus the 2,9-substituent penalty seen in BCP", 0.47),
    "T2T":    ("triazine N", 0.271),
    "PO-T2T": ("Ph3P=O", 0.408),
    "DPEPO":  ("Ph3P=O", 0.408),
    "TSPO1":  ("Ph3P=O", 0.408),
    "Alq3":   ("Liq (metal quinolinolate O)", 0.398),
    "CBP":    ("PhCz", 0.279),
    "mCP":    ("PhCz", 0.279),
    "TCTA":   ("PhCz / TPA", 0.279),
    "NPB":    ("TPA", 0.259),
    "TPD":    ("TPA", 0.259),
    "TAPC":   ("TPA", 0.259),
}


def main():
    sv = json.load(open(os.path.join(ROOT, "runs", "common_survey_xtb.json")))
    fr = json.load(open(os.path.join(ROOT, "runs", "common_pbe0_fragment.json")))
    old = json.load(open(os.path.join(ROOT, "runs", "xtb_survey.json")))
    rows = []
    for name, (role, eb, site) in EARLIER.items():
        dq = round(2 * old[name]["q_Ag2_per_atom"] - old[name]["q_Ag1"], 3)
        rows.append(dict(name=name, role=role, E_b=eb, site=site, dq=dq, source="pbe0"))
    for name, v in sv.items():
        r = dict(name=name, role=ROLE.get(name, ""), dq=v["delta_q_Ag1_Ag2"])
        f = fr.get(name, {})
        if "_E_b_deepest" in f:
            lab = f["_deepest"]
            r.update(E_b=f["_E_b_deepest"], site=f"{lab} (near {f[lab]['near']})",
                     source="pbe0_frag",
                     other={k: x["E_b"] for k, x in f.items()
                            if not k.startswith("_") and "E_b" in x and k != lab})
        else:
            a, e = ANALOG[name]
            r.update(E_b=e, site=f"{v['deepest_site']} (GFN2)", source="analog: " + a)
        rows.append(r)
    for r in rows:
        r["veto"] = "fail" if r["dq"] >= 0.3 else ("ambiguous" if r["dq"] > 0.1 else "pass")
        r["gate_0.5"] = r["E_b"] >= 0.5
    rows.sort(key=lambda r: -r["E_b"])
    json.dump(rows, open(os.path.join(ROOT, "runs", "common_ranking.json"), "w"), indent=1)
    print(f"{'#':>2} {'material':<10} {'role':<24} {'E_b':>6}  {'dq':>6} veto  gate  source")
    for i, r in enumerate(rows, 1):
        print(f"{i:>2} {r['name']:<10} {r['role']:<24} {r['E_b']:6.3f}  {r['dq']:+.2f} "
              f"{r['veto']:<5} {'yes' if r['gate_0.5'] else 'no':<4}  {r['source']}")


if __name__ == "__main__":
    main()
