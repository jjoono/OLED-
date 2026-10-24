"""Descriptor table for the two screening slides.

    python scripts/152_slide_descriptors.py   ->  runs/slide_descriptors.csv

One row per organic computed at PBE0 (19), Cs2CO3 (n-type EIL, inorganic),
plus MoOx as a reference row.
Each column comes from one method only, so a bar chart of a column compares
like with like:

  E_b       PBE0-D3BJ/def2-SVP, Ag adatom at the deepest site (multi-site)
  E_elec    electronic part of E_b (PBE0, no D3)
  E_disp    dispersion part of E_b (D3BJ, whole molecule)
  q_Ag1     Mulliken charge on the single Ag adatom, PBE0
  dq_Ag1_2  total Ag charge with two atoms minus with one, GFN2
            (<= 0: charge handed back, Ag2 metallic; >= +0.3: veto)

The four pi-only materials still on motif values (CBP, TCTA, TPD, TAPC) are
left out: mCP and NPB, computed directly, stand for that class.
"""
import csv
import json
import os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
R = lambda f: json.load(open(os.path.join(ROOT, "runs", f)))

ROLE = {"HATCN": "HIL", "F4TCNQ": "HIL (p-dopant)", "F6-TCNNQ": "HIL (p-dopant)",
        "p-bPPhenB": "ETL", "Bphen": "ETL", "B3PyMPM": "ETL", "Liq": "EIL",
        "TPBi": "ETL", "BCP": "ETL", "NBphen": "ETL", "TmPyPB": "ETL", "3TPYMB": "ETL",
        "T2T": "ETL", "PO-T2T": "ETL", "Alq3": "ETL", "DPEPO": "host", "TSPO1": "host",
        "mCP": "host", "NPB": "HTL", "Cs2CO3": "EIL (inorganic)"}

# earlier full-molecule PBE0 values (binding_energies_pbe0 + multi-site
# updates); electronic/dispersion split from the 149 validation entries,
# eb_bsse_basis.json (F4TCNQ) and a D3 evaluation on Liq_Ag_pbe0.xyz;
# q(Ag1) from the same PBE0 runs (package G for HATCN and F4TCNQ,
# hidden_sites_pbe0.json for Liq)
EARLIER = {
    #            E_b    E_elec  E_disp   q_Ag1
    "HATCN":    (1.631, 1.456, 0.175, 0.572),
    "F4TCNQ":   (1.208, 1.147, 0.061, 0.545),
    "p-bPPhenB": (0.785, 0.589, 0.195, 0.148),
    "Bphen":    (0.621, 0.475, 0.134, -0.085),
    "B3PyMPM":  (0.416, 0.195, 0.221, -0.071),
    "Liq":      (0.398, 0.307, 0.091, -0.104),
    # PBE0 single point on Cs2CO3_Ag_pbe0.xyz (reproduces E_b 0.7006 exactly)
    "Cs2CO3":   (0.701, 0.609, 0.091, -0.112),
}
# MoOx reference: Mo3O9 cluster for E_b and q (PBE0), MoO3(010) slab for dq
# (GFN2, the same method as the organics' dq). On the small cluster the PBE0
# Ag1->Ag2 step is only +0.10; the growth shows by Ag4 (+1.22 e total).
MOOX = dict(name="MoOx (ref.)", role="oxide HIL", E_b=2.164, E_elec=2.039, E_disp=0.125,
            q_Ag1=0.722, source="Mo3O9 cluster (PBE0) / MoO3(010) slab (GFN2)")


def main():
    old = R("xtb_survey.json")
    sv = R("common_survey_xtb.json")
    fr = R("common_pbe0_fragment.json")
    slab = R("moo3_slab_xtb.json")["MoO3_010"]
    rows = []
    for n, (eb, el, di, q) in EARLIER.items():
        dq = 2 * old[n]["q_Ag2_per_atom"] - old[n]["q_Ag1"]
        rows.append(dict(name=n, role=ROLE[n], E_b=eb, E_elec=el, E_disp=di, q_Ag1=q,
                         dq_Ag1_2=round(dq, 3), source="full-molecule PBE0"))
    for n, row in fr.items():
        if n.startswith("_") or "_deepest" not in row:
            continue
        k = row["_deepest"]
        d = row[k]
        rows.append(dict(name=n, role=ROLE[n], E_b=round(d["E_b"], 3),
                         E_elec=round(d["elec"], 3), E_disp=round(-d["D3"], 3),
                         q_Ag1=round(d["q_mulliken"], 3),
                         dq_Ag1_2=round(sv[n]["delta_q_Ag1_Ag2"], 3),
                         source="fragment PBE0" + (" (Ag re-placed)" if "relax" in k else "")))
    rows.sort(key=lambda r: -r["E_b"])
    m = dict(MOOX, dq_Ag1_2=round(slab["2"]["q_total"] - slab["1"]["q_total"], 3))
    rows.append(m)
    out = os.path.join(ROOT, "runs", "slide_descriptors.csv")
    keys = ["name", "role", "E_b", "E_elec", "E_disp", "q_Ag1", "dq_Ag1_2", "source"]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in keys})
    print(f"{'name':<12} {'role':<15} {'E_b':>6} {'elec':>6} {'disp':>6} {'q_Ag1':>7} {'dq':>7}")
    for r in rows:
        print(f"{r['name']:<12} {r['role']:<15} {r['E_b']:6.3f} {r['E_elec']:6.3f} "
              f"{r['E_disp']:6.3f} {r['q_Ag1']:+7.3f} {r['dq_Ag1_2']:+7.3f}")


if __name__ == "__main__":
    main()
