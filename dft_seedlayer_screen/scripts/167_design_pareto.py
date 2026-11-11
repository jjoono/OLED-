"""Step 2 end: hard constraints + Pareto ranking of the design library.

    python scripts/167_design_pareto.py  ->  runs/design/pareto.json

Hard constraints and objectives from runs/design/spec.json
(notes/DESIGN_SPEC_CGL_SEED.md). E_b per molecule, best available source:
  1. PBE0 anchor (script 164) on the production scale (HATCN offset)
  2. Delta-ML GP prediction (script 163) with its sigma
  3. ridge prediction, kept alongside as the model-disagreement check
"Safe" E_b = estimate - sigma; a candidate meets the E_b target only if the
safe value does.

Objectives (all four, non-dominated sorting):
  E_b (max), crystallisation risk (min), GFN2 gap (max, transparency proxy),
  synthesis cost (min: homoleptic one-pot 0, statistical mixture 1).
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "runs", "design")


def jl(p):
    return {json.loads(l)["name"]: json.loads(l) for l in open(p) if l.strip()}


def fronts(F):
    """F: n x k, all objectives to MAXIMISE. Returns front index per row."""
    n = len(F)
    rank = np.full(n, -1)
    left = set(range(n))
    r = 0
    while left:
        cur = [i for i in left if not any(
            np.all(F[j] >= F[i]) and np.any(F[j] > F[i]) for j in left if j != i)]
        for i in cur:
            rank[i] = r
        left -= set(cur)
        r += 1
    return rank


def main():
    spec = json.load(open(os.path.join(D, "spec.json")))
    h = spec["hard"]
    lib = json.load(open(os.path.join(D, "library.json")))
    P = jl(os.path.join(D, "props.jsonl"))
    S = jl(os.path.join(D, "sites.jsonl"))
    pred = json.load(open(os.path.join(D, "dml_pred.json")))
    A = json.load(open(os.path.join(D, "pbe0_anchor.json"))) if os.path.exists(os.path.join(D, "pbe0_anchor.json")) else {}
    off = 1.631 - A["HATCN"]["result"]["E_b"] if "result" in A.get("HATCN", {}) else None
    td = json.load(open(os.path.join(D, "tddft.json"))) if os.path.exists(os.path.join(D, "tddft.json")) else {}
    rows = []
    for n, s in S.items():
        if "error" in s or n not in pred:
            continue
        L, p = lib[n], P[n]
        gp, sg = pred[n]["gp"]
        src = "GP"
        eb = gp
        if n in A and "result" in A[n] and off is not None:
            eb = A[n]["result"]["E_b"] + off if n != "HATCN" else 1.631
            sg = 0.10
            src = "PBE0 anchor"
        cls = L.get("synth_class") or 1
        fails = []
        if p["EA"] < h["EA_min_eV"]:
            fails.append(f"EA {p['EA']:.2f}")
        if s.get("dq", 1) > h["dq_max_e"]:
            fails.append(f"dq {s.get('dq')}")
        if p["gap_gfn2"] < h["gap_gfn2_min_eV"]:
            fails.append(f"gap {p['gap_gfn2']:.2f}")
        if not h["MW_range"][0] <= L["MW"] <= h["MW_range"][1]:
            fails.append(f"MW {L['MW']}")
        if L["HBD"] > h["HBD_max"]:
            fails.append("HBD")
        if cls > h["synth_class_max"]:
            fails.append("route")
        if L["SA"] > h["SA_max"]:
            fails.append(f"SA {L['SA']}")
        row = {"name": n, "family": L["family"], "E_b": round(eb, 3), "sigma": round(sg, 3),
               "E_b_safe": round(eb - sg, 3), "E_b_ridge": pred[n]["ridge"][0], "source": src,
               "EA": p["EA"], "gap": p["gap_gfn2"], "risk": p["cryst_risk"], "dq": s.get("dq"),
               "MW": L["MW"], "SA": L["SA"], "synth_class": cls, "route": L["route"],
               "fails": fails, "meets_Eb_target": bool(eb - sg >= spec["Eb_target_eV"])}
        if n in td:
            t = td[n]
            row["tddft"] = {"neutral_first_bright_eV": t.get("neutral", {}).get("first_bright_eV"),
                            "anion_vis_sum_f": t.get("anion", {}).get("vis_sum_f")}
        rows.append(row)
    ok = [r for r in rows if not r["fails"] and r["family"] != "reference"]
    F = np.array([[r["E_b"], -r["risk"], r["gap"], -(r["synth_class"] - 1)] for r in ok])
    rk = fronts(F) if len(ok) else []
    for r, k in zip(ok, rk):
        r["pareto_front"] = int(k)
    ok.sort(key=lambda r: (r["pareto_front"], -r["E_b"]))
    out = {"n_surveyed": len(rows), "n_pass_hard": len(ok), "HATCN": [r for r in rows if r["name"] == "HATCN"],
           "pass": ok, "fail_counts": {}}
    for r in rows:
        for f in r["fails"]:
            k = f.split()[0]
            out["fail_counts"][k] = out["fail_counts"].get(k, 0) + 1
    json.dump(out, open(os.path.join(D, "pareto.json"), "w"), indent=1)
    print(f"surveyed {len(rows)}, pass hard {len(ok)}; fails {out['fail_counts']}")
    for r in ok:
        if r["pareto_front"] > 1:
            continue
        print(f"F{r['pareto_front']} {r['name']:<34} Eb {r['E_b']:.2f}±{r['sigma']:.2f} ({r['source']}) "
              f"ridge {r['E_b_ridge']:.2f} EA {r['EA']:.2f} gap {r['gap']:.2f} risk {r['risk']:.0f} "
              f"cls {r['synth_class']} SA {r['SA']}")


if __name__ == "__main__":
    main()
