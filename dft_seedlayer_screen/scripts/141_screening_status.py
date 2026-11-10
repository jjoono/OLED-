"""The screening table as it can honestly be stated now.

    python scripts/141_screening_status.py

Script 127's table gave one E_b per candidate. Packages E and F showed that a
single starting site is not a sample of the landscape -- HATCN went 0.604 ->
1.631 eV, Mo3O8 0.875 -> 2.185 -- so most of that column is a lower bound, and
one entry sat on a spin-contaminated branch. This rewrites the table with what
each number is: how many sites were sampled, the deepest found, the charge on
silver where it was measured, and a status that says whether the row can be
compared with any other row.

Nothing here is recomputed; it collates runs/*.json.
"""
import os, json

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
R = os.path.join(ROOT, "runs")


def load(n):
    p = os.path.join(R, n)
    return json.load(open(p)) if os.path.exists(p) else {}


def main():
    fs = load("final_screen.json")
    ss = load("site_spread.json")
    cg = load("cluster_growth.json")
    ox = {t: load(f"oxide_sites_{t}.json") for t in ("Mo3O9", "Mo3O8")}
    xs = load("xtb_survey.json")
    hid = load("hidden_sites_pbe0.json")
    # sites rejected on spin (S^2 >> 0.75) never count as the deepest
    rejected = {"Mo3O8": {"termO_top"}}
    rows = []
    for t, v in fs.items():
        r = {"candidate": t, "E_b_stage1": v.get("E_b"), "E_d_inter": v.get("E_d_inter"),
             "sites": 1, "E_b_deepest": v.get("E_b"), "q_Ag1": None, "q_Ag4_per_atom": None}
        if t in ss:
            r["sites"] = len(ss[t]["sites"])
            r["E_b_deepest"] = ss[t]["deepest"]
            r["spread"] = ss[t]["spread"]
        if t in ox and ox[t]:
            good = {k: e for k, e in ox[t].items() if k not in rejected.get(t, ())}
            r["sites"] = len(good)
            r["E_b_deepest"] = max(good.values())
            r["spread"] = max(good.values()) - min(good.values())
        if t in cg:
            r["q_Ag1"] = cg[t]["q1"]
            if "4" in cg[t]["steps"]:
                r["q_Ag4_per_atom"] = cg[t]["steps"]["4"]["q_per_atom"]
        if t in xs:
            x = xs[t]
            r["xtb_sites"] = x["n_sites"]
            r["q_Ag2_per_atom_xtb"] = x.get("q_Ag2_per_atom")
            st = x["sites"].get("stage1", {}).get("E_b")
            r["xtb_hidden_deeper"] = (st is not None and x["E_b_deepest"] - st > 0.15)
        if r["sites"] >= 3:
            r["status"] = "multi-site PBE0: comparable"
        elif t in xs and not r["xtb_hidden_deeper"]:
            r["status"] = "PBE0 single site, confirmed deepest by xTB survey"
        elif t in xs and r["xtb_hidden_deeper"]:
            if t in hid:
                d = hid[t]["delta_deep_minus_stage1"]
                if d > 0.05:
                    deep = hid[t]["deepest_xtb"]
                    # a PBE0 single point at the GFN2 position: not relaxed at
                    # PBE0, so itself a lower bound on that site
                    r["E_b_deepest"] = hid[t][deep]["E_b_pbe0_at_xtb_geom"]
                    r["sites"] = 2
                r["status"] = (f"xTB found deeper site; PBE0 says {d:+.2f} eV -> "
                               + ("row is a lower bound" if d > 0.05 else "stage-1 site stands"))
            else:
                r["status"] = "xTB found a deeper site; PBE0 check pending"
        else:
            r["status"] = "single site: lower bound only"
        if t == "Mo3O8":
            r["status"] += "; stage-1 value was on a spin-contaminated branch"
        if t == "Cu4I4":
            r["status"] += "; E_d unconverged"
        rows.append(r)
    rows.sort(key=lambda r: (r["sites"] < 3, -(r["E_b_deepest"] or 0)))
    json.dump(rows, open(os.path.join(R, "screening_status.json"), "w"), indent=1)

    def f(x, fmt="{:.3f}"):
        return "—" if x is None else fmt.format(x)
    print("| candidate | E_b stage 1 | E_b deepest | sites | q(Ag₁) | q/Ag at Ag₄ | q/Ag at Ag₂ (xTB) | E_d inter | status |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['candidate']} | {f(r['E_b_stage1'])} | **{f(r['E_b_deepest'])}** | "
              f"{r['sites']} | {f(r['q_Ag1'], '{:+.3f}')} | {f(r['q_Ag4_per_atom'], '{:+.3f}')} | "
              f"{f(r.get('q_Ag2_per_atom_xtb'), '{:+.3f}')} | "
              f"{f(r['E_d_inter'])} | {r['status']} |")
    from collections import Counter
    c = Counter(r["status"].split(";")[0].split(" ->")[0] for r in rows)
    print()
    for k, v in c.most_common():
        print(f"  {v:>2}  {k}")


if __name__ == "__main__":
    main()
