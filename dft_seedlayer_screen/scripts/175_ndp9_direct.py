"""NDP-9 (Novaled p-dopant) E_b by the common-ranking protocol, replacing the
EA-correlation estimate of script 113 (1.17 +- 0.14 eV).

    QC_THREADS=8 python scripts/175_ndp9_direct.py   -> runs/ndp9_direct.json

Structure: the [3]radialene tris[cyano(4-cyano-2,3,5,6-tetrafluorophenyl)
methylene]cyclopropane, C30F12N6, as drawn for NDP-9 in patents (checked
against structures/NDP9.xyz: same molecule). Novaled has not disclosed NDP-9
officially; this is the structure the literature and patents assign to it.
GFN2 site survey (148) -> PBE0-D3/def2-SVP fragment (149) at the three
deepest distinct sites; sites with Ag-H < 2.4 A do not arise (no H).
"""
import importlib.util, json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "runs", "ndp9_direct.json")
os.environ.setdefault("FRAG_RCUT", "5.0")


def load(fn, name):
    s = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m


def main():
    from pathgeom import read_xyz
    out = json.load(open(OUT)) if os.path.exists(OUT) else {}
    if "survey" not in out:
        s148 = load("148_common_materials_survey.py", "s148")
        s148.COMMON = os.path.join(ROOT, "structures", "ndp9")
        t0 = time.time(); out["survey"] = s148.survey("NDP9"); out["survey"]["seconds"] = round(time.time() - t0)
        json.dump(out, open(OUT, "w"), indent=1)
        sv = out["survey"]
        print("survey", sv["deepest_site"], sv["sites"][sv["deepest_site"]], "dq", sv.get("delta_q_Ag1_Ag2"), flush=True)
    f149 = load("149_fragment_pbe0.py", "f149")
    S, X = read_xyz(os.path.join(ROOT, "structures", "ndp9", "NDP9.xyz")); X = np.asarray(X, float)
    sv = out["survey"]; out.setdefault("pbe0", {})
    for lab in sv["pbe0_sites"]:
        if lab in out["pbe0"]:
            continue
        t0 = time.time()
        try:
            r = f149.eb(S, X, np.array(sv["sites"][lab]["pos"]))
        except Exception as ex:
            r = {"error": repr(ex)[:200]}
        r.update(near=sv["sites"][lab]["near"], gfn2_E_b=sv["sites"][lab]["E_b"], seconds=round(time.time() - t0))
        out["pbe0"][lab] = r; json.dump(out, open(OUT, "w"), indent=1)
        print("pbe0", lab, r, flush=True)
    ok = {k: v["E_b"] for k, v in out["pbe0"].items() if "E_b" in v}
    if ok:
        out["E_b_deepest"] = max(ok.values()); json.dump(out, open(OUT, "w"), indent=1)
        print("E_b(NDP-9) =", out["E_b_deepest"])


if __name__ == "__main__":
    main()
