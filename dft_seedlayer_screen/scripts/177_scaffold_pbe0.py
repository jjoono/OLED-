"""PBE0 E_b for new scaffolds from script 176 (same protocol as 148/149/175).

    QC_THREADS=4 python scripts/177_scaffold_pbe0.py NAME ...   -> runs/design/scaffold_pbe0.json

EA does not predict E_b across families (script 167: r = -0.09 inside HAT,
and GFN2 inverts HATCN vs F4TCNQ), so a new core needs its own PBE0 number
before D (which assumes every anchor traps) means anything.
Geometry: the GFN2 structure written by script 176 (structures/design/NAME.xyz).
"""
import importlib.util, json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "runs", "design", "scaffold_pbe0.json")
os.environ.setdefault("FRAG_RCUT", "5.0")


def load(fn, name):
    s = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m


def main(names):
    from pathgeom import read_xyz
    s148 = load("148_common_materials_survey.py", "s148")
    s148.COMMON = os.path.join(ROOT, "structures", "design")
    f149 = load("149_fragment_pbe0.py", "f149")
    for name in names:
        out = json.load(open(OUT)) if os.path.exists(OUT) else {}
        row = out.setdefault(name, {})
        if "survey" not in row:
            t0 = time.time(); row["survey"] = s148.survey(name); row["survey"]["seconds"] = round(time.time() - t0)
            json.dump(out, open(OUT, "w"), indent=1)
            sv = row["survey"]
            print(name, "survey deepest", sv["deepest_site"], sv["sites"][sv["deepest_site"]]["near"],
                  "GFN2", sv["sites"][sv["deepest_site"]]["E_b"], "dq", sv.get("delta_q_Ag1_Ag2"), flush=True)
        S, X = read_xyz(os.path.join(ROOT, "structures", "design", f"{name}.xyz")); X = np.asarray(X, float)
        sv = row["survey"]; row.setdefault("pbe0", {})
        for lab in sv["pbe0_sites"]:
            if lab in row["pbe0"]:
                continue
            ag = np.array(sv["sites"][lab]["pos"]); d = np.linalg.norm(X - ag, axis=1)
            t0 = time.time()
            try:
                r = f149.eb(S, X, ag)
            except Exception as ex:
                r = {"error": repr(ex)[:200]}
            hs = [d[i] for i, s in enumerate(S) if s == "H"]
            r.update(near=sv["sites"][lab]["near"], gfn2_E_b=sv["sites"][lab]["E_b"],
                     d_H_min=round(float(min(hs)), 3) if hs else None, seconds=round(time.time() - t0))
            row["pbe0"][lab] = r
            out = json.load(open(OUT)) if os.path.exists(OUT) else {}
            out[name] = row; json.dump(out, open(OUT, "w"), indent=1)
            print(name, "pbe0", lab, "E_b", r.get("E_b"), "elec", r.get("elec"), "H_min", r["d_H_min"], r["seconds"], "s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
