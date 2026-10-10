"""Check, at PBE0, the four candidates where GFN2-xTB found a site deeper than
the one stage 1 used.

    python scripts/144_hidden_sites_pbe0.py

The xTB survey (script 142) re-sampled every candidate and agreed with stage 1
on fifteen of the nineteen single-site rows. On four -- DMABN, Liq, Al4O6,
Cu4I4 -- it found a site 0.2-0.9 eV deeper at its own level. GFN2 overbinds
nitrogen chelates by about a factor of two, so its margin is not evidence; this
re-evaluates both sites at PBE0-D3/def2-SVP (PySCF, validated against Gaussian
to 5 meV) on the GFN2-relaxed silver positions.

Both complexes are single points at GFN2 geometries, so neither is at its PBE0
minimum and the unrelaxed error is shared; the difference between them is the
number that decides whether the stage-1 row should be replaced.
"""
import os, sys, json, importlib

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qclocal as q
from pathgeom import CANDIDATES, STRUCT, read_xyz, relaxed_file
from pyscf import lib

s142 = importlib.import_module("142_xtb_survey")
lib.num_threads(int(os.environ.get("QC_THREADS", "12")))
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "runs", "hidden_sites_pbe0.json")
FLAGGED = ["DMABN", "Liq", "Cu4I4", "Al4O6"]
EV = 27.211386


def geometries(tag):
    fn = {t: f for t, f, _, _ in CANDIDATES}[tag]
    S, X = read_xyz(os.path.join(STRUCT, relaxed_file(fn)))
    k = S.index("Ag")
    ss = [s for i, s in enumerate(S) if i != k]
    sx = np.delete(np.asarray(X, float), k, axis=0)
    deep = json.load(open(os.path.join(ROOT, "runs", "xtb_survey.json")))[tag]["deepest_site"]
    starts = {f"{kind}{i}": pos for kind, i, pos, _ in s142.s131.unique_sites(ss, sx)}
    starts.update(dict(s142.pockets(ss, sx)))
    starts["stage1"] = np.asarray(X[k], float)
    out = {}
    for lab in ("stage1", deep):
        r = s142.relax(ss + ["Ag"], np.vstack([sx, starts[lab]]), len(ss), 1)
        out[lab] = r[2][len(ss)]
    return ss, sx, deep, out


def main():
    done = json.load(open(OUT)) if os.path.exists(OUT) else {}
    e_ag = None
    for tag in FLAGGED:
        if tag in done:
            continue
        ss, sx, deep, pos = geometries(tag)
        os.makedirs(os.path.join(STRUCT, "hidden_sites"), exist_ok=True)
        if e_ag is None:
            e_ag = float(q.scf(q.mol(["Ag"], [[0, 0, 0]], 2))[1])
        e_mol = float(q.scf(q.mol(ss, sx, 1))[1])
        d3m = q.d3bj(ss, sx)
        row = {"deepest_xtb": deep}
        for lab, p in pos.items():
            S = ss + ["Ag"]
            X = np.vstack([sx, p])
            mf, e = q.scf(q.mol(S, X, 2))
            eb = (e_mol + e_ag - e) * EV - (q.d3bj(S, X) - d3m) * EV
            qa = float(q.mulliken(mf)[-1])
            row[lab] = {"E_b_pbe0_at_xtb_geom": round(eb, 4), "q_mulliken": round(qa, 4),
                        "converged": bool(mf.converged)}
            with open(os.path.join(STRUCT, "hidden_sites", f"{tag}_{lab}.xyz"), "w") as f:
                f.write(f"{len(S)}\n{tag} Ag at {lab}, GFN2-xTB relaxed, molecule frozen\n")
                for a, c in zip(S, X):
                    f.write(f"{a} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n")
        row["delta_deep_minus_stage1"] = round(
            row[deep]["E_b_pbe0_at_xtb_geom"] - row["stage1"]["E_b_pbe0_at_xtb_geom"], 4)
        done[tag] = row
        json.dump(done, open(OUT, "w"), indent=1)
        print(f"{tag:<8} stage1 {row['stage1']['E_b_pbe0_at_xtb_geom']:.3f}  "
              f"{deep} {row[deep]['E_b_pbe0_at_xtb_geom']:.3f}  "
              f"delta {row['delta_deep_minus_stage1']:+.3f} eV", flush=True)


if __name__ == "__main__":
    main()
