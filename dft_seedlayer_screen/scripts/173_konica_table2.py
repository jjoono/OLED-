"""Konica Minolta US 9,923,166 B2, Example 2 / Table 2, re-scored with our protocol.

    python scripts/173_konica_table2.py --build
    python scripts/173_konica_table2.py --survey No03 No04 ...   # GFN2 site survey + Ag2 charge
    QC_THREADS=4 python scripts/173_konica_table2.py --pbe0 No03 ...
    python scripts/173_konica_table2.py --summary
    -> structures/konica/, runs/konica/{survey,pbe0}_<No>.json, runs/konica/summary.json

The patent ranks ten N-heterocycle seed layers (25 nm) under 5 nm Ag
(0.1-0.2 nm/s) with the "effective action energy"
    dEef = n * dE / s
(dE: Ag-N interaction energy from Gaussian 03, method not stated; n: N atoms
whose Ag-N dihedral D < 10 deg; s: molecular surface area from Tencube/WM)
and reports the sheet resistance. Structures were read from Chemical
Formulas 70-71 of the patent PDF (USPTO); dibenzothiophene links are taken as
2,8 and carbazole links as 3,6 (the drawings fix which ring, the standard
positions are assumed). No.01 is anthracene (no N, Rs unmeasurable).

Our protocol is the one of the common-materials ranking
(notes/COMMON_OLED_MATERIALS_SEED_RANKING.md):
  build   SMILES -> ETKDG/MMFF best of 20 -> GFN2 relaxed (script 147)
  survey  GFN2: lone-pair, heteroatom, hollow and pocket sites, Ag relaxed on
          the frozen molecule; Ag2 at the deepest site for dq(Ag1->Ag2) (148)
  pbe0    PBE0/def2-SVP electronic part on a 5 A ring-completed fragment +
          whole-molecule D3BJ (149), at the three deepest distinct GFN2 sites.
          Sites with Ag-H < 2.4 A are flagged (GFN2 crowds Ag against H; 151).
"""
import os, sys, json, time, importlib
import importlib.util

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
STRUCT = os.path.join(ROOT, "structures", "konica")
RUNS = os.path.join(ROOT, "runs", "konica")
os.environ.setdefault("FRAG_RCUT", "5.0")   # the common ranking used 5 A


def _load(fn, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


PY = "c9ccccn9"                 # 2-pyridyl
PYPH = "c8ccccc8-c8ccccn8"      # 2-(pyridin-2-yl)phenyl


def cz36(r3, r6):               # carbazol-9-yl, 3- and 6-substituents ("" = H)
    a = f"(-{r3})" if r3 else ""
    b = f"(-{r6})" if r6 else ""
    return f"n3c4ccc{a}cc4c4cc{b}ccc43"


def dbt28(r2, r8):              # dibenzothiophene, 2- and 8-substituents
    return f"s1c2ccc(-{r2})cc2c2cc(-{r8})ccc21"


AR4 = "c2cc(-c3ccccc3)c(-c3ccccc3)c(-c3ccccc3)c2-c2ccccc2"   # 2,3,4,5-tetraphenylphenyl
DCARB = f"n2c3cccnc3c3cc(-{PYPH})ccc32"                       # delta-carbolin-5-yl, 8-PyPh

COMPOUNDS = {
    # name: (SMILES, patent n, dE kcal/mol, s A^2, dEef, Rs ohm/sq or None = unmeasurable)
    "No01": ("c1ccc2cc3ccccc3cc2c1", None, None, None, None, None),
    "No02": (f"c1cc(-{AR4})nc(-{AR4})c1", 1, -51.6, 902.7, -0.057, None),
    "No03": (dbt28(cz36(PYPH, PYPH), cz36(PYPH, "")), 3, -36.418, 981, -0.111, 963),
    "No04": (dbt28(cz36(PYPH, ""), PYPH), 2, -48.933, 666, -0.147, 915750),
    "No05": (f"c1cc2sc3ccccc3c2cc1-n1c2ccc(-{PYPH})cc2c2cc(-{PYPH})ccc21", 2, -49.419, 666, -0.148, 39120),
    "No06": (dbt28(cz36(PYPH, PYPH), PYPH), 3, -43.768, 821, -0.160, 769),
    "No07": (dbt28(cz36(PYPH, PYPH), PY), 3, -45.834, 739, -0.186, 3337),
    "No08": (f"c1cc(-{DCARB})cc(-{DCARB})c1", 4, -41.989, 729, -0.230, 278),
    "No09": (f"c1c(-{PY})c(-{PY})c(-{PY})c(-{PY})c1-{PY}", 2, -56.976, 479, -0.238, 122),
    "No10": (f"c1cc(-{cz36(PY, PY)})cc(-{cz36(PY, PY)})c1", 4, -45.01, 723, -0.249, 792),
    "No11": (f"c1c(-{PY})c(-{PY})c(-{PY})cc1-c1cc(-{PY})c(-{PY})c(-{PY})c1", 4, -57.302, 634, -0.362, 52),
}


def build(names):
    from rdkit import Chem
    s147 = _load("147_common_materials_build.py", "s147")
    os.makedirs(STRUCT, exist_ok=True)
    for name in names:
        smi = COMPOUNDS[name][0]
        p = os.path.join(STRUCT, f"{name}.xyz")
        if os.path.exists(p):
            continue
        t0 = time.time()
        m, cid = s147.embed(smi)
        S, X = s147.gfn2_opt(m, cid)
        conf = m.GetConformer(cid)
        for i, x in enumerate(X):
            conf.SetAtomPosition(i, x.tolist())
        Chem.MolToMolFile(m, os.path.join(STRUCT, f"{name}.mol"), confId=cid)
        with open(p, "w") as f:
            f.write(f"{len(S)}\n{name} Konica US9923166 Table 2; ETKDG/MMFF -> GFN2-xTB\n")
            for s, x in zip(S, X):
                f.write(f"{s} {x[0]:.6f} {x[1]:.6f} {x[2]:.6f}\n")
        print(f"{name} {len(S)} atoms {Chem.rdMolDescriptors.CalcMolFormula(m)} "
              f"{time.time()-t0:.0f}s", flush=True)


def survey(names):
    s148 = _load("148_common_materials_survey.py", "s148")
    s148.COMMON = STRUCT
    os.makedirs(RUNS, exist_ok=True)
    for name in names:
        out = os.path.join(RUNS, f"survey_{name}.json")
        if os.path.exists(out):
            continue
        t0 = time.time()
        r = s148.survey(name)
        r["seconds"] = round(time.time() - t0)
        json.dump(r, open(out, "w"), indent=1)
        print(f"{name} {r['n_atoms']} atoms {r['n_sites']} sites deepest {r['deepest_site']} "
              f"({r['sites'][r['deepest_site']]['near']}) GFN2 E_b "
              f"{r['sites'][r['deepest_site']]['E_b']:.3f} dq {r.get('delta_q_Ag1_Ag2', float('nan')):+.3f} "
              f"{r['seconds']}s", flush=True)


def pbe0(names):
    f149 = _load("149_fragment_pbe0.py", "f149")
    from pathgeom import read_xyz
    for name in names:
        sv_p = os.path.join(RUNS, f"survey_{name}.json")
        out = os.path.join(RUNS, f"pbe0_{name}.json")
        if not os.path.exists(sv_p):
            continue
        sv = json.load(open(sv_p))
        row = json.load(open(out)) if os.path.exists(out) else {}
        S, X = read_xyz(os.path.join(STRUCT, f"{name}.xyz"))
        X = np.asarray(X, float)
        for lab in sv["pbe0_sites"]:
            if lab in row:
                continue
            ag = np.array(sv["sites"][lab]["pos"])
            d = np.linalg.norm(X - ag, axis=1)
            t0 = time.time()
            try:
                r = f149.eb(S, X, ag)
            except Exception as ex:
                r = {"error": repr(ex)[:200]}
            r.update(seconds=round(time.time() - t0), near=sv["sites"][lab]["near"],
                     gfn2_E_b=sv["sites"][lab]["E_b"], R_cut=f149.R_CUT,
                     d_H_min=round(float(min(d[i] for i, s in enumerate(S) if s == "H")), 3),
                     n_N_within_3A=int(sum(1 for i, s in enumerate(S) if s == "N" and d[i] < 3.0)))
            row[lab] = r
            json.dump(row, open(out, "w"), indent=1)
            print(f"{name} {lab:<18} near {r['near']:<5} E_b {r.get('E_b', float('nan')):.3f} "
                  f"(elec {r.get('elec', float('nan')):.3f}) H_min {r['d_H_min']:.2f} "
                  f"frag {r.get('frag_atoms')} {r['seconds']}s", flush=True)


def open_sites(names):
    """Crowded GFN2 sites (Ag-H < 2.4 A) re-placed at the nearest open point.

    GFN2 packs Ag against ortho hydrogens, which PBE0 finds repulsive (No09
    lpN7: Ag-H 1.96 A, E_b 0.07 eV, electronic part negative). Script 151
    relaxed such sites at PBE-D3, 6-10 h each; for eleven compounds that is
    too slow, so the start point of 151 (push_out: no H within 2.3 A, no other
    heavy atom within 2.8 A, closest to the N anchors, between both N for a
    pocket) is scored directly at PBE0 as in 149. Without the relaxation this
    is a lower bound on the site's E_b.
    """
    f149 = _load("149_fragment_pbe0.py", "f149")
    s151 = _load("151_crowded_sites_pbe0.py", "s151")
    from pathgeom import read_xyz
    for name in names:
        out = os.path.join(RUNS, f"pbe0_{name}.json")
        if not os.path.exists(out):
            continue
        row = json.load(open(out))
        S, X = read_xyz(os.path.join(STRUCT, f"{name}.xyz"))
        X = np.asarray(X, float)
        for lab in [k for k, v in row.items() if "E_b" in v and v["d_H_min"] < 2.4]:
            key = lab + "_open"
            if key in row:
                continue
            sv = json.load(open(os.path.join(RUNS, f"survey_{name}.json")))
            ag0, het = s151.push_out(S, X, np.array(sv["sites"][lab]["pos"]))
            if ag0 is None:
                row[key] = {"buried": het}
                json.dump(row, open(out, "w"), indent=1)
                print(f"{name} {key} anchors {het} buried", flush=True)
                continue
            d = np.linalg.norm(X - ag0, axis=1)
            t0 = time.time()
            try:
                r = f149.eb(S, X, ag0)
            except Exception as ex:
                r = {"error": repr(ex)[:200]}
            r.update(seconds=round(time.time() - t0), anchors=het, R_cut=f149.R_CUT,
                     near=f"{S[int(d.argmin())]}{int(d.argmin())}", d_near=round(float(d.min()), 3),
                     d_H_min=round(float(min(d[i] for i, s in enumerate(S) if s == "H")), 3),
                     n_N_within_3A=int(sum(1 for i, s in enumerate(S) if s == "N" and d[i] < 3.0)),
                     pos=[round(float(v), 4) for v in ag0], from_site=lab)
            row = json.load(open(out))
            row[key] = r
            json.dump(row, open(out, "w"), indent=1)
            print(f"{name} {key:<22} anchors {het} E_b {r.get('E_b', float('nan')):.3f} "
                  f"(elec {r.get('elec', float('nan')):.3f}) H_min {r['d_H_min']:.2f} {r['seconds']}s",
                  flush=True)


def spearman(a, b):
    from scipy.stats import spearmanr
    return round(float(spearmanr(a, b).correlation), 3)


def summary():
    from rdkit import Chem
    from rdkit.Chem import rdFreeSASA
    rows = {}
    for name, (smi, n, dE, s, dEef, rs) in COMPOUNDS.items():
        sv_p = os.path.join(RUNS, f"survey_{name}.json")
        pb_p = os.path.join(RUNS, f"pbe0_{name}.json")
        if not (os.path.exists(sv_p) and os.path.exists(pb_p)):
            continue
        sv, pb = json.load(open(sv_p)), json.load(open(pb_p))
        # crowded GFN2 positions (Ag-H < 2.3 A) are artefacts: PBE0 finds them
        # repulsive; their re-placed "_open" versions stand in for them
        ok = {k: v for k, v in pb.items() if "E_b" in v and v["d_H_min"] >= 2.3}
        open_ = {k: v for k, v in ok.items() if v["d_H_min"] >= 2.3}
        best = max(ok, key=lambda k: ok[k]["E_b"]) if ok else None
        best_open = max(open_, key=lambda k: open_[k]["E_b"]) if open_ else None
        # SASA (probe 1.4 A) of the GFN2 structure, our counterpart of the patent's s
        m = Chem.MolFromMolFile(os.path.join(STRUCT, f"{name}.mol"), removeHs=False)
        radii = rdFreeSASA.classifyAtoms(m)
        sasa = rdFreeSASA.CalcSASA(m, radii)
        # GFN2 lone-pair N sites, a count of exposed N comparable to the patent's n
        lp = [v for k, v in sv["sites"].items() if k.startswith("lpN") and v["near_el"] == "N"]
        rows[name] = {"patent": {"n": n, "dE_kcal": dE, "s_A2": s, "dEef": dEef, "Rs": rs},
                      "E_b": ok[best]["E_b"] if best else None,
                      "elec": ok[best]["elec"] if best else None,
                      "disp": round(-ok[best]["D3"], 4) if best else None,
                      "site": ok[best]["near"] if best else None,
                      "d_H_min": ok[best]["d_H_min"] if best else None,
                      "n_N_at_site": ok[best]["n_N_within_3A"] if best else None,
                      "E_b_open_site": open_[best_open]["E_b"] if best_open else None,
                      "dq_Ag1_Ag2": sv.get("delta_q_Ag1_Ag2"),
                      "q_Ag1": sv.get("q_Ag1"),
                      "n_lpN_sites": len(lp),
                      "sasa_A2": round(float(sasa), 1),
                      "n_atoms": sv["n_atoms"]}
    # rank against Rs; "unmeasurable" = worse than every measured value
    meas = [k for k, r in rows.items() if r["patent"]["Rs"] is not None or k == "No02"]
    logrs = {k: (np.log10(rows[k]["patent"]["Rs"]) if rows[k]["patent"]["Rs"] else 9.0) for k in meas}
    corr = {}
    if len(meas) >= 4:
        y = [-logrs[k] for k in meas]       # higher = better seed
        corr["patent_dEef"] = spearman([-rows[k]["patent"]["dEef"] for k in meas], y)
        corr["our_E_b"] = spearman([rows[k]["E_b"] for k in meas], y)
        corr["our_elec"] = spearman([rows[k]["elec"] for k in meas], y)
        corr["our_E_b_open"] = spearman([rows[k]["E_b_open_site"] or 0 for k in meas], y)
        corr["our_E_b_x_nlp_per_sasa"] = spearman(
            [rows[k]["E_b"] * rows[k]["n_lpN_sites"] / rows[k]["sasa_A2"] for k in meas], y)
        corr["patent_n_over_s"] = spearman([rows[k]["patent"]["n"] / rows[k]["patent"]["s_A2"] for k in meas], y)
        corr["n_compounds"] = len(meas)
    out = {"rows": rows, "spearman_vs_minus_log_Rs": corr,
           "note": "Rs unmeasurable (No02) ranked worst. Spearman over the patent's compounds with N."}
    json.dump(out, open(os.path.join(RUNS, "summary.json"), "w"), indent=1)
    print(f"{'':5} {'dEef':>7} {'Rs':>8} | {'E_b':>6} {'elec':>6} {'disp':>6} {'site':>6} {'H_min':>5} "
          f"{'dq':>6} {'nlpN':>4} {'SASA':>6}")
    for k, r in rows.items():
        p = r["patent"]
        print(f"{k:5} {p['dEef'] if p['dEef'] is not None else float('nan'):7.3f} "
              f"{p['Rs'] if p['Rs'] else 'unmeas.':>8} | {r['E_b'] or float('nan'):6.3f} "
              f"{r['elec'] or float('nan'):6.3f} {r['disp'] or float('nan'):6.3f} {str(r['site']):>6} "
              f"{r['d_H_min'] or float('nan'):5.2f} {r['dq_Ag1_Ag2'] or float('nan'):+6.3f} "
              f"{r['n_lpN_sites']:4d} {r['sasa_A2']:6.0f}")
    print(json.dumps(corr, indent=1))


if __name__ == "__main__":
    args = sys.argv[1:]
    names = [a for a in args if not a.startswith("--")] or list(COMPOUNDS)
    if "--build" in args:
        build(names)
    if "--survey" in args:
        survey(names)
    if "--pbe0" in args:
        pbe0(names)
    if "--open" in args:
        open_sites(names)
    if "--summary" in args:
        summary()
