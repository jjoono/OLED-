"""NDP-9 E_b on the whole molecule (no fragment), at the three sites of 175.

    QC_THREADS=8 python scripts/178_ndp9_whole.py   -> runs/ndp9_direct.json ["whole"]

The fragment values of 175 move with the fragment (1.95 eV on 15 atoms,
2.42 on 23, 2.20 on 32): NDP-9's LUMO spans the whole radialene, so a cut
changes the acceptor. FRAG_RCUT = 50 A makes 149's fragment the full molecule.
"""
import importlib.util, json, os, sys, time
import numpy as np
os.environ["FRAG_RCUT"] = "50"
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "runs", "ndp9_direct.json")
s = importlib.util.spec_from_file_location("f149", os.path.join(HERE, "149_fragment_pbe0.py"))
f149 = importlib.util.module_from_spec(s); s.loader.exec_module(f149)
from pathgeom import read_xyz
S, X = read_xyz(os.path.join(ROOT, "structures", "ndp9", "NDP9.xyz")); X = np.asarray(X, float)
for lab in json.load(open(OUT))["survey"]["pbe0_sites"]:
    d = json.load(open(OUT)); w = d.setdefault("whole", {})
    if lab in w:
        continue
    t0 = time.time()
    try:
        r = f149.eb(S, X, np.array(d["survey"]["sites"][lab]["pos"]))
    except Exception as ex:
        r = {"error": repr(ex)[:200]}
    r["seconds"] = round(time.time() - t0)
    d = json.load(open(OUT)); d.setdefault("whole", {})[lab] = r; json.dump(d, open(OUT, "w"), indent=1)
    print("whole", lab, "E_b", r.get("E_b"), "elec", r.get("elec"), "atoms", r.get("frag_atoms"), r["seconds"], "s", flush=True)
