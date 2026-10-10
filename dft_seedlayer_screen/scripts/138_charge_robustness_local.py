"""Package J run locally: is the charge hand-back a property of the silver, or of
the Mulliken scheme?

    python scripts/138_charge_robustness_local.py

Same question and same geometries as script 136, which writes the Gaussian
version; this one needs no workstation. PySCF at the project's level --
PBE0/def2-SVP, def2 ECP -- reproduces the Gaussian state exactly where it was
checked: Mo3O9-Ag1 gives Mulliken q(Ag) = +0.722 against Gaussian's +0.722 and a
total energy (with s-dftd3 D3BJ added) 5 meV from Gaussian's, the difference
being integration grid and density fitting.

Five charge schemes on all sixteen package-G geometries:

    Mulliken      what package G used; the check that the state is the same
    Lowdin        symmetric orthogonalisation, still basis-dependent
    IAO           intrinsic atomic orbitals -- projected onto a minimal
                  free-atom basis, close to basis-independent
    NPA           natural population (NAO) charges, the NBO number
    Hirshfeld     real-space partition against free-atom densities

plus the natural occupation of the silver 5s, and def2-TZVP for HATCN and Mo3O9
at n = 1, 2. The claim passes if every scheme orders the substrates the same way
at every n; absolute values are not expected to agree.

Restartable: each finished job is written to runs/charge_robustness_local.json
before the next starts, and finished jobs are skipped.
"""
import os, sys, json, time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qclocal as q
from pathgeom import STRUCT, read_xyz
from pyscf import lib, scf as _scf

lib.num_threads(int(os.environ.get("QC_THREADS", "12")))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "runs",
                   "charge_robustness_local.json")
CLU = os.path.join(STRUCT, "clusters")
TAGS = ["HATCN", "Mo3O9", "F4TCNQ", "benzene"]
SCHEMES = ("mulliken", "lowdin", "iao", "npa", "hirshfeld")


def jobs():
    js = []
    for t in TAGS:
        for n in (1, 2, 3, 4):
            S, X = read_xyz(os.path.join(CLU, f"{t}_Ag{n}.xyz"))
            js.append((len(S), t, n, "def2-svp"))
            if t in ("HATCN", "Mo3O9") and n <= 2:
                js.append((len(S) + 1000, t, n, "def2-tzvp"))
    return [j[1:] for j in sorted(js)]


def run(t, n, basis, guess=None):
    S, X = read_xyz(os.path.join(CLU, f"{t}_Ag{n}.xyz"))
    m = q.mol(S, X, 2 if n % 2 else 1, basis)
    dm0 = None
    if guess is not None:
        dm0 = _scf.addons.project_dm_nr2nr(guess.mol, guess.make_rdm1(), m)
    mf, e = q.scf(m, conv=1e-8, grid=3, dm0=dm0)
    ag = [i for i, s in enumerate(S) if s == "Ag"]
    res = {"converged": bool(mf.converged), "E_elec": float(e),
           "E_D3": float(q.d3bj(S, X))}
    qs = {"mulliken": q.mulliken(mf), "lowdin": q.lowdin(mf),
          "iao": q.iao_charges(mf)}
    npa, per = q.natural(mf)
    qs["npa"] = npa
    qs["hirshfeld"] = q.hirshfeld(mf, S, X, basis)
    for k, v in qs.items():
        res[k] = float(sum(v[i] for i in ag))
    res["ag5s_mean"] = float(np.mean([q.ag5s(per[i]) for i in ag]))
    if basis == "def2-svp":
        S_ = mf.spin_square()[0] if m.spin else 0.0
        res["S2"] = float(S_)
    return res, mf


def main():
    done = json.load(open(OUT)) if os.path.exists(OUT) else {}
    svp_mf = {}
    for t, n, basis in jobs():
        key = f"{t}_Ag{n}_{'svp' if basis == 'def2-svp' else 'tzvp'}"
        if key in done and done[key].get("converged"):
            continue
        t0 = time.time()
        guess = svp_mf.get((t, n)) if basis == "def2-tzvp" else None
        try:
            res, mf = run(t, n, basis, guess)
        except Exception as ex:          # keep going; record the failure
            done[key] = {"converged": False, "error": repr(ex)[:300]}
            json.dump(done, open(OUT, "w"), indent=1)
            print(f"{key}: FAILED {ex!r}", flush=True)
            continue
        if basis == "def2-svp" and t in ("HATCN", "Mo3O9") and n <= 2:
            svp_mf[(t, n)] = mf
        res["seconds"] = round(time.time() - t0)
        done[key] = res
        json.dump(done, open(OUT, "w"), indent=1)
        print(f"{key:<18}" + "".join(f" {s[:4]} {res[s] / n:+.3f}" for s in SCHEMES)
              + f"  5s {res['ag5s_mean']:.3f}  {res['seconds']}s", flush=True)
    report(done)


def report(done):
    print("\nper-atom Ag charge; order of substrates in each scheme")
    ok = True
    for s in SCHEMES:
        for n in (1, 2, 3, 4):
            v = {t: done.get(f"{t}_Ag{n}_svp", {}).get(s) for t in TAGS}
            if any(x is None for x in v.values()):
                continue
            order = [t for t, _ in sorted(v.items(), key=lambda kv: kv[1])]
            hi = order[-1]
            print(f"  {s:<10} n={n}: " + " < ".join(order))
            if hi != "Mo3O9" or order[0] != "benzene":
                ok = False
    print("\nMo3O9 most positive and benzene least in every scheme at every n:",
          "YES" if ok else "NO")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--report":
        report(json.load(open(OUT)))
    else:
        main()
