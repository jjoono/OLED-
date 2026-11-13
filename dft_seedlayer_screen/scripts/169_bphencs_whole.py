"""Ag on Bphen(-)Cs(+): the whole-molecule PBE0 single points the Gaussian
package (bphencs_jobs*) was meant to deliver, run in the cloud.

    QC_THREADS=4 python scripts/169_bphencs_whole.py   (resumable)
    -> runs/bphencs_whole.json, SCF chkfiles in runs/bphencs_whole_ckpt/ (git-ignored)

Script 155 placed Ag (PBE-D3 relaxation on a 4 A fragment that keeps the Cs)
and scored it on that 23-33 atom fragment. Here the same Ag positions are
scored on the whole 43-atom Bphen-Cs molecule, so the fragment cut drops out
and the numbers are on the footing of the workstation Gaussian results.

    E_b = E(BphenCs, doublet) + E(Ag) - min_spin E(BphenCs+Ag)     (PBE0-D3BJ/def2-SVP)

D3(BJ) is on the whole molecule in both 155 and here, so only the electronic
part changes. Complex spin: singlet (closed-shell RKS, the radical pair
couples) and triplet, the lower one taken, as in 155. Ag2 at the Cs site is a
doublet (odd electron count).

SCF: converged PBE density as the start, 40 DIIS cycles, then second-order
SCF (the acceptor/radical complexes slosh under DIIS: scripts 157, 164).
"""
import importlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
XYZ = os.path.join(ROOT, "structures", "common", "BphenCs_opt.xyz")
FRAG = os.path.join(ROOT, "runs", "bphencs_fragment.json")
OUT = os.path.join(ROOT, "runs", "bphencs_whole.json")
CK = os.path.join(ROOT, "runs", "bphencs_whole_ckpt")
EV = 27.211386


def pbe_guess(m):
    from pyscf import dft
    mf = (dft.UKS(m) if m.spin else dft.RKS(m)).density_fit()
    mf.xc = "PBE"
    mf.grids.level = 2
    mf.conv_tol = 1e-6
    mf.max_cycle = 150
    mf.level_shift = 0.3
    mf.damp = 0.3
    mf.kernel()
    return mf.make_rdm1()


def scf_ck(m, tag):
    from pyscf import dft
    os.makedirs(CK, exist_ok=True)
    ck = os.path.join(CK, tag + ".chk")
    mf = (dft.UKS(m) if m.spin else dft.RKS(m)).density_fit()
    mf.xc = "PBE0"
    mf.grids.level = 4
    mf.conv_tol = 1e-9
    mf.level_shift = 0.2
    dm0 = mf.from_chk(ck) if os.path.exists(ck) else pbe_guess(m)
    mf.chkfile = ck
    mf.max_cycle = 40
    e = mf.kernel(dm0)
    if not mf.converged:
        mf = mf.newton()
        mf.chkfile = ck
        mf.max_cycle = 60
        e = mf.kernel(mf.make_rdm1())
    return mf, e


def load():
    return json.load(open(OUT)) if os.path.exists(OUT) else {}


def put(key, val):
    d = load()
    d[key] = val
    json.dump(d, open(OUT, "w"), indent=1)


def main():
    from pyscf import lib
    lib.num_threads(int(os.environ.get("QC_THREADS", "4")))
    q = importlib.import_module("qclocal")
    from pathgeom import read_xyz
    S, X = read_xyz(XYZ)
    X = np.asarray(X, float)
    frag = json.load(open(FRAG))
    e_ag = frag["E_ag"]
    icx = S.index("Cs")

    def job(key, Sx, Xx, mult):
        if key in load():
            return load()[key]
        t0 = time.time()
        mf, e = scf_ck(q.mol(Sx, Xx, mult), key)
        qq = q.mulliken(mf)
        r = {"E": float(e), "conv": bool(mf.converged), "sec": round(time.time() - t0),
             "q_Ag": [round(float(v), 4) for v, s in zip(qq, Sx) if s == "Ag"],
             "q_Cs": round(float(qq[icx]), 4)}
        if mf.converged and m_spin(mf):
            r["S2"] = round(float(mf.spin_square()[0]), 4)
        put(key, r)
        print(key, {k: v for k, v in r.items()}, flush=True)
        return r

    def m_spin(mf):
        return getattr(mf.mol, "spin", 0) > 0

    mol = job("mol", S, X, 2)
    # deepest first, then the second spin of each, then Ag2
    order = [("Cs4", 1), ("hollow8", 1), ("hollow3", 3), ("Cs4", 3), ("hollow8", 3), ("hollow3", 1)]
    if os.environ.get("SKIP_SPIN_CHECK"):
        # the fragment runs (155) already fixed the spin ordering with margins
        # of 0.18-0.86 eV, and whole vs fragment agree to 0.01 eV: run Ag2
        # first, the second-spin checks later
        order = order[:3]
    for site, mult in order:
        ag = np.array(frag["sites"][site]["pos"])
        job(f"{site}_m{mult}", S + ["Ag"], np.vstack([X, ag]), mult)
    a2 = frag["ag2"]
    r2 = job("ag2_Cs4", S + ["Ag", "Ag"], np.vstack([X, np.array(a2["pos"])]), 2)

    # summary
    d = load()
    summ = {}
    for site in ("Cs4", "hollow8", "hollow3"):
        es = {m: d[f"{site}_m{m}"]["E"] for m in (1, 3) if f"{site}_m{m}" in d and d[f"{site}_m{m}"]["conv"]}
        if not es:
            continue
        m_best = min(es, key=es.get)
        ag = np.array(frag["sites"][site]["pos"])
        d3 = (q.d3bj(S + ["Ag"], np.vstack([X, ag])) - q.d3bj(S, X)) * EV
        elec = (mol["E"] + e_ag - es[m_best]) * EV
        summ[site] = {"E_b": round(elec - d3, 4), "elec": round(elec, 4), "disp": round(-d3, 4),
                      "state": "singlet" if m_best == 1 else "triplet",
                      "singlet_minus_triplet_eV": round((es[1] - es[3]) * EV, 4) if len(es) == 2 else None,
                      "q_Ag": d[f"{site}_m{m_best}"]["q_Ag"][0],
                      "fragment_155": frag["sites"][site]["E_b"],
                      "whole_minus_fragment": round(elec - d3 - frag["sites"][site]["E_b"], 4)}
    if "ag2_Cs4" in d and "Cs4" in summ:
        q2 = sum(d["ag2_Cs4"]["q_Ag"])
        summ["ag2"] = {"q_Ag2_total": round(q2, 4), "q_Ag1": summ["Cs4"]["q_Ag"],
                       "dq_Ag1_Ag2": round(q2 - summ["Cs4"]["q_Ag"], 4)}
    put("summary", summ)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
