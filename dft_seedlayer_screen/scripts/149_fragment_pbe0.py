"""PBE0 binding energy of Ag on a large molecule, from a fragment.

    python scripts/149_fragment_pbe0.py --validate
    python scripts/149_fragment_pbe0.py              # common materials, survey sites

A PBE0/def2-SVP single point on a 80-110 atom molecule with Ag takes hours
locally; the binding is local. So the electronic part of E_b is computed on a
fragment and the dispersion part on the whole molecule:

    E_b = [E(frag) + E(Ag) - E(frag+Ag)]_PBE0, no D3
        + [D3(mol+Ag) - D3(mol)]_full molecule

The fragment keeps every heavy atom within R_CUT of Ag, completes any ring
system it touches (a cut aromatic ring is not the same chemistry), keeps every
substituent chain hanging off what is kept (nitriles, P=O, methyls, amine N),
and caps each remaining cut bond -- always a bond into a further ring -- with H
along the bond. --validate checks the construction
against E_b values already computed on whole molecules at PBE0-D3/def2-SVP.

Ag positions come from the GFN2 survey and are not re-relaxed at PBE0. On the
three cases where both exist (DMABN, Liq, Cu4I4 stage-1 sites) the single point
at the GFN2 position was within 0.03 eV of the PBE0-relaxed E_b.
"""
import os, sys, json, time

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qclocal as q
from pathgeom import read_xyz, STRUCT
from pyscf import lib

lib.num_threads(int(os.environ.get("QC_THREADS", "12")))
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "runs", "common_pbe0_fragment.json")
EV = 27.211386
R_CUT = float(os.environ.get("FRAG_RCUT", "7.0"))
CAP = {"C": 1.09, "N": 1.01, "O": 0.97, "P": 1.42, "Si": 1.48, "Al": 1.58, "B": 1.19}
q.Z.update({"Si": 14, "Al": 13})


def connectivity(S, X):
    xyz = f"{len(S)}\n\n" + "\n".join(f"{s} {x[0]} {x[1]} {x[2]}" for s, x in zip(S, X))
    m = Chem.MolFromXYZBlock(xyz)
    rdDetermineBonds.DetermineConnectivity(m)
    Chem.FastFindRings(m)
    return m


def fragment(S, X, ag):
    m = connectivity(S, X)
    n = len(S)
    heavy = [i for i in range(n) if S[i] != "H"]
    keep = {i for i in heavy if np.linalg.norm(X[i] - ag) < R_CUT}
    # complete ring systems (fused rings as one unit)
    rings = [set(r) for r in m.GetRingInfo().AtomRings()]
    systems = []
    for r in rings:
        merged = [s for s in systems if s & r]
        for s in merged:
            systems.remove(s)
            r |= s
        systems.append(r)
    changed = True
    while changed:
        changed = False
        for s in systems:
            if s & keep and not s <= keep:
                keep |= s
                changed = True
    # Substituent chains are never cut: every non-ring heavy atom reachable from
    # a kept atom without passing through a ring atom is kept (C#N, P=O and its
    # P, CH3, the N of a triarylamine), because cutting a nitrile off an
    # acceptor changes the acceptor. Rings beyond such a chain are cut.
    ring_atoms = set().union(*systems) if systems else set()
    frontier = list(keep)
    while frontier:
        i = frontier.pop()
        for nb in m.GetAtomWithIdx(i).GetNeighbors():
            j = nb.GetIdx()
            if j in keep or S[j] == "H" or j in ring_atoms:
                continue
            keep.add(j)
            frontier.append(j)
    for i in list(keep):
        for nb in m.GetAtomWithIdx(i).GetNeighbors():
            if S[nb.GetIdx()] == "H":
                keep.add(nb.GetIdx())
    FS, FX = [], []
    for i in sorted(keep):
        FS.append(S[i])
        FX.append(X[i])
    ncap = 0
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        for a, c in ((i, j), (j, i)):
            if a in keep and c not in keep:
                u = X[c] - X[a]
                u /= np.linalg.norm(u)
                FS.append("H")
                FX.append(X[a] + CAP.get(S[a], 1.09) * u)
                ncap += 1
    return FS, np.array(FX), ncap


def eb(S, X, ag, mult_mol=1):
    FS, FX, ncap = fragment(S, X, ag)
    e_cx = q.scf(q.mol(FS + ["Ag"], np.vstack([FX, ag]), 2 if mult_mol == 1 else 1))
    e_fr = q.scf(q.mol(FS, FX, mult_mol))
    e_ag = q.scf(q.mol(["Ag"], [[0, 0, 0]], 2))
    for mf, _ in (e_cx, e_fr, e_ag):
        if not mf.converged:
            raise RuntimeError("SCF not converged")
    elec = (e_fr[1] + e_ag[1] - e_cx[1]) * EV
    d3 = (q.d3bj(S + ["Ag"], np.vstack([X, ag])) - q.d3bj(S, X)) * EV
    qag = float(q.mulliken(e_cx[0])[-1])
    return {"E_b": round(elec + d3, 4), "elec": round(elec, 4), "D3": round(d3, 4),
            "q_mulliken": round(qag, 4), "frag_atoms": len(FS), "caps": ncap,
            "mol_atoms": len(S)}


def validate():
    cases = [("HATCN aza", "clusters/HATCN_Ag1.xyz", 1.631),
             ("Bphen chelate", "Bphen_Ag_pbe0.xyz", 0.621),
             ("B3PyMPM", "B3PyMPM_Ag_pbe0.xyz", 0.416),
             ("p-bPPhenB", "pbPPhenB_Ag_chelate_pbe0.xyz", 0.785)]
    done = json.load(open(OUT)) if os.path.exists(OUT) else {}
    val = done.get("_validation", {})
    for lab, fn, ref in cases:
        if lab in val:
            continue
        p = os.path.join(STRUCT, fn)
        if not os.path.exists(p):
            print(f"{lab}: {fn} missing")
            continue
        S, X = read_xyz(p)
        X = np.asarray(X, float)
        k = S.index("Ag")
        ag = X[k]
        S2 = [s for i, s in enumerate(S) if i != k]
        X2 = np.delete(X, k, axis=0)
        t0 = time.time()
        r = eb(S2, X2, ag)
        r["reference_full_pbe0"] = ref
        r["error"] = round(r["E_b"] - ref, 4)
        r["seconds"] = round(time.time() - t0)
        val[lab] = r
        done["_validation"] = val
        json.dump(done, open(OUT, "w"), indent=1)
        print(f"{lab:<14} fragment {r['frag_atoms']}/{r['mol_atoms']} atoms  E_b {r['E_b']:.3f} "
              f"(ref {ref:.3f}, error {r['error']:+.3f})  {r['seconds']}s", flush=True)


def run(names):
    sv = json.load(open(os.path.join(ROOT, "runs", "common_survey_xtb.json")))
    for name in names:
        done = json.load(open(OUT)) if os.path.exists(OUT) else {}
        if name in done or name not in sv:
            continue
        S, X = read_xyz(os.path.join(ROOT, "structures", "common", f"{name}.xyz"))
        X = np.asarray(X, float)
        row = {}
        # the GFN2-deepest site, plus the next one only if it binds a different
        # element: GFN2's overbinding depends on the element (nitrile N 4-5x,
        # carbon 1.2x), so its order is trusted within an element, not across
        picks = sv[name]["pbe0_sites"][:1]
        first_el = sv[name]["sites"][picks[0]]["near_el"]
        for lab in sv[name]["pbe0_sites"][1:]:
            if sv[name]["sites"][lab]["near_el"] != first_el:
                picks.append(lab)
                break
        for lab in picks:
            ag = np.array(sv[name]["sites"][lab]["pos"])
            t0 = time.time()
            try:
                r = eb(S, X, ag)
            except Exception as ex:
                r = {"error": repr(ex)[:200]}
            r["seconds"] = round(time.time() - t0)
            r["near"] = sv[name]["sites"][lab]["near"]
            r["R_cut"] = R_CUT
            row[lab] = r
            print(f"{name:<9} {lab:<22} near {r['near']:<5} E_b "
                  f"{r.get('E_b', float('nan')):.3f}  frag {r.get('frag_atoms')}  {r['seconds']}s",
                  flush=True)
        good = {k: v["E_b"] for k, v in row.items() if "E_b" in v}
        if good:
            row["_deepest"] = max(good, key=good.get)
            row["_E_b_deepest"] = good[row["_deepest"]]
        done = json.load(open(OUT)) if os.path.exists(OUT) else {}
        done[name] = row
        json.dump(done, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    if "--validate" in sys.argv:
        validate()
    else:
        names = [a for a in sys.argv[1:] if not a.startswith("--")]
        row_rcut = R_CUT
        if not names:
            names = sorted(json.load(open(os.path.join(ROOT, "runs", "common_survey_xtb.json"))))
        run(names)
