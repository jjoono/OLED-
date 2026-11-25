"""Electrical figures of merit for a HIL/seed candidate vs HATCN.

    QC_THREADS=4 python scripts/179_electrical.py [names...]   -> runs/design/electrical.json

Role in the inverted TEOLED (HTL / HIL=seed / Ag anode) and in a CGL: the HIL
takes an electron from the HTL HOMO (charge generation) and carries it on its
LUMO to the Ag. So the numbers that matter are
  EA      electron affinity, vertical and adiabatic: depth of the LUMO, i.e.
          the driving force for HTL -> HIL electron transfer;
  dG_CT   IP(NPB) - EA(HIL), gas phase: the charge-generation driving force
          relative to HATCN (film polarisation shifts both by about the same);
  lambda  inner-sphere reorganisation energy for electron hopping (four-point
          scheme): enters the Marcus rate as exp(-lambda/4kT);
  PR      participation ratio of the anion's spin density over heavy atoms:
          how delocalised the transported electron is.
Geometries: GFN2-xTB (neutral and anion, or cation for NPB) because geomeTRIC
could not be installed here; energies: PBE0/def2-SVP (DF). The same protocol
for every molecule, HATCN included, so differences are the result. Transfer
integrals need film packing (MD) and are not computed.
"""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "runs", "design", "electrical.json")
GEO = os.path.join(ROOT, "structures", "electrical")
HA = 27.211386
SRC = {"NPB": os.path.join(ROOT, "structures", "common", "NPB.xyz")}


def read(path):
    L = open(path).read().split("\n")
    n = int(L[0]); S, X = [], []
    for l in L[2:2 + n]:
        p = l.split(); S.append(p[0]); X.append([float(v) for v in p[1:4]])
    return S, np.array(X)


def write(path, S, X, title):
    with open(path, "w") as f:
        f.write(f"{len(S)}\n{title}\n" + "".join(f"{s} {x[0]:.6f} {x[1]:.6f} {x[2]:.6f}\n" for s, x in zip(S, X)))


def gfn2_opt(S, X, charge):
    from ase import Atoms
    from ase.optimize import BFGS
    from tblite.ase import TBLite
    a = Atoms(S, positions=X)
    a.calc = TBLite(method="GFN2-xTB", charge=charge, multiplicity=1 + abs(charge) % 2, verbosity=0,
                    max_iterations=500)
    BFGS(a, logfile=None, maxstep=0.15).run(fmax=0.03, steps=800)
    return a.positions.copy(), float(a.get_potential_energy())


def pbe0(S, X, charge):
    from pyscf import gto, dft
    m = gto.Mole(); m.atom = [(s, tuple(x)) for s, x in zip(S, X)]; m.unit = "Angstrom"
    m.basis = "def2-svp"; m.charge = charge; m.spin = abs(charge) % 2; m.verbose = 0
    m.max_memory = int(os.environ.get("QC_MEM_MB", "4000")); m.build()
    mf = (dft.UKS(m) if m.spin else dft.RKS(m)).density_fit()
    mf.xc = "pbe0"; mf.grids.level = 3; mf.conv_tol = 1e-8; mf.max_cycle = 80; mf.level_shift = 0.1
    mf.kernel()
    if not mf.converged:
        mf = mf.newton(); mf.kernel(mf.make_rdm1())
    pr = None
    if m.spin:
        from pyscf.scf import uhf
        pa, pb = mf.mulliken_pop(verbose=0)[0] if False else (None, None)
        dm = mf.make_rdm1(); s = m.intor("int1e_ovlp")
        spin_ao = np.einsum("ij,ji->i", dm[0] - dm[1], s)
        q = np.zeros(m.natm)
        for i, (a, *_r) in enumerate(m.ao_labels(fmt=False)):
            q[a] += spin_ao[i]
        heavy = [i for i, sym in enumerate(S) if sym != "H"]
        w = np.abs(q[heavy]); w = w / w.sum()
        pr = float(1 / np.sum(w ** 2))
    return float(mf.e_tot), bool(mf.converged), pr


def main(names):
    from pyscf import lib
    lib.num_threads(int(os.environ.get("QC_THREADS", "4")))
    os.makedirs(GEO, exist_ok=True)
    for name in names:
        out = json.load(open(OUT)) if os.path.exists(OUT) else {}
        if "lambda_eV" in out.get(name, {}) and not os.environ.get("REDO"):
            continue
        ion = +1 if name == "NPB" else -1
        S, XN0 = read(SRC.get(name, os.path.join(ROOT, "structures", "design", f"{name}.xyz")))
        t0 = time.time()
        gN = os.path.join(GEO, f"{name}_0.xyz"); gI = os.path.join(GEO, f"{name}_{'cat' if ion > 0 else 'an'}.xyz")
        prev = out.get(name, {})
        if (os.path.exists(gN) and os.path.exists(gI) and "E_N@N" in prev and "gfn2_adiabatic_eV" in prev
                and not os.environ.get("REDO")):
            # resume after an interrupted run: same geometries, keep finished single points
            XN, XI = read(gN)[1], read(gI)[1]
            r = {"gfn2_adiabatic_eV": prev["gfn2_adiabatic_eV"]}
            e = {t: prev[f"E_{t}"] for t in ("N@N", "I@N", "I@I", "N@I") if f"E_{t}" in prev}
            print(name, "resume with", sorted(e), flush=True)
        else:
            XN, XI, r, e = None, None, None, {}
        if XN is None:
            XN, eN_x = gfn2_opt(S, XN0, 0)
        # a random 0.03 A kick before the ion relaxation: from an exactly
        # symmetric start (HATCN, D3h, degenerate LUMO) BFGS cannot leave the
        # symmetric point and misses the Jahn-Teller distortion, which
        # underestimates lambda (first HATCN run: 0.034 eV)
            kick = float(os.environ.get("ION_KICK", "0.03"))
            XI, eI_x = gfn2_opt(S, XN + np.random.default_rng(1).normal(0, kick, XN.shape), ion)
            write(gN, S, XN, f"{name} neutral GFN2")
            write(gI, S, XI, f"{name} ion GFN2")
            r = {"gfn2_adiabatic_eV": round(abs(eI_x - eN_x), 3)}
        for tag, X, ch in (("N@N", XN, 0), ("I@N", XN, ion), ("I@I", XI, ion), ("N@I", XI, 0)):
            if tag in e:
                continue
            e[tag], conv, pr = pbe0(S, X, ch)
            r[f"conv_{tag}"] = conv
            if pr is not None and tag == "I@I":
                r["spin_PR_heavy_atoms"] = round(pr, 1)
            out = json.load(open(OUT)) if os.path.exists(OUT) else {}
            out.setdefault(name, {}).update({**r, f"E_{tag}": e[tag]})
            json.dump(out, open(OUT, "w"), indent=1)
        sign = -1 if ion < 0 else 1   # EA = E(N) - E(A);  IP = E(C) - E(N)
        r["vertical_eV"] = round(sign * (e["N@N"] - e["I@N"]) * HA * (-1 if ion > 0 else 1) * -1, 3)
        if ion < 0:
            r["EA_vert_eV"] = round((e["N@N"] - e["I@N"]) * HA, 3)
            r["EA_ad_eV"] = round((e["N@N"] - e["I@I"]) * HA, 3)
        else:
            r["IP_vert_eV"] = round((e["I@N"] - e["N@N"]) * HA, 3)
            r["IP_ad_eV"] = round((e["I@I"] - e["N@N"]) * HA, 3)
        r.pop("vertical_eV")
        r["lambda_eV"] = round(((e["I@N"] - e["I@I"]) + (e["N@I"] - e["N@N"])) * HA, 3)
        r["seconds"] = round(time.time() - t0)
        out = json.load(open(OUT)) if os.path.exists(OUT) else {}
        out.setdefault(name, {}).update(r)
        json.dump(out, open(OUT, "w"), indent=1)
        print(name, {k: v for k, v in r.items() if not k.startswith("conv")}, flush=True)
    out = json.load(open(OUT)) if os.path.exists(OUT) else {}
    if "NPB" in out and "IP_ad_eV" in out["NPB"]:
        ip = out["NPB"]["IP_ad_eV"]
        for n, v in out.items():
            if "EA_ad_eV" in v:
                v["dG_CT_vs_NPB_gas_eV"] = round(ip - v["EA_ad_eV"], 3)
        json.dump(out, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1:] or ["HATCN", "Hz-NCN-Hz", "NPB"])
