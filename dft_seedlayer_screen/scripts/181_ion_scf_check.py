"""Make sure the ion single points of 179/180 are the lowest SCF solution.

    QC_THREADS=4 python scripts/181_ion_scf_check.py TAG name [name...]
    TAG = T30 (geometries of 180, results in electrical_T30.json)
        | 300K (geometries of 179, results in electrical.json)

Why: for HATCN (180 geometries) the default-guess UKS anion at the neutral
geometry converged to a solution 0.24 eV above the one reached from the
relaxed anion's density (I@N -1340.37569 vs -1340.38462 Eh). That inflated
the anion relaxation (0.32 instead of 0.077 eV) and lambda (0.434 instead of
0.192 eV), and deflated EA_vert. Molecules with (near-)degenerate LUMOs
(HATCN, and the two-unit heptazine dimers) are exposed to this.

Protocol, per molecule: ion SCF at both geometries from the default guess,
then each again from the other geometry's converged density; keep the lower
energy at each geometry. Neutral (closed shell) values are kept from 179/180.
Writes E_ion@N_min, E_ion@I_min and the corrected EA/IP, lambda, spin PR
into the same json under "scf_checked".
"""
import gc
import importlib.util
import json
import types
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
GEO = os.path.join(ROOT, "structures", "electrical")
HA = 27.211386
_s = importlib.util.spec_from_file_location("e179", os.path.join(HERE, "179_electrical.py"))
e179 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(e179)


def scf(S, X, ch, dm0=None):
    from pyscf import gto, dft
    m = gto.Mole()
    m.atom = [(s, tuple(x)) for s, x in zip(S, X)]
    m.unit = "Angstrom"; m.basis = "def2-svp"; m.charge = ch; m.spin = 1; m.verbose = 0
    m.max_memory = int(os.environ.get("QC_MEM_MB", "3500")); m.build()
    mf = dft.UKS(m).density_fit()
    mf.xc = "pbe0"; mf.grids.level = 3; mf.conv_tol = 1e-8; mf.max_cycle = 100
    if dm0 is None:
        mf.level_shift = 0.1
    mf.kernel(dm0)
    if not mf.converged:
        mf = mf.newton(); mf.kernel(mf.make_rdm1())
    # keep only what is needed and drop the density-fitting tensor (3-4 GB scratch file each)
    out = types.SimpleNamespace(e_tot=mf.e_tot, dm=mf.make_rdm1(), pr=spin_pr(mf, S))
    try:
        mf.with_df.reset(); mf.with_df._cderi_to_save = None   # the temp file closes and is removed
    except AttributeError:
        pass
    del mf; gc.collect()
    return out


def spin_pr(mf, S):
    dm = mf.make_rdm1(); s = mf.mol.intor("int1e_ovlp")
    spin_ao = np.einsum("ij,ji->i", dm[0] - dm[1], s)
    q = np.zeros(mf.mol.natm)
    for i, (a, *_r) in enumerate(mf.mol.ao_labels(fmt=False)):
        q[a] += spin_ao[i]
    w = np.abs(q[[i for i, x in enumerate(S) if x != "H"]]); w = w / w.sum()
    return float(1 / np.sum(w ** 2))


def main(tag, names):
    from pyscf import lib
    lib.num_threads(int(os.environ.get("QC_THREADS", "4")))
    out_p = os.path.join(ROOT, "runs", "design", "electrical_T30.json" if tag == "T30" else "electrical.json")
    sfx = "_T30" if tag == "T30" else ""
    for name in names:
        row = json.load(open(out_p))[name]
        if "scf_checked" in row:
            continue
        ion = +1 if name == "NPB" else -1
        S, XN = e179.read(os.path.join(GEO, f"{name}{sfx}_0.xyz"))
        XI = e179.read(os.path.join(GEO, f"{name}{sfx}_{'cat' if ion > 0 else 'an'}.xyz"))[1]
        aN = scf(S, XN, ion); aI = scf(S, XI, ion)
        bN = scf(S, XN, ion, aI.dm); bI = scf(S, XI, ion, aN.dm)
        eIN = min(aN.e_tot, bN.e_tot); mI = aI if aI.e_tot <= bI.e_tot else bI
        eII = mI.e_tot
        eNN, eNI = row["E_N@N"], row["E_N@I"]
        r = {"E_ion@N": {"default": aN.e_tot, "cross": bN.e_tot},
             "E_ion@I": {"default": aI.e_tot, "cross": bI.e_tot},
             "lambda_eV": round(((eIN - eII) + (eNI - eNN)) * HA, 3),
             "ion_relax_eV": round((eIN - eII) * HA, 3), "neutral_relax_eV": round((eNI - eNN) * HA, 3),
             "spin_PR_heavy_atoms": round(mI.pr, 1)}
        if ion < 0:
            r["EA_vert_eV"] = round((eNN - eIN) * HA, 3); r["EA_ad_eV"] = round((eNN - eII) * HA, 3)
        else:
            r["IP_vert_eV"] = round((eIN - eNN) * HA, 3); r["IP_ad_eV"] = round((eII - eNN) * HA, 3)
        d = json.load(open(out_p)); d[name]["scf_checked"] = r; json.dump(d, open(out_p, "w"), indent=1)
        print(name, tag, "scf_checked", {k: v for k, v in r.items() if not k.startswith("E_")},
              "shift@N", round((aN.e_tot - bN.e_tot) * HA, 3), "shift@I", round((aI.e_tot - bI.e_tot) * HA, 3),
              flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
