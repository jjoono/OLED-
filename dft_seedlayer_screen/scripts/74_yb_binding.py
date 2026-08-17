"""Yb-substrate / Ag-Yb binding: the Yb half of scripts/71, rescued.

WHY A SEPARATE SCRIPT. scripts/71 cannot do Yb at all, and the reason is not
SCF convergence (the risk we planned for) but something simpler: psi4's
bundled def2 basis sets stop before the lanthanides, so any Yb geometry dies
at BasisSetNotFound before an SCF is ever attempted.

THE FIX, and why it costs nothing in comparability. def2-SVP and def2-ECP ARE
defined for Yb (Weigend et al.); only psi4's shipped copy omits them. They are
pulled from the Basis Set Exchange into def2svp_yb.gbs (regenerated here if
missing) and used for every element, so the level of theory is IDENTICAL to
the rest of the screen: PBE-D3BJ/def2-SVP, def2-ECP on Ag and Yb,
counterpoise-corrected rigid scan.

VALIDATED BEFORE USE: the BSE file reproduces psi4's bundled def2-SVP on
Mg2 at 3.6 A to 0.03 uHa (-399.751763 both ways), so numbers from this script
sit on the same scale as Mg2 = 0.216 eV and Ag/PhCN = 0.29 eV from 71/30.

Yb is [Xe]4f14 6s2 -- a FULL 4f shell, so the def2 small-core ECP28 (4f in
valence) is chemically the right choice and RKS is tried first; UKS singlet
and triplet follow as fallbacks. Separate checkpoint file from 71 so the two
processes can never clobber each other's JSON.
"""
import importlib.util
import os
import sys

import basis_set_exchange as bse
import psi4

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(os.path.join(HERE, ".."))
RUNS = os.path.join(BASE, "runs")
GBS = os.path.join(BASE, "def2svp_yb.gbs")

if not os.path.exists(GBS):
    open(GBS, "w").write(bse.get_basis(
        "def2-svp", elements=["H", "C", "N", "Mg", "Ag", "Yb"], fmt="psi4"))

# import 71's machinery (its run loop is guarded by __main__)
spec = importlib.util.spec_from_file_location(
    "s71", os.path.join(HERE, "71_dopant_substrate_Eb.py"))
s71 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s71)

# override what 71 set at import time: own output file, own basis, fewer
# threads (scripts/71 may still be running on the same 4-core box)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_yb.out"), False)
psi4.set_num_threads(3)
psi4.set_memory("8 GB")
s71.base_opts["basis"] = "def2svp_yb"
# SAD is psi4's default guess and 71 sets it explicitly, but the SAD-FIT
# auxiliary basis is undefined for Yb -- that, not the SCF itself, is what
# killed the first Yb attempt (error names SAD-FIT, not BASIS). GWH needs no
# per-element fit basis and converges Ag-Yb at 3.0 A on the first try
# (E = -1306.525174 Ha). The guess cannot change a converged energy, only
# which solution is reached, so this stays comparable to the SAD-guessed rows.
s71.base_opts["guess"] = "gwh"

from _ckpt import Checkpoint                                    # noqa: E402
s71.ck = Checkpoint(os.path.join(RUNS, "yb_binding_eV.json"), label="Yb binding")

SYSTEMS = [
    # Ag-Yb dimer first: 2 atoms, answers "does Yb wet silver" immediately
    ("AgYb",    lambda r: s71.dimer("Ag", "Yb", r), [2], 2, 1,
     [2.6, 2.8, 3.0, 3.2, 3.5, 3.8]),
    ("Yb_PhCN", lambda r: s71.phcn_complex("Yb", r), [1, 3], 1, 1,
     [2.3, 2.5, 2.7, 3.0, 3.3, 3.7]),
]

if __name__ == "__main__":
    print(f"basis file: {GBS}", flush=True)
    for name, build, cmults, multA, multB, rs in SYSTEMS:
        s71.run_system(name, build, cmults, multA, multB, rs)
        s71.clear_scratch()
    print("YB DONE", flush=True)
