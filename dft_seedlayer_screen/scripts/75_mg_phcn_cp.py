"""Finish Mg/PhCN: counterpoise fragments only, with the robust guess.

scripts/71 computed the whole Mg-PhCN rigid scan (minimum at r = 3.00 A) and
then could not converge the two counterpoise fragments: with guess = SAD both
the ghost-substrate and ghost-metal SCFs threw, so the run stalled one step
short of E_b. scripts/74 had already shown GWH to be the robust guess on this
box, so that is the only thing changed here. A guess cannot move a converged
energy, so the CP numbers stay comparable with every other row of the screen.

Reuses 71's machinery and its checkpoint file (71 itself is stopped, so there
is no second writer). Only Mg_PhCN is requested; its scan points are already
checkpointed and get skipped, so this is two SCFs, not nine.
"""
import importlib.util
import os

import psi4

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.abspath(os.path.join(HERE, "..", "runs"))

spec = importlib.util.spec_from_file_location(
    "s71", os.path.join(HERE, "71_dopant_substrate_Eb.py"))
s71 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s71)

psi4.core.set_output_file(os.path.join(RUNS, "psi4_mg_cp.out"), False)
psi4.set_num_threads(3)
psi4.set_memory("6 GB")
s71.base_opts["guess"] = "gwh"

if __name__ == "__main__":
    row = [s for s in s71.SYSTEMS if s[0] == "Mg_PhCN"][0]
    s71.run_system(*row)
    print("MG CP DONE", flush=True)
