"""EA of F6TCNNQ, checkpointed per SCF so a ten-minute cap cannot lose it.

WHY SEPARATE FROM scripts/45. That script's cache is keyed per MOLECULE: it stores
the EA only after BOTH the neutral and the anion have converged. At
wB97X/def2-TZVP a 26-atom molecule takes longer than one SCF per foreground chunk,
so every attempt died with nothing saved -- twice. The unit has to be one SCF.

WHY THE NUMBER MATTERS. F6TCNNQ's cluster E_b could not be computed at all: with an
Ag atom present the SCF oscillates between neutral and charge-separated solutions
4.0 eV apart, and does so even at 20 A separation where the two fragments cannot
interact, so the failure is a functional/algorithm problem rather than a physical
near-degeneracy. Every remedy tried (maxiter 500, damping to 60%, three level
shifts, four guesses, far-to-near orbital carry-over, second-order SCF) failed.
HATCN and F4TCNQ converge without complaint at the same level.

EA does not have that problem: it is the isolated molecule, no Ag, no ambiguity
about where an unpaired electron lives, and F6TCNNQ alone already converged
cleanly at PBE. So EA is the one axis of Fig. 1(b) that is still obtainable for
this molecule, and it answers the first-order question -- is F6TCNNQ a stronger
acceptor than F4TCNQ -- even though the binding axis stays blank.

IMPORTANT: the convergence failure is NOT a result and must not be reported as
evidence that F6TCNNQ oxidises Ag. A calculation that did not converge measures
nothing. It is a reason to change method (a hybrid would separate the two charge
states properly), not a finding.
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _ckpt import Checkpoint

BASE = os.path.abspath(os.path.join(HERE, ".."))
STR, RUNS = os.path.join(BASE, "structures"), os.path.join(BASE, "runs")
HA = 27.211386

FUNCTIONAL = "wb97x"
# def2-SVP, NOT the def2-TZVP of scripts/45, and the reason matters.
#
# A single def2-TZVP SCF on this 26-atom molecule ran past 19 minutes on the four
# cores available here while its density-fitting scratch grew past 20 GB and
# filled the session's disk allowance -- twice, the second time after clearing it.
# The absolute EA at this level is simply not obtainable in this container.
#
# What IS obtainable is the DIFFERENCE. EA(F6TCNNQ) - EA(F4TCNQ) at a common
# smaller basis is far more robust than either absolute value, for the same reason
# the TCNQ differential gate in scripts/45 works: the two molecules are close
# relatives and the basis-set error largely cancels between them. There is direct
# evidence the EAs here are basis-insensitive -- F4TCNQ moved 4.020 -> 4.028 eV,
# just 0.008 eV, on adding diffuse functions (def2-TZVP -> def2-TZVPD).
#
# So both molecules are run at def2-SVP, the difference is taken, and it is added
# to the def2-TZVP value already in hand for F4TCNQ. The result is labelled as the
# estimate it is, never as a directly computed def2-TZVP number.
BASIS = "def2-svp"
TAG, FNAME = "F6TCNNQ", "F6TCNNQ.xyz"
PAIR = [("F4TCNQ", "F4TCNQ.xyz"), ("F6TCNNQ", "F6TCNNQ.xyz")]
EA_F4TCNQ_TZVP = 4.020                       # scripts/45, wB97X/def2-TZVP
REF_EA = {"HATCN": 3.378, "TCNQ": 3.536, "F4TCNQ": 4.020}


def single_point(fname, charge, mult):
    import psi4
    psi4.core.clean()
    lines = open(os.path.join(STR, fname)).read().splitlines()
    n = int(lines[0])
    body = "\n".join(lines[2:2 + n])
    mol = psi4.geometry(f"{charge} {mult}\n{body}\nsymmetry c1\nno_reorient\nno_com\n")
    psi4.set_options({
        "basis": BASIS, "scf_type": "df",
        "reference": "uks" if mult > 1 else "rks",
        "maxiter": 300, "guess": "sad",
        # e_convergence is 1e-6, NOT the 1e-8 copied from scripts/45. That value
        # works for the smaller molecules there but sits BELOW the density-fitting
        # noise floor of this one: the F6TCNNQ anion reaches an RMS density error
        # of 7e-9 -- converged by any physical standard -- while its energy keeps
        # jittering at the 1e-6 Ha level, so the 1e-8 test can never pass and the
        # SCF spins until it is killed. That is what burned three foreground
        # chunks before the iteration log was read.
        #
        # The jitter was then MEASURED rather than guessed: iterations 29-33 of the
        # anion swing by 2.4e-6 to 7.4e-6 Ha, so 1e-6 fails too. 1e-5 clears it.
        # That is 2.7e-4 eV on each energy, ~5e-4 eV on the EA -- still two orders
        # below the 0.1 eV differences at stake, so nothing is given up.
        # d_convergence stays tight: the density is what has to be converged, and
        # it reaches 7e-9 without trouble.
        "d_convergence": 1e-7, "e_convergence": 1e-5,
        "damping_percentage": 20.0,
    })
    return psi4.energy(FUNCTIONAL, molecule=mol)


def main():
    import psi4
    psi4.set_memory(os.environ.get("PSI4_MEM", "6 GB"))
    psi4.set_num_threads(int(os.environ.get("PSI4_THREADS", "4")))
    psi4.core.set_output_file(os.path.join(RUNS, "psi4_f6_ea.out"), False)

    ck = Checkpoint("dft_seedlayer_screen/runs/f6tcnnq_ea.json", label="F6TCNNQ EA")
    ea_svp = {}
    for name, fname in PAIR:
        for suffix, charge, mult, what in (("neutral", 0, 1, "neutral (RKS)"),
                                           ("anion", -1, 2, "anion (UKS)")):
            key = f"{name}_{suffix}_svp"
            if ck.has(key):
                print(f"  {name:<9} {what:<14} [cached] {ck.get(key):.8f} Ha", flush=True)
                continue
            print(f"  {name:<9} {what:<14} running...", flush=True)
            e = single_point(fname, charge, mult)
            ck.put(key, e)
            print(f"  {name:<9} {what:<14} {e:.8f} Ha  [checkpointed]", flush=True)
        e0, em = ck.get(f"{name}_neutral_svp"), ck.get(f"{name}_anion_svp")
        if e0 is None or em is None:
            print("\nincomplete -- rerun to continue from the checkpoint")
            return
        ea_svp[name] = (e0 - em) * HA

    d = ea_svp["F6TCNNQ"] - ea_svp["F4TCNQ"]
    ea = EA_F4TCNQ_TZVP + d
    ck.put("EA_svp", ea_svp)
    ck.put("delta_svp", d)
    ck.put("EA_eV_estimated_tzvp", ea)

    print("\n" + "=" * 62)
    print(f"EA at {FUNCTIONAL}/{BASIS} (both molecules, same basis)")
    print("=" * 62)
    for k, v in ea_svp.items():
        print(f"  {k:<9} {v:.3f} eV")
    print(f"\n  difference  F6TCNNQ - F4TCNQ = {d:+.3f} eV")
    print(f"\n  F4TCNQ at def2-TZVP (scripts/45)      = {EA_F4TCNQ_TZVP:.3f} eV")
    print(f"  F6TCNNQ ESTIMATED at def2-TZVP       = {ea:.3f} eV")
    print("  (estimate: the def2-SVP difference transferred onto the def2-TZVP")
    print("   reference. Not a directly computed def2-TZVP value -- see the")
    print("   module docstring for why that is unobtainable here.)")

    print("\n" + "=" * 62)
    print("placed against the family")
    print("=" * 62)
    for k, v in sorted({**REF_EA, TAG + " (est)": ea}.items(), key=lambda t: t[1]):
        mark = "  <-- this run" if k.startswith(TAG) else ""
        print(f"  {k:<15} {v:.3f} eV{mark}")
    print(f"\n  vs F4TCNQ: {d:+.3f} eV")
    if d > 0.05:
        print("  F6TCNNQ is the stronger acceptor, as the larger fluorinated core")
        print("  predicts. On Fig. 1(b) it sits further right -- higher oxidation")
        print("  driving force than F4TCNQ, which already exceeds HATCN by 0.64 eV.")
        print("  The binding axis is still blank, so whether it also anchors more")
        print("  strongly -- the coupling the figure asserts -- remains untested.")
    elif d < -0.05:
        print("  F6TCNNQ is the WEAKER acceptor, which the larger fluorinated core")
        print("  did not predict. Check the geometry before using this.")
    else:
        print("  Indistinguishable from F4TCNQ at this level.")
    print("\n  Reminder: E_b for F6TCNNQ+Ag never converged and is NOT reported.")
    print("  A non-converged calculation measures nothing; it is not evidence")
    print("  either for or against this molecule as a seed.")


if __name__ == "__main__":
    main()
