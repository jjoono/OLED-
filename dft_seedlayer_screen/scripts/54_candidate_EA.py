"""EA of the low-crystallisation candidates, as the validated anchor-strength proxy.

WHY EA AND NOT E_b. scripts/53 established, over three end-on fragments spanning
donor to acceptor substitution, that EA tracks the counterpoise binding energy with
r = 0.962. That makes EA a usable stand-in, and it is the only one available here:
E_b needs an Ag atom, which for the strongest acceptors makes the SCF oscillate
between neutral and charge-separated solutions and never converge (F6TCNNQ,
scripts/52). EA is the isolated molecule -- no Ag, no ambiguity, and a fraction of
the cost.

WHAT IS BEING RANKED. The sp3 and spiro architectures that scripts/53 found can
keep nitrile anchors while halving the crystallisation risk:

    TCPM    tetraphenylmethane + 4 benzonitrile arms   risk 35, 4x CN
    TCPSi   the silicon analogue                        risk 35, 4x CN
    SBF2CN  spirobifluorene dicarbonitrile              risk 32, 2x CN

HATCN is computed at the SAME basis as the reference point. Its def2-TZVP EA is
already known (3.378 eV) but a def2-SVP number is needed here, because a ranking
is only meaningful within one basis.

READ THE SLOPE BEFORE READING THE RANKING. The fitted proxy slope is 0.055 eV of
E_b per eV of EA. Even a 2 eV spread in EA is worth about 0.1 eV of binding. So
this ranks the candidates against each other; it does NOT predict that any of them
approaches HATCN's 1.03 eV, which comes mostly from having six nitriles and from
the inter-molecular pocket, not from core electronics. Treat the output as "which
of these is the better acceptor", not as a binding-energy prediction.

ORDER is cheapest-first so that a container rollback or a ten-minute cap leaves
the most useful subset done: the HATCN reference is worthless without candidates,
but a candidate is worthless without the reference.
"""
import os, sys, glob, shutil, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _ckpt import Checkpoint

BASE = os.path.abspath(os.path.join(HERE, ".."))
STR, RUNS = os.path.join(BASE, "structures"), os.path.join(BASE, "runs")
HA = 27.211386

FUNCTIONAL, BASIS = "wb97x", "def2-svp"      # the level scripts/53 validated at

# cheapest first; atom counts include hydrogens
ORDER = [
    ("HATCN",  30, "reference, 6x CN, planar (risk 92)"),
    ("SBF2CN", 43, "spiro, 2x CN (risk 32)"),
    ("TCPM",   49, "tetrahedral C, 4x CN (risk 35)"),
    ("TCPSi",  49, "tetrahedral Si, 4x CN (risk 35)"),
]
PROXY_SLOPE = 0.055          # eV(E_b) per eV(EA), scripts/53
PROXY_R = 0.962
FRAG_ANCHOR = (-0.462, 0.292)   # (EA, E_b) of benzonitrile, the proxy's anchor


def clear_scratch():
    """psi4's DF scratch from a killed run is not cleaned up, and 22 GB of it once
    filled this session's disk allowance mid-calculation (scripts/52). Clearing is
    safe here because nothing else is running."""
    for pat in ("/tmp/psi.*", "/tmp/psi4_*.npy", "/tmp/dfh.*"):
        for p in glob.glob(pat):
            try:
                shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
            except OSError:
                pass


def run_xtb(tag):
    """GFN2 geometry, matching every other structure in the project."""
    wd = os.path.join(RUNS, tag)
    out = os.path.join(wd, "xtbopt.xyz")
    if os.path.exists(out):
        return out
    os.makedirs(wd, exist_ok=True)
    src = os.path.join(STR, f"{tag}.xyz")
    if not os.path.exists(src):
        return None
    shutil.copy(src, os.path.join(wd, "in.xyz"))
    with open(os.path.join(wd, "xtb.log"), "w") as log:
        subprocess.run(["xtb", "in.xyz", "--gfn", "2", "--opt", "tight"],
                       cwd=wd, stdout=log, stderr=subprocess.STDOUT)
    return out if os.path.exists(out) else None


def single_point(geom_path, charge, mult):
    import psi4
    psi4.core.clean()
    lines = open(geom_path).read().splitlines()
    n = int(lines[0])
    body = "\n".join(lines[2:2 + n])
    mol = psi4.geometry(f"{charge} {mult}\n{body}\n"
                        "symmetry c1\nno_reorient\nno_com\n")
    psi4.set_options({
        "basis": BASIS, "scf_type": "df",
        "reference": "uks" if mult > 1 else "rks",
        "maxiter": 300, "guess": "sad",
        # 1e-5, set from the DF noise floor measured in scripts/52. Tighter cannot
        # be reached on molecules this size and simply hangs the run.
        "d_convergence": 1e-7, "e_convergence": 1e-5,
        "damping_percentage": 20.0,
    })
    return psi4.energy(FUNCTIONAL, molecule=mol)


def main():
    import psi4
    psi4.set_memory(os.environ.get("PSI4_MEM", "5 GB"))
    psi4.set_num_threads(int(os.environ.get("PSI4_THREADS", "4")))
    psi4.core.set_output_file(os.path.join(RUNS, "psi4_cand_ea.out"), False)
    clear_scratch()

    ck = Checkpoint("dft_seedlayer_screen/runs/candidate_EA.json",
                    label="candidate EA")
    ea = {}
    for tag, natoms, note in ORDER:
        geom = run_xtb(tag)
        if geom is None:
            print(f"  {tag:<8} geometry missing -- skipped", flush=True)
            continue
        for suffix, chg, mult in (("neutral", 0, 1), ("anion", -1, 2)):
            key = f"{tag}_{suffix}"
            if ck.has(key):
                print(f"  {tag:<8} {suffix:<8} [cached] {ck.get(key):.6f}", flush=True)
                continue
            print(f"  {tag:<8} {suffix:<8} running ({natoms} atoms)...", flush=True)
            e = single_point(geom, chg, mult)
            ck.put(key, e)
            print(f"  {tag:<8} {suffix:<8} {e:.6f} Ha  [checkpointed]", flush=True)
            clear_scratch()
        e0, em = ck.get(f"{tag}_neutral"), ck.get(f"{tag}_anion")
        if e0 is not None and em is not None:
            ea[tag] = (e0 - em) * HA
            ck.put(f"{tag}_EA_eV", ea[tag])

    if "HATCN" not in ea:
        print("\nHATCN reference not yet done -- rerun to continue. A candidate EA")
        print("means nothing without the reference at the same basis.")
        return
    if len(ea) < 2:
        print("\nonly the reference is done -- rerun to continue")
        return

    print("\n" + "=" * 74)
    print(f"EA at {FUNCTIONAL}/{BASIS} (vertical, gas phase) -- anchor-strength proxy")
    print("=" * 74)
    print(f"{'candidate':<9}{'EA (eV)':>9}{'vs HATCN':>10}{'est dE_b':>10}   note")
    print("-" * 74)
    notes = {t: n for t, _, n in ORDER}
    for t in sorted(ea, key=lambda k: -ea[k]):
        d = ea[t] - ea["HATCN"]
        print(f"{t:<9}{ea[t]:>9.3f}{d:>10.3f}{d*PROXY_SLOPE:>10.3f}   {notes[t]}")

    print(f"\nHOW TO READ THIS (proxy slope {PROXY_SLOPE} eV(E_b)/eV(EA), "
          f"r = {PROXY_R})")
    best = max(ea, key=lambda k: ea[k])
    if best != "HATCN":
        print(f"  {best} is the stronger acceptor of the set, but the est dE_b")
        print("  column is what matters: core electronics move the binding by")
        print("  only tenths of an eV at most.")
    else:
        print("  HATCN remains the strongest acceptor of the set. The sp3 cores buy")
        print("  glass formation at some cost in acceptor strength.")
    print("  NONE of these numbers predicts HATCN-strength binding. HATCN's")
    print("  1.03 eV comes from six nitriles plus the inter-molecular pocket")
    print("  (+0.53 eV, scripts/41), not from core electronics -- and the pocket")
    print("  is exactly what a tetrahedral molecule may fail to form, since its")
    print("  four nitriles point in four different directions. That is the")
    print("  question a monolayer slab has to answer, not this script.")


if __name__ == "__main__":
    main()
