"""Does F4TCNQ stay transparent on Ag? TD-DFT of the neutral and the radical anion.

WHY THIS MATTERS. F4TCNQ is the only material in the screen whose cluster E_b
(0.97 eV) is close enough to HATCN's (1.03 eV) that the monolayer correction could
flip the ranking. But anchoring is not the only axis: the whole point of the HATCN
seed is that it costs nothing optically (scripts/32: +0.00 %p absorption against
neat Ag, vs +2.7 %p for 2 nm Au and +11 %p for 2 nm Al).

F4TCNQ cannot be assumed to be in the same position. Its electron affinity
(~5.2 eV) is well above the Ag work function (~4.3 eV), so electron transfer from
Ag to F4TCNQ is strongly downhill and the interfacial molecules become F4TCNQ
radical anions. That is not a side effect -- it is exactly the mechanism that
makes F4TCNQ a p-dopant. And the anion is coloured: experiment puts bands at
410 nm and ~766 nm (CrystEngComm 18, 8906 (2016); Sci. Rep. 6, 28510 (2016)),
both inside the visible.

WHAT THIS SCRIPT DOES. Computes the vertical excitation spectrum of

    F4TCNQ neutral / F4TCNQ radical anion
    HATCN  neutral / HATCN  radical anion

at the same level, so the two materials can be compared internally. Absolute
positions from a vacuum TD-DFT calculation are not expected to match solid-state
spectra; what is being tested is whether HATCN's anion is as strongly coloured as
F4TCNQ's, since both accept charge from Ag.

Functional: a range-separated hybrid. Anion excitations of an acceptor have
substantial charge-transfer character and a global hybrid (let alone a GGA)
underestimates them badly.
"""
import os, json
import numpy as np
import psi4
from psi4.driver.procrouting.response.scf_response import tdscf_excitations

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STR, RUNS = os.path.join(BASE, "structures"), os.path.join(BASE, "runs")
NM_PER_EV = 1239.841984

CASES = [
    ("F4TCNQ", "F4TCNQ.xyz", 0, 1),
    ("F4TCNQ_anion", "F4TCNQ.xyz", -1, 2),
    ("HATCN", "HATCN.xyz", 0, 1),
    ("HATCN_anion", "HATCN.xyz", -1, 2),
]
FUNCTIONAL = "wb97x"        # range-separated: needed for acceptor CT states
BASIS = "def2-svp"
N_STATES = 12


def geometry(fname, charge, mult):
    lines = open(os.path.join(STR, fname)).read().splitlines()
    n = int(lines[0])
    body = "\n".join(lines[2:2 + n])
    return psi4.geometry(f"{charge} {mult}\n{body}\nsymmetry c1\nno_reorient\nno_com\n")


def run(tag, fname, charge, mult):
    print(f"\n=== {tag}  (charge {charge}, mult {mult}) ===", flush=True)
    psi4.core.clean()
    mol = geometry(fname, charge, mult)
    psi4.set_options({
        "basis": BASIS,
        "scf_type": "df",
        "reference": "uks" if mult > 1 else "rks",
        "maxiter": 300,
        "guess": "sad",
        # TDSCF reuses the SCF's JK object; without this psi4 tears it down after
        # the ground state and tdscf_excitations dies with
        # "JK object is not initialized, please set option SAVE_JK to True".
        "save_jk": True,
    })
    e, wfn = psi4.energy(FUNCTIONAL, molecule=mol, return_wfn=True)
    # tdscf_excitations lives in the procrouting submodule, not at psi4 top level,
    # and it RETURNS the results rather than only stashing them on the wavefunction.
    res = tdscf_excitations(wfn, states=N_STATES, tda=False)

    rows = []
    for i, r in enumerate(res):
        ev = float(r["EXCITATION ENERGY"]) * 27.211386
        f = float(r.get("OSCILLATOR STRENGTH (LEN)",
                        r.get("OSCILLATOR STRENGTH (VEL)", 0.0)))
        rows.append({"state": i + 1, "eV": ev, "nm": NM_PER_EV / ev, "f": f})
    return rows


def report(all_rows):
    print("\n" + "=" * 74)
    print("Vertical excitations with oscillator strength f > 0.01, 300-900 nm")
    print(f"{'system':<16}{'lambda (nm)':>13}{'E (eV)':>9}{'f':>9}   band")
    print("-" * 74)
    for tag, rows in all_rows.items():
        shown = 0
        for r in rows:
            if r["f"] < 0.01 or not (300 <= r["nm"] <= 900):
                continue
            band = ("VISIBLE" if 400 <= r["nm"] <= 700 else
                    "NIR" if r["nm"] > 700 else "UV")
            print(f"{tag:<16}{r['nm']:>13.0f}{r['eV']:>9.2f}{r['f']:>9.3f}   {band}")
            shown += 1
        if shown == 0:
            print(f"{tag:<16}{'-- none in range --':>31}")

    print("\n" + "=" * 74)
    print("Visible-range oscillator strength (400-700 nm), summed:")
    tot = {}
    for tag, rows in all_rows.items():
        s = sum(r["f"] for r in rows if 400 <= r["nm"] <= 700)
        tot[tag] = s
        print(f"  {tag:<16}{s:>8.3f}")
    print("\nThe number that decides the seed question is the ANION comparison:")
    a, b = tot.get("HATCN_anion", 0), tot.get("F4TCNQ_anion", 0)
    print(f"  HATCN(-)  {a:.3f}   vs   F4TCNQ(-)  {b:.3f}")
    if b > 2 * max(a, 1e-6):
        print("  -> F4TCNQ's anion is far more strongly coloured in the visible.")
        print("     Since Ag necessarily reduces F4TCNQ at the interface, this is an")
        print("     intrinsic optical cost that HATCN does not pay, and it stands")
        print("     whatever the anchoring comparison turns out to be.")
    elif a > 2 * max(b, 1e-6):
        print("  -> HATCN's anion is the more coloured of the two, which would")
        print("     undercut the transparency argument for HATCN. Check carefully.")
    else:
        print("  -> Comparable. Transparency does not separate them at this level;")
        print("     decide on anchoring and on the measured stack absorption.")
    print("\nCaveat: vacuum TD-DFT vertical excitations. Absolute positions will not")
    print("match solid-state spectra (experiment puts the F4TCNQ anion at 410 and")
    print("766 nm); the internal HATCN-vs-F4TCNQ comparison is the useful output.")


if __name__ == "__main__":
    os.makedirs(RUNS, exist_ok=True)
    psi4.set_memory("3 GB")
    psi4.set_num_threads(int(os.environ.get("PSI4_THREADS", "2")))
    psi4.core.set_output_file(os.path.join(RUNS, "psi4_tddft.out"), False)

    all_rows = {}
    for tag, fname, chg, mult in CASES:
        try:
            all_rows[tag] = run(tag, fname, chg, mult)
            print(f"  {len(all_rows[tag])} states", flush=True)
        except Exception as ex:
            print(f"  FAILED: {type(ex).__name__}: {ex}", flush=True)
            all_rows[tag] = []
    json.dump(all_rows, open(os.path.join(RUNS, "tddft_absorption.json"), "w"), indent=2)
    report(all_rows)
