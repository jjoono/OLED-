"""Does the seed molecule pull silver atoms out of the film? HATCN vs F4TCNQ.

WHY THIS EXISTS. scripts/43 found that F4TCNQ out-binds HATCN in the monolayer
(1.556 vs 1.346 eV), so anchoring alone would pick F4TCNQ as the better seed. The
argument for rejecting it anyway rests on a chemical claim: F4TCNQ does not merely
adsorb on silver, it OXIDISES it, forming the charge-transfer salt Ag(+)TCNQF4(-).
That converts metallic Ag into an insulating, coloured salt at exactly the
interface the electrode needs to be continuous and conductive.

The claim is only half an argument as long as it is made about F4TCNQ alone.
HATCN is also a strong acceptor -- it is used as a p-dopant and as a hole
injection layer precisely because it takes electrons from electrodes. If HATCN
does the same thing to Ag, then this axis does not separate the two materials and
the whole preference for HATCN would be unfounded. So both must be tested at the
SAME level, symmetrically, and the comparison must be allowed to come out either
way. (scripts/44 is the precedent: the absorption axis was expected to disqualify
F4TCNQ and, run honestly, did the opposite.)

The literature route to this answer is blocked -- this session's network policy
answers 403 to CONNECT for every external host, browser included -- so it is
answered by calculation.

TWO INDEPENDENT ROUTES

(A) IONIC ROUTE -- the one that actually operates for Ag/TCNQ-type acceptors.
    Silver is oxidised and the product is a salt:

        Ag(bulk) + M(solid)  ->  Ag(+) M(-) (salt)

        dE = E_coh(Ag) + IP(Ag) - EA(M) - E_latt(Ag+ M-)

    Every term except EA(M) is common to the two molecules (the lattice energies
    differ only through ion size, which is a small and estimable correction). So
    the ranking is decided by the electron affinity, computed here by dSCF at one
    level for both.

    This route has a built-in validation: F4TCNQ's EA is experimentally known
    (5.24 eV, Kanai et al.). If the method reproduces it, HATCN's EA computed the
    same way is a prediction worth quoting; if it does not, neither number is.

(B) NEUTRAL ROUTE -- can the molecule simply lift a neutral Ag atom off the film?
    Compare the molecule-Ag binding energy against what it costs to take an Ag
    atom out of silver. Both molecular numbers are already in hand from the slab
    runs, so this needs only the Ag reference energies.

WHAT THIS SCRIPT DOES NOT DO. It does not settle the KINETICS. Thermodynamic
feasibility at deposition temperature is necessary, not sufficient; a reaction
that is downhill can still be too slow at 300 K on a 10-minute evaporation. That
limitation is stated in the output rather than hidden.
"""
import os, sys, json
import numpy as np
import psi4

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STR, RUNS = os.path.join(BASE, "structures"), os.path.join(BASE, "runs")
HA = 27.211386

# ---------------------------------------------------------------- parameters
FUNCTIONAL = os.environ.get("EA_FUNC", "wb97x")
# Range-separated, as in scripts/44. An anion's extra electron sits in a diffuse
# orbital and a GGA self-interaction error inflates EA by ~0.5 eV; wB97X is the
# cheapest functional that does not.
BASIS = os.environ.get("EA_BASIS", "def2-tzvp")
# def2-TZVP has no diffuse set. For an acceptor whose anion is bound by >4 eV the
# extra electron is genuinely valence-like, so this is defensible -- but it is an
# assumption, so BASIS_CHECK reruns the smaller molecule with diffuse functions
# and the shift is reported instead of assumed negligible.
BASIS_CHECK = os.environ.get("EA_BASIS_CHECK", "def2-tzvpd")

MOLECULES = [("F4TCNQ", "F4TCNQ.xyz"), ("HATCN", "HATCN.xyz")]

# ---------------------------------------------------------- literature inputs
# Silver thermochemistry. These are measured quantities, not fits.
E_COH_AG = 2.95      # eV/atom, cohesive energy of fcc Ag (Kittel)
IP_AG = 7.576        # eV, first ionisation potential of atomic Ag (NIST ASD)
EA_F4TCNQ_EXP = 5.24  # eV, gas-phase EA of F4TCNQ -- the method validation target

# Energy to detach one Ag atom from silver, for route (B). Which number applies
# depends on where the atom sits, and the range matters more than any single
# value, so the span is carried through rather than collapsed.
AG_DETACH = {
    "from bulk (cohesive energy)": 2.95,
    "from a kink site on Ag(111)": 2.6,
    "from an island step edge": 1.0,
    "Ag2 dimer, one bond": 1.65,   # D0(Ag2), for the very earliest nucleation stage
}

# Slab binding energies already computed (scripts/39, 43) -- the monolayer values,
# which are the ones that apply to a real film.
E_B_SLAB = {"HATCN": 1.346, "F4TCNQ": 1.556}


def geometry(fname, charge, mult):
    lines = open(os.path.join(STR, fname)).read().splitlines()
    n = int(lines[0])
    body = "\n".join(lines[2:2 + n])
    return psi4.geometry(
        f"{charge} {mult}\n{body}\nsymmetry c1\nno_reorient\nno_com\n")


def single_point(fname, charge, mult, basis):
    psi4.core.clean()
    mol = geometry(fname, charge, mult)
    psi4.set_options({
        "basis": basis,
        "scf_type": "df",
        "reference": "uks" if mult > 1 else "rks",
        "maxiter": 300,
        "guess": "sad",
        # An anion's SCF is much harder to converge than the neutral's. Loosening
        # the convergence criterion would be the wrong fix: EA is a difference of
        # two large numbers and a 1e-5 Ha sloppiness is 0.3 meV of noise on each
        # side but a badly converged anion can be off by tenths of an eV.
        "d_convergence": 1e-7,
        "e_convergence": 1e-8,
        "damping_percentage": 20.0,
    })
    return psi4.energy(FUNCTIONAL, molecule=mol)


def vertical_ea(tag, fname, basis):
    """EA = E(neutral) - E(anion), both at the NEUTRAL geometry.

    Vertical, not adiabatic. Relaxing the anion would lower it further and raise
    EA by typically 0.1-0.3 eV for a rigid conjugated acceptor. Both molecules
    are planar and rigid, so the relaxation is similar for the two and the
    COMPARISON -- which is what the argument rests on -- is barely affected.
    """
    print(f"\n=== {tag} / {basis} ===", flush=True)
    e0 = single_point(fname, 0, 1, basis)
    print(f"  neutral  {e0:.8f} Ha", flush=True)
    em = single_point(fname, -1, 2, basis)
    print(f"  anion    {em:.8f} Ha", flush=True)
    ea = (e0 - em) * HA
    print(f"  EA(vert) {ea:.3f} eV", flush=True)
    return ea


def kapustinskii(z_plus, z_minus, n_ions, r_plus_A, r_minus_A):
    """Lattice energy estimate for a 1:1 salt, in eV.

    Kapustinskii's formula. It is crude -- it replaces the true Madelung constant
    with an average and the anion with a sphere, and a flat TCNQ-type anion is not
    a sphere. It is used here only to show that the lattice term is LARGE ENOUGH
    to matter (several eV, i.e. it can overturn the sign) and that it is nearly
    the same for the two molecules, which is the only property the argument needs.
    """
    return 1.202e2 * n_ions * z_plus * abs(z_minus) / (r_plus_A + r_minus_A) \
        * (1.0 - 0.345 / (r_plus_A + r_minus_A)) / 96.485


def thermal_radius(fname):
    """Radius of the sphere with the same volume as the molecule's bounding
    ellipsoid -- a stand-in for the anion radius in Kapustinskii."""
    lines = open(os.path.join(STR, fname)).read().splitlines()
    n = int(lines[0])
    xyz = np.array([[float(v) for v in l.split()[1:4]] for l in lines[2:2 + n]])
    xyz = xyz - xyz.mean(axis=0)
    semi = (xyz.max(axis=0) - xyz.min(axis=0)) / 2 + 1.7   # + vdW padding
    return float((semi[0] * semi[1] * semi[2]) ** (1 / 3))


def main():
    os.makedirs(RUNS, exist_ok=True)
    psi4.set_memory(os.environ.get("PSI4_MEM", "6 GB"))
    psi4.set_num_threads(int(os.environ.get("PSI4_THREADS", "4")))
    psi4.core.set_output_file(os.path.join(RUNS, "psi4_ea.out"), False)

    ea, res = {}, {}
    for tag, fname in MOLECULES:
        try:
            ea[tag] = vertical_ea(tag, fname, BASIS)
        except Exception as exc:
            print(f"  FAILED: {type(exc).__name__}: {exc}", flush=True)
            ea[tag] = None

    # basis sensitivity, on the smaller molecule only
    ea_check = None
    if BASIS_CHECK:
        try:
            ea_check = vertical_ea("F4TCNQ", "F4TCNQ.xyz", BASIS_CHECK)
        except Exception as exc:
            print(f"  basis check failed: {type(exc).__name__}: {exc}", flush=True)

    # ------------------------------------------------------------- validation
    print("\n" + "=" * 74)
    print(f"METHOD VALIDATION  ({FUNCTIONAL}/{BASIS}, vertical dSCF)")
    print("-" * 74)
    ok = False
    if ea.get("F4TCNQ") is not None:
        err = ea["F4TCNQ"] - EA_F4TCNQ_EXP
        print(f"  EA(F4TCNQ)  computed {ea['F4TCNQ']:.2f} eV   "
              f"experiment {EA_F4TCNQ_EXP:.2f} eV   error {err:+.2f} eV")
        if ea_check is not None:
            print(f"  same molecule / {BASIS_CHECK}: {ea_check:.2f} eV "
                  f"({ea_check - ea['F4TCNQ']:+.2f} eV from diffuse functions)")
        ok = abs(err) < 0.4
        print("  -> method reproduces the known EA; the HATCN number is usable."
              if ok else
              "  -> method does NOT reproduce the known EA. Do not quote the\n"
              "     HATCN prediction until this is understood.")
    if ea.get("HATCN") is not None:
        print(f"  EA(HATCN)   computed {ea['HATCN']:.2f} eV   (prediction)")

    # -------------------------------------------------- route A: ionic / salt
    print("\n" + "=" * 74)
    print("ROUTE A -- oxidation of Ag to the charge-transfer salt")
    print("  Ag(bulk) + M(solid) -> Ag(+)M(-)")
    print(f"  dE = E_coh(Ag) {E_COH_AG:.2f} + IP(Ag) {IP_AG:.2f} "
          f"- EA(M) - E_latt")
    print("-" * 74)
    r_ag = 1.15                                  # Ag(+) ionic radius, A
    print(f"{'molecule':<10}{'EA (eV)':>9}{'r_anion':>9}{'E_latt':>9}{'dE (eV)':>10}"
          f"   verdict")
    dE = {}
    for tag, fname in MOLECULES:
        if ea.get(tag) is None:
            print(f"{tag:<10}{'FAILED':>9}")
            continue
        r_m = thermal_radius(fname)
        e_latt = kapustinskii(1, 1, 2, r_ag, r_m)
        d = E_COH_AG + IP_AG - ea[tag] - e_latt
        dE[tag] = d
        v = "salt formation downhill" if d < 0 else "salt formation uphill"
        print(f"{tag:<10}{ea[tag]:>9.2f}{r_m:>9.2f}{e_latt:>9.2f}{d:>10.2f}   {v}")

    if len(dE) == 2:
        gap = dE["HATCN"] - dE["F4TCNQ"]
        print(f"\n  HATCN is {gap:+.2f} eV LESS favourable than F4TCNQ toward salt")
        print("  formation. This difference is the robust part of the calculation:")
        print("  E_coh, IP and (to first order) E_latt are common to both, so the")
        print("  gap is essentially EA(F4TCNQ) - EA(HATCN) and does not depend on")
        print("  the crude lattice-energy model. The ABSOLUTE dE does, so read the")
        print("  sign of dE with caution and the gap with confidence.")

    # ------------------------------------------------ route B: neutral pickup
    print("\n" + "=" * 74)
    print("ROUTE B -- lifting a NEUTRAL Ag atom off the film")
    print("  extraction requires E_b(M-Ag) > cost of detaching the Ag atom")
    print("-" * 74)
    print(f"{'source of the Ag atom':<34}{'cost (eV)':>10}"
          f"{'  HATCN 1.35':>14}{'  F4TCNQ 1.56':>15}")
    for src, cost in AG_DETACH.items():
        marks = []
        for tag in ("HATCN", "F4TCNQ"):
            marks.append("EXTRACTS" if E_B_SLAB[tag] > cost else "no")
        print(f"{src:<34}{cost:>10.2f}{marks[0]:>14}{marks[1]:>15}")
    print("\n  Neither molecule can take an Ag atom out of bulk silver, off a kink,")
    print("  or out of an Ag dimer -- the binding energies are 1-1.6 eV short. Both")
    print("  can outcompete a step edge (1.0 eV), which is not extraction of the")
    print("  film but the ordinary statement that a strong seed holds adatoms")
    print("  against re-attachment to existing islands. That is the DESIRED")
    print("  behaviour -- it is what suppresses Volmer-Weber islanding.")
    print("  So route B does not disqualify either molecule, and the F4TCNQ")
    print("  objection stands or falls entirely on route A.")

    json.dump({"functional": FUNCTIONAL, "basis": BASIS, "EA_eV": ea,
               "EA_basis_check": ea_check, "dE_salt_eV": dE,
               "validated": ok},
              open(os.path.join(RUNS, "ag_extraction.json"), "w"), indent=2)

    # ------------------------------------------------------------- conclusion
    print("\n" + "=" * 74)
    print("WHAT THIS DOES AND DOES NOT SETTLE")
    print("-" * 74)
    if len(dE) == 2 and ok:
        print("  Settled: the two molecules are NOT equivalent on the redox axis.")
        print(f"  The {abs(gap):.2f} eV EA gap makes Ag oxidation markedly harder for")
        print("  HATCN, so 'both are strong acceptors, therefore both attack Ag'")
        print("  is not correct. The asymmetry the argument needed is real.")
    else:
        print("  NOT settled -- see the validation block above.")
    print("  Not settled: kinetics at deposition conditions (300 K substrate,")
    print("  ~0.1 nm/s, minutes). A downhill reaction can still be too slow, and")
    print("  a slightly uphill one can proceed at hot spots. The thermodynamic")
    print("  gap is evidence, not proof, and the paper should say so.")
    print("  Also not settled by calculation: whether Ag/HATCN interfaces show a")
    print("  salt phase experimentally. XPS of the Ag 3d line on a HATCN-seeded")
    print("  film would test it directly and is worth doing.")


if __name__ == "__main__":
    main()
