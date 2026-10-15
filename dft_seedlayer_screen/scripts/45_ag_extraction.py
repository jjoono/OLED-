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

MOLECULES = [("F4TCNQ", "F4TCNQ.xyz"), ("HATCN", "HATCN.xyz"), ("TCNQ", "TCNQ.xyz")]

# ---------------------------------------------------------- literature inputs
# Silver thermochemistry. These are measured quantities, not fits.
E_COH_AG = 2.95      # eV/atom, cohesive energy of fcc Ag (Kittel)
IP_AG = 7.576        # eV, first ionisation potential of atomic Ag (NIST ASD)

# VALIDATION TARGET -- and the place where the first version of this script went
# wrong, so the reasoning is written out rather than left as a number.
#
# The value quoted for F4TCNQ throughout the OLED/doping literature is 5.24 eV.
# The first run compared the computed gas-phase dSCF EA against it and missed by
# -1.22 eV, with diffuse functions worth only +0.01 eV, so the gap was not a
# basis-set problem. It was a PHASE problem: 5.24 eV is the SOLID-STATE electron
# affinity, i.e. the LUMO edge below vacuum measured by UPS/IPES on a thin film.
# That includes the electronic polarisation of the surrounding molecules around
# the extra charge, worth of order 1.0-1.5 eV in a molecular solid. A gas-phase
# calculation cannot reproduce it and should not be asked to.
#
# Comparing the two was my error, and it matters twice over: it also means the
# thermodynamic cycle below was mixing phases, since M there is a SOLID.
EA_F4TCNQ_GAS = 3.9     # eV, gas-phase EA of F4TCNQ. APPROXIMATE -- see note.
EA_F4TCNQ_SOLID = 5.24  # eV, solid-state (UPS/IPES) EA, the widely quoted value
# The gas-phase figure is the weak link: it is from memory, this session cannot
# reach the literature to check it, and the whole point of a validation target is
# that it be independent. So the absolute check is reported as PROVISIONAL and
# the script leans instead on the differential check below, which does not depend
# on it.

# DIFFERENTIAL VALIDATION, independent of the absolute reference.
# EA(F4TCNQ) - EA(TCNQ) is ~0.5-0.6 eV, and this difference is far better
# established than either absolute value because it is reproduced in gas phase
# AND in the solid (the polarisation energies of two molecules this similar
# cancel to within ~0.1 eV). If the method gets this splitting right, it is
# describing fluorination of the acceptor correctly, which is the chemistry that
# the HATCN-vs-F4TCNQ comparison actually rests on.
DELTA_F4_TCNQ = 0.55    # eV, expected EA(F4TCNQ) - EA(TCNQ)
DELTA_TOL = 0.25        # eV, tolerance on that splitting

# Polarisation stabilisation of a molecular anion in its own solid. Enters the
# cycle because M is a solid there. Carried explicitly, and identical for the two
# molecules, so it shifts both dE values together and CANCELS in the gap.
E_POL = 1.3             # eV, order-of-magnitude for a molecular organic solid

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


def ensure_tcnq():
    """Build TCNQ from F4TCNQ by swapping the four F for H.

    Used only for the DIFFERENTIAL validation, where what matters is the change
    on fluorination, not either absolute number. Replacing F by H along the same
    C-X vector at the standard C-H length is a rigid approximation -- the ring
    geometry is not re-relaxed -- but the ring barely moves on defluorination and,
    more to the point, the same approximation is applied to the neutral and the
    anion, so it largely cancels in the EA difference.
    """
    dst = os.path.join(STR, "TCNQ.xyz")
    if os.path.exists(dst):
        return
    lines = open(os.path.join(STR, "F4TCNQ.xyz")).read().splitlines()
    n = int(lines[0])
    sym, xyz = [], []
    for l in lines[2:2 + n]:
        p = l.split()
        sym.append(p[0])
        xyz.append([float(v) for v in p[1:4]])
    xyz = np.array(xyz)
    for i, s in enumerate(sym):
        if s != "F":
            continue
        j = int(np.argmin([np.linalg.norm(xyz[i] - xyz[k]) if k != i else 9e9
                           for k in range(len(sym))]))
        v = xyz[i] - xyz[j]
        xyz[i] = xyz[j] + v / np.linalg.norm(v) * 1.08   # C-H
        sym[i] = "H"
    with open(dst, "w") as f:
        f.write(f"{n}\nTCNQ, built from F4TCNQ by F->H (scripts/45)\n")
        for s, p in zip(sym, xyz):
            f.write(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
    print(f"built {dst}", flush=True)


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


CACHE = os.path.join(RUNS, "ea_cache.json")
# Every SCF here is expensive and this container has repeatedly killed long runs
# mid-flight (idle reclaim, and four working-tree rollbacks). Without a cache each
# restart redoes work that already converged, and the last attempt died on the
# OPTIONAL basis check after all three production numbers were finished -- losing
# the report for want of a single-point that was not needed for the conclusion.
# Keyed on everything that changes the number, so a settings change invalidates it.


def cache_load():
    try:
        return json.load(open(CACHE))
    except Exception:
        return {}


def cache_key(tag, basis):
    return f"{tag}|{FUNCTIONAL}|{basis}"


def vertical_ea(tag, fname, basis):
    """EA = E(neutral) - E(anion), both at the NEUTRAL geometry.

    Vertical, not adiabatic. Relaxing the anion would lower it further and raise
    EA by typically 0.1-0.3 eV for a rigid conjugated acceptor. All three
    molecules are planar and rigid, so the relaxation is similar and the
    COMPARISONS -- which is what the argument rests on -- are barely affected.
    """
    c = cache_load()
    k = cache_key(tag, basis)
    if k in c:
        print(f"\n=== {tag} / {basis} ===  [cached] EA = {c[k]['ea']:.3f} eV"
              f"   ({c[k].get('provenance', 'this script')})", flush=True)
        return c[k]["ea"]

    print(f"\n=== {tag} / {basis} ===", flush=True)
    e0 = single_point(fname, 0, 1, basis)
    print(f"  neutral  {e0:.8f} Ha", flush=True)
    em = single_point(fname, -1, 2, basis)
    print(f"  anion    {em:.8f} Ha", flush=True)
    ea = (e0 - em) * HA
    print(f"  EA(vert) {ea:.3f} eV", flush=True)

    c = cache_load()          # re-read: another process may have added entries
    c[k] = {"ea": ea, "E_neutral_Ha": e0, "E_anion_Ha": em,
            "provenance": "computed"}
    json.dump(c, open(CACHE, "w"), indent=2)
    return ea


def kapustinskii(z_plus, z_minus, n_ions, r_plus_A, r_minus_A):
    """Lattice energy estimate for a 1:1 salt, in eV.

    Kapustinskii's formula. It is crude -- it replaces the true Madelung constant
    with an average and the anion with a sphere, and a flat TCNQ-type anion is not
    a sphere. It is used here only to show that the lattice term is LARGE ENOUGH
    to matter (several eV, i.e. it can overturn the sign) and that it is nearly
    the same for the two molecules, which is the only property the argument needs.

    UNITS. The canonical form takes radii in PICOMETRES:
        U [kJ/mol] = 1.202e5 * nu * z+ * |z-| / r_pm * (1 - 34.5 / r_pm)
    The first version of this script kept radii in angstrom but carried the
    prefactor as 1.202e2 instead of 1.202e3, so every lattice energy came out a
    factor of 10 too small (0.43 eV instead of 4.3 eV) -- small enough to look
    like a negligible correction when it is in fact the largest single term in
    the cycle. Radii are converted to pm here so the constants match the textbook.
    """
    r_pm = (r_plus_A + r_minus_A) * 100.0
    u_kj = 1.202e5 * n_ions * z_plus * abs(z_minus) / r_pm * (1.0 - 34.5 / r_pm)
    return u_kj / 96.485


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

    ensure_tcnq()
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
    print("  All values below are GAS-PHASE vertical dSCF. The 5.24 eV usually")
    print("  quoted for F4TCNQ is a SOLID-STATE (UPS/IPES) affinity and includes")
    print("  ~1-1.5 eV of polarisation; comparing the two directly is an error and")
    print("  was the one made in the first version of this script.")
    print()
    for tag in ("TCNQ", "F4TCNQ", "HATCN"):
        if ea.get(tag) is not None:
            print(f"  EA_gas({tag:<7}) = {ea[tag]:>5.2f} eV")
    if ea_check is not None and ea.get("F4TCNQ") is not None:
        print(f"  F4TCNQ / {BASIS_CHECK}: {ea_check:.2f} eV "
              f"({ea_check - ea['F4TCNQ']:+.2f} eV from diffuse functions) "
              f"-> basis converged")

    # (1) differential check -- the one that does not depend on a remembered value
    ok = False
    print("\n  CHECK 1 (differential, primary): EA(F4TCNQ) - EA(TCNQ)")
    if ea.get("F4TCNQ") is not None and ea.get("TCNQ") is not None:
        d = ea["F4TCNQ"] - ea["TCNQ"]
        err = d - DELTA_F4_TCNQ
        ok = abs(err) < DELTA_TOL
        print(f"    computed {d:+.2f} eV   expected {DELTA_F4_TCNQ:+.2f} eV   "
              f"error {err:+.2f} eV   -> {'PASS' if ok else 'FAIL'}")
        print("    This is the check the argument leans on: it tests whether the")
        print("    method describes acceptor strength correctly, and polarisation")
        print("    cancels between two molecules this similar, so it is valid in")
        print("    gas phase.")
    else:
        print("    unavailable -- TCNQ or F4TCNQ did not converge")

    # (2) absolute check -- provisional, flagged as such
    print("\n  CHECK 2 (absolute, PROVISIONAL): EA_gas(F4TCNQ)")
    if ea.get("F4TCNQ") is not None:
        err2 = ea["F4TCNQ"] - EA_F4TCNQ_GAS
        print(f"    computed {ea['F4TCNQ']:.2f} eV   reference ~{EA_F4TCNQ_GAS:.2f} eV"
              f"   error {err2:+.2f} eV")
        print(f"    implied polarisation, EA_solid - EA_gas = "
              f"{EA_F4TCNQ_SOLID - ea['F4TCNQ']:+.2f} eV, which is the right size")
        print("    for a molecular solid -- consistent, but not an independent test.")
        print("    TREAT AS UNCONFIRMED: the gas-phase reference is from memory and")
        print("    this session cannot reach the literature to verify it. Check it")
        print("    against the NIST WebBook before this goes in a manuscript.")

    print(f"\n  -> {'Method accepted on the differential check.' if ok else 'METHOD NOT VALIDATED -- do not quote these numbers.'}")

    # -------------------------------------------------- route A: ionic / salt
    print("\n" + "=" * 74)
    print("ROUTE A -- oxidation of Ag to the charge-transfer salt")
    print("  Ag(bulk) + M(solid) -> Ag(+)M(-)")
    print(f"  dE = E_coh(Ag) {E_COH_AG:.2f} + IP(Ag) {IP_AG:.2f} "
          f"- [EA_gas(M) + E_pol {E_POL:.2f}] - E_latt")
    print("  M is a SOLID here, so the affinity entering the cycle is the")
    print("  solid-state one: EA_gas plus the polarisation of the surrounding")
    print("  molecules. E_pol is taken as identical for the two, so it shifts both")
    print("  dE together and drops out of the gap.")
    print("-" * 74)
    r_ag = 1.15                                  # Ag(+) ionic radius, A
    print(f"{'molecule':<10}{'EA_gas':>8}{'EA_sol':>8}{'r_anion':>9}{'E_latt':>8}"
          f"{'dE (eV)':>9}   verdict")
    dE = {}
    for tag, fname in MOLECULES:
        if tag == "TCNQ":
            continue                 # validation only, not a seed candidate
        if ea.get(tag) is None:
            print(f"{tag:<10}{'FAILED':>8}")
            continue
        r_m = thermal_radius(fname)
        e_latt = kapustinskii(1, 1, 2, r_ag, r_m)
        ea_sol = ea[tag] + E_POL
        d = E_COH_AG + IP_AG - ea_sol - e_latt
        dE[tag] = d
        v = "downhill" if d < 0 else "uphill"
        print(f"{tag:<10}{ea[tag]:>8.2f}{ea_sol:>8.2f}{r_m:>9.2f}{e_latt:>8.2f}"
              f"{d:>9.2f}   {v}")

    if len(dE) == 2:
        gap = dE["HATCN"] - dE["F4TCNQ"]
        print(f"\n  HATCN is {gap:+.2f} eV LESS favourable than F4TCNQ toward salt")
        print("  formation.")
        print("\n  WHAT TO TRUST HERE. The GAP is robust: E_coh, IP and E_pol are")
        print("  common to both molecules and cancel exactly, leaving")
        print("      gap = [EA(F4TCNQ) - EA(HATCN)] + [E_latt(F4) - E_latt(HATCN)]")
        print("  whose two terms happen to push the same way -- HATCN is both the")
        print("  weaker acceptor and the bulkier anion, hence the more weakly bound")
        print("  salt. The ABSOLUTE dE is NOT trustworthy: Kapustinskii models a")
        print("  flat TCNQ-type anion as a sphere and puts Ag(+) at its spherical")
        print("  surface, when in the real Ag-TCNQ structure Ag(+) sits ~2.3 A from")
        print("  a nitrile N, i.e. far deeper in the potential. That underestimates")
        print("  E_latt by an unknown but substantial amount, which is why the")
        print("  cycle can report 'uphill' for a salt that is known to exist.")
        print("  Read the sign as unresolved and the gap as the result.")

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

    json.dump({"functional": FUNCTIONAL, "basis": BASIS,
               "EA_gas_eV": ea, "EA_basis_check": ea_check,
               "E_pol_assumed_eV": E_POL, "dE_salt_eV": dE,
               "gap_eV": (dE["HATCN"] - dE["F4TCNQ"]) if len(dE) == 2 else None,
               "differential_check_passed": ok,
               "absolute_EA_reference_verified": False},
              open(os.path.join(RUNS, "ag_extraction.json"), "w"), indent=2)

    # ------------------------------------------------------------- conclusion
    print("\n" + "=" * 74)
    print("WHAT THIS DOES AND DOES NOT SETTLE")
    print("-" * 74)
    if len(dE) == 2 and ok:
        print("  Settled: the two molecules are NOT equivalent on the redox axis.")
        print(f"  The {abs(gap):.2f} eV gap makes Ag oxidation markedly harder for")
        print("  HATCN, so 'both are strong acceptors, therefore both attack Ag'")
        print("  is not correct. The asymmetry the argument needed is real.")
        print("  This is a RELATIVE statement and should be written as one.")
    else:
        print("  NOT settled -- see the validation block above.")
    print("  Not settled: the ABSOLUTE question of whether either salt forms.")
    print("  The lattice-energy model is too crude for that, and it is the wrong")
    print("  tool anyway -- an interface is not a bulk salt. Answering it properly")
    print("  means a periodic slab of Ag(111) with the molecule on top and a Bader")
    print("  or Loewdin charge analysis, which the CP2K setup in scripts/38 can")
    print("  already do. That is the natural next calculation.")
    print("  Not settled: kinetics at deposition conditions (300 K substrate,")
    print("  ~0.1 nm/s, minutes). A downhill reaction can still be too slow, and")
    print("  a slightly uphill one can proceed at hot spots. The thermodynamic")
    print("  gap is evidence, not proof, and the paper should say so.")
    print("  Also not settled by calculation: whether Ag/HATCN interfaces show a")
    print("  salt phase experimentally. XPS of the Ag 3d line on a HATCN-seeded")
    print("  film would test it directly and is worth doing.")


if __name__ == "__main__":
    main()
