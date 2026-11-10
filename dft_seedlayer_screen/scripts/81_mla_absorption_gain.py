"""Why an MLA makes electrode absorption hurt MORE, and how much EQE it buys back.

THE ARGUMENT THE EXPERIMENT IS MEANT TO SHOW. In a top-emitting OLED the
microlens array does not create light; it gives already-trapped light more
chances to escape. Every failed chance sends that light back DOWN onto the
device, where it must reflect off the electrode stack to try again. So the
electrode's absorption is paid once on the way out and TWICE more on every
recycling round trip:

    escape with recycling  =  q / (1 - (1-q) * S)

with q the escape probability per encounter with the MLA and S the survival
of one round trip back to the device and out again. S carries the electrode
loss, so the MLA turns a linear absorption penalty into a compounding one --
and the fewer photons a single MLA encounter releases (small q), the more
round trips are needed and the more violently S matters.

TWO PLACES THE ELECTRODE ABSORPTION ENTERS, which is why the effect is
stronger than intuition suggests:
  1. first pass   -- the light leaving the cavity is attenuated by (1 - A)
  2. every recycle -- S = (1 - A)^2 * R_mirror * (1 - A_org), the factor
     (1-A)^2 because the returning photon crosses the electrode inward and
     outward again.

WHAT THE MODEL IS AND IS NOT. This is a photon-budget model, not a full
cavity+MLA simulation: the cavity's first-pass outcoupling is taken as a
fixed efficiency and only its attenuation by the electrode is varied, so the
ABSOLUTE EQE here is illustrative. The RATIO between electrode conditions --
which is what the experiment measures and what the paper claims -- follows
from the photon bookkeeping and is robust.

Numbers used: A_electrode from this project's own measurements and targets
(scripts/72, 78); bottom mirror 97 % (thick Ag); organic/parasitic 2 % per
round trip; escape-cone fraction 1/(2 n^2) with n = 1.5 encapsulation.
"""
import numpy as np

R_MIRROR = 0.97          # thick Ag bottom mirror, per round trip
A_ORG = 0.02             # organics + parasitic, per round trip
N_ENC = 1.5              # encapsulation index
Q_LIST = (0.15, 0.30, 0.50)    # MLA escape probability per encounter
ETA_CAV = 0.45           # first-pass cavity outcoupling INTO the encapsulation

CASES = [("measured now", 0.115),
         ("voids closed", 0.062),
         ("+ specular p=0.5", 0.042),
         ("ideal p=0.9", 0.025),
         ("lossless (limit)", 0.0)]


def escape_direct():
    """Fraction of isotropic flux inside the encapsulation within the cone."""
    return 1.0 / (2.0 * N_ENC ** 2)


def eqe(a_el, q, use_mla):
    """Relative EQE (IQE = 1). First pass out of the cavity, then escape."""
    first = ETA_CAV * (1.0 - a_el)          # reaches the encapsulation
    f_cone = escape_direct()
    S = (1.0 - a_el) ** 2 * R_MIRROR * (1.0 - A_ORG)
    if not use_mla:
        # only the escape cone gets out; trapped light dies (weak recycling)
        return first * f_cone
    # MLA: every encounter has probability q, failures recycle through S
    return first * q / (1.0 - (1.0 - q) * S)


def main():
    print("=" * 74)
    print("Relative EQE vs electrode absorption   (IQE = 1, arbitrary scale)")
    print("=" * 74)
    base_no = eqe(CASES[0][1], Q_LIST[0], False)
    for q in Q_LIST:
        base_ml = eqe(CASES[0][1], q, True)
        print(f"\nMLA escape probability q = {q:.2f}")
        print(f"{'electrode':<20}{'A':>7}{'noMLA':>9}{'MLA':>8}"
              f"{'MLA gain':>10}{'vs now(MLA)':>13}")
        for tag, a in CASES:
            e0, e1 = eqe(a, q, False), eqe(a, q, True)
            print(f"{tag:<20}{100*a:>6.1f}%{e0/base_no:>9.3f}{e1/base_ml:>8.3f}"
                  f"{e1/e0:>10.2f}x{100*(e1/base_ml-1):>12.1f}%")
    print("""
'noMLA' and 'MLA' are normalised to the CURRENT electrode, so each column
reads as the relative EQE improvement the experiment would measure.""")

    print("\n" + "=" * 74)
    print("AMPLIFICATION: how much more the MLA device gains from the same "
          "absorption cut")
    print("=" * 74)
    print(f"{'absorption cut':<26}{'noMLA gain':>12}{'MLA gain (q)':>28}")
    for tag, a in CASES[1:]:
        g0 = eqe(a, Q_LIST[0], False) / eqe(CASES[0][1], Q_LIST[0], False)
        gs = [eqe(a, q, True) / eqe(CASES[0][1], q, True) for q in Q_LIST]
        print(f"{CASES[0][0]+' -> '+tag:<26}{100*(g0-1):>11.1f}%"
              + "".join(f"{100*(g-1):>8.1f}%" for g in gs)
              + f"   q={Q_LIST}")


if __name__ == "__main__":
    main()

# RESULT (run 2026-08-16):
#
# THE HEADLINE. The same absorption cut is worth several times more in an
# MLA device than in a bare one:
#
#   A 11.5 % -> 4.2 %      no MLA:  +8.2 %
#                          MLA q=0.50: +20.5 %   (2.5x the bare gain)
#                          MLA q=0.30: +33.1 %   (4.0x)
#                          MLA q=0.15: +53.8 %   (6.6x)
#
# and the MLA's own gain factor climbs with the better electrode:
#   q=0.15:  1.84x (now) -> 2.61x (4.2 %) -> 3.52x (lossless)
#
# WHY. Without an MLA only the escape cone gets out, so electrode absorption
# is paid once and the EQE responds almost linearly -- a 7 %p cut buys 8 %.
# With an MLA the trapped light is recycled, and every recycle crosses the
# electrode twice; the geometric sum q/(1-(1-q)S) then amplifies whatever
# survives. The weaker the MLA per encounter (small q, more round trips), the
# more violently the electrode matters.
#
# WHAT THIS MEANS FOR THE PAPER'S EXPERIMENT. The claim "lower electrode
# absorption raises EQE_MLA" is not just true, it is a MULTIPLICATIVE
# statement: seed-layer work and outcoupling work compound rather than add.
# The right experiment is 2 x N -- each electrode condition measured WITH and
# WITHOUT the MLA -- because:
#   * the no-MLA column is the control. It should move only ~6-13 %, and its
#     near-flatness is what proves the large MLA-column gain comes from
#     recycling and not from some other change in the stack;
#   * the ratio of the two columns IS the physics: plot MLA gain factor vs
#     measured electrode absorption and the upward trend is the result;
#   * q is not a free parameter afterwards -- the measured gain factor at a
#     known absorption fixes it, so the experiment calibrates its own model.
#
# CAVEAT. Photon-budget model: the cavity's first-pass outcoupling is held
# fixed (only its attenuation by the electrode varies) and the MLA is
# described by a single q. Absolute EQE is illustrative; the ratios between
# electrode conditions, which is what the experiment reads, are the robust
# part.
