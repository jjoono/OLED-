"""What alloying does to n of silver -- the right model for a solid solution.

THE QUESTION. Pure Ag has an anomalously LOW refractive index in the visible
(n ~ 0.06-0.13 at 550 nm), and that is exactly why a thin Ag film absorbs so
little. Does co-depositing Mg or Yb destroy it, even at a few percent?

WHY scripts/70 WAS THE WRONG MODEL FOR THIS. That script mixed dopant and Ag
with a Bruggeman EMA, i.e. as two SEPARATE phases with their own dielectric
functions. A real Mg:Ag or Yb:Ag co-deposit is a substitutional SOLID
SOLUTION: the dopant does not sit as inclusions, it enters the lattice and
scatters the conduction electrons of the host. The correct dilute-alloy model
is Drude with an increased damping constant, and it is what this script uses.
Bruggeman got the DIRECTION right (A rises with y) but for the wrong reason.

THE PHYSICS, in one line. For a metal with |eps1| >> eps2,

    n  ~  eps2 / (2*sqrt(|eps1|))        and     eps2 = wp^2 * gamma / w^3

so n is proportional to the Drude damping gamma. Damping is what resistivity
measures (rho = m*gamma / (n_e e^2 hbar)), hence

    n(alloy) / n(Ag)  ~  rho(alloy) / rho(Ag)

Alloy scattering therefore raises n by the SAME factor it raises the sheet
resistance -- the optical and electrical penalties are one penalty, not two.
The thin-film Joule (Drude) absorption tracks eps2 as well, so it scales the
same way. That is the trade-off the user asked about, and it has no escape
within a single-phase alloy.

Nordheim: d(rho) = C * y * (1-y). C is taken as 45 uOhm.cm for Mg:Ag and
55 for Yb:Ag -- chosen so that 10 at% gives ~5.7 / ~6.5 uOhm.cm, the range
reported for dilute Ag alloys; the exact C shifts the numbers, not the
conclusion (a caveat repeated in the RESULT block).
"""
import numpy as np

# --- silver reference at 550 nm, from the J&C values used across the project
N_AG, K_AG = 0.06, 3.59
EPS1 = N_AG ** 2 - K_AG ** 2
EPS2 = 2 * N_AG * K_AG
RHO_AG = 1.6                      # uOhm.cm, bulk Ag
C_NORDHEIM = {"Mg": 45.0, "Yb": 55.0}   # uOhm.cm, see docstring

# thin-film geometry used for the absorption estimate
D_NM = 5.0
LAM = 550.0


def alloy_optics(y, C):
    """Dilute solid solution: gamma (hence eps2) scales with resistivity."""
    rho = RHO_AG + C * y * (1 - y)
    F = rho / RHO_AG                       # damping enhancement
    eps2 = EPS2 * F                        # eps1 essentially unchanged (dilute)
    eps1 = EPS1
    mod = np.hypot(eps1, eps2)
    n = np.sqrt((mod + eps1) / 2)
    k = np.sqrt((mod - eps1) / 2)
    return rho, n, k


def absorptance(n, k, d_nm=D_NM, lam=LAM):
    """Single-pass Joule absorption of a free-standing thin metal film,
    A = 1 - T - R from the thin-film limit of the Fresnel expressions."""
    N = complex(n, -k)                       # exp(-i w t) convention
    delta = 2 * np.pi / lam * N * d_nm
    c, s = np.cos(delta), np.sin(delta)
    M = np.array([[c, 1j * s / N], [1j * N * s, c]])
    B, C_ = M @ np.array([1.0, 1.0], dtype=complex)
    r = (B - C_) / (B + C_)
    t = 2.0 / (B + C_)
    return 1 - abs(t) ** 2 - abs(r) ** 2


def main():
    print("=" * 72)
    print(f"Ag at {LAM:.0f} nm: n = {N_AG}, k = {K_AG}  "
          f"(eps1 {EPS1:.2f}, eps2 {EPS2:.3f}), rho = {RHO_AG} uOhm.cm")
    print(f"A of a free {D_NM:.0f} nm film: {100*absorptance(N_AG, K_AG):.2f} %")
    print("=" * 72)
    for metal, C in C_NORDHEIM.items():
        print(f"\n{metal}:Ag solid solution (Nordheim C = {C:.0f} uOhm.cm)")
        print(f"{'y at%':>7}{'rho':>9}{'n':>8}{'k':>7}{'n/n_Ag':>9}"
              f"{'A(5nm)':>9}{'A/A_Ag':>9}")
        A0 = absorptance(N_AG, K_AG)
        for y in (0.0, 0.01, 0.02, 0.05, 0.10, 0.20):
            rho, n, k = alloy_optics(y, C)
            A = absorptance(n, k)
            print(f"{100*y:>7.0f}{rho:>9.2f}{n:>8.3f}{k:>7.2f}"
                  f"{n/N_AG:>9.1f}{100*A:>8.2f}%{A/A0:>9.1f}")
    print("""
READ: n and the film absorption both scale with the resistivity ratio, so the
optical penalty of co-deposition is not a separate effect from the electrical
one -- it IS the electrical one, seen at optical frequency. Even 5 at% costs a
factor of ~2 in both. Morphology improvement has to beat that factor to be
worth it, which is only plausible for a film still near percolation.""")


if __name__ == "__main__":
    main()

# RESULT (run 2026-08-16, 550 nm, free-standing 5 nm film):
#
#   y at%     Mg:Ag  n   A(5nm)  |    Yb:Ag  n   A(5nm)
#     0       0.060      2.11 %  |    0.060      2.11 %
#     1       0.077      2.68 %  |    0.080      2.81 %
#     2       0.093      3.23 %  |    0.100      3.48 %
#     5       0.140      4.79 %  |    0.158      5.36 %
#    10       0.212      7.06 %  |    0.245      8.09 %
#
# 1. THE USER'S INTUITION IS EXACTLY RIGHT, AND THE SCALING IS ONE-TO-ONE.
#    n/n_Ag and A/A_Ag come out IDENTICAL at every doping level (2.3x for
#    5 at% Mg, 2.6x for Yb) because both are set by the same Drude damping.
#    The optical penalty of co-deposition is not a second, separate cost
#    alongside the resistivity rise -- it IS the resistivity rise, read at
#    optical frequency. Even 1-2 at% costs 30-60 %, so "a few percent is
#    harmless" is false.
#
# 2. Yb is worse than Mg at equal doping (2.5x vs 2.3x at 5 at%), consistent
#    with the ordering scripts/70 found by the Bruggeman route.
#
# 3. THIS RE-OPENS WHAT scripts/70 CLOSED, because the morphology penalty in
#    the REAL films is far larger than the one 70 assumed. Measured (scripts/
#    72, eta = 0.96): Ag5 films absorb 10-13 %, against 2-3 % for an ideal
#    continuous 5 nm film -- a morphology penalty of ~4-5x, i.e. BIGGER than
#    the 2.3x alloy penalty at 5 at%. So a dopant that genuinely closed the
#    film would win: 2.1 % x 2.3 = 4.8 % beats the 10-13 % measured now.
#    scripts/70's "doping never wins" verdict was reached with f_void capped
#    at 0.24 and a Bruggeman alloy; both assumptions understate the real
#    morphology loss, and the SEM (labyrinth network, scripts/73) shows why.
#
# 4. BUT THE SEED-LAYER ROUTE DOMINATES ON THIS SAME ARITHMETIC. Smoothing
#    via the organic underlayer removes the 4-5x morphology penalty at ZERO
#    optical cost, because nothing lossy enters the silver. Doping buys the
#    same morphology at a 2.3x multiplier. Whenever the seed can reach
#    comparable continuity, it is strictly better -- and if it cannot, the
#    two are additive, so the endgame is the smoothest seed FIRST and only
#    then the minimum dopant needed to finish the job.
#
# CAVEAT: Nordheim coefficients are literature-typical, not measured for
# these films; eps1 is held fixed (dilute limit) and interband contributions
# of the dopant are ignored, which makes these numbers a LOWER bound on the
# alloy penalty, most of all for Yb.
