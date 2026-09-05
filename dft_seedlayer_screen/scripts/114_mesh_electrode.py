"""A metal mesh over an organic that carries the current sideways -- does it work?

The idea: instead of a continuous 8 nm silver film, put down metal lines and let
the HIL (HATCN, or an NDP-9-doped transporter) conduct laterally between them.
Light through the openings never meets silver, so the one-pass absorption
should drop. Two questions, answered separately, because they get different
answers.

OPTICS. A mesh with fill fraction f loses roughly f*A_line + (1-f)*A_open per
pass. Thick lines reflect rather than absorb, so A_line is the ~3% of bulk
silver, not much worse than the thin film's 2.5% -- the gain is entirely the
(1-f) that sees only organic.

ELECTRICS. Between lines the current has to travel sideways through the organic,
half a pitch each way. A strip of pitch p, fed from both edges, drawing a uniform
J out of its face, drops

    dV = J * Rs_HIL * (p/2)^2 / 2

from edge to centre. That is the constraint, and it is brutal, because Rs of an
organic is measured in gigaohms per square where silver is measured in ohms.
The pitch it allows is then set against the optical regimes: below lambda/n the
mesh is a wire-grid polariser, below a few lambda it is a diffraction grating,
and neither is an electrode any more.
"""
import math

L, N_ORG = 550.0, 1.8            # nm
J = 100.0                        # A/m^2  = 10 mA/cm^2, a display-relevant current
DV_MAX = 0.10                    # V, tolerable edge-to-centre drop
A_THIN, A_OPEN = 0.0245, 0.0056  # one-pass: Ag 8 nm measured; HATCN 5 nm alone
A_THICK = 0.03                   # bulk silver mirror, per hit

# lateral conductivity, S/cm, literature order of magnitude
HIL = [
    ("HATCN neat, 5 nm",              1e-6,   5.0),
    ("HATCN neat, 50 nm",             1e-6,  50.0),
    ("NDP-9:HTL 5 wt%, 30 nm",        1e-4,  30.0),
    ("NDP-9:HTL 10 wt%, 100 nm",      1e-3, 100.0),
    ("best reported p-doped, 100 nm", 1e-2, 100.0),
    ("Ag 8 nm on HATCN (reference)",  1.0 / (7.28e-6), 8.0),
]


def rs(sigma_S_cm, d_nm):
    return 1.0 / (sigma_S_cm * 100.0 * d_nm * 1e-9)          # ohm/sq


def max_pitch(rs_):
    return 2.0 * math.sqrt(2.0 * DV_MAX / (J * rs_)) * 1e6     # um


def main():
    print(f"tolerated drop {DV_MAX} V at {J/10:.0f} mA/cm2\n")
    print(f"{'lateral conductor':<32} {'sigma':>8} {'Rs':>11} {'max pitch':>11}   regime")
    print("-" * 84)
    for name, s, d in HIL:
        r = rs(s, d)
        p = max_pitch(r)
        p_nm = p * 1000
        if p_nm < L / N_ORG / 2:
            reg = "wire-grid polariser -- not an electrode"
        elif p_nm < 3 * L:
            reg = "diffraction grating"
        elif p < 20:
            reg = "lithography on top of the organics"
        else:
            reg = "shadow mask feasible"
        print(f"{name:<32} {s:>8.0e} {r:>9.1e} {p:>9.2f} um   {reg}")

    print("\nShadow-mask lines are 20 um pitch at best. The Rs that allows that:")
    r_need = 2 * DV_MAX / (J * (10e-6)**2)
    print(f"  Rs <= {r_need:.0f} ohm/sq, i.e. sigma >= {1/(r_need*100*50e-9):.1e} S/cm at 50 nm")
    print("  -- roughly the conductivity of the silver film the mesh was meant to")
    print("  replace. No organic is within six orders of magnitude of it.\n")

    print("optics, per pass, if the electrics somehow worked:")
    print(f"{'fill f':>7} {'thin lines':>11} {'thick lines':>12}   vs continuous {A_THIN*100:.2f}%")
    for f in (0.05, 0.10, 0.20, 0.30):
        a1 = f * A_THIN + (1 - f) * A_OPEN
        a2 = f * A_THICK + (1 - f) * A_OPEN
        print(f"{f:>7.2f} {a1*100:>10.2f}% {a2*100:>11.2f}%")
    print("\n  A real gain -- three to four times -- and the reason the idea is")
    print("  attractive. It is not reachable, because the pitch that makes the")
    print("  organic conduct is the pitch that turns the mesh into a grating.")


if __name__ == "__main__":
    main()
