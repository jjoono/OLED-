"""What is left to win, and can an alloy win any of it?

Absorption follows resistivity, so the way to cut one-pass loss is to cut
resistivity, and the first thing to do is see where it currently goes. At 8 nm
on HATCN the 7.28 uOhm cm splits three ways, and the three respond to entirely
different interventions -- one of which alloying cannot touch at all.

THE SURFACE TERM IS OFF LIMITS TO ALLOYING. Fuchs-Sondheimer gives

    rho_film = rho_0 * (1 + 0.375 (1-p) l / d)

and adding impurities raises rho_0 while shortening l by the same factor,
because rho_0 * l is a material constant. The surface contribution
0.375 (1-p) rho_0 l / d therefore does not move. An alloy pays the impurity
term in full and gets nothing back from the largest single contribution.

So an alloy is only worth it if it buys something else: a smaller
grain-boundary term, a higher specularity p, or a film that closes thinner. This
puts a number on how much it has to buy to break even, per dopant.

Residual resistivities are dilute-limit literature values in uOhm cm per at.%,
good to maybe 30%. The ORDERING is what the argument rests on and that is solid:
Linde's rule makes the increment scale with the square of the valence
difference, so monovalent Cu and Au sit far below divalent Mg and trivalent Al.
"""
import numpy as np

RHO_BULK, MFP = 1.59, 52.0
D0, RS0 = 8.0, 9.1
P_SPEC = 0.0

# uOhm cm per at.%, dilute limit in a silver host
DOPANTS = [
    ("Cu", 0.10, "monovalent, same column"),
    ("Au", 0.36, "monovalent, near-perfect size match"),
    ("Pd", 0.44, ""),
    ("Zn", 0.63, "divalent"),
    ("Mg", 2.00, "divalent, large size mismatch"),
    ("Al", 2.60, "trivalent"),
    ("Ge", 5.00, "tetravalent"),
]


def surface_term(d, p=P_SPEC):
    """0.375 (1-p) rho_0 l / d -- invariant under alloying, since rho_0*l is."""
    return 0.375 * (1 - p) * RHO_BULK * MFP / d


def main():
    rho = RS0 * D0 * 0.1
    surf = surface_term(D0)
    gb = rho - RHO_BULK - surf
    print(f"HATCN 5 / Ag {D0:.0f} nm, rho = {rho:.2f} uOhm cm\n")
    for name, v, note in (("bulk silver", RHO_BULK, "irreducible"),
                          ("surface (Fuchs, p=0)", surf, "thickness, or specularity p"),
                          ("grain boundaries", gb, "deposition rate, cooling, annealing")):
        print(f"  {name:<22} {v:5.2f}  {v/rho*100:4.1f} %   lever: {note}")

    print(f"\nAlloying cannot touch the surface term. It can at best erase the")
    print(f"grain-boundary term, {gb:.2f} uOhm cm. Break-even doping is where the")
    print(f"impurity term equals that:\n")
    print(f"  {'dopant':<7} {'uOhm cm/at%':>12} {'break-even':>12}   note")
    for name, c, note in DOPANTS:
        print(f"  {name:<7} {c:>12.2f} {gb/c:>10.1f} at%   {note}")

    print(f"\nAt 5 at%, a typical smoothing level:\n")
    print(f"  {'dopant':<7} {'rho if gb -> 0':>15} {'vs now':>9} {'A change':>10}")
    for name, c, _ in DOPANTS:
        r = RHO_BULK + surf + 5.0 * c
        print(f"  {name:<7} {r:>15.2f} {r/rho:>8.2f}x {(r/rho-1)*100:>+9.0f} %")

    print("\nOnly Cu and Au come out ahead, and only because their impurity cost is")
    print("small enough to be paid out of the grain-boundary savings. Mg and Al,")
    print("the usual smoothing choices, cost several times what they could save.")

    print("\n" + "=" * 68)
    print("THE LEVERS, RANKED BY WHAT THEY ARE WORTH HERE")
    print("=" * 68)
    A0 = 2.60                                    # device one-pass A at 8 nm, %
    cases = [
        ("as built", rho, D0),
        ("grain boundaries removed", RHO_BULK + surf, D0),
        ("specularity p = 0.5", RHO_BULK + surface_term(D0, 0.5) + gb, D0),
        ("both", RHO_BULK + surface_term(D0, 0.5), D0),
        ("both, + 5 at% Cu", RHO_BULK + surface_term(D0, 0.5) + 0.5, D0),
        ("both, at 7 nm", RHO_BULK + surface_term(7.0, 0.5), 7.0),
    ]
    print(f"  {'':<26} {'rho':>6} {'Rs':>7} {'one-pass A':>12}")
    for label, r, d in cases:
        rs = r / (d * 0.1)
        # A scales as eps2 * d, and eps2 as rho
        a = A0 * (r / rho) * (d / D0)
        print(f"  {label:<26} {r:>6.2f} {rs:>6.1f}  {a:>11.2f} %")

    print("\n  Absorption goes as rho*d and sheet resistance as rho/d, so thickness")
    print("  trades one against the other while resistivity improves both. That is")
    print("  why the work belongs in the film quality, not in the thickness -- 7 to")
    print("  8 nm was already the optimum and still is.")
    print("\n  And none of this touches the architecture. Removing the trapped modes")
    print("  with an index-matched substrate took the light lost at the electrode")
    print("  from 44.7 % to 7.8 %. Every row above is a refinement inside that.")


if __name__ == "__main__":
    main()
