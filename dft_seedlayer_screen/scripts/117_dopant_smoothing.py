"""One percent of another metal, for smoothness, at no cost in absorption?

Two questions are tangled here and they pull in opposite directions.

WHAT A DOPANT COSTS is settled by Linde's rule (script 107): the impurity
resistivity per atomic percent scales with the square of the valence
difference. At 1 at% copper costs 0.10 uOhm cm against a film at 7.28; aluminium
costs 2.6. So "one percent at no cost" is possible -- for the monovalent
dopants, and only those.

WHAT A DOPANT DOES is the problem. A dopant smooths a film by pinning adatoms:
less surface mobility, more nuclei, smaller grains, less height variation
between them. But smaller grains mean more grain boundaries, and at 8 nm on
HATCN the grain boundaries are the ONLY part of the absorption film quality can
still reach (script 115). Roughness at the nanometre scale costs nothing --
photons cannot see it, and the electrons already see p = 0 through the organic
underneath. So a dopant that buys smoothness by pinning spends the one currency
that matters to buy the one that does not.

The route that is not self-defeating is a SURFACTANT rather than a solute: a
species that floats on the growth front, lowers the step-edge barrier so
arriving atoms fill in laterally instead of stacking, and does not incorporate.
Sb, Bi, Pb and oxygen do this on silver in the epitaxy literature. The catch
is where a floating species ends up -- at the top interface -- and a monolayer
of bismuth there is itself lossy. Oxygen is the one that is not.

What one percent of anything CAN legitimately buy is a continuous film at a
thickness where pure silver is still islands. That is worth a great deal on a
weak seed. On HATCN, which closes at 5 nm, there is nothing left for it to buy.
"""
RHO_BULK, MFP, D = 1.59, 52.0, 8.0
RHO, GB, A0 = 7.28, 1.81, 2.60           # 8 nm on HATCN: total, grain-boundary, one-pass %
COST = {"Cu": 0.10, "Au": 0.36, "Zn": 0.63, "Mg": 2.00, "Al": 2.60}   # uOhm cm / at%


def A_of(rho):
    return A0 * rho / RHO


def main():
    print(f"8 nm on HATCN: rho {RHO:.2f} uOhm cm, of which grain boundaries {GB:.2f}, "
          f"one-pass A {A0:.2f} %\n")
    print("cost of 1 at% as a solute, before it does anything to the film:")
    for el, c in COST.items():
        print(f"  {el:<3} +{c:.2f} uOhm cm  ->  A {A_of(RHO + c):.2f} %  ({c/RHO*100:+.0f} %)")

    print("\nwhat pinning does: grain size D -> D/x, grain-boundary term x-fold")
    print(f"  {'grain size':>12} {'gb term':>8} {'rho':>6} {'A, no dopant':>13} "
          f"{'A, +1% Cu':>10} {'A, +1% Al':>10}")
    for x in (1.0, 1.5, 2.0, 3.0):
        gb = GB * x
        rho = RHO - GB + gb
        print(f"  {'D' if x == 1 else f'D/{x:g}':>12} {gb:>8.2f} {rho:>6.2f} "
              f"{A_of(rho):>12.2f} % {A_of(rho + 0.10):>9.2f} % {A_of(rho + 2.6):>9.2f} %")
    print("\n  Halving the grain size to get a smoother film raises absorption by a")
    print("  quarter, and the dopant's own cost sits on top. The smoothest film is")
    print("  the most absorbing one.")

    print("\nthe other direction -- a surfactant that grows grains instead:")
    for x in (1.0, 0.5, 0.25, 0.0):
        gb = GB * x
        print(f"  gb x {x:<4g} rho {RHO - GB + gb:.2f}  A {A_of(RHO - GB + gb):.2f} %")
    print("  ceiling 1.95 %, reached only with no boundaries at all. And a floating")
    print("  species ends at the Ag/CPL interface: a monolayer of Bi or Sb there is")
    print("  a lossy film in its own right; oxygen is the one that is not.")


if __name__ == "__main__":
    main()
