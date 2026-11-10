"""How smooth does the silver have to be, and how smooth can evaporated silver get?

Three questions hide in "why is it not smooth even on HATCN".

1. IS THERE LOCALISED-PLASMON ABSORPTION? Compare what the optics measured with
   what the resistivity predicts (script 96). Whatever the DC route cannot see
   -- scattering, island resonances, anything not electron damping -- shows up
   as the difference. It is there below percolation and gone by 7 nm.

2. WHAT DOES ROUGHNESS COST, AND TO WHOM? Photons and electrons see it on
   different scales. For a photon at 550 nm a nanometre of roughness is nothing
   (script 106). For an electron the yardstick is the Fermi wavelength, 0.52 nm
   in silver, and Soffer's expression for the specularity,

       p = exp[-(4 pi sigma / lambda_F)^2]

   goes to zero for any roughness a thermal evaporator can produce. This
   retracts the "p = 0.5" row in script 107: that row asked for an interface
   flat to a tenth of an angstrom, which is an epitaxial film, not an
   evaporated one. The specularity lever is not available.

3. WHY IS IT ROUGH AT ALL? Silver on an organic is Volmer-Weber whatever the
   seed: gamma_Ag is several times gamma_organic, so islands are the
   thermodynamic answer and a seed only raises how many there are. After they
   coalesce the film is polycrystalline with grains the size of the island
   spacing, and its roughness is the height spread of those grains. Smoothness
   has to be bought kinetically -- by denying the adatoms the mobility to find
   the thermodynamic answer -- and every way of doing that has a price in
   resistivity, which is a price in absorption.

Written against the standard library only.
"""
import math

LAMBDA_F = 0.52          # nm, silver Fermi wavelength

# HATCN 5 / Ag d at 550 nm, from script 96: measured 1-T-R, and the absorptance
# a flat film with the measured resistivity would have
SERIES = [(4, 14.59, 11.29), (5, 11.17, 8.28), (6, 9.72, 9.20),
          (7, 8.09, 8.18), (8, 7.32, 7.32), (10, 7.04, 8.52), (12, 6.27, 7.61)]


def main():
    print("1. localised-plasmon / island excess: measured minus DC-predicted\n")
    print(f"{'d':>4} {'A meas':>8} {'A from Rs':>10} {'excess':>8}")
    for d, am, adc in SERIES:
        ex = am - adc
        tag = ("  <- below/at percolation: island resonances live here" if ex > 1
               else ("  <- gone" if abs(ex) < 0.6 else
                     "  <- negative: DC over-counts grain boundaries (script 96)"))
        print(f"{d:>4} {am:>7.2f}% {adc:>9.2f}% {ex:>+7.2f}%p{tag}")
    print("\n   On HATCN the excess is 3 %p at 4-5 nm and zero from 7 nm on. The")
    print("   localised-plasmon problem is solved by thickness, and 7-8 nm is")
    print("   already past it. What remains at 8 nm is electron damping -- grain")
    print("   boundaries and surface scattering -- not a plasmon.\n")

    print("2. electron specularity vs roughness (Soffer)\n")
    print(f"{'RMS (nm)':>9} {'p':>8}")
    for s in (0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 1.00):
        p = math.exp(-(4 * math.pi * s / LAMBDA_F)**2)
        print(f"{s:>9.2f} {p:>8.3f}")
    print("\n   p = 0.5 needs 0.035 nm RMS -- a tenth of an atomic step. An")
    print("   evaporated film at 0.5-1 nm RMS has p = 0 to all practical purposes.")
    print("   The 'p = 0.5' case in script 107 is withdrawn. Of the 7.28 uOhm cm")
    print("   at 8 nm, the surface term 3.88 is fixed by thickness alone; the only")
    print("   part film quality can still touch is the grain-boundary 1.81.\n")

    rho, gb, A0 = 7.28, 1.81, 2.60
    print("3. the ceiling for a continuous evaporated film at 8 nm on HATCN\n")
    print(f"   grain boundaries removed entirely: rho {rho:.2f} -> {rho-gb:.2f}, "
          f"one-pass A {A0:.2f}% -> {A0*(rho-gb)/rho:.2f}%")
    print("   That is the most any smoothing recipe can deliver at this thickness.")
    print("   Everything larger -- 44.7% -> 7.8% of light lost at the electrode --")
    print("   comes from removing the trapped modes, not from the film.")


if __name__ == "__main__":
    main()
