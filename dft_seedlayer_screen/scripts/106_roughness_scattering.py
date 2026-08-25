"""Does the granular texture on HATCN scatter enough to matter?

The SEM shows the continuous film on HATCN is not smooth -- it has a rice-grain
texture. Two different questions follow, and they have opposite answers.

FAR-FIELD SCATTERING is negligible. Rayleigh-Rice gives the scattered fraction
off a rough mirror as 1 - exp[-(4 pi sigma n cos(theta) / lambda)^2], and for an
8 nm film sigma cannot be more than a nanometre or two before the film stops
existing. That puts scattering in the tenths of a percent.

ROUGHNESS-INDUCED PLASMON COUPLING is not negligible, and it is the reason the
answer is not simply "no". A flat metal surface cannot absorb normally incident
light into a surface plasmon: the plasmon carries more in-plane momentum than
the photon has. Roughness supplies the missing momentum, so a textured film
opens an absorption channel a flat one does not have. That energy ends as heat
in the silver, so it is absorption, not scattering.

This matters for how the DC route was described earlier in the project. Sheet
resistance measures electron damping, and the plasmon dissipates through that
same damping -- but the COUPLING is geometry, and a flat-film transfer matrix
does not know about it. So the DC-derived one-pass absorption is a lower bound
in two senses, not one: it excludes scattering, and it excludes the extra
plasmon coupling that roughness buys.

The measurement can tell which is happening, and for these two films it does.
"""
import numpy as np

L = 550.0
N_MEAS, N_DEV = 1.0, 1.80        # medium above the film: air, then the organic
# 8 nm, 550 nm, from the T/R inversion and the Fuchs floor
A_MEAS, A_FLOOR = 0.0732, 0.0583


def tis(sigma, n, theta_deg=0.0):
    """Total integrated scatter off a rough reflector, Rayleigh-Rice."""
    return 1.0 - np.exp(-(4 * np.pi * sigma * n * np.cos(np.deg2rad(theta_deg)) / L)**2)


def sigma_for(frac, n, theta_deg=0.0):
    """RMS roughness that would scatter this fraction."""
    x = np.sqrt(-np.log(1.0 - frac))
    return x * L / (4 * np.pi * n * np.cos(np.deg2rad(theta_deg)))


def main():
    print(f"Rayleigh-Rice scatter at {L:.0f} nm\n")
    print(f"{'RMS (nm)':>9} {'in air':>10} {'in organic n=1.8':>18} "
          f"{'organic, 60 deg':>17}")
    for s in (0.5, 1.0, 1.5, 2.0, 3.0, 5.0):
        print(f"{s:>9.1f} {tis(s, N_MEAS)*100:>9.2f}% "
              f"{tis(s, N_DEV)*100:>17.2f}% {tis(s, N_DEV, 60)*100:>16.2f}%")

    print("\nNote which way the index goes. Embedding the metal in a high-index")
    print("medium makes roughness scattering WORSE, by n^2 -- the wavelength in")
    print("the medium is shorter, so the same bumps are optically larger. An")
    print("earlier suggestion in this project that an index-matched sandwich")
    print("would suppress scattering was wrong for a metal surface. It holds for")
    print("a rough interface between two dielectrics of similar index, where the")
    print("CONTRAST is what scatters; a metal keeps its contrast whatever it is")
    print("embedded in.\n")

    excess = A_MEAS - A_FLOOR
    s_needed = sigma_for(excess, N_MEAS)
    print(f"HATCN 5 / Ag 8 nm: measured A {A_MEAS*100:.2f}%, Fuchs floor "
          f"{A_FLOOR*100:.2f}%, excess {excess*100:.2f} %p")
    print(f"For scattering alone to account for that excess would take")
    print(f"RMS = {s_needed:.2f} nm on an 8 nm film. That is most of the film")
    print("thickness; a film that rough would show holes, and the SEM shows none.")
    print("So on HATCN the excess is genuinely absorbed, not scattered.\n")

    print("The MoOx film gives the opposite answer, and by the same arithmetic in")
    print("reverse: it is 1.19x worse than HATCN electrically and 3.32x worse")
    print("optically, and holes big enough to see are exactly what scatters")
    print("without impeding current. Two films, two mechanisms, and the pair of")
    print("SEM images plus a four-point probe distinguishes them.\n")

    print("WHAT WOULD SETTLE THE PLASMON QUESTION")
    print("  AFM gives sigma. With sigma in hand the scattered fraction above is")
    print("  fixed, and whatever excess absorption remains after subtracting the")
    print("  grain-boundary term is the roughness-plasmon channel. The three")
    print("  contributions are then separated with no free parameters:")
    print("    Fuchs surface term      from thickness")
    print("    grain-boundary term     from resistivity")
    print("    scattering              from AFM sigma")
    print("    remainder               roughness-assisted plasmon coupling")


if __name__ == "__main__":
    main()
