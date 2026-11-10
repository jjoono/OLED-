"""What grain size does the measured resistivity require, and does SEM agree?

The 8 nm film on HATCN carries 7.28 uOhm cm. Fuchs-Sondheimer surface
scattering accounts for 5.47 of that, leaving 1.81 uOhm cm that has to come from
grain boundaries. Mayadas-Shatzkes turns that number into a grain size -- given
a reflection coefficient R, which is the one thing transport cannot separate.

So the prediction is a curve, not a number, and SEM is what collapses it. A
measured grain size picks out R, and if the R it picks is outside the 0.3-0.5
that silver films are normally reported at, something in the picture is wrong:
either the bright features are surface topography rather than grains, or the
thickness is not what the crystal monitor said.

The commonly quoted form rho/rho0 = 1 + 1.5 alpha is a small-alpha expansion and
alpha here is around 0.7, so the full expression is solved numerically instead;
the difference is not cosmetic at this alpha and is printed for comparison.
"""
import numpy as np
from scipy.optimize import brentq

RHO_BULK = 1.59      # uOhm cm, Ag at 300 K
MFP = 52.0           # nm
D_AG = 8.0           # nm
RS = 9.1             # ohm/sq, measured on HATCN 5 / Ag 8
P_SPEC = 0.0         # Fuchs specularity, worst case


def ms_ratio(a):
    """Mayadas-Shatzkes rho_grain / rho_bulk, full expression."""
    return 1.0 / (3.0 * (1.0 / 3.0 - a / 2.0 + a**2 - a**3 * np.log(1.0 + 1.0 / a)))


def main():
    rho = RS * D_AG * 0.1                      # 1 ohm nm = 0.1 uOhm cm
    rho_surf = RHO_BULK * (1.0 + 0.375 * (1.0 - P_SPEC) * MFP / D_AG)
    d_gb = rho - rho_surf
    target = 1.0 + d_gb / RHO_BULK

    print(f"HATCN 5 / Ag {D_AG:.0f} nm,  Rs = {RS} ohm/sq")
    print(f"  rho measured        {rho:6.2f} uOhm cm")
    print(f"  Fuchs surface term  {rho_surf:6.2f}")
    print(f"  grain-boundary term {d_gb:6.2f}   -> rho_gb/rho_bulk = {target:.3f}\n")

    a_full = brentq(lambda a: ms_ratio(a) - target, 1e-6, 50.0)
    a_lin = (target - 1.0) / 1.5
    print(f"  alpha, full Mayadas-Shatzkes  {a_full:.3f}")
    print(f"  alpha, 1+1.5a expansion       {a_lin:.3f}   "
          f"({(a_lin/a_full-1)*100:+.0f}% -- the expansion is not valid here)\n")

    print(f"  grain size D required, as a function of the reflection coefficient")
    print(f"  {'R':>5} {'D (nm)':>9}   {'':<3}")
    for R in (0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6):
        D = MFP * (R / (1 - R)) / a_full
        note = "  <- silver is usually reported here" if R in (0.4, 0.45) else ""
        print(f"  {R:>5.2f} {D:>9.1f}{note}")

    print("\n  read the SEM the other way round: a measured D implies")
    print(f"  {'D (nm)':>7} {'R':>7}")
    for D in (10, 15, 20, 25, 30, 40, 50, 70):
        x = a_full * D / MFP
        R = x / (1 + x)
        print(f"  {D:>7.0f} {R:>7.2f}")

    print("\n  a film this thin cannot have grains much larger than a few times its")
    print(f"  thickness in the growth direction, but lateral grains of 20-50 nm on")
    print(f"  an {D_AG:.0f} nm film are ordinary. Anything below ~10 nm would mean the")
    print("  bright features are not grains.")


if __name__ == "__main__":
    main()
