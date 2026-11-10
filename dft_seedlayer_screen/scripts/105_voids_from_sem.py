"""MoOx/Ag 8 nm has holes in it, and that changes what the numbers mean.

SEM at 65 kX shows the 8 nm film on MoOx riddled with voids -- dark irregular
patches roughly 50-150 nm across, covering a substantial fraction of the field.
The film on HATCN at the same thickness has none.

TWO THINGS FOLLOW.

First, a correction. This project has been saying MoOx "closes" at 7 nm, on the
strength of the knee in rho = rho0 + C/d and the step in the optical
absorptance. Both of those mark PERCOLATION -- the thickness at which a
connected path first spans the film -- and a percolating network can be, and
here is, full of holes. Closure in the sense of complete coverage has not
happened at 8 nm on MoOx. The transport and optical knees were read correctly;
the word put on them was wrong.

Second, and more useful: voids separate the electrical and optical penalties,
and they separate them by a lot. Current flows around a hole with a modest
detour. Light at 550 nm meets a scatterer comparable to its own wavelength. So
if the excess optical loss on MoOx were absorption, it would have to show up in
the resistivity too; if it is scattering, it need not. This script puts numbers
on both sides of that.
"""
import numpy as np
from scipy.optimize import brentq

RHO_BULK, MFP = 1.59, 52.0
D = 8.0
RS = {"HATCN": 9.1, "MoOx": 10.8}          # ohm/sq at 8 nm
EPS2 = {"HATCN": 1.525, "MoOx": 5.062}     # from the T/R inversion at 550 nm
EPS2_BULK = 0.316


def sigma_ratio_2d(f):
    """Conductivity of a sheet with an area fraction f of insulating holes.

    Two-dimensional Maxwell-Garnett for circular voids: sigma_eff/sigma =
    (1-f)/(1+f). Exact for dilute, and low by a few percent as f grows toward
    the 2D percolation threshold near 0.5.
    """
    return (1.0 - f) / (1.0 + f)


def main():
    print(f"At {D:.0f} nm, both seeds, 550 nm\n")
    print(f"{'':<8} {'Rs':>6} {'rho app':>9} {'rho/bulk':>9} {'eps2':>7} "
          f"{'eps2/bulk':>10}")
    for k in ("HATCN", "MoOx"):
        rho = RS[k] * D * 0.1
        print(f"{k:<8} {RS[k]:>6.1f} {rho:>9.2f} {rho/RHO_BULK:>9.2f} "
              f"{EPS2[k]:>7.3f} {EPS2_BULK and EPS2[k]/EPS2_BULK:>10.2f}")

    r_elec = (RS['MoOx'] / RS['HATCN'])
    r_opt = EPS2['MoOx'] / EPS2['HATCN']
    print(f"\nMoOx / HATCN:   electrical {r_elec:.2f}x     optical {r_opt:.2f}x")
    print(f"The optical penalty is {r_opt/r_elec:.1f} times the electrical one.\n")

    print("If the extra optical loss on MoOx were absorption, it would be extra")
    print("electron damping, and the resistivity would have to rise by the same")
    print(f"factor. To match eps2 the sheet resistance would have to be")
    print(f"{RS['HATCN'] * r_opt:.1f} ohm/sq; it is {RS['MoOx']:.1f}. So most of what")
    print("1-T-R counts as absorption on MoOx is light thrown sideways by the")
    print("holes, not light absorbed by the silver.\n")

    f = brentq(lambda x: sigma_ratio_2d(x) - RS["HATCN"] / RS["MoOx"], 1e-6, 0.49)
    print(f"Void fraction that explains the resistivity alone: {f*100:.0f} %")
    print("Compare that against the SEM: measure the dark-area fraction. If the")
    print("image shows appreciably more void than this, the silver between the")
    print("holes is BETTER than the sheet resistance suggests, because the")
    print("measurement has been charging it for the detour.\n")

    for fs in (0.10, 0.15, 0.20, 0.25, 0.30):
        rho_app = RS["MoOx"] * D * 0.1
        rho_true = rho_app * sigma_ratio_2d(fs)
        rho_surf = RHO_BULK * (1.0 + 0.375 * MFP / D)
        print(f"  SEM void {fs*100:>3.0f} %  ->  true rho {rho_true:5.2f} uOhm cm"
              f"   (HATCN {RS['HATCN']*D*0.1:.2f}, Fuchs floor {rho_surf:.2f})")

    # The Fuchs floor bounds this from above: correcting for more void than the
    # film can support drives the "true" resistivity below what a perfectly
    # specular film of this thickness could ever have.
    rho_surf = RHO_BULK * (1.0 + 0.375 * MFP / D)
    fmax = brentq(lambda x: RS["MoOx"] * D * 0.1 * sigma_ratio_2d(x) - rho_surf,
                  1e-6, 0.49)
    print(f"\nUPPER BOUND. Past {fmax*100:.0f} % void the corrected resistivity falls")
    print(f"below the Fuchs floor of {rho_surf:.2f} uOhm cm, which no film of this")
    print("thickness can be under. So the SEM should show somewhere between 9 and")
    print(f"{fmax*100:.0f} % dark area. That is a real prediction and worth checking against")
    print("the image.")

    print("\nCONSEQUENCE FOR THE EXTRACTED CONSTANTS. The MoOx inversion treated")
    print("the film as one homogeneous layer. It is not -- it is silver plus")
    print("holes. Its n and k are effective-medium values for that mixture and")
    print("must not be quoted as optical constants of silver. The HATCN films,")
    print("which the SEM shows continuous, are not affected.")


if __name__ == "__main__":
    main()
