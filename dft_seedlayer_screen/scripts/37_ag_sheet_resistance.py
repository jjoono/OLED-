"""Sheet resistance of an IDEAL continuous Ag film, 1-15 nm.

"Ideal" still is not bulk: once the thickness drops below the bulk mean free path
(lambda = 53 nm for Ag at 300 K) the carriers scatter off the two surfaces, and in
a real polycrystalline film off grain boundaries as well. Both are intrinsic to a
thin film, not defects, so the naive rho_bulk/d is a floor that no film reaches.

Models
------
Fuchs-Sondheimer (surface scattering, specularity p in [0,1]):

    rho_0/rho_FS = 1 - (3/(2k))(1-p) * INT_1^inf (1/t^3 - 1/t^5)
                                        (1 - e^-kt)/(1 - p e^-kt) dt,  k = d/lambda

    p = 1 : fully specular, no surface contribution
    p = 0 : fully diffuse, the usual pessimistic bound for an evaporated film

Mayadas-Shatzkes (grain-boundary scattering, reflection coefficient Rgb):

    rho_0/rho_MS = 3[1/3 - a/2 + a^2 - a^3 ln(1 + 1/a)],  a = (lambda/D) Rgb/(1-Rgb)

    D = lateral grain size. The textbook default D = d (columnar grains) is
    WRONG for this system: the SEM series on HATCN shows lateral grain features
    of ~50-150 nm at 12 nm Ag and ~100-300 nm at 25 nm, i.e. D >> d. Using D = d
    overestimates the resistivity of the 20 and 25 nm films by ~2x against the
    measured values. D is therefore an explicit parameter here, with the
    SEM-informed D = 100 nm as the default.

Combined with Matthiessen: rho = rho_0 + (rho_FS - rho_0) + (rho_MS - rho_0).

WHAT THIS DOES NOT COVER: below the percolation threshold the film is not
continuous and the sheet resistance diverges -- no continuum model applies.
kMC (script 36) puts percolation for Ag on HATCN near 8 nm, so the 1-7 nm rows
here are the continuous-film limit, not a prediction for a real deposit.
"""
import numpy as np
from scipy import integrate

RHO_BULK = 1.587e-8      # Ohm.m, Ag at 300 K
LAMBDA = 53.0            # nm, bulk mean free path of Ag at 300 K


def fuchs_sondheimer(d_nm, p=0.5, lam=LAMBDA):
    k = d_nm / lam

    def integrand(t):
        e = np.exp(-k * t)
        return (1.0 / t ** 3 - 1.0 / t ** 5) * (1.0 - e) / (1.0 - p * e)

    I, _ = integrate.quad(integrand, 1.0, np.inf, limit=200)
    ratio = 1.0 - (3.0 / (2.0 * k)) * (1.0 - p) * I      # rho_0 / rho_FS
    return RHO_BULK / max(ratio, 1e-9)


def mayadas_shatzkes(d_nm, Rgb=0.3, D_nm=100.0, lam=LAMBDA):
    D = d_nm if D_nm is None else D_nm
    a = (lam / D) * Rgb / (1.0 - Rgb)
    f = 3.0 * (1.0 / 3.0 - a / 2.0 + a ** 2 - a ** 3 * np.log(1.0 + 1.0 / a))
    return RHO_BULK / max(f, 1e-9)


def combined(d_nm, p=0.5, Rgb=0.3, D_nm=100.0):
    rho_fs = fuchs_sondheimer(d_nm, p)
    rho_ms = mayadas_shatzkes(d_nm, Rgb, D_nm)
    return RHO_BULK + (rho_fs - RHO_BULK) + (rho_ms - RHO_BULK)


def Rs(rho, d_nm):
    """Ohm per square from resistivity (Ohm.m) and thickness (nm)."""
    return rho / (d_nm * 1e-9)


if __name__ == "__main__":
    ds = np.arange(1, 16)
    print("Sheet resistance of a CONTINUOUS Ag film (Ohm/sq), 300 K")
    print(f"rho_bulk = {RHO_BULK*1e8:.3f} uOhm.cm,  lambda = {LAMBDA:.0f} nm,"
          f"  grain size D = 100 nm (from SEM)\n")
    print(f"{'d':>4}{'floor':>9}{'p=1':>9}{'p=0.5':>9}{'p=0':>9}"
          f"{'+GB':>9}{'+GB':>9}   {'':<8}")
    print(f"{'nm':>4}{'rho_b/d':>9}{'ideal':>9}{'FS':>9}{'FS':>9}"
          f"{'p=0.5':>9}{'p=0':>9}   {'':<8}")
    print("-" * 66)
    for d in ds:
        r_floor = Rs(RHO_BULK, d)
        r_p1 = Rs(fuchs_sondheimer(d, 1.0), d)
        r_p5 = Rs(fuchs_sondheimer(d, 0.5), d)
        r_p0 = Rs(fuchs_sondheimer(d, 0.0), d)
        r_g5 = Rs(combined(d, 0.5, 0.3, 100.0), d)
        r_g0 = Rs(combined(d, 0.0, 0.3, 100.0), d)
        flag = "  below percolation" if d < 8 else ""
        print(f"{d:>4}{r_floor:>9.2f}{r_p1:>9.2f}{r_p5:>9.2f}{r_p0:>9.2f}"
              f"{r_g5:>9.2f}{r_g0:>9.2f}{flag}")

    print("\nvalidation against Park & Suh 2018 (HATCN 7 nm / Ag), measured 4-point:")
    print(f"  {'d':>4}{'measured':>10}{'FS+GB p=0.5':>14}{'ratio':>8}")
    for d, meas in [(15, 10.7), (20, 1.4), (25, 0.9)]:
        pred = Rs(combined(d, 0.5, 0.3, 100.0), d)
        print(f"  {d:>4}{meas:>10.2f}{pred:>14.2f}{meas/pred:>8.1f}x")
    print("""
  The 20 and 25 nm films land within ~10% of the continuous-film model, so the
  SEM-informed grain size is the right one. The 15 nm film is 4x above it: that
  film was still morphology-limited (voids, incomplete percolation), which is
  exactly what its SEM showed. Sheet resistance is therefore the sharpest probe
  of percolation available -- it departs from the continuum model by a factor of
  several the moment the film stops being continuous, while transmittance barely
  moves.""")
