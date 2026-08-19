"""Optical constants of ITO / IZO derived from transport, not assumed.

A TCO's visible loss is free-carrier (Drude) absorption, and the same carriers
set the sheet resistance, so n and k follow from (N, mu, m*, eps_inf) with no
free optical parameter:

    wp^2 = N e^2 / (eps0 m*)          gamma = e / (m* mu)
    eps  = eps_inf - wp^2/(w^2 + i gamma w)
    rho  = 1/(N e mu)                 Rs = rho/d

That gives the FLOOR.  Real films add sub-gap absorption from oxygen vacancies
and disorder, which room-temperature deposition cannot anneal out; that excess
is quoted separately as a range rather than derived.
"""
import numpy as np

E, EPS0, ME, HBAR = 1.602176634e-19, 8.8541878128e-12, 9.1093837015e-31, 1.054571817e-34

CASES = [
    # label,                    N (cm^-3), mu (cm^2/Vs), m*/me, eps_inf, d (nm)
    ("ITO, annealed 300 C",       8.0e20,  35.0, 0.40, 4.0, 150),
    ("ITO, RT sputtered",         5.0e20,  20.0, 0.40, 4.0, 100),
    ("IZO, RT sputtered",         4.0e20,  40.0, 0.30, 3.8, 100),
]
# excess k at 550 nm from defect / sub-gap absorption, beyond the Drude floor
DEFECT_K = {"ITO, annealed 300 C": (0.000, 0.008),
            "ITO, RT sputtered":   (0.005, 0.030),
            "IZO, RT sputtered":   (0.003, 0.020)}


def drude(N_cm3, mu_cm2, mstar, eps_inf, lam_nm):
    N = N_cm3*1e6
    m = mstar*ME
    wp2 = N*E**2/(EPS0*m)                       # rad^2/s^2
    gam = E/(m*mu_cm2*1e-4)                     # rad/s
    w = 2*np.pi*2.99792458e17/lam_nm            # rad/s
    eps = eps_inf - wp2/(w**2 + 1j*gam*w)
    nk = np.sqrt(eps)
    rho = 1.0/(N*E*mu_cm2*1e-4)*100             # Ohm cm
    return nk, rho, HBAR*np.sqrt(wp2)/E, HBAR*gam/E


print("Derived from transport -- no optical parameter is fitted\n")
print(f"{'':<24}{'hw_p':>7}{'h_gam':>7}{'rho':>9}{'Rs':>8} | "
      f"{'n550':>7}{'k550(Drude)':>13}{'k550(+defects)':>16}")
print("-"*100)
for lab, N, mu, ms, ei, d in CASES:
    nk, rho, hwp, hg = drude(N, mu, ms, ei, 550.0)
    Rs = rho/(d*1e-7)
    lo, hi = DEFECT_K[lab]
    print(f"{lab:<24}{hwp:6.2f}eV{hg:6.3f}eV{rho*1e6:8.1f}u{Rs:7.1f} | "
          f"{nk.real:7.3f}{nk.imag:13.4f}{f'{nk.imag+lo:.3f} - {nk.imag+hi:.3f}':>16}")

print("\nDispersion of the Drude part (k grows as lam^3 in the transparent window):")
print(f"{'':<24}" + "".join(f"{l:>9}" for l in [450, 550, 650, 750]))
for lab, N, mu, ms, ei, d in CASES:
    row = f"{lab:<24}"
    for l in [450, 550, 650, 750]:
        nk, *_ = drude(N, mu, ms, ei, l)
        row += f"{nk.imag:9.4f}"
    print(row)
