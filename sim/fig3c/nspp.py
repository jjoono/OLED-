"""Surface-plasmon index at a metal / uniaxial-dielectric interface.

For a uniaxial dielectric whose optic axis is normal to the interface, the TM
wave obeys  k_z^2 = eps_o k0^2 - (eps_o/eps_e) k_x^2.  Matching eps_o/kappa_d
= -eps_m/kappa_m at the interface gives

    (k_SPP/k0)^2 = eps_m (eps_o - eps_m) / (eps_o - eps_m^2/eps_e)

which reduces to the familiar eps_m eps_o/(eps_m + eps_o) when eps_e = eps_o.
"""
import numpy as np


def n_spp(n_metal, n_o, n_e):
    em = np.asarray(n_metal, dtype=complex) ** 2
    eo = np.asarray(n_o, dtype=complex) ** 2
    ee = np.asarray(n_e, dtype=complex) ** 2
    return np.sqrt(em * (eo - em) / (eo - em ** 2 / ee)).real


def n_e_threshold(n_metal, n_o, n_target):
    """The n_e at which n_SPP falls to n_target (bisection)."""
    lo, hi = 0.5, float(np.real(n_o))
    if n_spp(n_metal, n_o, hi) < n_target:
        return float('nan')
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if n_spp(n_metal, n_o, mid) > n_target:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)
