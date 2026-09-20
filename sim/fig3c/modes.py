"""Complex TM mode index of the stack, found as a pole of the cavity factor.

A guided or surface mode is a pole of 1/(1 - r_b r_t), i.e. a root of

    f(k_x) = 1 - r_b(k_x) r_t(k_x) = 0

in the complex k_x plane.  The round-trip phase makes f independent of where
the reference plane sits, so the emitter plane may be used directly.

Once the pole is known, the field of that mode decays away from the metal as
exp(-kappa z) with, in a uniaxial layer of optic axis along z,

    kappa = (n_o/n_e) sqrt(k_x^2 - n_e^2 k0^2)

which is the quantity that governs how fast the dipole decouples from the mode
as the layer between it and the metal is thickened.  Note that it is the index
of *that* layer that enters, not the index of the layer the dipole sits in.
"""
import numpy as np
import cps2


def _f(neff, S):
    k0 = 2 * np.pi / S.lam
    kx = np.atleast_1d(np.asarray(neff, dtype=complex)) * k0
    rb, _ = cps2.stack_rt(S.bot_no, S.bot_ne, S.bot_d, kx, k0, 'p')
    rt, _ = cps2.stack_rt(S.top_no, S.top_ne, S.top_d, kx, k0, 'p')
    return (1.0 - rb * rt)[0]


def tm_pole(S, lo, hi, n_im=0.08, rounds=5, n=161):
    """Complex TM mode index: the minimum of |1 - r_b r_t| found by zooming a
    grid.  A secant search fails here because the branch enforcement Im(k_z)>=0
    breaks analyticity, so the root is located rather than iterated to."""
    a, b, c, d = lo, hi, 1e-5, n_im
    for _ in range(rounds):
        re = np.linspace(a, b, n)
        im = np.linspace(c, d, n)
        R, I = np.meshgrid(re, im, indexing='ij')
        Z = (R + 1j * I).ravel()
        k0 = 2 * np.pi / S.lam
        kx = Z * k0
        rb, _ = cps2.stack_rt(S.bot_no, S.bot_ne, S.bot_d, kx, k0, 'p')
        rt, _ = cps2.stack_rt(S.top_no, S.top_ne, S.top_d, kx, k0, 'p')
        F = np.abs(1.0 - rb * rt).reshape(R.shape)
        F[~np.isfinite(F)] = np.inf
        i0, j0 = np.unravel_index(np.argmin(F), F.shape)
        dr, di = (b - a) / (n - 1), (d - c) / (n - 1)
        a, b = re[i0] - 3 * dr, re[i0] + 3 * dr
        c, d = max(1e-7, im[j0] - 3 * di), im[j0] + 3 * di
    return re[i0] + 1j * im[j0]


def kappa(neff, n_o, n_e, lam=550.0):
    """Field decay constant of a mode of index `neff` inside a uniaxial layer,
    returned as the 1/e power-coupling length 1/(2 Re kappa) in nm."""
    k0 = 2 * np.pi / lam
    kap = (n_o / n_e) * np.sqrt((neff * k0) ** 2 - (n_e * k0) ** 2 + 0j)
    kap = kap if kap.real > 0 else -kap
    return kap, 1.0 / (2 * kap.real)
