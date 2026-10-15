"""CompleteEASE Gen-Osc oscillators with ANALYTIC eps1.

Analytic (not truncated-KK) matters here: Einf is written into the .mod, and a
truncated KK integral would silently dump the missing high-energy weight into
Einf, so my Einf would not mean what CompleteEASE's Einf means.
"""
import numpy as np
from scipy.special import dawsn

HB_EVFS = 0.6582119569
HB_EVS  = 6.582119569e-16
EPS0    = 8.8541878128e-12


def drude_rt(E, rho, tau):
    """CompleteEASE Drude(RT): rho in Ohm.cm, tau in fs."""
    Br = HB_EVFS / tau
    A  = HB_EVS**2 / (EPS0 * (rho / 100.0) * (tau * 1e-15))       # eV^2
    den = E**2 + Br**2
    return -A / den, A * Br / (E * den)


def gaussian(E, A, Br, En):
    """CompleteEASE Gaussian; eps1 via the Dawson function (exact KK partner)."""
    s = Br / (2.0 * np.sqrt(np.log(2.0)))
    e2 = A * (np.exp(-((E - En) / s)**2) - np.exp(-((E + En) / s)**2))
    e1 = (2.0 * A / np.sqrt(np.pi)) * (dawsn((E + En) / s) - dawsn((E - En) / s))
    return e1, e2


def tauc_lorentz(E, A, C, E0, Eg):
    """Jellison-Modine Tauc-Lorentz, analytic eps1 (APL 69, 371 + erratum)."""
    E = np.asarray(E, float)
    e2 = np.where(E > Eg,
                  (A * E0 * C * (E - Eg)**2) /
                  ((E**2 - E0**2)**2 + C**2 * E**2) / np.maximum(E, 1e-12),
                  0.0)

    aln  = (Eg**2 - E0**2) * E**2 + Eg**2 * C**2 - E0**2 * (E0**2 + 3.0 * Eg**2)
    aatn = (E**2 - E0**2) * (E0**2 + Eg**2) + Eg**2 * C**2
    alpha = np.sqrt(max(4.0 * E0**2 - C**2, 1e-30))
    gam2  = E0**2 - C**2 / 2.0
    zeta4 = (E**2 - gam2)**2 + alpha**2 * C**2 / 4.0

    t1 = (A * C * aln) / (2.0 * np.pi * zeta4 * alpha * E0) * \
         np.log((E0**2 + Eg**2 + alpha * Eg) / (E0**2 + Eg**2 - alpha * Eg))
    t2 = -(A * aatn) / (np.pi * zeta4 * E0) * \
         (np.pi - np.arctan((2.0 * Eg + alpha) / C) + np.arctan((alpha - 2.0 * Eg) / C))
    t3 = (4.0 * A * E0 * Eg * (E**2 - gam2)) / (np.pi * zeta4 * alpha) * \
         (np.pi / 2.0 + np.arctan(2.0 * (gam2 - Eg**2) / (alpha * C)))
    d  = np.abs(E - Eg)
    t4 = -(A * E0 * C * (E**2 + Eg**2)) / (np.pi * zeta4 * np.maximum(E, 1e-12)) * \
         np.log(np.maximum(d, 1e-12) / (E + Eg))
    t5 = (2.0 * A * E0 * C * Eg) / (np.pi * zeta4) * \
         np.log(np.maximum(d, 1e-12) * (E + Eg) /
                np.sqrt((E0**2 - Eg**2)**2 + Eg**2 * C**2))
    return t1 + t2 + t3 + t4 + t5, e2


def _kk_ref(E_out, e2_fun, emax=400.0, n=400001):
    """brute-force principal-value KK on a very wide fine grid, for validation"""
    Eg_ = np.linspace(1e-4, emax, n)
    e2 = e2_fun(Eg_)
    dE = Eg_[1] - Eg_[0]
    out = np.empty_like(E_out)
    for i, Ei in enumerate(E_out):
        den = Eg_**2 - Ei**2
        j = np.argmin(np.abs(den))
        num = Eg_ * e2
        g = Ei * np.interp(Ei, Eg_, e2)
        f = (num - g) / np.where(np.abs(den) < 1e-30, np.nan, den)
        f[j] = 0.0
        a, b = Eg_[0], Eg_[-1]
        Ian = (1.0 / (2 * Ei)) * (np.log(abs((b - Ei) / (b + Ei))) -
                                  np.log(abs((a - Ei) / (a + Ei))))
        out[i] = (2 / np.pi) * (np.nansum(f) * dE + g * Ian)
    return out


if __name__ == '__main__':
    Et = np.array([0.8, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 6.45])
    print('Tauc-Lorentz  A=165.5 C=1.055 E0=4.0 Eg=3.337')
    p = (165.546636633031, 1.0550595369907967, 4.0, 3.3367948707895163)
    e1a, _ = tauc_lorentz(Et, *p)
    e1n = _kk_ref(Et, lambda x: tauc_lorentz(x, *p)[1])
    for E, a, b in zip(Et, e1a, e1n):
        print('  E=%5.2f  analytic %10.5f   numeric %10.5f   diff %+.2e' % (E, a, b, a - b))
    print('  max |diff| = %.2e' % np.abs(e1a - e1n).max())

    print('\nGaussian  A=0.072 Br=0.825 En=1.911')
    q = (0.07201885654941534, 0.8249493158431357, 1.911146995799902)
    g1a, _ = gaussian(Et, *q)
    g1n = _kk_ref(Et, lambda x: gaussian(x, *q)[1], emax=60.0, n=300001)
    print('  max |diff| = %.2e' % np.abs(g1a - g1n).max())
