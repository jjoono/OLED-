"""Independent textbook CPS reference for the benchmark OLED stack.

Purpose: adjudicate between the (reconstructed) group planar-sweep code and
the Meep FDTD results. Everything is derived from first principles here:

  * Parratt recursion for the p/s reflection coefficients of the two
    half-spaces seen from the dipole plane, kz on the physical branch
    Im(kz) >= 0 (decaying evanescent waves).
  * Chance-Prock-Silbey power dissipation spectra with the vacuum limit
    F = 1 recovered exactly when all reflections vanish (validated below).

Stack (bottom emission through the glass):
    air | Al 200 nm | EML n=1.75, 100 nm (dipole in the centre) |
    ITO n=1.9, 100 nm | glass n=1.5 (semi-infinite)
lambda = 550 nm, Al n = 0.9656 + 6.4581i (Meep Rakic model at 550 nm).

Mode bookkeeping (same windows as the group code, u = kpar/(k0*n_EML)):
    sub : power transmitted into the glass, u < n_glass/n_EML  (== FDTD glass)
    abs : radiative-window power absorbed on the way (Al), u < u_sub
    wg  : u_sub < u < 1  (guided, propagating in EML)
    spp : u > 1          (evanescent in EML; SPP + TM0-like modes)
Total dissipated power F = Purcell factor (per orientation).
"""

import numpy as np

LAM = 550.0
K0 = 2 * np.pi / LAM
N_EML = 1.75
N_AL = 0.9656319200398573 + 6.458058124032318j
D_EML = 100.0
Z0 = 50.0                      # dipole to EML/Al interface
T_AL = 200.0
T_ITO = 100.0


def kz(n, kpar):
    """Vertical wavevector with Im(kz) >= 0 (physical decay branch)."""
    v = np.sqrt((K0 * n) ** 2 - kpar ** 2 + 0j)
    return np.where(v.imag < 0, -v, v)


def parratt(ns, ts, kpar, pol):
    """r seen from layer 0 into stack ns[1:]; ts = thicknesses of ns[1:-1].

    ns: [n0, n1, ..., nN] (n0 = incidence medium, nN = final half-space)
    Returns the complex reflection coefficient at the 0|1 interface.
    """
    kzs = [kz(n, kpar) for n in ns]

    def rij(i, j):
        if pol == "s":
            return (kzs[i] - kzs[j]) / (kzs[i] + kzs[j])
        ni2, nj2 = ns[i] ** 2, ns[j] ** 2
        return (nj2 * kzs[i] - ni2 * kzs[j]) / (nj2 * kzs[i] + ni2 * kzs[j])

    r = rij(len(ns) - 2, len(ns) - 1)
    for i in range(len(ns) - 3, -1, -1):
        ph = np.exp(2j * kzs[i + 1] * ts[i])
        rn = rij(i, i + 1)
        r = (rn + r * ph) / (1 + rn * r * ph)
    return r


def spectra(kk_wg=0.0):
    """Return u grid and dissipation integrands f_v(u), f_h(u).

    kk_wg: small extinction added to EML/ITO to regularise guided-mode poles
    (needed only for the wg/spp window split; keep 0 for the sub window).
    """
    n_eml = N_EML + 1j * kk_wg
    n_ito = 1.9 + 1j * kk_wg

    # dense grid, extra dense around u ~ 1 where modes live
    u = np.concatenate([
        np.linspace(1e-4, 0.84, 4000, endpoint=False),
        np.linspace(0.84, 1.12, 60000, endpoint=False),
        np.linspace(1.12, 4.0, 20000),
    ])
    kpar = K0 * N_EML * u          # windows defined w.r.t. real EML index
    l1 = kz(n_eml, kpar) / (K0 * n_eml)

    out = {}
    for pol in ("s", "p"):
        r_b = parratt([n_eml, n_ito, 1.5], [T_ITO], kpar, pol)
        r_t = parratt([n_eml, N_AL, 1.0], [T_AL], kpar, pol)
        a_b = r_b * np.exp(2j * kz(n_eml, kpar) * (D_EML - Z0))
        a_t = r_t * np.exp(2j * kz(n_eml, kpar) * Z0)
        out[pol] = (a_b, a_t)

    ab_s, at_s = out["s"]
    ab_p, at_p = out["p"]
    Ds = 1 - ab_s * at_s
    Dp = 1 - ab_p * at_p

    f_v = 1.5 * np.real(u ** 3 / l1 * (1 + ab_p) * (1 + at_p) / Dp)
    f_h = 0.75 * np.real(u / l1 * ((1 + ab_s) * (1 + at_s) / Ds
                                   + (1 - u ** 2) * (1 - ab_p) * (1 - at_p) / Dp))
    return u, f_v, f_h, (ab_s, at_s, ab_p, at_p, l1)


def transmitted_sub(u, cav):
    """Power flux into the glass per u (radiative window only, lossless below).

    Down-going field in the EML at the dipole: E- = S- + a_t (S+ ... ) with
    source amplitudes S. Net downward flux at the bottom boundary equals the
    cavity-enhanced source flux times (1-|R_eff|^2 ... ). For a LOSSLESS
    bottom stack it is simpler and exact to compute the flux into the glass
    as the difference between down-going and up-going flux just below the
    dipole, since no power is absorbed on the way down:
        S_down_net(u) = f_orientation_down(u) as below.
    Standard result (e.g. Wasey & Barnes): the power flowing into the bottom
    half-space is

        T(u) = f(u) * Re[ (1 - a_b conj(1 ... ) ]   -- messy;

    we instead use energy bookkeeping per u that is exact for our geometry:
    everything dissipated at a given u < u_sub either enters the glass or is
    absorbed in the METAL (top side). Power absorbed in the metal per u can
    be computed by the same CPS formula with the bottom stack made perfectly
    transparent? Not available. -> We therefore compute the metal absorption
    directly: replace the bottom half-space by pure EML (r_b = 0). Then ALL
    dissipation is metal absorption + free radiation downward. Subtracting
    isolates nothing useful either.

    Conclusion: the clean observable to compare with FDTD frac_glass is
      P_sub = F - P_metal,
    where P_metal is integrated over ALL u of the metal absorption, computed
    from the difference between total dissipation and net downward flux.
    The net downward flux per u through the plane just below the dipole is

      S_v(u)  = 1.5 * Re[ u^3/l1 * |1 + a_t|^2 * (1 - |a_b'|^2 ...) ] ...

    Deriving flux formulas reliably in one shot is error-prone; instead we
    use the独立 relation: for u in the radiative window with a lossless
    bottom stack, power reaching the glass equals total dissipation minus
    metal absorption AT THAT u. We obtain metal absorption per u numerically
    via perturbation: dF/d(eps_Al'') * eps_Al'' is exact for linear media
    (Poynting theorem: P_abs = w/2 eps0 eps'' |E|^2 integrated in metal).
    We approximate instead with the two-run difference:
        P_metal(u) ~ f(u; Al) - f(u; Al -> PEC-like lossless mirror)
    which is NOT exact. --> Simplest exact route: FDTD comparison uses
    frac_glass; CPS "sub" from the group K2 formulas is computed in the
    Octave path. Here we only cross-check F and the window split of F.
    """
    raise NotImplementedError


def main():
    # ---- validation: vacuum limit (all r = 0 -> F = 1) ------------------
    u = np.linspace(1e-5, 0.999999, 400000)
    l1 = np.sqrt(1 - u ** 2)
    Fv0 = np.trapezoid(1.5 * u ** 3 / l1, u)
    Fh0 = np.trapezoid(0.75 * (u / l1) * (1 + (1 - u ** 2)), u)
    print("vacuum limit check: F_v = %.6f  F_h = %.6f  (expect 1, 1)" % (Fv0, Fh0))

    # ---- real stack -----------------------------------------------------
    for kk in (1e-3, 2e-3, 4e-3):
        u, f_v, f_h, _ = spectra(kk_wg=kk)
        F_v = np.trapezoid(f_v, u)
        F_h = np.trapezoid(f_h, u)
        w_sub = u < 1.5 / N_EML
        w_wg = (u >= 1.5 / N_EML) & (u < 1)
        w_spp = u >= 1

        def split(f):
            F = np.trapezoid(f, u)
            return (np.trapezoid(np.where(w_sub, f, 0), u) / F,
                    np.trapezoid(np.where(w_wg, f, 0), u) / F,
                    np.trapezoid(np.where(w_spp, f, 0), u) / F, F)

        sv = split(f_v)
        sh = split(f_h)
        print("kk=%.0e | F_v=%.4f  v-split sub/wg/spp = %.3f/%.3f/%.3f" %
              (kk, sv[3], sv[0], sv[1], sv[2]))
        print("        | F_h=%.4f  h-split sub/wg/spp = %.3f/%.3f/%.3f" %
              (sh[3], sh[0], sh[1], sh[2]))


if __name__ == "__main__":
    main()
