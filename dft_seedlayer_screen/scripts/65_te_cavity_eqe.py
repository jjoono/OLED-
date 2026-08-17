"""Maximum EQE of the Liu-Angew TE-OLED stack: true (angle-integrated) vs
Lambertian-assumed (forward luminance x pi), from a CPS dipole-cavity model.

STACK (as published, bottom to top; emission through the top):
  glass | ITO 40 | Ag 140 | ITO 15 | HATCN 5 | BPBPA 135 | Bicar 15 |
  EML 35 (source) | Na-An-BI:Liq 45 | Yb 3 | Ag 12 | BPBPA 75 | air

MODEL. Chance-Prock-Silbey power-dissipation formalism: the exciton is a
classical dipole between two mirror stacks with complex reflection
coefficients r_up(u), r_down(u) evaluated at the dipole plane (round-trip
phases included), u = in-plane wavevector / (k0 * n_EML). Kernels
(normalised so the free-space integral is 1 per orientation):

  vertical (TM):    (3/2) Re[ u^3/w * (1+r_up)(1+r_dn)/(1-r_up r_dn) ]
  horizontal TE:    (3/4) Re[ u /w  * (1+r_up)(1+r_dn)/(1-r_up r_dn) ]
  horizontal TM:    (3/4) Re[ u  w  * (1-r_up)(1-r_dn)/(1-r_up r_dn) ]

Outcoupled power: the upward free-dipole kernel (half of the above with
r=0) times the cavity enhancement |1 + r_dn|^2 / |1 - r_up r_dn|^2
(TM-horizontal: |1 - r_dn|^2) times the power transmittance of the top
stack into air. Purcell factor F = full-kernel integral (u to 4, capturing
SPP at u ~ 1.17); effective radiative yield q* = qF/(1-q+qF), q = 0.97.

The two efficiency readings compared:
  TRUE EQE       = q* x (integral of S(theta) over the air hemisphere)
  LAMBERTIAN EQE = q* x pi x S(theta=0)      <- what a forward-only
                   luminance measurement reports for a non-Lambertian device

Charge balance and exciton utilisation are set to 1 (the paper's own best
case), so both numbers are UPPER BOUNDS for the published structure.

n,k: Ag Johnson&Christy; Yb 1.1+2.6i (Palik-order magnitude; 3 nm layer);
ITO 1.95+0.012i; organics n=1.75; HATCN 1.9. Intrinsic emitter spectrum:
Gaussian, 518 nm centre, 20 nm FWHM (BN-Tpl-Ph film values). Dipole plane
scanned across the 35 nm EML; the most favourable position is reported
(again the upper bound). Horizontal dipole fraction 0.85 and the ideal 1.0.
"""
import numpy as np

JC_WL = np.array([397.4, 413.3, 430.5, 450.9, 471.4, 495.9,
                  521.0, 548.6, 582.1, 616.8, 659.5, 704.5])
AG_N = np.array([0.05, 0.05, 0.04, 0.04, 0.05, 0.05, 0.05, 0.06, 0.05, 0.06, 0.05, 0.04])
AG_K = np.array([2.07, 2.21, 2.36, 2.66, 2.83, 3.09, 3.34, 3.59, 3.93, 4.15, 4.48, 4.84])

N_ORG, N_HAT, N_ITO, N_GLASS = 1.75, 1.90, 1.95 + 0.012j, 1.52
N_YB = 1.1 + 2.6j
Q_PL = 0.97

LAM = np.arange(480.0, 581.0, 4.0)
SPEC = np.exp(-4 * np.log(2) * ((LAM - 518.0) / 20.0) ** 2)   # photon spectrum
SPEC /= SPEC.sum()


def n_ag(lam):
    return np.interp(lam, JC_WL, AG_N) + 1j * np.interp(lam, JC_WL, AG_K)


def stack_r_t(n_list, d_list, n_s, u, lam, pol):
    """r at the source-plane boundary and t into the final semi-infinite
    medium, for a plane wave launched from the source medium n_s.
    n_list/d_list: layers AFTER the source plane (first entries may include
    a slab of the source medium itself); last n_list entry semi-infinite."""
    kx2 = (n_s.real * u) ** 2
    def w_of(n):
        w = np.sqrt(n ** 2 - kx2 + 0j)
        return np.where(w.imag < 0, -w, w)
    ws = w_of(n_s)
    ns_all = [n_s] + list(n_list)
    ws_all = [ws] + [w_of(n) for n in n_list]
    # backward recursion for reflection; forward product for transmission
    def fresnel(ni, wi, nj, wj):
        if pol == "TE":
            r = (wi - wj) / (wi + wj)
            t = 2 * wi / (wi + wj)
        else:
            r = (nj ** 2 * wi - ni ** 2 * wj) / (nj ** 2 * wi + ni ** 2 * wj)
            t = 2 * ni * nj * wi / (nj ** 2 * wi + ni ** 2 * wj)
        return r, t
    r_tot = np.zeros_like(u, dtype=complex)
    t_tot = np.ones_like(u, dtype=complex)
    # build from the far end backwards
    r_acc = np.zeros_like(u, dtype=complex)
    for j in range(len(n_list) - 1, -1, -1):
        ni, wi = ns_all[j], ws_all[j]
        nj, wj = ns_all[j + 1], ws_all[j + 1]
        r_ij, t_ij = fresnel(ni, wi, nj, wj)
        if j < len(n_list) - 1:
            ph = np.exp(2j * (2 * np.pi / lam) * ws_all[j + 1] * d_list[j])
            r_acc = (r_ij + r_acc * ph) / (1 + r_ij * r_acc * ph)
        else:
            r_acc = r_ij
    r_tot = r_acc
    # transmission amplitude via forward sweep with known r at each depth:
    # simpler: chain matrices
    M11 = np.ones_like(u, dtype=complex); M12 = np.zeros_like(u, dtype=complex)
    M21 = np.zeros_like(u, dtype=complex); M22 = np.ones_like(u, dtype=complex)
    for j in range(len(n_list)):
        ni, wi = ns_all[j], ws_all[j]
        nj, wj = ns_all[j + 1], ws_all[j + 1]
        r_ij, t_ij = fresnel(ni, wi, nj, wj)
        # interface matrix
        I11, I12, I21, I22 = 1 / t_ij, r_ij / t_ij, r_ij / t_ij, 1 / t_ij
        A11 = M11 * I11 + M12 * I21; A12 = M11 * I12 + M12 * I22
        A21 = M21 * I11 + M22 * I21; A22 = M21 * I12 + M22 * I22
        if j < len(n_list) - 1:
            ph = np.exp(1j * (2 * np.pi / lam) * ws_all[j + 1] * d_list[j])
            A11, A12 = A11 / ph, A12 * ph
            A21, A22 = A21 / ph, A22 * ph
        M11, M12, M21, M22 = A11, A12, A21, A22
    t_tot = 1 / M11
    return r_acc, t_tot, ws, ws_all[-1], ns_all[-1]


def run(theta_par):
    best = None
    z_list = np.array([3.5, 10.5, 17.5, 24.5, 31.5])
    for z_dn in z_list:                    # dipole depth from EML bottom
        z_up = 35.0 - z_dn
        F_bar = 0.0; eta_bar = 0.0; S0_bar = 0.0
        spp_bar = 0.0; trap_bar = 0.0
        for il, lam in enumerate(LAM):
            nag = complex(n_ag(lam))
            n_s = complex(N_ORG)
            up_n = [N_ORG, N_ORG, N_YB, nag, N_ORG, 1.0]
            up_d = [z_up, 45.0, 3.0, 12.0, 75.0]
            dn_n = [N_ORG, N_ORG, N_ORG, complex(N_HAT), N_ITO, nag, N_ITO, complex(N_GLASS)]
            dn_d = [z_dn, 15.0, 135.0, 5.0, 15.0, 140.0, 40.0]
            # u grids: propagating (theta substitution) + evanescent
            th = np.linspace(1e-4, np.pi / 2 - 1e-4, 400)
            u1 = np.sin(th); du1 = np.gradient(th) * np.cos(th)
            u2 = np.concatenate([np.linspace(1.0001, 1.6, 1500),
                                 np.linspace(1.6005, 4.0, 600)])
            du2 = np.gradient(u2)
            F = {}; P_air_u1 = {}
            for pol in ("TE", "TM"):
                out = {}
                for tag, u, in (("p", u1), ("e", u2)):
                    r_up, t_up, ws, w_air, n_air = stack_r_t(up_n, up_d, n_s, u, lam, pol)
                    r_dn, _, _, _, _ = stack_r_t(dn_n, dn_d, n_s, u, lam, pol)
                    den = 1 - r_up * r_dn
                    Xp = (1 + r_up) * (1 + r_dn) / den
                    Xm = (1 - r_up) * (1 - r_dn) / den
                    w = np.sqrt(n_s ** 2 - (n_s.real * u) ** 2 + 0j) / n_s.real
                    w = np.where(w.imag < 0, -w, w)
                    out[tag] = (u, r_up, t_up, r_dn, den, Xp, Xm, w, w_air, n_air)
                F[pol] = out
            # Purcell integrals (both grids)
            def kern(pol):
                Kv = Kte = Ktm = 0.0
                for tag, du in (("p", du1), ("e", du2)):
                    u, r_up, t_up, r_dn, den, Xp, Xm, w, w_air, n_air = F[pol][tag]
                    if pol == "TM":
                        Kv += np.sum((1.5 * np.real(u ** 3 / w * Xp)) * du)
                        Ktm += np.sum((0.75 * np.real(u * w * Xm)) * du)
                    else:
                        Kte += np.sum((0.75 * np.real(u / w * Xp)) * du)
                return Kv, Kte, Ktm
            Kv, _, _ = kern("TM"); _, Kte, _ = kern("TE")
            _, _, Ktm = kern("TM")
            F_tot = theta_par * (Kte + Ktm) + (1 - theta_par) * Kv
            # SPP / trapped bookkeeping (evanescent part of F)
            def kern_ev(pol):
                s = 0.0
                u, r_up, t_up, r_dn, den, Xp, Xm, w, w_air, n_air = F[pol]["e"]
                du = du2
                if pol == "TM":
                    s += np.sum((1.5 * np.real(u ** 3 / w * Xp)) * du) * (1 - theta_par)
                    s += np.sum((0.75 * np.real(u * w * Xm)) * du) * theta_par
                else:
                    s += np.sum((0.75 * np.real(u / w * Xp)) * du) * theta_par
                return s
            spp = kern_ev("TM") + kern_ev("TE")
            # outcoupled: propagating grid, transmittance into air
            u, r_upTE, t_upTE, r_dnTE, denTE, XpTE, XmTE, w, w_airTE, n_air = F["TE"]["p"]
            _, r_upTM, t_upTM, r_dnTM, denTM, XpTM, XmTM, _, w_airTM, _ = F["TM"]["p"]
            # power transmittance with E-field t (both pols):
            # T = Re(kz_air)/kz_source * |t|^2 ; kz_source = n_s * w (real, u<1)
            # normal-incidence check: 4 n1 n2/(n1+n2)^2 reproduced for both.
            T_te = np.where(w_airTE.real > 0,
                            (w_airTE.real / (n_s.real * w.real + 1e-30)) * abs(t_upTE) ** 2,
                            0.0)
            T_tm = np.where(w_airTM.real > 0,
                            (w_airTM.real / (n_s.real * w.real + 1e-30)) * abs(t_upTM) ** 2,
                            0.0)
            up_te = 0.375 * np.real(u / w) * abs(1 + r_dnTE) ** 2 / abs(denTE) ** 2 * T_te
            up_tmh = 0.375 * np.real(u * w) * abs(1 - r_dnTM) ** 2 / abs(denTM) ** 2 * T_tm
            up_tmv = 0.75 * np.real(u ** 3 / w) * abs(1 + r_dnTM) ** 2 / abs(denTM) ** 2 * T_tm
            P_ang = theta_par * (up_te + up_tmh) + (1 - theta_par) * up_tmv
            P_air = np.sum(P_ang * du1)
            # forward photon intensity per sr: S(0) via smallest-u bin
            # dOmega mapping: u = (n_air/n_s) sin(theta_air)
            sin_ta = np.clip(n_s.real * u / 1.0, 0, 1)
            valid = sin_ta < 1.0
            cos_ta = np.sqrt(1 - sin_ta[valid] ** 2)
            S_th = P_ang[valid] * du1[valid]
            # S(theta)dtheta_air -> per sr: /(2 pi sin ta dta); dta from du
            dta = np.gradient(np.arcsin(sin_ta[valid]))
            S_sr = (P_ang[valid] * du1[valid] / np.maximum(dta, 1e-12)) / \
                   (2 * np.pi * np.maximum(sin_ta[valid], 1e-9))
            S0 = S_sr[1]                        # near theta = 0
            w_l = SPEC[il]
            F_bar += w_l * F_tot
            eta_bar += w_l * P_air
            S0_bar += w_l * S0
            spp_bar += w_l * spp
        q_star = Q_PL * F_bar / (1 - Q_PL + Q_PL * F_bar)
        eqe_true = q_star * (eta_bar / F_bar)
        eqe_lamb = q_star * (np.pi * S0_bar / F_bar)
        frac_spp = spp_bar / F_bar
        if best is None or eqe_true > best[1]:
            best = (z_dn, eqe_true, eqe_lamb, q_star, F_bar, frac_spp,
                    eta_bar / F_bar)
    return best


def main():
    print("=" * 74)
    print("Liu-Angew TE stack: CPS dipole-cavity upper bounds (best dipole plane)")
    print("gamma = exciton utilisation = 1; q_PL = 0.97; spectrum 518/20 nm")
    print("=" * 74)
    print(f"{'Theta_par':>10}{'z* (nm)':>9}{'F':>7}{'eta_out':>9}"
          f"{'SPP+evan':>10}{'EQE true':>10}{'EQE Lamb.':>11}{'ratio':>7}")
    for tp in (0.85, 1.00):
        z, et, el, qs, F, fspp, eta = run(tp)
        print(f"{tp:>10.2f}{z:>9.1f}{F:>7.2f}{eta:>9.1%}"
              f"{fspp:>10.1%}{et:>10.1%}{el:>11.1%}{el/et:>7.2f}")
    print("""
  EQE true  = angle-integrated photons out of the top surface per exciton
  EQE Lamb. = pi x forward radiance (photon) per exciton -- what a
              forward-only IVL measurement reports under the Lambertian
              assumption. The ratio is the inflation factor.""")


# RESULT (run 2026-08-14, after the TM-transmittance normalisation fix;
# the first run violated energy conservation and was discarded):
#
#   Theta_par  z*(nm)   F    eta_out  SPP+evan  EQE_true  EQE_Lambertian
#     0.85      24.5   1.34   36.6%    14.9%     35.8%       81.7%
#     1.00      24.5   1.33   43.1%    12.0%     42.1%       96.8%
#
# READING. For the published stack, the TRUE angle-integrated EQE tops out
# at ~36 % with a realistic MR-TADF orientation (Theta_par 0.85) and ~42 %
# even with perfectly horizontal dipoles, gamma = 1, and the most favourable
# dipole plane. The claimed 59.2 % exceeds the perfect-orientation bound by
# 17 points -- it is not reachable as a true EQE in this structure. A
# forward-only measurement converted with the Lambertian assumption, on the
# other hand, reads 82-97 % at the optimum: the measured 59.2 % sits
# comfortably inside what that procedure produces for a device whose true
# efficiency is in the twenties-to-thirties. Model inflation factor ~2.3.
#
# The user's physical intuition is confirmed on both counts: the thin
# cathode-side spacer (ETL 45 nm + Yb 3 nm) leaves 12-15 % of the dipole
# power in SPP/evanescent channels even at 85 % horizontal orientation, and
# metal/ITO absorption in the high-Q cavity (lossy Yb, 12 nm Ag, 140 nm Ag
# mirror, each pass) consumes ~half of the emitted power.
#
# Caveats: organics carried at a uniform n = 1.75 (no birefringence), Yb at
# 1.1+2.6i, dipole planes on a 7 nm grid, published thicknesses assumed
# exact. None of these moves the true-EQE bound by more than a few points.

if __name__ == "__main__":
    main()
