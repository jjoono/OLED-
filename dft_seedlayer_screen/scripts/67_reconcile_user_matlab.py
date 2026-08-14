"""Reconciling the user's MATLAB result (60 % @ optimised red TE stack)
with the 44 % green bound: replicate their structure in the same CPS model
and switch the differing assumptions one at a time.

USER'S STRUCTURE (their MATLAB, emission through the thin Ag + capping):
  air | Ag 100 (cathode mirror) | ETL 228.8 | EML 25 (z0 12.5) | TCTA 10 |
  HTL 252.4 | Ag 12 | cap 129 (n = 2.3!) | air
  red phosphor (~600 nm), Theta_par 0.95, eta_rad 1, no Yb.
  Their reported optimum: "Totalpower" 0.5994 with lambda restricted to
  590-600 nm -- their own header warns this window renormalisation
  distorts absolute EQE.

LADDER (each row changes ONE thing from the row above):
  A  their stack, lambda window 590-600 nm only     <- should reproduce ~0.60
  B  full red spectrum (600 nm, 45 nm FWHM)         <- window effect
  C  + Yb 3 nm under the thin Ag                    <- wetting-layer cost
  D  green (518 nm / 20 nm) with green-optimised    <- Ag loss at 520
     re-optimised spacers (coarse rescan)
  E  D with Theta_par 0.85                          <- orientation
"""
import importlib.util
import os
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "te65", os.path.join(_HERE, "65_te_cavity_eqe.py"))
te65 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(te65)
stack_r_t, n_ag = te65.stack_r_t, te65.n_ag

N_ORG = 1.75
N_CAP = 2.3


def compute(etl, htl, cap, lam_grid, spec_w, theta_par, yb_nm=0.0,
            z0=12.5, eml=25.0, spacer=10.0):
    """User's architecture: emission side = TCTA spacer/HTL/(Yb)/Ag12/cap/air;
    back side = ETL/Ag100/air. Source in EML at z0 from the back side."""
    n_s = complex(N_ORG)
    F_bar = P_bar = S0_bar = 0.0
    th = np.linspace(1e-4, np.pi / 2 - 1e-4, 350)
    u1 = np.sin(th); du1 = np.gradient(th) * np.cos(th)
    u2 = np.concatenate([np.linspace(1.0001, 1.6, 900),
                         np.linspace(1.6005, 3.5, 300)])
    du2 = np.gradient(u2)
    for il, lam in enumerate(lam_grid):
        nag = complex(n_ag(lam))
        up_n = [N_ORG, N_ORG, N_ORG]
        up_d = [eml - z0, spacer, htl]
        if yb_nm > 0:
            up_n += [te65.N_YB]; up_d += [yb_nm]
        up_n += [nag, complex(N_CAP), 1.0]
        up_d += [12.0, cap]
        dn_n = [N_ORG, N_ORG, nag, 1.0]
        dn_d = [z0, etl, 100.0]
        F_tot = P_air = S0 = 0.0
        for pol in ("TE", "TM"):
            for u, du, tag in ((u1, du1, "p"), (u2, du2, "e")):
                r_up, t_up, ws, w_air, n_air = stack_r_t(up_n, up_d, n_s, u, lam, pol)
                r_dn, _, _, _, _ = stack_r_t(dn_n, dn_d, n_s, u, lam, pol)
                den = 1 - r_up * r_dn
                Xp = (1 + r_up) * (1 + r_dn) / den
                Xm = (1 - r_up) * (1 - r_dn) / den
                w = np.sqrt(n_s ** 2 - (n_s.real * u) ** 2 + 0j) / n_s.real
                w = np.where(w.imag < 0, -w, w)
                if pol == "TM":
                    F_tot += (1 - theta_par) * np.sum(1.5 * np.real(u ** 3 / w * Xp) * du)
                    F_tot += theta_par * np.sum(0.75 * np.real(u * w * Xm) * du)
                else:
                    F_tot += theta_par * np.sum(0.75 * np.real(u / w * Xp) * du)
                if tag == "p":
                    T = np.where(w_air.real > 0,
                                 (w_air.real / (n_s.real * w.real + 1e-30))
                                 * abs(t_up) ** 2, 0.0)
                    if pol == "TE":
                        pa = theta_par * 0.375 * np.real(u / w) \
                            * abs(1 + r_dn) ** 2 / abs(den) ** 2 * T
                    else:
                        pa = (theta_par * 0.375 * np.real(u * w)
                              * abs(1 - r_dn) ** 2 / abs(den) ** 2
                              + (1 - theta_par) * 0.75 * np.real(u ** 3 / w)
                              * abs(1 + r_dn) ** 2 / abs(den) ** 2) * T
                    P_air += np.sum(pa * du)
                    sin_ta = np.clip(n_s.real * u, 0, 1)
                    v = sin_ta < 1.0
                    dta = np.gradient(np.arcsin(sin_ta[v]))
                    S_sr = (pa[v] * du[v] / np.maximum(dta, 1e-12)) / \
                           (2 * np.pi * np.maximum(sin_ta[v], 1e-9))
                    S0 += S_sr[1]
        wl = spec_w[il]
        F_bar += wl * F_tot; P_bar += wl * P_air; S0_bar += wl * S0
    q = 1.0  # eta_rad = 1 as in the user's run
    return P_bar / F_bar, np.pi * S0_bar / F_bar, F_bar


def gauss(lg, c, fw):
    s = np.exp(-4 * np.log(2) * ((lg - c) / fw) ** 2)
    return s / s.sum()


def main():
    ETL, HTL, CAP = 228.85, 252.45, 129.04

    rows = []
    # A: their lambda window
    lgA = np.arange(590.0, 601.0, 2.0); wA = gauss(lgA, 598.0, 25.0)
    a = compute(ETL, HTL, CAP, lgA, wA, 0.95)
    rows.append(("A window 590-600, red, no Yb, Tpar .95", a))
    # B: full red spectrum
    lgB = np.arange(540.0, 701.0, 5.0); wB = gauss(lgB, 600.0, 45.0)
    b = compute(ETL, HTL, CAP, lgB, wB, 0.95)
    rows.append(("B full red spectrum (600/45)", b))
    # C: + Yb 3 nm
    c = compute(ETL, HTL, CAP, lgB, wB, 0.95, yb_nm=3.0)
    rows.append(("C + Yb 3 nm wetting layer", c))
    # D: green, re-optimised spacers (coarse)
    lgD = te65.LAM; wD = te65.SPEC
    bestD = None
    for etl in np.arange(40, 261, 20):
        for htl in np.arange(40, 261, 20):
            for cap in np.arange(40, 161, 20):
                r = compute(etl, htl, cap, lgD[::2], wD[::2] / wD[::2].sum(),
                            0.95, yb_nm=3.0)
                if bestD is None or r[0] > bestD[0][0]:
                    bestD = (r, etl, htl, cap)
    (rD, etlD, htlD, capD) = bestD
    rD = compute(etlD, htlD, capD, lgD, wD, 0.95, yb_nm=3.0)
    rows.append((f"D green 518/20, Yb, reopt ({etlD:.0f}/{htlD:.0f}/{capD:.0f})", rD))
    # E: orientation 0.85
    e = compute(etlD, htlD, capD, lgD, wD, 0.85, yb_nm=3.0)
    rows.append(("E same, Theta_par 0.85", e))

    print("=" * 74)
    print("Assumption ladder: user's 60 % result -> Liu-paper conditions")
    print("(eta_rad = 1 everywhere, cap n = 2.3, user's optimised thicknesses)")
    print("=" * 74)
    print(f"{'step':<46}{'true EQE':>10}{'Lamb.read':>11}{'F':>6}")
    for tag, (eta, lamb, F) in rows:
        print(f"{tag:<46}{eta:>10.1%}{lamb:>11.1%}{F:>6.2f}")


if __name__ == "__main__":
    main()
