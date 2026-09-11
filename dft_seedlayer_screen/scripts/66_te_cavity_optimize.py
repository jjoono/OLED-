"""What if the Liu-Angew TE stack could use ETL / HTL / capping up to 300 nm?

Same CPS dipole-cavity model as scripts/65 (imported from there), but the
three optical spacers become free variables:

    HTL (BPBPA, published 135 nm)      : 40 - 300 nm
    ETL (Na-An-BI:Liq, published 45 nm): 30 - 300 nm
    capping (BPBPA, published 75 nm)   : 0 - 300 nm

Everything else fixed as published (HATCN 5, Bicar 15, EML 35, Yb 3, Ag 12,
mirror side unchanged). Two-stage search: coarse grid at a single mid-EML
dipole plane and a thinned spectrum, then refinement around the coarse
optimum with the full spectrum and a 5-plane dipole scan (same settings as
scripts/65, so numbers are directly comparable).

Physical expectation being tested: a thicker cathode-side spacer pulls the
emitter away from the Yb/Ag surface (kills SPP), and moving to a higher
cavity order with an optimised capping layer can trade a little mode
confinement for less metal absorption per photon.

Electrical caveat, stated up front: 300 nm transport layers cost drive
voltage, so power efficiency would NOT scale with the EQE found here even
where the optics improve. This is an optics-only bound.
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
N_ORG, N_HAT, N_ITO, N_GLASS, N_YB, Q_PL = (te65.N_ORG, te65.N_HAT,
                                            te65.N_ITO, te65.N_GLASS,
                                            te65.N_YB, te65.Q_PL)
THETA_PAR = 0.85


def point(htl, etl, cap, z_dn, lam_grid, spec_w, n_u_ev=900):
    """Return (F, P_air, S0) spectrum-averaged for one geometry."""
    z_up = 35.0 - z_dn
    F_bar = P_bar = S0_bar = 0.0
    th = np.linspace(1e-4, np.pi / 2 - 1e-4, 300)
    u1 = np.sin(th); du1 = np.gradient(th) * np.cos(th)
    u2 = np.concatenate([np.linspace(1.0001, 1.6, int(n_u_ev * 0.75)),
                         np.linspace(1.6005, 3.5, int(n_u_ev * 0.25))])
    du2 = np.gradient(u2)
    n_s = complex(N_ORG)
    for il, lam in enumerate(lam_grid):
        nag = complex(n_ag(lam))
        up_n = [N_ORG, N_ORG, N_YB, nag, N_ORG, 1.0]
        up_d = [z_up, etl, 3.0, 12.0, cap]
        dn_n = [N_ORG, N_ORG, N_ORG, complex(N_HAT), N_ITO, nag, N_ITO,
                complex(N_GLASS)]
        dn_d = [z_dn, 15.0, htl, 5.0, 15.0, 140.0, 40.0]
        F_tot = 0.0; P_air = 0.0; S0 = 0.0
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
                    F_tot += (1 - THETA_PAR) * np.sum(1.5 * np.real(u ** 3 / w * Xp) * du)
                    F_tot += THETA_PAR * np.sum(0.75 * np.real(u * w * Xm) * du)
                else:
                    F_tot += THETA_PAR * np.sum(0.75 * np.real(u / w * Xp) * du)
                if tag == "p":
                    T = np.where(w_air.real > 0,
                                 (w_air.real / (n_s.real * w.real + 1e-30))
                                 * abs(t_up) ** 2, 0.0)
                    if pol == "TE":
                        pa = THETA_PAR * 0.375 * np.real(u / w) \
                            * abs(1 + r_dn) ** 2 / abs(den) ** 2 * T
                    else:
                        pa = (THETA_PAR * 0.375 * np.real(u * w)
                              * abs(1 - r_dn) ** 2 / abs(den) ** 2
                              + (1 - THETA_PAR) * 0.75 * np.real(u ** 3 / w)
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
    return F_bar, P_bar, S0_bar


def eqe_from(F, P, S0):
    q = Q_PL * F / (1 - Q_PL + Q_PL * F)
    return q * P / F, q * np.pi * S0 / F


def main():
    lam_c = np.arange(484.0, 581.0, 8.0)
    w_c = np.exp(-4 * np.log(2) * ((lam_c - 518.0) / 20.0) ** 2); w_c /= w_c.sum()
    lam_f = te65.LAM; w_f = te65.SPEC

    print("Stage 1: coarse grid (dipole mid-EML, thinned spectrum)")
    best = None
    for htl in np.arange(40, 301, 20):
        for etl in np.arange(30, 301, 30):
            for cap in np.arange(0, 301, 40):
                F, P, S0 = point(htl, etl, cap, 17.5, lam_c, w_c, n_u_ev=500)
                et, _ = eqe_from(F, P, S0)
                if best is None or et > best[0]:
                    best = (et, htl, etl, cap)
    et, htl0, etl0, cap0 = best
    print(f"  coarse optimum: HTL {htl0:.0f}, ETL {etl0:.0f}, cap {cap0:.0f}"
          f"  -> true EQE ~{et:.1%}")

    print("Stage 2: refine +-25 nm, full spectrum, 5 dipole planes")
    best2 = None
    for htl in np.arange(max(40, htl0 - 25), htl0 + 26, 10):
        for etl in np.arange(max(30, etl0 - 25), etl0 + 26, 10):
            for cap in np.arange(max(0, cap0 - 30), cap0 + 31, 15):
                for z in (3.5, 10.5, 17.5, 24.5, 31.5):
                    F, P, S0 = point(htl, etl, cap, z, lam_f, w_f)
                    et, el = eqe_from(F, P, S0)
                    if best2 is None or et > best2[0]:
                        best2 = (et, el, htl, etl, cap, z, F, P / F)
    et, el, htl, etl, cap, z, F, eta = best2
    print("=" * 70)
    print(f"OPTIMISED (Theta_par = {THETA_PAR}):")
    print(f"  HTL {htl:.0f} nm | ETL {etl:.0f} nm | cap {cap:.0f} nm | "
          f"dipole z {z:.1f} nm | F {F:.2f}")
    print(f"  eta_out {eta:.1%}   TRUE EQE {et:.1%}   Lambertian-read {el:.1%}")
    print(f"  (published geometry, same model: true 35.8%, Lambertian 81.7%)")


# RESULT (run 2026-08-14, 3 min):
#
#   coarse optimum:  HTL 140 / ETL 60 / cap 80   (true EQE ~42 %)
#   refined optimum: HTL 145 / ETL 55 / cap 80, dipole z 24.5, F 1.52
#       eta_out 45.0 %   TRUE EQE 44.1 %   Lambertian-read 26.3 %
#   published stack, same model: true 35.8 %, Lambertian-read 81.7 %
#
# THREE READINGS.
#
# 1. Full 300 nm freedom on all three spacers buys ~8 points of true EQE
#    (35.8 -> 44.1 %). The optimiser REFUSES the thick-layer route: the
#    optimum sits a few nm from the published design (145/55/80 vs
#    135/45/75). Moving to a higher cavity order with 100-300 nm spacers
#    loses more to per-pass metal absorption than it gains from SPP
#    suppression, so "thicker ETL kills SPP" does not pay in this
#    two-Ag-mirror architecture. The published thicknesses are essentially
#    optimal; the claimed NUMBER, not the design, is the anomaly. Even the
#    44.1 % bound (and ~50 % at a perfect Theta_par = 1) stays far below
#    the claimed 59.2 %.
#
# 2. The metric inversion is the striking part: at the true-EQE optimum
#    the forward-Lambertian reading COLLAPSES to 26.3 % (sub-Lambertian,
#    wide angular lobe), while the published geometry reads 81.7 % on the
#    same instrument logic while genuinely delivering 35.8 %. A forward-
#    only IVL metric therefore actively mis-ranks TE designs: it rewards
#    piling photons at normal incidence and penalises the geometry that
#    actually emits the most light. Chasing the IVL number and chasing
#    real efficiency lead to DIFFERENT devices.
#
# 3. For the editor letter: no spacer configuration reachable with these
#    materials supports 59.2 % as an angle-integrated EQE (optics-only
#    upper bounds, gamma = 1). Electrical caveat unchanged: 300 nm
#    transport layers would also cost drive voltage, so even the 44 %
#    optimum would sacrifice power efficiency.

if __name__ == "__main__":
    main()
