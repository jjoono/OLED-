"""Alloy-doped ultrathin Ag electrode: does Mg/Yb co-deposition ever WIN on
the absorption-first criterion (min visible A, subject to Rs <= 100 Ohm/sq)?

THE TRADE-OFF UNDER TEST. Co-depositing a few percent Mg or Yb with Ag is
known to smooth ultrathin films (suppressed island growth -> less void ->
less plasmonic/void absorption), but the alloy itself is optically lossier
than pure Ag (the dopant raises k and kills silver's uniquely low Im(eps)).
The user predicts a U-shaped A(y): absorption first drops as morphology
improves, then rises as alloy loss takes over. This script quantifies where
(and whether) the minimum beats the undoped film.

MODEL (every constant sourced or flagged as an assumption):

1. OPTICS. Alloy dielectric = 3D Bruggeman mix of bulk Ag (Johnson & Christy
   1972 n,k arrays, copied from scripts/62_eta_from_EMA.py) with the dopant
   metal at volume fraction y. Dopant optical constants are APPROXIMATIONS,
   hardcoded from literature bulk values:
     - Mg: n ~ 0.4, k ~ 4.4 at 520 nm (Mg is a decent free-electron
       reflector; e.g. Hagemann/Palik compilations put k ~ 4-5 mid-visible).
       k is scaled linearly with wavelength (k ~ lambda, the Drude limit
       omega*tau >> 1); n held constant. Crude but adequate at y <= 0.25.
     - Yb: n ~ 1.1, k ~ 2.6, weak dispersion (both held constant). Yb is a
       LOSSY metal in the visible (divalent rare earth, low k/n ratio),
       which is why it is expected to hurt more than Mg.
   The granular film is then a second Bruggeman mix of that alloy with void
   (eps = 1) at void fraction f_void, i.e. the same effective-medium
   description that fit our measured 5 nm film in scripts/62 (metal fraction
   0.88 there -> f_void ~ 0.12 at d = 5 nm).

2. MORPHOLOGY COUPLING (model assumption, to be pinned by the kMC track):
       f_void(d, y) = f0(d) * s(y)
   f0(d) = 0.12 * (7.5/d - 0.5), clipped to [0, 0.35]: reproduces the
   measured f_void = 0.12 at d = 5 nm, grows for thinner films, vanishes by
   d ~ 15 nm (continuous film). s(y) ramps linearly from 1 at y = 0 down to
   s_min at y >= y_sat. Three smoothing scenarios are carried in parallel
   because the true s(y) is unknown until the kMC agent reports:
       strong: s_min = 0.2, y_sat = 0.10
       medium: s_min = 0.5, y_sat = 0.10
       weak:   s_min = 0.8, y_sat = 0.15

3. STACK. air / film(d) / HATCN 15 nm / glass (incoherent), film-side
   incidence, using stack_TR copied from scripts/62 (same n - ik sign
   convention, hence the np.conj calls). HATCN n = 1.95, k = 2e-4; glass
   n = 1.52 (as in 62). Figure of merit A = 1 - T - R averaged 430-700 nm.

4. SHEET RESISTANCE.  Rs = rho(y) / (d - d_dead(y)).
   rho(y) = rho_Ag + C * y * (1 - y)   (Nordheim rule for dilute alloys),
   rho_Ag = 1.6 uOhm.cm (bulk Ag, standard). C = 48.9 uOhm.cm chosen so
   rho(0.10) ~ 6 uOhm.cm -- literature-typical for Mg:Ag / dilute Ag alloy
   films. CAVEAT: reported thin-film alloy resistivities scatter by 2-3x
   with deposition conditions; C is an order-of-magnitude anchor, not a
   measured constant. d_dead(y) = 2.5 nm * s(y): the dead/percolation-
   deficit thickness shrinks as the film smooths (model assumption tied to
   the same s(y) as the optics).

5. OPTIMIZATION. Grid d = 3..10 nm (0.25 steps), y = 0..0.25 (0.01 steps),
   for {Mg, Yb} x {strong, medium, weak}: minimize A subject to
   Rs <= 100 Ohm/sq; also report the unconstrained minimum and the A(y)
   curve at the optimal d (the U-shape test).

Output: printed tables + dft_seedlayer_screen/runs/alloy_optics.json.
"""
import json
import os
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
F_OUT = os.path.join(_HERE, "..", "runs", "alloy_optics.json")

# ---- constants copied from scripts/62_eta_from_EMA.py -----------------------
JC_WL = np.array([397.4, 413.3, 430.5, 450.9, 471.4, 495.9,
                  521.0, 548.6, 582.1, 616.8, 659.5, 704.5])
AG_N = np.array([0.05, 0.05, 0.04, 0.04, 0.05, 0.05, 0.05, 0.06, 0.05, 0.06, 0.05, 0.04])
AG_K = np.array([2.07, 2.21, 2.36, 2.66, 2.83, 3.09, 3.34, 3.59, 3.93, 4.15, 4.48, 4.84])
N_GLASS, N_HATCN, K_HATCN = 1.52, 1.95, 2.0e-4
D_HATCN = 15.0                      # nm, as in the 62 fit's preferred branch

WL = np.arange(430.0, 700.01, 10.0)  # visible average window for A

# ---- dopant optical constants (APPROXIMATE, see header) ---------------------
def eps_dopant(metal, wl):
    if metal == "Mg":
        n = np.full_like(wl, 0.40)          # ~const across visible (approx.)
        k = 4.4 * wl / 520.0                # Drude-like k ~ lambda, 4.4 @ 520
    elif metal == "Yb":
        n = np.full_like(wl, 1.10)          # weak dispersion (approx.)
        k = np.full_like(wl, 2.60)
    else:
        raise ValueError(metal)
    return (n + 1j * k) ** 2

# ---- morphology coupling -----------------------------------------------------
def f0_void(d):
    return float(np.clip(0.12 * (7.5 / d - 0.5), 0.0, 0.35))

def smooth_factor(y, s_min, y_sat):
    if y >= y_sat:
        return s_min
    return 1.0 - (1.0 - s_min) * (y / y_sat)

SCENARIOS = {"strong": (0.2, 0.10), "medium": (0.5, 0.10), "weak": (0.8, 0.15)}

# ---- sheet resistance ---------------------------------------------------------
RHO_AG = 1.6          # uOhm.cm, bulk Ag
C_NORD = (6.0 - RHO_AG) / (0.10 * 0.90)   # = 48.9 uOhm.cm -> rho(0.10) = 6

def sheet_res(d, y, s):
    rho = RHO_AG + C_NORD * y * (1.0 - y)   # uOhm.cm
    d_eff = d - 2.5 * s                     # nm
    if d_eff <= 0:
        return np.inf
    return rho * 1e-8 / (d_eff * 1e-9)      # Ohm/sq (1 uOhm.cm = 1e-8 Ohm.m)

# ---- EMA + TMM machinery (bruggeman generalized; coherent/stack_TR from 62) ---
def bruggeman2(eps1, eps2, f1):
    """General two-phase 3D Bruggeman (spherical inclusions), phase-1
    fraction f1. Quadratic 2e^2 - Be - eps1*eps2 = 0 with
    B = eps1*(3f1-1) + eps2*(2-3f1); absorbing root (Im >= 0)."""
    B = eps1 * (3 * f1 - 1) + eps2 * (2 - 3 * f1)
    e = (B + np.sqrt(B * B + 8 * eps1 * eps2)) / 4
    bad = e.imag < 0
    if np.any(bad):
        e2 = (B - np.sqrt(B * B + 8 * eps1 * eps2)) / 4
        e = np.where(bad, e2, e)
    return e

def coherent(ns, ds, lam):
    M = np.eye(2, dtype=complex)
    for N, d in zip(ns[1:-1], ds[1:-1]):
        delta = 2 * np.pi / lam * N * d
        c, s = np.cos(delta), np.sin(delta)
        M = M @ np.array([[c, 1j * s / N], [1j * N * s, c]])
    n0, nsub = ns[0], ns[-1]
    B, C = M @ np.array([1.0, nsub], dtype=complex)
    return (float(4 * n0.real * nsub.real / abs(n0 * B + C) ** 2),
            float(abs((n0 * B - C) / (n0 * B + C)) ** 2))

def stack_TR(layers, lam):
    """Film-side (air) incidence; incoherent glass slab behind."""
    ns = [complex(1.0)] + [l[0] for l in layers] + [complex(N_GLASS)]
    ds = [0.0] + [l[1] for l in layers] + [0.0]
    T_in, R_front = coherent(ns, ds, lam)
    T_back, R_back = coherent(ns[::-1], ds[::-1], lam)
    r = (N_GLASS - 1.0) / (N_GLASS + 1.0)
    Rr, Tr = r ** 2, 1 - r ** 2
    den = 1 - R_back * Rr
    return T_in * Tr / den, R_front + T_in * Rr * T_back / den

# ---- absorption of one design --------------------------------------------------
EPS_AG = (np.interp(WL, JC_WL, AG_N) + 1j * np.interp(WL, JC_WL, AG_K)) ** 2
HAT = np.conj(np.full_like(WL, N_HATCN + 1j * K_HATCN, dtype=complex))
_EPS_DOP = {}   # cache per metal
_cache = {}

def absorption(metal, d, y, s_min, y_sat):
    key = (metal, round(d, 3), round(y, 3), s_min, y_sat)
    if key in _cache:
        return _cache[key]
    if metal not in _EPS_DOP:
        _EPS_DOP[metal] = eps_dopant(metal, WL)
    s = smooth_factor(y, s_min, y_sat)
    fv = f0_void(d) * s
    eps_alloy = EPS_AG if y == 0 else bruggeman2(_EPS_DOP[metal], EPS_AG, y)
    eps_film = bruggeman2(eps_alloy, np.ones_like(WL, dtype=complex), 1.0 - fv)
    n_eff = np.conj(np.sqrt(eps_film))     # n - ik convention of coherent()
    A = 0.0
    for i in range(len(WL)):
        T, R = stack_TR([(n_eff[i], d), (HAT[i], D_HATCN)], WL[i])
        A += 1.0 - T - R
    A /= len(WL)
    _cache[key] = A
    return A

# ---- optimization ---------------------------------------------------------------
def main():
    D_GRID = np.arange(3.0, 10.01, 0.25)
    Y_GRID = np.arange(0.0, 0.2501, 0.01)
    results = {}
    for metal in ("Mg", "Yb"):
        for scen, (s_min, y_sat) in SCENARIOS.items():
            best_c = None   # constrained
            best_u = None   # unconstrained
            base_c = None   # y = 0, constrained (its own best d)
            for d in D_GRID:
                for y in Y_GRID:
                    A = absorption(metal, d, y, s_min, y_sat)
                    Rs = sheet_res(d, y, smooth_factor(y, s_min, y_sat))
                    rec = (A, float(y), float(d), Rs)
                    if best_u is None or A < best_u[0]:
                        best_u = rec
                    if Rs <= 100.0:
                        if best_c is None or A < best_c[0]:
                            best_c = rec
                        if y == 0.0 and (base_c is None or A < base_c[0]):
                            base_c = rec
            A_c, y_c, d_c, Rs_c = best_c
            A_b, _, d_b, Rs_b = base_c
            curve = [{"y": float(y),
                      "A": absorption(metal, d_c, float(y), s_min, y_sat),
                      "Rs": sheet_res(d_c, float(y),
                                      smooth_factor(float(y), s_min, y_sat))}
                     for y in Y_GRID]
            results[f"{metal}_{scen}"] = {
                "metal": metal, "scenario": scen,
                "s_min": s_min, "y_sat": y_sat,
                "constrained_opt": {"y": y_c, "d_nm": d_c, "A": A_c,
                                    "Rs_ohm_sq": Rs_c},
                "unconstrained_opt": {"y": best_u[1], "d_nm": best_u[2],
                                      "A": best_u[0], "Rs_ohm_sq": best_u[3]},
                "baseline_y0": {"d_nm": d_b, "A": A_b, "Rs_ohm_sq": Rs_b},
                "doping_wins": bool(A_c < A_b - 1e-6),
                "A_vs_y_at_dopt": curve,
            }
            print(f"\n== {metal}  |  smoothing {scen} "
                  f"(s_min={s_min}, y_sat={y_sat}) ==")
            print(f"  constrained opt : y*={y_c:.2f}  d*={d_c:.2f} nm  "
                  f"A={100*A_c:.2f}%  Rs={Rs_c:.1f} Ohm/sq")
            print(f"  unconstrained   : y ={best_u[1]:.2f}  d ={best_u[2]:.2f} nm  "
                  f"A={100*best_u[0]:.2f}%  Rs={best_u[3]:.1f} Ohm/sq")
            print(f"  y=0 baseline    : d ={d_b:.2f} nm  A={100*A_b:.2f}%  "
                  f"Rs={Rs_b:.1f} Ohm/sq")
            print(f"  doping wins on A? {'YES' if A_c < A_b - 1e-6 else 'no'}"
                  f"  (delta A = {100*(A_c - A_b):+.2f} %p)")
            ys = [c["y"] for c in curve]; As = [100 * c["A"] for c in curve]
            print(f"  A(y) at d*={d_c:.2f} nm:")
            for j in range(0, len(ys), 5):
                print("    " + "  ".join(f"y={ys[k]:.2f}:{As[k]:5.2f}%"
                                         for k in range(j, min(j + 5, len(ys)))))
    os.makedirs(os.path.dirname(F_OUT), exist_ok=True)
    with open(F_OUT, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nwrote {os.path.normpath(F_OUT)}")


# RESULT (run 2026-08-15, grid d=3..10 nm x y=0..0.25, both metals, 3 scenarios)
#
# DOES DOPING EVER WIN ON THE ABSORPTION-FIRST CRITERION?  No -- in all six
# cases (Mg/Yb x strong/medium/weak smoothing) the constrained optimum is
# y* = 0, d* = 3.0 nm, A = 0.51 %, Rs = 32 Ohm/sq, identical to the y = 0
# baseline. The unconstrained and constrained optima coincide because the
# Rs <= 100 Ohm/sq constraint NEVER BINDS: even at d = 3 nm with the full
# 2.5 nm dead layer, bulk-rho Ag gives Rs = 32 Ohm/sq. With the constraint
# slack, the thinnest film wins on A and doping has nothing to buy back.
#
# WHY THE PREDICTED U-SHAPE DOES NOT APPEAR. A(y) at d* is monotonically
# INCREASING for every metal and scenario (Mg strong: 0.51 % -> 1.94 % over
# y = 0 -> 0.25; Yb strong: 0.51 % -> 4.10 %). Two reasons, both properties
# of this model rather than of nature:
#   1. In a 3D Bruggeman film the void-driven "plasmonic" loss only becomes
#      large near the percolation anomaly (f_void >~ 0.5 here: A jumps from
#      ~0.8 % at f_void = 0.35 to 4-8 % at f_void = 0.5-0.6). Our calibrated
#      f0(d) tops out at 0.24 (d = 3 nm), pinned by the measured f_void =
#      0.12 at d = 5 nm from scripts/62 -- so the film never sits close
#      enough to percolation for smoothing to pay. In 0.06 < f_void < 0.24
#      the Bruggeman loss is in fact locally DECREASING in f_void, so
#      smoothing can even raise A slightly (strong scenario A(y=0.10) = 1.18 %
#      vs weak 0.89 % at d = 3, Mg).
#   2. The alloy loss is strictly increasing in y from the first grid point,
#      and it is never offset. Sensitivity check: had f0(3 nm) been ~0.5
#      (i.e. films near percolation), strong smoothing at y = 0.10 would cut
#      A by ~3 %p against an alloy penalty of ~0.7 %p and doping would win
#      decisively -- the verdict is entirely controlled by how close the
#      real 3-5 nm film is to percolation, which is exactly what the kMC
#      track must pin down (along with s(y)).
#
# Mg vs Yb: as expected Yb loses HARDER -- its A(y) slope is ~2x Mg's
# (Yb: +14 %p/unit-y vs Mg: +7 %p/unit-y near y = 0, strong scenario), because
# Yb's k/n ratio makes the alloy far lossier. If any smoothing benefit ever
# materializes (near-percolation films), Mg is the dopant of choice.
#
# SCENARIO SENSITIVITY: the no-win verdict is INSENSITIVE to the smoothing
# scenario (identical optimum in all three) -- but only because the Rs
# constraint is slack and f_void stays below percolation. Both hinge on model
# assumptions flagged in the header: bulk rho at 0.5 nm effective thickness
# (a 3 nm pure-Ag film is optimistic there; real percolation failure at
# d ~ 3 nm would push d* up, re-tighten Rs, and could re-open the door for
# doping), and the f0(d) extrapolation below 5 nm. Treat the d* = 3 nm corner
# as "as thin as the model allows", not a process recommendation.

if __name__ == "__main__":
    main()
