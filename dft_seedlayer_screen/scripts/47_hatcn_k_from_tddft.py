"""A first-principles k(lambda) for neutral HATCN, to replace the k = 0 assumption.

WHY. scripts/32 assumed k = 0 for HATCN and scripts/46 showed what that assumption
is worth (nothing -- it returns the input). The literature route to a measured k is
closed: this session's network policy answers 403 to CONNECT for every external
host, so no reported n,k can be retrieved here. What CAN be done is compute the
spectrum and convert it to k, which is what this does.

THE INPUT. scripts/44 already ran TD-DFT (wB97X/def2-SVP, TDA) on neutral HATCN.
The result is the useful one:

    lowest excited state   323.4 nm (3.83 eV), DARK (f = 0.0000)
    lowest bright state    316.2 nm (3.92 eV), f = 0.0133
    main band              266.5 nm (4.65 eV), f = 0.73

Twelve states were requested and ALL of them lie above 3.8 eV. Because TD-DFT
returns the LOWEST n states, the absence of a visible transition is a real result
and not a truncation artefact -- had a visible state existed it would have been
returned first. Neutral HATCN has no electronic transition in 400-700 nm.

SO WHY IS k NOT ZERO. Two reasons, and only the first is captured by a bulk k:

  1. The absorption edge is not a step. Vibronic and inhomogeneous broadening put
     an exponential (Urbach) tail below the lowest transition, and that tail
     reaches into the blue. This is what the script computes.

  2. Interfacial charge-transfer absorption. Where HATCN meets a donor (the HTL)
     or the Ag electrode, CT states form that DO absorb in the visible -- that is
     the mechanism that makes HATCN work as a HIL in the first place. This scales
     with interface AREA, not with layer thickness, so it is not a bulk k at all
     and this script does not and cannot bound it. It is the reason the claim must
     ultimately be checked by measuring the real stack, not a neat film.

CAVEATS ON THE POSITION OF THE EDGE. Vacuum, vertical, TDA. TDA overestimates
excitation energies by ~0.1-0.2 eV and solid-state polarisation red-shifts them by
~0.1-0.3 eV, both pushing the edge to lower energy. EDGE_SHIFT applies that
correction pessimistically -- i.e. in the direction that makes the transparency
claim harder -- so the resulting k is an upper bound rather than a best estimate.
"""
import os, json
import numpy as np

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")

# physical constants (SI)
E_CHG, EPS0, M_E, C0 = 1.602176634e-19, 8.8541878128e-12, 9.1093837015e-31, 2.99792458e8
HC_EV_NM = 1239.841984

# HATCN film. MW of C18N12 = 384.3 g/mol. Amorphous organic films of this kind sit
# near 1.6 g/cm3; the result scales linearly in density, so it is reported.
MW_HATCN = 384.28          # g/mol
RHO = 1.6                  # g/cm3
N_DENS = RHO * 1e3 / (MW_HATCN * 1e-3) * 6.02214076e23   # molecules / m^3

FWHM_EV = 0.35             # vibronic + inhomogeneous broadening of a solid film
URBACH_EV = 0.10           # Urbach energy; 50-150 meV typical for amorphous organics
EDGE_SHIFT = -0.35         # eV, TDA overestimate + solid-state red shift (pessimistic)

WL = np.arange(380.0, 781.0, 2.0)


def load_states(tag="HATCN"):
    d = json.load(open(os.path.join(RUNS, "tddft_absorption.json")))
    rows = d.get(tag) or []
    return np.array([(r["eV"], r["f"]) for r in rows])


def k_spectrum(states, wl=WL):
    """Convert oscillator strengths to k(lambda) via Gaussian bands + Urbach tail.

    Integrated absorption coefficient of a transition, in SI with nu in Hz:
        int(alpha) d(nu) = N f e^2 / (4 eps0 m_e c)
    Each state is spread as a normalised Gaussian in energy carrying that integral.
    """
    ev = HC_EV_NM / wl
    alpha = np.zeros_like(ev)                       # m^-1
    sigma = FWHM_EV / (2 * np.sqrt(2 * np.log(2)))

    for e0, f in states:
        if f <= 0:
            continue                                # dark state carries no intensity
        e0 = e0 + EDGE_SHIFT
        integ = N_DENS * f * E_CHG ** 2 / (4 * EPS0 * M_E * C0)     # m^-1 * Hz
        # Change the integration variable from nu to E. E = h*nu, so d(nu) = d(E)/h
        # and int(alpha)dE = h * int(alpha)d(nu) -- MULTIPLY by h. Dividing instead
        # inflates every band by ~1/h^2 ~ 1e29 and produced k ~ 1e25, which is what
        # made the error visible; a subtler unit slip here would not have been.
        integ_ev = integ * 4.135667696e-15          # m^-1 * eV
        alpha += integ_ev * np.exp(-0.5 * ((ev - e0) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))

    # Urbach tail below the lowest BRIGHT transition. A Gaussian tail decays far too
    # fast to represent a real band edge, and using it alone would understate k in
    # the blue -- exactly the region the transparency claim is weakest.
    bright = [e + EDGE_SHIFT for e, f in states if f > 1e-3]
    if bright:
        e_edge = min(bright)
        a_edge = np.interp(e_edge, ev[::-1], alpha[::-1])
        tail = a_edge * np.exp((ev - e_edge) / URBACH_EV)
        alpha = np.where(ev < e_edge, np.maximum(alpha, tail), alpha)

    k = alpha * (wl * 1e-9) / (4 * np.pi)
    return k, alpha


def main():
    states = load_states("HATCN")
    print(f"HATCN neutral: {len(states)} TD-DFT states, "
          f"lowest {HC_EV_NM/states[:,0].max():.0f}-{HC_EV_NM/states[:,0].min():.0f} nm")
    print(f"states in 400-700 nm: "
          f"{int(((HC_EV_NM/states[:,0] >= 400) & (HC_EV_NM/states[:,0] <= 700)).sum())}")
    print(f"film: rho = {RHO} g/cm3 -> N = {N_DENS:.3e} m^-3")
    print(f"edge shift applied: {EDGE_SHIFT:+.2f} eV (pessimistic)\n")

    k, alpha = k_spectrum(states)
    print(f"{'lambda (nm)':>12}{'k':>12}{'alpha (cm^-1)':>16}")
    for target in (400, 450, 500, 550, 600, 650, 700):
        i = int(np.argmin(abs(WL - target)))
        print(f"{WL[i]:>12.0f}{k[i]:>12.2e}{alpha[i]/100:>16.2f}")

    kmax_vis = k[(WL >= 400) & (WL <= 700)].max()
    print(f"\nmax k over 400-700 nm: {kmax_vis:.2e}")

    # compare against the bounds from scripts/46
    bounds = {3: 0.03, 5: 0.01, 10: 0.01, 30: 0.003}
    print("\nagainst the tolerance bounds of scripts/46 (<0.5 %p absorption):")
    print(f"{'HATCN (nm)':>12}{'k allowed':>12}{'k computed':>13}{'margin':>10}")
    for d, kb in bounds.items():
        print(f"{d:>12}{kb:>12.3f}{kmax_vis:>13.2e}{kb/kmax_vis:>9.0f}x")

    json.dump({"wl_nm": WL.tolist(), "k": k.tolist(),
               "alpha_cm-1": (alpha / 100).tolist(),
               "rho_g_cm3": RHO, "fwhm_eV": FWHM_EV, "urbach_eV": URBACH_EV,
               "edge_shift_eV": EDGE_SHIFT, "k_max_visible": float(kmax_vis)},
              open(os.path.join(RUNS, "hatcn_k_tddft.json"), "w"), indent=2)

    print("\nWHAT THIS DOES AND DOES NOT SUPPORT")
    print("  Supports: neutral HATCN has no visible electronic transition, and the")
    print("  band-edge tail leaves k far below the level at which a 3-5 nm layer")
    print("  would absorb measurably. The transparency claim survives without")
    print("  needing k = 0, which is the point.")
    print("  Does NOT support: any statement about the HATCN/donor or HATCN/Ag")
    print("  interface. CT absorption there is a real visible-range loss, is not a")
    print("  bulk k, and is not bounded by this calculation.")
    print("  Sensitivity: k scales linearly with film density and exponentially")
    print("  with the Urbach energy. Doubling URBACH_EV to 0.2 eV is the single")
    print("  most damaging assumption -- rerun with it before quoting a margin.")
    print("\n  MEASURE IT: ellipsometry or transmission on the 30 nm HATCN film")
    print("  already deposited settles this directly, and 30 nm is thick enough")
    print("  that a k of 1e-3 is detectable.")


if __name__ == "__main__":
    main()
