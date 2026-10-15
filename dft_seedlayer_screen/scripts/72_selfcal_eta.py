"""Self-calibration of the R geometry factor from the HATCN2_Ag25 sample,
then absorption of the new HATCN(2-30)/Ag(3-5) series, vs the SEM impression.

WHY THIS WORKS NOW. The earlier eta estimate (scripts/62) had to lean on a
Bruggeman model of a granular film - the honest answer was "0.70, uncertainty
reaching 0.9". This new series contains a 25 nm Ag sample: thick enough to be
bulk-like (percolation long past, J&C n,k valid, roughness a minor
perturbation), so its absolute R(lambda) is COMPUTABLE. Comparing computed
R_abs with the file value R_file = eta * R_abs / T_base pins eta from
measured data, replacing the model-derived number.

DATA SOURCE / PRECISION. Values are read off the user's Origin plot at anchor
wavelengths (no CSVs yet) - carry +-2 %p on every input. CSVs will sharpen
numbers but not verdicts.

MEASUREMENT MODEL (same session geometry as scripts/61-63):
    T_file = T_abs / T_base            (clean relative transmittance)
    R_file = eta * R_abs / T_base      (transmission baseline reused in the
                                        6/12-degree reflection accessory)
"""
import numpy as np

JC_WL = np.array([397.4, 413.3, 430.5, 450.9, 471.4, 495.9,
                  521.0, 548.6, 582.1, 616.8, 659.5, 704.5, 756.0, 800.0])
AG_N = np.array([0.05, 0.05, 0.04, 0.04, 0.05, 0.05, 0.05, 0.06, 0.05, 0.06,
                 0.05, 0.04, 0.03, 0.04])
AG_K = np.array([2.07, 2.21, 2.36, 2.66, 2.83, 3.09, 3.34, 3.59, 3.93, 4.15,
                 4.48, 4.84, 5.28, 5.63])
N_GLASS, N_HATCN = 1.52, 1.95

# ---- plot-read inputs (percent), +-2 %p ----------------------------------
WL = np.array([450.0, 550.0, 650.0, 750.0])
AG25_R_FILE = np.array([65.0, 82.0, 89.0, 93.5])
SAMPLES = {  # name: (d_hatcn, T_rel(%), R_file(%))
    "HATCN2_Ag3":  (2.0,  np.array([84, 77, 71, 67]), np.array([14, 16, 19, 22])),
    "HATCN30_Ag3": (30.0, np.array([87, 80, 74, 69]), np.array([14, 16, 19, 22])),
    "HATCN2_Ag5":  (2.0,  np.array([81, 77, 74, 71]), np.array([14, 17, 20, 24])),
    "HATCN30_Ag5": (30.0, np.array([85, 82, 78, 74]), np.array([15, 18, 22, 26])),
}


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


def glass_stack(layers, lam):
    """air | layers | glass slab (incoherent back)."""
    ns = [complex(1.0)] + [l[0] for l in layers] + [complex(N_GLASS)]
    ds = [0.0] + [l[1] for l in layers] + [0.0]
    T_in, R_f = coherent(ns, ds, lam)
    T_b, R_b = coherent(ns[::-1], ds[::-1], lam)
    r = (N_GLASS - 1) / (N_GLASS + 1)
    Rr, Tr = r * r, 1 - r * r
    den = 1 - R_b * Rr
    return T_in * Tr / den, R_f + T_in * Rr * T_b / den


def main():
    nag = np.interp(WL, JC_WL, AG_N) + 1j * np.interp(WL, JC_WL, AG_K)
    hat = N_HATCN + 2e-4j

    # baselines: glass/HATCN(d) relative-T references
    Tb2 = np.array([glass_stack([(np.conj(hat), 2.0)], w)[0] for w in WL])
    Tb30 = np.array([glass_stack([(np.conj(hat), 30.0)], w)[0] for w in WL])

    # 1. computable absolute R of the 25 nm sample (film-side incidence)
    R25 = np.array([glass_stack([(np.conj(nag[i]), 25.0),
                                 (np.conj(hat), 2.0)], WL[i])[1]
                    for i in range(len(WL))])
    eta = AG25_R_FILE / 100 * Tb2 / R25
    print("SELF-CALIBRATION from HATCN2_Ag25  (R_theory: TMM, J&C bulk):")
    for i, w in enumerate(WL):
        print(f"  {w:4.0f} nm  R_theory {100*R25[i]:5.1f}%  R_file "
              f"{AG25_R_FILE[i]:5.1f}%  -> eta {eta[i]:.3f}")
    e0 = float(np.mean(eta[1:]))       # 550-750; 450 nm carries HATCN edge risk
    print(f"  eta = {e0:.2f}  (mean 550-750; scripts/62 EMA said 0.70 with "
          f"honest ceiling ~0.9 - the DATA now picks the ceiling)")

    # 2. absolute quantities for the thin samples
    print("\nABSORPTION of the thin series (A = 1 - T_rel*T_base - "
          "R_file*T_base/eta):")
    print(f"{'sample':<14}{'A@450':>8}{'A@550':>8}{'A@650':>8}{'A@750':>8}"
          f"{'  R_abs@750':>11}")
    for name, (dh, T_rel, R_file) in SAMPLES.items():
        Tb = Tb2 if dh < 10 else Tb30
        A = 1 - T_rel / 100 * Tb - R_file / 100 * Tb / e0
        Rab = R_file[-1] / 100 * Tb[-1] / e0
        print(f"{name:<14}" + "".join(f"{100*a:7.1f}%" for a in A)
              + f"{100*Rab:10.1f}%")
    print("""
NOTE +-2 %p on every plot-read input propagates to ~+-3 %p on A.
Ideal continuous films at these thicknesses absorb 1-2 % (Ag3) / 2-3 % (Ag5).
""")


if __name__ == "__main__":
    main()

# RESULT (run 2026-08-15, plot-read inputs +-2 %p):
#
# 1. eta SELF-CALIBRATED FROM DATA: 0.96 +- 0.01, FLAT across 550-750 nm
#    (0.960/0.965/0.961). Flatness was exactly the self-test the EMA
#    estimate failed (15 % ripple) - the Ag25 anchor passes it. This
#    SUPERSEDES the scripts/62 model estimate (0.70): the reflection
#    accessory loses almost nothing vs the open transmission path, and
#    R_file ~ R_abs/T_base within ~4 %. Caveat: if the 25 nm film has
#    pinholes its real R is below the TMM value and eta -> 1.0; so quote
#    eta = 0.96-1.00. All earlier absorption numbers move to the HIGH
#    branch: sample 5 (HATCN/Ag~5) A ~ 9-10 % visible-average (was
#    "3-9 %"); sample 4 (holey 3-4 nm) A ~ 16-17 %.
#
# 2. NEW SERIES ABSORPTION (eta 0.96): every thin film absorbs 10-19 %,
#    vs 1-3 % for ideal continuous films. Spectral shapes are diagnostic:
#    Ag3 films: A RISES to the red (14 -> 17.5 % for HATCN2_Ag3) =
#    plasmonic tail of a barely-percolated network. Ag5 films: A FALLS to
#    the red (13 -> 10.5-12 %) = past percolation, residual LSPR toward
#    green. HATCN30_Ag5 is the best film in the red (A 10.5 %).
#
# 3. SEM CONTRADICTION: all nine SEMs look featureless at 75 kX, but
#    10-17 % absorption CANNOT come from a smooth continuous film. The
#    morphology that absorbs (few-10-nm voids/grain boundaries) is below
#    this magnification/contrast. SEM at this setting is NOT a valid
#    smoothness check for the absorption criterion - need 200 kX+ InLens,
#    AFM roughness, and per-sample Rs (the conductivity axis would also
#    separate Ag3-barely-percolated from Ag5-percolated immediately).
#
# 4. INTERNAL CHECK THAT FLAGS A T PROBLEM: measured T of HATCN2_Ag25 at
#    800 nm reads ~10 % relative (~9 % absolute); TMM says a continuous
#    25 nm Ag film transmits 2-3 %. Either the plot was misread (+-2 %p
#    does not cover this) or there is a stray-light/pinhole floor in the
#    T channel at low transmittance. Check against the CSV.
#
# ACTIONS: (i) get the CSVs to replace plot-read values, (ii) measure Rs
# for all nine samples, (iii) re-image 2-3 samples at 200 kX / AFM,
# (iv) propagate eta = 0.96-1.00 through scripts/62/63 quoted numbers.
