"""Is the measured T/R on HATCN/Ag usable, and what does it say? A data audit.

THE FILES. 24T.csv and 24R.csv, HATCN(nominally 30 nm, likely 15-20) plus Ag
(nominally 5 nm) on soda lime glass, 300-800 nm. Transmittance is RELATIVE, with
HATCN-on-glass as the baseline; reflectance is presumed absolute.

WHAT THIS CHECKS, in order:
  1. Where the data stops being data. Soda lime glass cuts off in the near UV, so
     the baseline transmittance goes to zero and the ratio explodes -- the file
     reaches 580 % at 318 nm. The usable edge is found from the point-to-point
     scatter rather than assumed.
  2. Whether T and R can be combined. They cannot as they stand: T is relative and
     R is absolute, so 1 - T - R is not the absorptance and comes out NEGATIVE
     below ~470 nm in this file. That is a normalisation mismatch, not a
     measurement fault, and it is fixed by converting T to absolute first.
  3. What the numbers imply for the film, once the normalisation is consistent:
     the absorptance spectrum, and the Ag thickness that best reproduces both T
     and R through the project's TMM.

The Ag thickness fit is the useful output. The nominal 5 nm is a quartz-crystal
reading; what the optics sees is an effective thickness, and for a film this thin
the two need not agree.
"""
import os, sys
import numpy as np

# The upload directory is wiped by rollbacks (it took the only copy of these
# measurements with it once); the repo copy in data/ is authoritative now.
_HERE = os.path.dirname(os.path.abspath(__file__))
F_T = os.path.join(_HERE, "..", "data", "24T.csv")
F_R = os.path.join(_HERE, "..", "data", "24R.csv")

JC_WL = np.array([397.4, 413.3, 430.5, 450.9, 471.4, 495.9,
                  521.0, 548.6, 582.1, 616.8, 659.5, 704.5])
AG_N = np.array([0.05, 0.05, 0.04, 0.04, 0.05, 0.05, 0.05, 0.06, 0.05, 0.06, 0.05, 0.04])
AG_K = np.array([2.07, 2.21, 2.36, 2.66, 2.83, 3.09, 3.34, 3.59, 3.93, 4.15, 4.48, 4.84])
N_GLASS, N_AIR, N_HATCN, K_HATCN = 1.52, 1.0, 1.95, 2.0e-4

# The instrument swaps lamps near 720 nm and the data steps there. It is visible
# in R -- 23.19 at 720, 23.62 at 718, 22.93 at 714 -- a ~0.7 %p excursion against
# ~0.05 %p scatter either side. Excluded from the fit rather than smoothed over,
# because a lamp step is a systematic offset between two halves of the spectrum,
# not noise that averages away.
LAMP_LO, LAMP_HI = 700.0, 735.0
FIT_LO, FIT_HI = 430.0, 800.0


def fit_mask(wl):
    return (wl >= FIT_LO) & (wl <= FIT_HI) & ~((wl > LAMP_LO) & (wl < LAMP_HI))


def load(path, col):
    wl, y = [], []
    for line in open(path):
        p = line.strip().split(",")
        if len(p) < 2:
            continue
        try:
            wl.append(float(p[0])); y.append(float(p[1]))
        except ValueError:
            continue
    o = np.argsort(wl)
    return np.array(wl)[o], np.array(y)[o]


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


def stack_TR(layers, lam, n_ag_l):
    """Coherent film stack on a thick incoherent glass substrate, air both sides."""
    ns = [complex(N_AIR)] + [l[0] for l in layers] + [complex(N_GLASS)]
    ds = [0.0] + [l[1] for l in layers] + [0.0]
    T_in, R_front = coherent(ns, ds, lam)
    R_back_stack = coherent(ns[::-1], ds[::-1], lam)[1]
    r = (N_GLASS - N_AIR) / (N_GLASS + N_AIR)
    R_rear, T_rear = r ** 2, 1 - r ** 2
    denom = 1 - R_back_stack * R_rear
    T = T_in * T_rear / denom
    R = R_front + T_in * R_rear * T_in / denom * 0 + \
        (T_in * R_rear * R_back_stack * 0)          # front-surface R dominates
    # proper incoherent R: front reflection plus light that returns through stack
    R = R_front + (T_in * R_rear * coherent(ns[::-1], ds[::-1], lam)[0]) / denom
    return T, R


def main():
    wlT, T = load(F_T, 1)
    wlR, R = load(F_R, 1)
    assert np.allclose(wlT, wlR), "wavelength grids differ"
    wl = wlT

    # ---- 1. where does the data stop being data ----
    print("=" * 72)
    print("1. USABLE RANGE  (point-to-point scatter, soda lime glass cutoff)")
    print("=" * 72)
    print(f"{'band (nm)':>14}{'T noise':>10}{'R noise':>10}{'T mean':>10}   verdict")
    bands = [(735, 800), (700, 735), (600, 700), (500, 600), (450, 500),
             (420, 450), (400, 420), (380, 400), (350, 380), (300, 350)]
    for lo, hi in bands:
        m = (wl >= lo) & (wl < hi)
        if m.sum() < 3:
            continue
        # scatter = RMS of the second difference, insensitive to real slope
        nT = float(np.std(np.diff(T[m], 2)) / np.sqrt(6))
        nR = float(np.std(np.diff(R[m], 2)) / np.sqrt(6))
        v = ("excellent" if nT < 0.15 else "good" if nT < 0.4 else
             "noisy" if nT < 1.5 else "UNUSABLE")
        if lo == 700:
            v += "  <-- LAMP CHANGEOVER, excluded from fit"
        print(f"{lo:>6}-{hi:<7}{nT:>10.3f}{nR:>10.3f}{T[m].mean():>10.2f}   {v}")
    print("""
  Soda lime glass absorbs in the near UV, so the HATCN-on-glass baseline goes to
  zero there and the ratio blows up: the file reaches 580 % at 318 nm and 340 %
  at 324 nm. Those are not measurements. The user's instinct to discard below
  400 nm is right, and the scatter puts the honest edge at about 430 nm.""")

    # ---- 2. can T and R be combined ----
    print("\n" + "=" * 72)
    print("2. NORMALISATION CHECK  (is 1 - T - R a valid absorptance?)")
    print("=" * 72)
    print(f"{'lambda':>8}{'T rel':>9}{'R':>8}{'100-T-R':>10}   ")
    for target in (750, 700, 650, 600, 550, 500, 460, 440, 420):
        i = int(np.argmin(abs(wl - target)))
        a = 100 - T[i] - R[i]
        flag = "  <-- NEGATIVE, unphysical" if a < 0 else ""
        print(f"{wl[i]:>8.0f}{T[i]:>9.2f}{R[i]:>8.2f}{a:>10.2f}{flag}")
    print("""
  Negative absorptance below ~460 nm is the giveaway. It is NOT a bad
  measurement: T is relative to HATCN-on-glass while R is absolute, so the two
  are on different denominators and 1 - T - R mixes them. Multiplying T by the
  baseline's own transmittance fixes it.""")

    # ---- 3. absolute T, absorptance, and the Ag thickness that fits ----
    n_ag_arr = np.conj(np.interp(wl, JC_WL, AG_N) + 1j * np.interp(wl, JC_WL, AG_K))
    hat = np.conj(np.full_like(wl, N_HATCN + 1j * K_HATCN, dtype=complex))

    print("\n" + "=" * 72)
    print("3. WITH T PUT ON AN ABSOLUTE FOOTING")
    print("=" * 72)
    for d_hatcn in (15.0, 20.0, 30.0):
        Tb = np.array([stack_TR([(hat[i], d_hatcn)], wl[i], n_ag_arr[i])[0]
                       for i in range(len(wl))])
        m = (wl >= 430) & (wl <= 800)
        Tabs = T[m] / 100 * Tb[m]
        A = 1 - Tabs - R[m] / 100
        print(f"  HATCN {d_hatcn:.0f} nm baseline: T_base(550) = "
              f"{100*Tb[int(np.argmin(abs(wl-550)))]:.1f} %, "
              f"absorptance = {100*A.min():.1f} to {100*A.max():.1f} %")

    print("\n" + "=" * 72)
    print("4. WHAT Ag THICKNESS FITS BOTH T AND R")
    print("=" * 72)
    m = fit_mask(wl)
    best = None
    print(f"{'d_HATCN':>9}{'d_Ag':>7}{'RMS T':>9}{'RMS R':>9}{'RMS tot':>10}")
    for d_h in (15.0, 20.0, 30.0):
        Tb = np.array([stack_TR([(hat[i], d_h)], wl[i], n_ag_arr[i])[0]
                       for i in range(len(wl))])
        for d_ag in np.arange(2.0, 9.01, 0.25):
            Tm, Rm = np.zeros(len(wl)), np.zeros(len(wl))
            for i in range(len(wl)):
                Tm[i], Rm[i] = stack_TR([(n_ag_arr[i], d_ag), (hat[i], d_h)],
                                        wl[i], n_ag_arr[i])
            Trel_model = Tm / Tb * 100
            eT = float(np.sqrt(np.mean((Trel_model[m] - T[m]) ** 2)))
            eR = float(np.sqrt(np.mean((Rm[m] * 100 - R[m]) ** 2)))
            tot = np.hypot(eT, eR)
            if best is None or tot < best[0]:
                best = (tot, d_h, d_ag, eT, eR)
        # print the best for this HATCN thickness
        bb = None
        for d_ag in np.arange(2.0, 9.01, 0.25):
            Tm, Rm = np.zeros(len(wl)), np.zeros(len(wl))
            for i in range(len(wl)):
                Tm[i], Rm[i] = stack_TR([(n_ag_arr[i], d_ag), (hat[i], d_h)],
                                        wl[i], n_ag_arr[i])
            eT = float(np.sqrt(np.mean((Tm[m] / Tb[m] * 100 - T[m]) ** 2)))
            eR = float(np.sqrt(np.mean((Rm[m] * 100 - R[m]) ** 2)))
            if bb is None or np.hypot(eT, eR) < bb[0]:
                bb = (np.hypot(eT, eR), d_ag, eT, eR)
        print(f"{d_h:>9.0f}{bb[1]:>7.2f}{bb[2]:>9.2f}{bb[3]:>9.2f}{bb[0]:>10.2f}")
    print(f"\n  best overall: HATCN {best[1]:.0f} nm, Ag {best[2]:.2f} nm "
          f"(RMS {best[3]:.2f} %p in T, {best[4]:.2f} %p in R)")


if __name__ == "__main__":
    main()
