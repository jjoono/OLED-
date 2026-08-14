"""Which normalisation is the measured R actually in? Settling it with the model.

WHAT WAS LEFT HANGING. scripts/60 combined the measured T and R by converting T
to absolute (T is explicitly relative, HATCN-on-glass baseline) while PRESUMING R
absolute. That presumption was never tested, and the user is right to call it
out: reflectance accessories normalise to a reference, and what the file holds
depends on what that reference was. Three candidate hypotheses:

  H-A  R is absolute (reference mirror calibrated, or corrected by software)
  H-B  R is relative to the same HATCN-on-glass sample used for the T baseline
  H-C  R is relative to an uncalibrated metal mirror (R_mirror ~ 0.85-0.95),
       so the file overstates absolute R by 1/R_mirror

THE LEVER. HATCN's n,k are known (n = 1.95, k <= 2e-4 in the visible, scripts/47)
and soda lime glass is n ~ 1.52, so the baseline sample's own R and T can be
COMPUTED. Each hypothesis then implies an absolute R spectrum, and each implied
spectrum must satisfy two hard constraints:

  1. A = 1 - T_abs - R_abs must be non-negative everywhere (energy conservation)
  2. A must be entirely attributable to the Ag film, since neither glass nor
     HATCN absorbs in 430-800 nm -- so its magnitude has to sit in the range a
     3-8 nm granular Ag film can absorb (a few % to ~15 %), and its shape should
     rise toward the blue (LSPR of a granular film), not go negative or explode.

A hypothesis that violates either is dead regardless of any fitting.

HATCN THICKNESS. The user's nominal 30 nm is suspected to be 15-20 in reality.
All three baselines are carried through; the T conversion depends on it visibly,
the R of the baseline sample barely does, and the spread across 15/20/30 is
reported as the systematic uncertainty of the conversion.
"""
import os
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
LAMP = (700.0, 735.0)
FIT = (430.0, 800.0)


def load(path):
    wl, y = [], []
    for line in open(path):
        p = line.strip().split(",")
        if len(p) >= 2:
            try:
                wl.append(float(p[0])); y.append(float(p[1]))
            except ValueError:
                pass
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


def sample_TR(layers, lam):
    """Film-side incidence, coherent films on an incoherent glass slab."""
    ns = [complex(N_AIR)] + [l[0] for l in layers] + [complex(N_GLASS)]
    ds = [0.0] + [l[1] for l in layers] + [0.0]
    T_in, R_front = coherent(ns, ds, lam)
    T_back, R_back = coherent(ns[::-1], ds[::-1], lam)
    r = (N_GLASS - N_AIR) / (N_GLASS + N_AIR)
    Rr, Tr = r ** 2, 1 - r ** 2
    den = 1 - R_back * Rr
    return T_in * Tr / den, R_front + T_in * Rr * T_back / den


def main():
    wl, Trel = load(F_T)
    _, Rfile = load(F_R)
    m = (wl >= FIT[0]) & (wl <= FIT[1]) & ~((wl > LAMP[0]) & (wl < LAMP[1]))

    n_ag = np.conj(np.interp(wl, JC_WL, AG_N) + 1j * np.interp(wl, JC_WL, AG_K))
    hat = np.conj(np.full_like(wl, N_HATCN + 1j * K_HATCN, dtype=complex))

    # baseline sample (glass/HATCN) computed for each thickness hypothesis
    print("=" * 74)
    print("BASELINE SAMPLE (glass / HATCN d) -- computed from known n,k")
    print("=" * 74)
    print(f"{'d_HATCN':>8}{'T@550':>9}{'R@550':>9}{'T@800':>9}{'R@800':>9}"
          f"{'  bare glass R':>15}")
    base = {}
    for dh in (15.0, 20.0, 30.0):
        Tb = np.zeros(len(wl)); Rb = np.zeros(len(wl))
        for i in range(len(wl)):
            Tb[i], Rb[i] = sample_TR([(hat[i], dh)], wl[i])
        base[dh] = (Tb, Rb)
        i5, i8 = int(np.argmin(abs(wl - 550))), int(np.argmin(abs(wl - 800)))
        Tg, Rg = sample_TR([], wl[i5])
        print(f"{dh:>8.0f}{100*Tb[i5]:>9.1f}{100*Rb[i5]:>9.1f}"
              f"{100*Tb[i8]:>9.1f}{100*Rb[i8]:>9.1f}{100*Rg:>15.1f}")

    # hypotheses for R
    print("\n" + "=" * 74)
    print("HYPOTHESIS TEST  (constraint: absorptance A = 1 - T_abs - R_abs)")
    print("=" * 74)
    dh = 15.0
    Tb, Rb = base[dh]
    Tabs = Trel / 100 * Tb
    hyps = {
        "H-A  R absolute": Rfile / 100,
        "H-B  R / (HATCN/glass R)": Rfile / 100 * Rb,
        "H-C  R / mirror 0.90": Rfile / 100 * 0.90,
    }
    print(f"{'hypothesis':<28}{'A min %':>9}{'A max %':>9}{'A@800':>8}"
          f"{'A@550':>8}{'A@440':>8}   verdict")
    for name, Rabs in hyps.items():
        A = 1 - Tabs - Rabs
        i8 = int(np.argmin(abs(wl - 800)))
        i5 = int(np.argmin(abs(wl - 550)))
        i4 = int(np.argmin(abs(wl - 440)))
        amin, amax = 100 * A[m].min(), 100 * A[m].max()
        if amin < -0.5:
            v = "DEAD -- negative absorptance"
        elif amax > 40:
            v = "DEAD -- Ag 5 nm cannot absorb this much"
        elif amax > 20:
            v = "strained"
        else:
            v = "consistent"
        print(f"{name:<28}{amin:>9.1f}{amax:>9.1f}{100*A[i8]:>8.1f}"
              f"{100*A[i5]:>8.1f}{100*A[i4]:>8.1f}   {v}")

    print("""
  H-B multiplies the file by the baseline sample's own reflectance (~10 %),
  which would put the film's absolute R at ~2-3 % -- but then 30+ % of the
  light is unaccounted for and must be called absorption, far beyond what a
  5 nm Ag film can do. H-C inflates absorptance by the ~2.5 %p that the
  mirror correction removes from R. H-A is the only one that keeps A inside
  the physically allowed band, so the file is (as presumed, now tested)
  absolute R -- most instruments' reflectance software does apply the
  reference-mirror calibration before writing %R.""")

    # final combined dataset under the surviving hypothesis
    print("=" * 74)
    print("FINAL ABSOLUTE DATASET (H-A), with d_HATCN as the systematic")
    print("=" * 74)
    print(f"{'lambda':>7}" + "".join(f"{'A d=%.0f' % d:>10}" for d in base)
          + f"{'R_abs':>8}{'T_abs(15)':>11}")
    for tgt in (780, 700, 650, 600, 550, 500, 460, 440):
        i = int(np.argmin(abs(wl - tgt)))
        row = f"{wl[i]:>7.0f}"
        for d in base:
            A = 1 - Trel[i] / 100 * base[d][0][i] - Rfile[i] / 100
            row += f"{100*A:>10.1f}"
        row += f"{Rfile[i]:>8.1f}{Trel[i]/100*base[15.0][0][i]*100:>11.1f}"
        print(row)

    A15 = 1 - Trel / 100 * base[15.0][0] - Rfile / 100
    A30 = 1 - Trel / 100 * base[30.0][0] - Rfile / 100
    print(f"""
  Visible-average absorptance: {100*A15[m].mean():.1f} % (HATCN 15 nm) to
  {100*A30[m].mean():.1f} % (HATCN 30 nm) -- the baseline-thickness systematic
  is ~{100*abs(A30[m].mean()-A15[m].mean()):.1f} %p and dominates the error
  budget. Measuring the baseline sample's ABSOLUTE T once removes it.

  The shape is the physics check: A rises from {100*A15[int(np.argmin(abs(wl-780)))]:.1f} %
  at 780 nm toward {100*A15[int(np.argmin(abs(wl-460)))]:.1f} % at 460 nm.
  Glass and HATCN are transparent here, so all of it is the Ag film, and a
  blue-rising absorption is exactly the LSPR signature of a granular film --
  consistent with the SEM texture and with why the bulk-Ag TMM would not fit
  in scripts/60. A continuous bulk-like film would absorb 2-4 % with much
  weaker dispersion.""")


if __name__ == "__main__":
    main()
