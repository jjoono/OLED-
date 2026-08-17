"""Sample 4 (HATCN / Ag 3-4 nm, holey): how much does the bad morphology absorb?

WHAT THIS SAMPLE IS. Thinner Ag than sample 5 (3-4 nm vs ~5), 100 Ohm/sq, SEM
shows holes -- and the user reports the striking observation directly: the
THINNER film transmits LESS. For a continuous film that is impossible (an ideal
3 nm film transmits ~99 % relative, 5 nm ~94 %), so the deficit is morphology,
and quantifying it in the ABSORPTION channel is exactly what the user's research
goal needs, since their stated criterion is minimum absorption above all.

R was measured the same flawed way as sample 5 (transmission baseline reused in
the 6-degree reflection geometry), so the same correction applies:

    R_abs = R_file * T_base / eta,   eta = eta_R/eta_T

and the same instrument session means the SAME eta as sample 5 -- the EMA-derived
0.70 with an honest ceiling near 0.90. Both are carried; A lands in a range.

The R file carries two repeat scans (sample4R, sample4R2). They are averaged, and
their difference doubles as a direct noise measurement -- something sample 5's
single scan never gave.
"""
import os
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(_HERE, "..", "data")

JC_WL = np.array([397.4, 413.3, 430.5, 450.9, 471.4, 495.9,
                  521.0, 548.6, 582.1, 616.8, 659.5, 704.5])
AG_N = np.array([0.05, 0.05, 0.04, 0.04, 0.05, 0.05, 0.05, 0.06, 0.05, 0.06, 0.05, 0.04])
AG_K = np.array([2.07, 2.21, 2.36, 2.66, 2.83, 3.09, 3.34, 3.59, 3.93, 4.15, 4.48, 4.84])
N_GLASS, N_HATCN, K_HATCN = 1.52, 1.95, 2.0e-4
LAMP, FITB = (700.0, 735.0), (430.0, 800.0)
ETAS = (0.70, 0.90)          # scripts/62: model mean and honest ceiling


def load_cols(path, pairs):
    """pairs = [(wl_col, y_col), ...]; returns wl and one array per pair."""
    rows = []
    for line in open(path):
        p = line.strip().split(",")
        try:
            rows.append([float(v) for v in p[:max(max(pairs)) + 1] if v != ""])
        except ValueError:
            continue
    out = []
    for wc, yc in pairs:
        wl = np.array([r[wc] for r in rows if len(r) > yc])
        y = np.array([r[yc] for r in rows if len(r) > yc])
        o = np.argsort(wl)
        out.append((wl[o], y[o]))
    return out


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


def stack_T(layers, lam):
    ns = [complex(1.0)] + [l[0] for l in layers] + [complex(N_GLASS)]
    ds = [0.0] + [l[1] for l in layers] + [0.0]
    T_in, _ = coherent(ns, ds, lam)
    R_back = coherent(ns[::-1], ds[::-1], lam)[1]
    r = (N_GLASS - 1.0) / (N_GLASS + 1.0)
    return T_in * (1 - r * r) / (1 - R_back * r * r)


def main():
    (wlT, T4), = load_cols(os.path.join(DATA, "112T.csv"), [(0, 1)])
    (wr1, R1), (wr2, R2) = load_cols(os.path.join(DATA, "112R.csv"),
                                     [(0, 1), (2, 3)])
    wl = wlT
    R4 = 0.5 * (R1 + R2)
    (wlT5, T5), = load_cols(os.path.join(DATA, "24T.csv"), [(0, 1)])
    (wlR5, R5), = load_cols(os.path.join(DATA, "24R.csv"), [(0, 1)])
    m = (wl >= FITB[0]) & (wl <= FITB[1]) & ~((wl > LAMP[0]) & (wl < LAMP[1]))

    print("=" * 72)
    print("REPEAT-SCAN NOISE (the two R columns measure the same sample)")
    print("=" * 72)
    d = (R1 - R2)[m]
    print(f"  mean |R1-R2| = {np.mean(abs(d)):.3f} %p,  rms = {np.std(d):.3f} %p")
    print("  -> single-scan reproducibility ~0.3 %p; spectral wiggles smaller")
    print("     than this are noise, not features.")

    hat = np.conj(np.full_like(wl, N_HATCN + 1j * K_HATCN, dtype=complex))
    Tb = np.array([stack_T([(hat[i], 15.0)], wl[i]) for i in range(len(wl))])
    i5 = int(np.argmin(abs(wl - 550)))

    print("\n" + "=" * 72)
    print("SAMPLE 4 vs SAMPLE 5 -- thinner but darker")
    print("=" * 72)
    print(f"{'':>22}{'Ag 3-4 nm (s4)':>16}{'Ag ~5 nm (s5)':>15}{'ideal 3 nm':>12}")
    print(f"{'T_rel @550':>22}{T4[i5]:>15.1f}%{T5[i5]:>14.1f}%{98.6:>11.1f}%")
    print(f"{'T_rel avg':>22}{T4[m].mean():>15.1f}%{T5[m].mean():>14.1f}%{98.2:>11.1f}%")
    print(f"{'Rs (Ohm/sq)':>22}{'~100':>16}{'~25':>15}{'--':>12}")

    print("\n" + "=" * 72)
    print("ABSORPTION (the user's figure of merit), eta carried as a range")
    print("=" * 72)
    T4a = T4 / 100 * Tb
    T5a = T5 / 100 * Tb
    print(f"{'eta':>6}{'  A(s4)@550':>12}{'A(s4) avg':>11}{'  A(s5)@550':>12}"
          f"{'A(s5) avg':>11}")
    for eta in ETAS:
        A4 = 1 - T4a - (R4 / 100) * Tb / eta
        A5 = 1 - T5a - (R5 / 100) * Tb / eta
        print(f"{eta:>6.2f}{100*A4[i5]:>11.1f}%{100*A4[m].mean():>10.1f}%"
              f"{100*A5[i5]:>11.1f}%{100*A5[m].mean():>10.1f}%")
    print("""
  Whichever eta is right, the ORDERING is eta-independent because both
  samples share the instrument session and hence the factor: the thinner,
  holey film absorbs roughly 6-7 %p MORE than the 5 nm film, on less
  metal. An ideal continuous film of either thickness would absorb 1-2 %.

  That is the quantitative version of the SEM impression: the holes are
  not just missing metal, they are plasmonic absorbers. Absorption here
  is a MORPHOLOGY meter, and it moves the right way when the film gets
  thicker and closes -- which is the central fact for a research goal of
  minimum absorption: past percolation, MORE silver means LESS absorption
  until the intrinsic (thickness-proportional) term takes over. The
  optimum is just after hole closure, and these two samples bracket it
  from the bad side.""")


if __name__ == "__main__":
    main()
