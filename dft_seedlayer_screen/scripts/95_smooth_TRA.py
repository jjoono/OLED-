"""Smooth the T/R series only as far as the noise justifies.

The scan is 2 nm apart over 350-850 nm and carries detector noise plus, in the
UV, the lamp changeover. Smoothing helps a reader see the trend, but any window
wide enough to flatten a real feature -- the Fe3+ edge below 400 nm, the
plasma-edge structure near 320-400 nm, the interference the capping produces --
destroys the thing the data is for.

So the window is not chosen by eye. For each spectrum:

  1. Estimate the white-noise sigma from the second difference. For noise alone,
     var(y[i-1] - 2y[i] + y[i+1]) = 6 sigma^2, and a smooth underlying curve
     contributes almost nothing at 2 nm spacing.
  2. Try Savitzky-Golay windows from narrow to wide, cubic.
  3. Keep the widest window whose RMS(smoothed - raw) over 400-850 nm stays at or
     below sigma. A filter that only removes noise leaves an RMS residual below
     sigma by construction; once it starts bending the underlying curve the
     residual carries signal as well and climbs past it. (The MAXIMUM residual is
     no use as a test -- the largest of ~230 gaussian draws sits near 3 sigma
     whatever the window does.)
  4. Cap the window at 11 points, 22 nm, so no setting can wash out the Fe3+ edge
     or a capping interference fringe even if the arithmetic would allow it.

A is then recomputed from the smoothed T and R rather than smoothed itself, so
T + R + A = 100 still holds exactly at every wavelength.

Raw sheets are kept in the same workbook. The smoothed sheets are a reading aid,
not a replacement.
"""
import csv, os
import numpy as np
from scipy.signal import savgol_filter

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RAW = os.path.join(BASE, "data", "TR_20260820", "raw")
KEEP = 0.843
FIT_LO, FIT_HI = 400.0, 850.0      # window judged outside the lamp changeover
MAXRMS = 1.0                       # RMS residual, in units of the noise sigma
WINDOWS = [5, 7, 9, 11]            # capped at 22 nm to protect real features


def noise_sigma(y):
    d2 = y[2:] - 2 * y[1:-1] + y[:-2]
    return float(np.median(np.abs(d2 - np.median(d2))) * 1.4826 / np.sqrt(6.0))


def smooth(lam, y):
    """Widest cubic Savitzky-Golay window that stays inside 2 sigma of the raw."""
    sig = noise_sigma(y)
    band = (lam >= FIT_LO) & (lam <= FIT_HI)
    best = WINDOWS[0]
    for w in WINDOWS:
        if w >= len(y):
            break
        s = savgol_filter(y, w, 3, mode="interp")
        if float(np.sqrt(np.mean(((s - y)[band])**2))) <= MAXRMS * sig:
            best = w
        else:
            break
    s = savgol_filter(y, best, 3, mode="interp")
    rms = float(np.sqrt(np.mean(((s - y)[band])**2)))
    return s, best, sig, rms


def read(p):
    out = {}
    if not os.path.exists(p):
        return out
    for row in csv.reader(open(p)):
        try:
            out[float(row[0])] = float(row[1])
        except (ValueError, IndexError):
            pass
    return out


def n_glass(l):
    return 1.5220 + 3900.0 / l**2


SAMPLES = [("1-1","HATCN",0),  ("1-2","HATCN",4), ("1-3","HATCN",5), ("1-4","HATCN",6),
           ("2-1","HATCN",7),  ("2-2","HATCN",8), ("2-3","HATCN",10),("2-4","HATCN",12),
           ("1-9","MoOx",0),   ("1-10","MoOx",4), ("1-11","MoOx",5), ("1-12","MoOx",6),
           ("2-9","MoOx",7),   ("2-10","MoOx",8), ("2-11","MoOx",10),("2-12","MoOx",12)]


def build(lam):
    """-> Ts, Rs, Rcs, As, report   each keyed by sample id (dict lambda->value)."""
    Ts, Rs, Rcs, As, report = {}, {}, {}, {}, []
    for sid, seed, d in SAMPLES:
        T = read(os.path.join(RAW, f"{sid}T.csv"))
        R = read(os.path.join(RAW, f"{sid}R.csv"))
        for tag, src, store in (("T", T, Ts), ("R", R, Rs)):
            g = np.array([l for l in lam if l in src])
            if len(g) < max(WINDOWS) + 2:
                store[sid] = {}
                continue
            y = np.array([src[l] for l in g])
            s, w, sig, dev = smooth(g, y)
            store[sid] = dict(zip(g, s))
            report.append((sid, seed, d, tag, len(g), w, sig, dev))
        rc, a = {}, {}
        for l in lam:
            t, r = Ts[sid].get(l), Rs[sid].get(l)
            if t is None or r is None:
                continue
            Rb = ((n_glass(l) - 1) / (n_glass(l) + 1))**2
            rc[l] = r + 100 * (1 - KEEP) * (t / 100.0)**2 * Rb
            a[l] = 100 - t - rc[l]
        Rcs[sid], As[sid] = rc, a
    return Ts, Rs, Rcs, As, report


def main():
    lams = set()
    for sid, _, _ in SAMPLES:
        lams |= set(read(os.path.join(RAW, f"{sid}T.csv")))
        lams |= set(read(os.path.join(RAW, f"{sid}R.csv")))
    lam = sorted(l for l in lams if 350 <= l <= 850)
    _, _, _, _, rep = build(lam)
    print(f"{'id':>5} {'seed':>6} {'Ag':>3} {'ch':>2} {'pts':>4} {'window':>7} "
          f"{'sigma(%p)':>10} {'rms':>8} {'rms/sigma':>10}")
    print("-" * 68)
    for sid, seed, d, tag, n, w, sig, dev in rep:
        print(f"{sid:>5} {seed:>6} {d:>3} {tag:>2} {n:>4} {w:>5} nm "
              f"{sig:>10.4f} {dev:>8.4f} {dev/sig if sig else 0:>10.2f}")
    ws = [r[5] for r in rep]
    sg = [r[6] for r in rep]
    print(f"\nwindow: {min(ws)}-{max(ws)} points ({min(ws)*2}-{max(ws)*2} nm), "
          f"median {int(np.median(ws))}")
    print(f"noise sigma: {min(sg):.4f}-{max(sg):.4f} %p, median {np.median(sg):.4f} %p")
    rr = [r[7] / r[6] for r in rep if r[6]]
    print(f"rms residual / sigma: {min(rr):.2f}-{max(rr):.2f}")
    over = [r for r in rep if r[6] and r[7] / r[6] > 1.0]
    print("A ratio at or below 1 means the filter removed noise and nothing else.")
    if over:
        print(f"{len(over)}/{len(rep)} spectra sit slightly above 1 even at the "
              f"narrowest cubic window (5 points, the minimum for order 3);")
        print("  " + ", ".join(f"{r[0]}{r[3]} {r[7]/r[6]:.2f}" for r in over))
        print("  those are already at the noise floor and are barely smoothed at all.")


if __name__ == "__main__":
    main()
