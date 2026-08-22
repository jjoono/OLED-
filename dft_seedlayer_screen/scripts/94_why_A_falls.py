"""Why the measured absorptance falls as the silver gets thicker.

A = 1 - T - R is the absorptance of the whole glass/seed/Ag/air structure under
external illumination, and across the HATCN series it falls monotonically,
14.59% at 4 nm down to 6.27% at 12 nm.  That is not the thick film absorbing
less per pass than a thin one at the same quality -- it is the thin films being
far worse silver.

Two things change together with thickness.  Geometry: a thicker metal reflects
more, so less light is admitted (1-R drops 85.7% -> 60.5%).  Film quality: the
Fuchs-Sondheimer and grain-boundary terms collapse as the film closes, and the
inverted eps2 falls 6.02 -> 1.04 against a bulk 0.32.

Geometry alone does not explain the ordering.  Held at bulk optical constants,
the external absorptance of this structure RISES over 4-12 nm and peaks near
10 nm.  The measured fall is therefore the size effect, and it is strong enough
to overturn the geometric trend.

The device asks a different question again.  There the electrode sits between
organic (n ~ 1.8) and capping (n ~ 2.1), the front reflection is largely gone,
and the one-pass absorption is set by eps2 x d rather than by how much light
gets in.  That quantity bottoms out at 8 nm, and 7-12 nm are all within 0.15 %p.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
m = importlib.import_module("91_nk_library_series")

L = 550.0
D = [4, 5, 6, 7, 8, 10, 12]
IDS = {4: "1-2", 5: "1-3", 6: "1-4", 7: "2-1", 8: "2-2", 10: "2-3", 12: "2-4"}


def measured(sid):
    T = m.read_csv(os.path.join(m.RAW, f"{sid}T.csv"))
    R = m.read_csv(os.path.join(m.RAW, f"{sid}R.csv"))
    return T[L], R[L]                     # read_csv already returns fractions


def main():
    ag = m.AG_IDEAL(L)
    e2_bulk = 2 * ag.real * ag.imag

    print(f"HATCN 5 nm / Ag d, at {L:.0f} nm.  Bulk silver (McPeak): "
          f"n={ag.real:.4f}, k={ag.imag:.4f}, eps2={e2_bulk:.3f}\n")
    print(f"{'d':>4} {'T':>7} {'R_corr':>7} {'A':>7} | {'1-R':>7} {'A/(1-R)':>8} | "
          f"{'n':>6} {'k':>6} {'eps2':>6} {'/bulk':>6} | {'A bulk nk':>9} | {'device':>7}")
    print("-" * 100)
    rows = []
    for d in D:
        T, R = measured(IDS[d])
        Rb = ((m.n_glass(L) - 1) / (m.n_glass(L) + 1))**2
        Rc = R + (1 - m.KEEP) * (T**2) * Rb          # for reporting only
        A = 1 - T - Rc
        # invert against the RAW R: observable() already carries the 0.843
        # back-surface factor, so feeding it a corrected R counts it twice
        nk, _ = m.invert(L, T, R, "HATCN", float(d))
        e2 = 2 * nk.real * nk.imag
        Tb, Rbk = m.observable(ag, L, "HATCN", float(d))   # same stack, bulk silver
        dev, _ = m.best_device_A(nk, L, float(d))
        rows.append((d, T, Rc, A, nk, e2, 1 - Tb - Rbk, dev))
        print(f"{d:>4} {T*100:>6.2f}% {Rc*100:>6.2f}% {A*100:>6.2f}% | "
              f"{(1-Rc)*100:>6.2f}% {A/(1-Rc)*100:>7.2f}% | "
              f"{nk.real:>6.3f} {nk.imag:>6.3f} {e2:>6.3f} {e2/e2_bulk:>5.1f}x | "
              f"{(1-Tb-Rbk)*100:>8.2f}% | {dev*100:>6.2f}%")

    a0, a1 = rows[0], rows[-1]
    print(f"\nMeasured A                 {a0[3]*100:5.2f}% -> {a1[3]*100:5.2f}%   "
          f"(4 -> 12 nm)   falls")
    print(f"Same stack, bulk n and k   {a0[6]*100:5.2f}% -> {a1[6]*100:5.2f}%   "
          f"               RISES")
    print("Geometry alone predicts the opposite ordering. What actually drives the")
    print(f"fall is film quality: eps2 goes {a0[5]:.2f} -> {a1[5]:.2f}, "
          f"{a0[5]/e2_bulk:.0f}x -> {a1[5]/e2_bulk:.0f}x bulk.")
    print("Even the 12 nm film is still 3x worse than bulk silver.\n")

    grid = np.arange(0.5, 40.01, 0.5)
    Aext = np.array([1 - sum(m.observable(ag, L, "HATCN", float(d))) for d in grid])
    i = int(np.argmax(Aext))
    print(f"With bulk silver held fixed, external A peaks at d = {grid[i]:.1f} nm "
          f"({Aext[i]*100:.2f}%) and falls to {Aext[-1]*100:.2f}% at 40 nm --")
    print("the familiar metal-film absorptance maximum, which the real series never")
    print("reaches because the size effect swamps it.\n")

    dev = [r[7] for r in rows]
    j = int(np.argmin(dev))
    print("Device one-pass absorption (organic 1.8 / Ag / CPL 2.1, CPL optimised):")
    print("  " + "   ".join(f"{r[0]}nm {r[7]*100:.2f}%" for r in rows))
    print(f"  minimum at {D[j]} nm ({dev[j]*100:.2f}%); 7-12 nm all within "
          f"{(max(dev[3:])-min(dev[3:]))*100:.2f} %p.")
    print("  This is the number that matters for the device, and it is NOT the")
    print("  measured A -- different question, different ordering below 7 nm.")


if __name__ == "__main__":
    main()
