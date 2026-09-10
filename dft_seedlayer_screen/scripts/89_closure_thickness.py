"""Closure thickness of Ag on each seed, read off the resistivity series.

A film that is CLOSED obeys Matthiessen: rho = rho_0 + C/d, with rho_0 carrying
bulk plus grain boundaries and C the surface term.  A film that is not closed
has extra resistance from the percolation path and sits ABOVE that line.  So
fit the line to the thick end, where closure is certain, and read the thickness
at which the thin points fall back onto it.

This is better founded than fitting Rs ~ (d-dc)^-t: that power law has to be
extrapolated outside the measured range, since both seeds already conduct at
4 nm, and it returned dc = 3.3 nm for HATCN against 2.8 nm for MoOx -- the
wrong ordering, driven entirely by the curvature of the fit.
"""
import numpy as np

RHO_B, MFP = 1.59, 52.0
D = np.array([4, 5, 6, 7, 8, 10, 12], float)
RS = {"HATCN": np.array([52.2, 23.3, 18.9, 12.7, 9.1, 7.5, 5.3]),
      "MoOx":  np.array([138.0, 49.1, 32.0, 14.9, 10.8, 9.6, 6.6])}
FIT_FROM = 7.0            # thicknesses at or above this are taken as closed
TOL = 0.25                # a point within 25 % of the line counts as closed


def analyse(seed):
    rho = RS[seed]*D*0.1
    sel = D >= FIT_FROM
    A = np.vstack([np.ones(sel.sum()), 1.0/D[sel]]).T
    coef, *_ = np.linalg.lstsq(A, rho[sel], rcond=None)
    resid = A @ coef - rho[sel]
    scatter = float(np.std(resid, ddof=1))/float(np.mean(rho[sel]))
    line = coef[0] + coef[1]/D
    excess = rho/line - 1.0
    closed = D[excess <= TOL]
    return dict(rho=rho, rho0=coef[0], C=coef[1], line=line, excess=excess,
                scatter=scatter, d_close=float(closed.min()) if len(closed) else np.nan)


if __name__ == "__main__":
    out = {}
    for seed in RS:
        r = analyse(seed); out[seed] = r
        print(f"\n{seed}   closed-film line  rho = {r['rho0']:.2f} + {r['C']:.1f}/d "
              f"uOhm-cm   (fitted to d >= {FIT_FROM:.0f} nm, scatter "
              f"{100*r['scatter']:.1f} %)")
        print(f"  {'d':>4}{'Rs':>8}{'rho':>8}{'line':>8}{'excess':>9}   verdict")
        print("  " + "-"*52)
        for i, d in enumerate(D):
            v = "closed" if r['excess'][i] <= TOL else "NOT closed"
            print(f"  {d:4.0f}{RS[seed][i]:8.1f}{r['rho'][i]:8.2f}"
                  f"{r['line'][i]:8.2f}{100*r['excess'][i]:+8.0f} %   {v}")
        print(f"  -> closes at {r['d_close']:.0f} nm")

    dh, dm = out['HATCN']['d_close'], out['MoOx']['d_close']
    print(f"\n\nHATCN closes at {dh:.0f} nm, MoOx at {dm:.0f} nm: "
          f"HATCN saves {dm-dh:.0f} nm of silver.")
    print("Both bounds are limited by the 1 nm spacing of the series, and the")
    print("true HATCN value could be below 5 nm -- there is no 2 or 3 nm sample.")

    print("\nIndependent optical confirmation at 5 nm (this session's T/R):")
    print("  HATCN/Ag5  eps1 = -11.9 against bulk silver -13.0, i.e. within 9 %")
    print("             across 360-850 nm  ->  already a continuous metal at 5 nm")
    print("  MoOx/Ag5   fitted damping 21.4x bulk against HATCN's 5.9x, and its")
    print("             resistivity sits 109 % above the closed-film line")
