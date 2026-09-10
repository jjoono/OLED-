"""Does the ellipsometric n,k reproduce the measured T and R? Three checks.

Ellipsometry on the 7 nm film returns n = 0.08, k = 3.5 at 550 nm. The T/R
inversion on the same sample returns n = 0.282, k = 3.479. Those are not small
differences: eps2 = 2nk is 0.56 against 1.96, a factor of 3.5 in absorption, and
a device model fed one or the other will disagree by that factor.

Ellipsometry and 1-T-R do not measure the same thing, so a disagreement is
expected in one direction. Ellipsometry reads psi and delta, the RATIO
r_p / r_s of the specularly reflected beam. Anything that removes light from the
specular direction without changing that ratio -- roughness scattering, haze,
a wedge -- is invisible to it, and any absolute-intensity error cancels. So the
n and k it returns describe the specular film and can legitimately sit below
what 1-T-R implies.

But it cannot sit below the physics. Three checks, in increasing severity:

  1. Feed the ellipsometric n,k through the same transfer matrix and predict T
     and R for the actual measured stack. Compare against the measurement.
  2. Compare eps2 against the resistivity. Sheet resistance is a DC measurement
     that cannot scatter, and grain boundaries load it while staying invisible
     at optical frequency, so the DC value is an upper bound on optical damping.
  3. Compare eps2 against the Fuchs-Sondheimer surface-scattering floor with
     p = 0. That floor is not a fit; it is what electrons bouncing diffusely off
     two surfaces cost, and no film of that thickness can be under it.

Check 3 is the one that decides, because it does not depend on the resistivity
measurement or on the T/R measurement at all -- only on the thickness.
"""
import csv, math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m = __import__("118_zns_cap_compare")

L = 550.0
D_AG, D_SEED = 7.0, 5.0
RS = 12.7                       # Ohm/sq, measured on HATCN 5 / Ag 7
RHO_BULK, MFP = 1.59, 52.0
KEEP = 0.843                    # back-surface collection factor
ELL = complex(0.08, 3.50)       # ellipsometry, this work
INV = complex(0.2820, 3.4791)   # T/R inversion, same sample
MCPEAK = complex(0.0438, 3.6101)


def n_glass(l):
    return 1.5220 + 3900.0 / l**2


def stack_TR(nk, lam):
    """T and R of glass / HATCN 5 / Ag d / air, illuminated through the glass,
    with the substrate's incoherent back surface handled as in the measurement."""
    hl, hn, hk = m.load("HATCN")
    ns = complex(m.interp(lam, hl, hn), m.interp(lam, hl, hk))
    ng = complex(n_glass(lam), 0)
    Rb = ((n_glass(lam) - 1) / (n_glass(lam) + 1))**2
    front = [ng, ns, nk, complex(1, 0)]
    d_f = [0.0, D_SEED, D_AG, 0.0]
    R1, T1 = m.tmm(front, d_f, lam, 0.0)
    back = [complex(1, 0), nk, ns, ng]
    d_b = [0.0, D_AG, D_SEED, 0.0]
    R1b, T1b = m.tmm(back, d_b, lam, 0.0)
    den = 1 - Rb * R1b
    T = T1 * (1 - Rb) / den
    R = Rb + (1 - Rb)**2 * R1 / den
    return T, R


def measured():
    raw = os.path.join(m.BASE, "data", "TR_20260820", "raw")
    out = {}
    for tag, fn in (("T", "2-1T.csv"), ("R", "2-1R.csv")):
        p = os.path.join(raw, fn)
        if not os.path.exists(p):
            return None
        for row in csv.reader(open(p)):
            try:
                if abs(float(row[0]) - L) < 0.6:
                    out[tag] = float(row[1]) / 100.0
            except (ValueError, IndexError):
                pass
    return out


def main():
    print(f"HATCN 5 / Ag {D_AG:.0f} nm at {L:.0f} nm\n")
    print(f"{'source':<22} {'n':>7} {'k':>7} {'eps2':>7} {'x bulk':>8}")
    for tag, nk in (("ellipsometry", ELL), ("T/R inversion", INV),
                    ("bulk silver (McPeak)", MCPEAK)):
        e2 = 2 * nk.real * nk.imag
        print(f"{tag:<22} {nk.real:>7.3f} {nk.imag:>7.3f} {e2:>7.3f} "
              f"{e2/(2*MCPEAK.real*MCPEAK.imag):>7.2f}x")

    print("\n--- check 1: predicted vs measured T and R ---\n")
    meas = measured()
    print(f"{'source':<22} {'T':>8} {'R':>8} {'1-T-R':>8}")
    for tag, nk in (("ellipsometry", ELL), ("T/R inversion", INV)):
        T, R = stack_TR(nk, L)
        print(f"{tag:<22} {T*100:>7.2f}% {R*100:>7.2f}% {(1-T-R)*100:>7.2f}%")
    if meas:
        Rc = meas["R"] + (1 - KEEP) * meas["T"]**2 * \
             ((n_glass(L) - 1) / (n_glass(L) + 1))**2
        print(f"{'MEASURED':<22} {meas['T']*100:>7.2f}% {Rc*100:>7.2f}% "
              f"{(1-meas['T']-Rc)*100:>7.2f}%")
        Te, Re = stack_TR(ELL, L)
        print(f"\n  ellipsometry misses T by {(Te-meas['T'])*100:+.2f} %p and "
              f"R by {(Re-Rc)*100:+.2f} %p")
        print(f"  and puts the absorptance at {(1-Te-Re)*100:.2f}% against a "
              f"measured {(1-meas['T']-Rc)*100:.2f}%")

    print("\n--- check 2: against the resistivity ---\n")
    rho = RS * D_AG * 0.1
    e2_dc = 2 * MCPEAK.real * MCPEAK.imag * rho / RHO_BULK
    print(f"  Rs {RS} Ohm/sq at {D_AG:.0f} nm -> rho {rho:.2f} uOhm cm, "
          f"{rho/RHO_BULK:.2f}x bulk")
    print(f"  DC-implied eps2 = {e2_dc:.3f}   (an upper bound: grain boundaries")
    print(f"  load the DC path but are invisible at 550 nm)")
    print(f"  ellipsometry {2*ELL.real*ELL.imag:.3f}, inversion "
          f"{2*INV.real*INV.imag:.3f}")

    print("\n--- check 3: against the Fuchs-Sondheimer floor ---\n")
    rho_surf = RHO_BULK * (1 + 0.375 * MFP / D_AG)
    e2_floor = 2 * MCPEAK.real * MCPEAK.imag * rho_surf / RHO_BULK
    n_floor = e2_floor / (2 * 3.5)
    print(f"  diffuse surfaces, p = 0, no grain boundaries at all:")
    print(f"    rho_floor = {rho_surf:.2f} uOhm cm = {rho_surf/RHO_BULK:.2f}x bulk")
    print(f"    eps2_floor = {e2_floor:.3f}, i.e. n >= {n_floor:.3f} at k = 3.5")
    for tag, nk in (("ellipsometry", ELL), ("T/R inversion", INV)):
        e2 = 2 * nk.real * nk.imag
        v = "BELOW THE FLOOR" if e2 < e2_floor else "above the floor, physical"
        print(f"    {tag:<16} eps2 {e2:.3f}  {v}")
    print(f"\n  A 7 nm film cannot damp less than its own surfaces do. Reaching")
    print(f"  n = {ELL.real:.2f} needs specularity p = "
          f"{1 - (2*ELL.real*ELL.imag/(2*MCPEAK.real*MCPEAK.imag) - 1)/(0.375*MFP/D_AG):.2f}, "
          f"which by Soffer")
    sig = 0.52 / (4 * math.pi) * math.sqrt(
        -math.log(max(1e-9, 1 - (2*ELL.real*ELL.imag /
                                 (2*MCPEAK.real*MCPEAK.imag) - 1) / (0.375*MFP/D_AG))))
    print(f"  is {sig:.3f} nm RMS on both faces -- an epitaxial film, not one")
    print(f"  evaporated onto an amorphous organic.")


if __name__ == "__main__":
    main()
