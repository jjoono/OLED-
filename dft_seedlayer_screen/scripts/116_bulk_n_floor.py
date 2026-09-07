"""Can a sub-10 nm silver film have bulk optical constants? What it would take.

The target: n ~ 0.04 -- bulk silver at 550 nm, eps2 = 0.32 -- in a film under
10 nm, giving one-pass absorption under 1 %. Three numbers decide it.

  1. The ideal film. With bulk constants held fixed, one-pass absorption in the
     device rises with thickness and crosses 1 % near 17 nm. So the target is
     consistent with thickness alone: an 8 nm film of bulk silver would sit at
     0.52 %.

  2. The floor with diffuse surfaces. Even a film with no grain boundaries at
     all carries the Fuchs-Sondheimer surface term, rho_bulk * 0.375 (1-p) l/d,
     and with p = 0 that alone puts an 8 nm film at 3.4x bulk damping and 1.77 %
     one-pass. The floor never dips below 1.5 % at any thickness. This is not a
     film-quality problem; it is the electrons hitting the surfaces.

  3. The specularity that would be needed. Solving for p such that the 8 nm
     film reaches 1.0 % gives p = 0.62, and Soffer's expression turns that into
     an interface roughness of 0.03 nm RMS -- on BOTH faces. The lower face is
     the seed layer's surface, and an amorphous organic film is rough at the
     scale of its own molecules, several angstroms. The seed forbids it before
     the silver is even deposited.

So the answer to "can it be done" splits cleanly. On a crystalline substrate by
epitaxy, films approaching this exist in the literature. As a top electrode on
an organic stack, at room temperature, no -- and no smoothing recipe changes
that, because the limit is set by the interface the silver is grown on.
"""
import math

RHO_BULK, MFP, LAMBDA_F = 1.59, 52.0, 0.52
# one-pass device absorption with bulk (McPeak) constants, from script 91
IDEAL = {7: 0.460, 8: 0.522, 10: 0.643, 12: 0.758, 15: 0.919, 20: 1.153}


def ratio(d, p):
    return 1.0 + 0.375 * (1.0 - p) * MFP / d


def main():
    print("1. ideal film, bulk constants:  one-pass A vs thickness")
    for d, a in IDEAL.items():
        print(f"   {d:>3} nm  {a:5.2f} %")
    print("   -> crosses 1 % near 17 nm; the target is compatible with d <= 15 nm\n")

    print("2. floor with diffuse surfaces (p = 0), no grain boundaries")
    print(f"   {'d':>4} {'rho/bulk':>9} {'eps2':>6} {'n':>6} {'one-pass A':>11}")
    for d, a in IDEAL.items():
        r = ratio(d, 0.0)
        e2 = 0.316 * r
        n = e2 / (2 * 3.53)                    # eps1 ~ -12.4 -> k ~ 3.53
        print(f"   {d:>4} {r:>9.2f} {e2:>6.2f} {n:>6.3f} {a*r:>10.2f} %")
    print("   -> never below 1.5 %. Measured 8 nm on HATCN: n = 0.217, A = 2.45 %,")
    print("      i.e. within 1.4x of this floor already.\n")

    d, a = 8, IDEAL[8]
    need = 1.0 / a                              # rho ratio that gives 1.0 %
    p = 1.0 - (need - 1.0) / (0.375 * MFP / d)
    sigma = LAMBDA_F / (4 * math.pi) * math.sqrt(-math.log(p))
    print(f"3. what 1.0 % at {d} nm would take")
    print(f"   rho/bulk <= {need:.2f}  ->  specularity p >= {p:.2f}")
    print(f"   Soffer: p = exp[-(4 pi sigma/lambda_F)^2]  ->  sigma <= {sigma:.3f} nm RMS")
    print("   on both interfaces. A silver atomic step is 0.24 nm. An amorphous")
    print("   organic surface is rougher than that by construction.")


if __name__ == "__main__":
    main()
