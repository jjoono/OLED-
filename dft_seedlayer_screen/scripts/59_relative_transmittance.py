"""How does an electrode transmit MORE than the bare substrate? Reproducing the DMD trick.

THE CLAIM BEING ANALYSED. Zhao et al. (L. J. Guo group), Nat. Commun. 11, 3367
(2020): Al2O3 / Cu-doped Ag 6.5 nm / ZnO on a polymer substrate gives 88.4 %
absolute transmittance against 88.1 % for the bare substrate -- relative
transmittance 100.3 %.

WHY THAT IS NOT A PARADOX. A bare substrate is not a perfect window: it reflects
at BOTH of its surfaces, roughly 4-6 % each for n = 1.5-1.65, so it starts around
88-92 % and the missing light is reflected, not absorbed. A dielectric-metal-
dielectric stack is a conductive ANTIREFLECTION coating: the two dielectrics are
chosen so their reflections interfere destructively with the metal's, killing the
front-surface reflection. If the reflection it removes exceeds the absorption the
metal adds, the coated substrate passes more light than the bare one. Nothing is
gained beyond unity absolute -- 88.4 % is still well under 100 % -- the ratio just
uses a denominator that was never 100 either.

WHY THIS MATTERS HERE. scripts/58 found HATCN already does a weak version of this:
at Ag 3 nm it transmits 0.79 %p MORE than no seed at all, because n = 1.95 sits
between glass and silver and partially index-matches. That is one dielectric on
one side. The Guo structure is the same physics carried to its conclusion with a
dielectric on BOTH sides and thicknesses tuned. So the question this script asks
is what a capping layer would buy on top of the measured HATCN/Ag stack.

SUBSTRATE HANDLING. The substrate is thick and its two surfaces do not interfere,
so it is treated incoherently: the coherent stack is solved into a semi-infinite
substrate, then combined with the back Fresnel surface allowing for multiple
incoherent bounces. Treating the substrate coherently would produce fringes that
no measurement sees and would make the comparison meaningless.

DIELECTRIC INDICES are nominal visible-range values (Al2O3 1.65, ZnO 2.00,
non-absorbing). They are adequate for showing the mechanism and reproducing the
~100 % figure, not for redesigning someone else's stack. Cu-doped Ag is
approximated by Ag; the doping shifts n,k somewhat but not the mechanism.
"""
import numpy as np

WL = np.arange(400.0, 701.0, 2.0)
JC_WL = np.array([397.4, 413.3, 430.5, 450.9, 471.4, 495.9,
                  521.0, 548.6, 582.1, 616.8, 659.5, 704.5])
AG_N = np.array([0.05, 0.05, 0.04, 0.04, 0.05, 0.05, 0.05, 0.06, 0.05, 0.06, 0.05, 0.04])
AG_K = np.array([2.07, 2.21, 2.36, 2.66, 2.83, 3.09, 3.34, 3.59, 3.93, 4.15, 4.48, 4.84])

N_AIR = 1.0
N_POLY = 1.65        # polymer substrate, Guo stack
N_GLASS = 1.52
N_AL2O3, N_ZNO = 1.65, 2.00
N_HATCN = 1.95
K_HATCN = 2.0e-4


def n_ag():
    return np.conj(np.interp(WL, JC_WL, AG_N) + 1j * np.interp(WL, JC_WL, AG_K))


def coherent(ns, ds, lam):
    """T, R of a coherent stack between ns[0] and ns[-1] (both semi-infinite)."""
    M = np.eye(2, dtype=complex)
    for N, d in zip(ns[1:-1], ds[1:-1]):
        delta = 2 * np.pi / lam * N * d
        c, s = np.cos(delta), np.sin(delta)
        M = M @ np.array([[c, 1j * s / N], [1j * N * s, c]])
    n0, nsub = ns[0], ns[-1]
    B, C = M @ np.array([1.0, nsub], dtype=complex)
    T = float(4 * n0.real * nsub.real / abs(n0 * B + C) ** 2)
    R = float(abs((n0 * B - C) / (n0 * B + C)) ** 2)
    return T, R


def coherent_rev(ns, ds, lam):
    """R seen from the substrate side -- needed for the incoherent bounce sum."""
    return coherent(ns[::-1], ds[::-1], lam)[1]


def stack_on_substrate(layers, n_sub, lam, i):
    """Coherent stack on the front, thick substrate behind, air behind that.

    `i` indexes the wavelength: each layer carries a full n(lambda) array and the
    matrix algebra needs the scalar at this wavelength. Passing the array through
    silently builds a 151-element "index" and matmul fails -- which is how this
    was caught, but a shape that happened to broadcast would not have been.
    """
    ns = [complex(N_AIR)] + [complex(l[0][i]) for l in layers] + [complex(n_sub)]
    ds = [0.0] + [l[1] for l in layers] + [0.0]
    T_in, _ = coherent(ns, ds, lam)
    R_back_of_stack = coherent_rev(ns, ds, lam)
    # substrate / air rear surface
    r = (n_sub - N_AIR) / (n_sub + N_AIR)
    R_rear, T_rear = r ** 2, 1 - r ** 2
    return T_in * T_rear / (1 - R_back_of_stack * R_rear)


def bare_substrate(n_sub, lam):
    r = (N_AIR - n_sub) / (N_AIR + n_sub)
    R = r ** 2
    return (1 - R) ** 2 / (1 - R ** 2)


def photopic(a):
    v = 1.019 * np.exp(-285.4 * ((WL / 1000.0) - 0.559) ** 2)
    return float((a * v).sum() / v.sum())


def evaluate(name, layers, n_sub):
    T = np.array([stack_on_substrate(layers, n_sub, l, i)
                  for i, l in enumerate(WL)])
    T0 = np.array([bare_substrate(n_sub, l) for l in WL])
    return name, photopic(T), photopic(T0), photopic(T) / photopic(T0) * 100


def const(n, k=0.0):
    return np.conj(np.full_like(WL, n + 1j * k, dtype=complex))


def main():
    ag = n_ag()

    print("=" * 74)
    print("PART 1 -- reproducing the Guo-type DMD on a polymer substrate")
    print("=" * 74)
    print(f"{'stack':<42}{'T abs':>9}{'T sub':>9}{'T rel':>9}")
    rows = [
        evaluate("bare polymer", [], N_POLY),
        evaluate("Ag 6.5 nm only", [(ag, 6.5)], N_POLY),
        evaluate("ZnO 40 / Ag 6.5 / Al2O3 40",
                 [(const(N_ZNO), 40.0), (ag, 6.5), (const(N_AL2O3), 40.0)], N_POLY),
    ]
    # scan the two dielectric thicknesses for the best relative T
    best = None
    for d_top in range(20, 71, 2):
        for d_bot in range(20, 71, 2):
            r = evaluate("", [(const(N_ZNO), float(d_top)), (ag, 6.5),
                              (const(N_AL2O3), float(d_bot))], N_POLY)
            if best is None or r[3] > best[0]:
                best = (r[3], d_top, d_bot, r[1])
    for name, T, T0, rel in rows:
        print(f"{name:<42}{100*T:>8.2f}%{100*T0:>8.2f}%{rel:>8.1f}%")
    print(f"{'optimised ZnO %d / Ag 6.5 / Al2O3 %d' % (best[1], best[2]):<42}"
          f"{100*best[3]:>8.2f}%{100*rows[0][2]:>8.2f}%{best[0]:>8.1f}%")
    print(f"\n  reported: 88.4 % absolute, 88.1 % substrate, 100.3 % relative")

    print("\n" + "=" * 74)
    print("PART 2 -- what a cap would buy on the MEASURED HATCN/Ag stack (glass)")
    print("=" * 74)
    hat = const(N_HATCN, K_HATCN)
    for d_ag in (3.0, 5.0):
        print(f"\n  Ag {d_ag:.0f} nm")
        print(f"    {'stack':<40}{'T abs':>9}{'T rel':>9}")
        cases = [
            ("bare glass", []),
            (f"HATCN 4 / Ag {d_ag:.0f}  (as measured)", [(hat, 4.0), (ag, d_ag)]),
        ]
        for dcap in (30.0, 40.0, 50.0, 60.0):
            cases.append((f"cap n=2.0 {dcap:.0f} / HATCN 4 / Ag {d_ag:.0f}",
                          [(const(2.0), dcap), (hat, 4.0), (ag, d_ag)]))
        best_c = None
        for nm, layers in cases:
            _, T, T0, rel = evaluate(nm, layers, N_GLASS)
            print(f"    {nm:<40}{100*T:>8.2f}%{rel:>8.1f}%")
            if nm.startswith("cap") and (best_c is None or rel > best_c[0]):
                best_c = (rel, nm)
        # finer optimisation over cap index and thickness
        bo = None
        for ncap in (1.8, 2.0, 2.2, 2.4):
            for dcap in range(20, 81, 2):
                _, T, T0, rel = evaluate("", [(const(ncap), float(dcap)),
                                              (hat, 4.0), (ag, d_ag)], N_GLASS)
                if bo is None or rel > bo[0]:
                    bo = (rel, ncap, dcap, T)
        print(f"    {'-> optimum: cap n=%.1f, %d nm' % (bo[1], bo[2]):<40}"
              f"{100*bo[3]:>8.2f}%{bo[0]:>8.1f}%")

    print("\n" + "=" * 74)
    print("WHAT THIS MEANS FOR THE PAPER")
    print("=" * 74)
    print("""
  Relative transmittance above 100 % is an antireflection result, not an
  amplification one, and it is available to any dielectric/metal/dielectric
  stack whose layers are tuned. It is not evidence that a particular metal or
  seed is special -- the Guo stack uses a Cu-doped Ag film 6.5 nm thick, which
  is thicker than the 3 nm measured here.

  So do not frame a relative-transmittance number as competing with theirs. The
  claim that survives comparison is the one scripts/58 quantifies: at the SAME
  Ag thickness an organic seed costs no absorption while a metal seed costs
  5-9 %p. That is a statement about the seed, and it is independent of how well
  the surrounding dielectrics happen to be tuned.

  The actionable part is in PART 2: the measured stack already has one
  index-matching layer underneath, and adding a tuned cap on top is the cheap
  route to a high relative transmittance in a device that needs a capping layer
  anyway.""")


if __name__ == "__main__":
    main()
